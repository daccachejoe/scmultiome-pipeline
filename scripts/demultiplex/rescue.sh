#!/bin/bash
# Demultiplex: pooled-genotype rescue for one library, then the combined per-cell calls.
# Branch: demultiplex (routes/00b_demultiplex.sh submits one job per library, after its
# souporcell job). Run from the project root; re-reads config/pipeline.config.
#
# souporcell leaves many low-RNA nuclei unassigned (13-38% of barcodes in
# prelim-long-data). This builds each donor's genotype from souporcell's singlets, pooling
# their ATAC and GEX reads, and calls every cell against those genotypes:
#   1  01_split_by_donor.py  split the BAMs by singlet donor and fold A/B (pysam, BAMs read in place)
#   -  bcftools mpileup/call on the per-donor pools (fold A, and all singlets), per chromosome
#   2  02_soft_sites.py      donor-informative SNPs from pooled allele fractions
#   -  vartrix --umi          per-cell GEX allele counts
#   3  03_count_fragments.py per-cell ATAC allele counts, once per fragment
#   4  04_calibrate.py       fit the ATAC weight on held-out fold-B singlets; call all cells
#   5  05_combine_calls.py   souporcell + pooled calls -> combined_clusters.tsv (thin-donor guard)
# No GEX BAM (e.g. healthy-human cntrl.2) -> ATAC-only: no GEX split, pools or vartrix.
# bcftools/samtools/vartrix run in the souporcell container, which needs its inputs inside
# the bind dir, so the pools, reference and GEX BAM are copied into a work dir under
# demux_tmp_dir; the ATAC BAM is only read by pysam, in place.
#
# Inputs:  <cellranger_outs>/{atac_possorted_bam.bam, gex_possorted_bam.bam (optional),
#          filtered_feature_bc_matrix/barcodes.tsv.gz}; souporcell clusters.tsv
# Outputs: output/demultiplex/<sample>/pooled/ (calls, held-out reports, soft sites, counts)
#          output/demultiplex/<sample>/combined_clusters.tsv
# Usage:   bash scripts/demultiplex/rescue.sh <sampleName> <cellranger_outs> <souporcell_clusters.tsv> \
#              <threads> [min_rna_margin] [gex_bam|none]
set -euo pipefail
set -a; source config/pipeline.config; set +a

sample=$1; outs=$2; clusters=$3; threads=$4; margin=${5:-0}
gex_bam=${6:-$outs/gex_possorted_bam.bam}
here=$(cd "$(dirname "$0")" && pwd)
out_dir=output/demultiplex/$sample
work=${demux_tmp_dir:-output/demultiplex}/demux-work-${project_prefix}-${sample}-pooled
atac_bam=$outs/atac_possorted_bam.bam
case "${species:-human}" in
    human) chroms="$(printf 'chr%s ' $(seq 1 22) X)" ;;
    mouse) chroms="$(printf 'chr%s ' $(seq 1 19) X)" ;;
    *) echo "unsupported species: $species" >&2; exit 1 ;;
esac

for f in "$atac_bam" "$atac_bam.bai" "$clusters" "$outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
         "$genome_fasta" "$genome_fasta.fai" "$souporcell_sif"; do
    [ -f "$f" ] || { echo "rescue $sample: missing input $f" >&2; exit 1; }
done
modalities="atac gex"
if [ "$gex_bam" = none ] || [ ! -f "$gex_bam" ]; then
    echo "rescue $sample: no GEX BAM ($gex_bam) -- ATAC-only rescue"
    modalities="atac"
fi
if [ -s "$out_dir/combined_clusters.tsv" ] && [ "${DEMUX_FORCE:-0}" != 1 ]; then
    echo "rescue $sample: $out_dir/combined_clusters.tsv exists; set DEMUX_FORCE=1 to rerun" >&2; exit 1
fi
module load "${singularity_module:-singularity}"
source "$personal_anaconda_path"
conda activate "${demux_env_name:-demux-env}"

mkdir -p "$work" "$out_dir/pooled"
clusters_abs=$(cd "$(dirname "$clusters")" && pwd)/$(basename "$clusters")
atac_abs=$(cd "$(dirname "$atac_bam")" && pwd)/$(basename "$atac_bam")
[ "$modalities" = "atac gex" ] && gex_abs=$(cd "$(dirname "$gex_bam")" && pwd)/$(basename "$gex_bam")
zcat "$outs/filtered_feature_bc_matrix/barcodes.tsv.gz" > "$work/barcodes.tsv"
cp "$clusters" "$work/clusters_gex.tsv"
cd "$work"

echo "== $(date) $sample: splitting BAMs by donor and fold ($modalities)"
bam_args="--bam atac=$atac_abs"
[ "$modalities" = "atac gex" ] && bam_args="$bam_args --bam gex=$gex_abs"
python "$here/01_split_by_donor.py" --clusters clusters_gex.tsv --out_dir . --species "${species:-human}" \
    --threads "$threads" $bam_args

cp "$souporcell_sif" souporcell.sif
cp "$genome_fasta" genome.fa
cp "$genome_fasta.fai" genome.fa.fai
work_abs=$(pwd)
run() { singularity exec --cleanenv --bind "$work_abs" souporcell.sif "$@"; }
donors=$(tail -n +2 singlet_folds.tsv | cut -f3 | sort -u)
echo "donors: $donors"
for fold in A B; do for d in $donors; do for m in $modalities; do
    run samtools cat -o pool_${fold}_${d}_${m}.bam $(for c in $chroms; do echo split/${fold}_${d}_${m}.$c.bam; done)
    run samtools index pool_${fold}_${d}_${m}.bam
done; done; done
rm -rf split

echo "== $(date) $sample: calling donor genotypes on the pools"
mkdir -p calls
call_chrom() {  # <set> <chrom> <bams...>
    local set=$1 chrom=$2; shift 2
    singularity exec --cleanenv --bind "$work_abs" souporcell.sif bash -c \
        "bcftools mpileup -f genome.fa -r $chrom -a FORMAT/AD,FORMAT/DP -q 30 -Q 20 -d 10000 -Ou $* 2>/dev/null |
         bcftools call -m -v -Oz -o calls/$set.$chrom.vcf.gz"
}
export -f call_chrom; export work_abs
for set in foldA all; do
    if [ $set = foldA ]; then bams=$(ls pool_A_*.bam); else bams=$(ls pool_*.bam); fi
    printf '%s\n' $chroms | xargs -P "$threads" -I{} bash -c "call_chrom $set {} $(echo $bams)"
    run bcftools concat -Oz -o calls/$set.vcf.gz $(for c in $chroms; do echo calls/$set.$c.vcf.gz; done)
    rm calls/$set.chr*.vcf.gz
    python "$here/02_soft_sites.py" calls/$set.vcf.gz soft_$set
done
rm -f pool_*.bam pool_*.bam.bai

if [ "$modalities" = "atac gex" ]; then
    echo "== $(date) $sample: counting GEX alleles per cell (vartrix, UMIs)"
    cp "$gex_abs" gex_possorted_bam.bam
    cp "$gex_abs.bai" gex_possorted_bam.bam.bai
    for set in foldA all; do
        mkdir -p vartrix_soft_${set}_gex
        run vartrix --bam gex_possorted_bam.bam --cell-barcodes barcodes.tsv --fasta genome.fa \
            --vcf soft_$set.vcf --out-matrix vartrix_soft_${set}_gex/alt.mtx --ref-matrix vartrix_soft_${set}_gex/ref.mtx \
            --scoring-method coverage --mapq 30 --threads "$threads" --umi
    done
    rm -f gex_possorted_bam.bam gex_possorted_bam.bam.bai
fi
rm -f souporcell.sif genome.fa genome.fa.fai

echo "== $(date) $sample: counting ATAC alleles per fragment"
python "$here/03_count_fragments.py" --run_dir . --bam "$atac_abs" --threads "$threads"
echo "== $(date) $sample: calibrating the ATAC weight on held-out singlets"
python "$here/04_calibrate.py" --run_dir . --clusters clusters_gex.tsv --min_rna_margin "$margin" \
    --min_posterior "${demux_min_posterior:-0.95}" --ambient "${demux_ambient:-0.10}" \
    --doublet_prior "${demux_doublet_prior:-0.05}"
cd - > /dev/null

tag=fragments; [ "$(awk "BEGIN{print ($margin > 0)}")" = 1 ] && tag=fragments_rnamargin$margin
cp -r "$work"/. "$out_dir/pooled/"
rm -rf "$work"
echo "== $(date) $sample: combining souporcell and pooled calls"
python "$here/05_combine_calls.py" --clusters "$clusters_abs" --out "$out_dir/combined_clusters.tsv" \
    --calls "$out_dir/pooled/calls_soft_all_calibrated_$tag.tsv" \
    --per_donor "$out_dir/pooled/heldout_per_donor_$tag.tsv" \
    --min_donor_precision "${demux_min_donor_precision:-0.98}"
echo "== $(date) $sample: done"
