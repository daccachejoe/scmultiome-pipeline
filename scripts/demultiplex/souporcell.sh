#!/bin/bash
# Demultiplex: souporcell on one pooled library's GEX BAM (no reference VCF; freebayes calls).
# Branch: demultiplex (routes/00b_demultiplex.sh submits one job per library). Run from the
# project root; re-reads config/pipeline.config so it doesn't depend on env vars reaching
# the job.
#
# GEX, not ATAC: souporcell on the ATAC BAM gave clusters unrelated to donors (PSO v0
# cross-tab ~50/50), so ATAC is only used in the pooled rescue (rescue.sh).
# The souporcell container only works when every input sits inside its bind directory
# (bind-mounting the originals failed), so the GEX BAM, reference and image are copied
# into a work dir under demux_tmp_dir, and the work dir is removed once the final outputs
# are copied out.
#
# Inputs:  <cellranger_outs>/{gex_possorted_bam.bam(.bai), filtered_feature_bc_matrix/barcodes.tsv.gz}
#          (or the GEX BAM given as the 5th argument), genome_fasta(.fai), souporcell_sif
# Outputs: output/demultiplex/<sample>/souporcell/{clusters.tsv, cluster_genotypes.vcf, ambient_rna.txt}
# Usage:   bash scripts/demultiplex/souporcell.sh <sampleName> <cellranger_outs> <n_donors> <threads> [gex_bam]
set -euo pipefail
set -a; source config/pipeline.config; set +a

sample=$1; outs=$2; k=$3; threads=$4
gex_bam=${5:-$outs/gex_possorted_bam.bam}
out_dir=output/demultiplex/$sample/souporcell
work=${demux_tmp_dir:-output/demultiplex}/demux-work-${project_prefix}-${sample}-souporcell

for f in "$gex_bam" "$gex_bam.bai" "$outs/filtered_feature_bc_matrix/barcodes.tsv.gz" "$genome_fasta" "$genome_fasta.fai" "$souporcell_sif"; do
    [ -f "$f" ] || { echo "souporcell $sample: missing input $f" >&2; exit 1; }
done
if [ -s "$out_dir/clusters.tsv" ] && [ "${DEMUX_FORCE:-0}" != 1 ]; then
    echo "souporcell $sample: $out_dir/clusters.tsv exists; set DEMUX_FORCE=1 to rerun" >&2; exit 1
fi
module load "${singularity_module:-singularity}"

echo "== $(date) souporcell $sample: k=$k, $threads threads, work dir $work"
mkdir -p "$work" "$out_dir"
cp "$souporcell_sif" "$work/souporcell.sif"
cp "$genome_fasta" "$work/genome.fa"
cp "$genome_fasta.fai" "$work/genome.fa.fai"
cp "$gex_bam" "$work/gex_possorted_bam.bam"
cp "$gex_bam.bai" "$work/gex_possorted_bam.bam.bai"
cp "$outs/filtered_feature_bc_matrix/barcodes.tsv.gz" "$work/barcodes.tsv.gz"

work_abs=$(cd "$work" && pwd)
(cd "$work" && singularity exec --cleanenv --bind "$work_abs" souporcell.sif souporcell_pipeline.py \
    -i gex_possorted_bam.bam -b barcodes.tsv.gz -t "$threads" -f genome.fa -k "$k" -o out_gex)

for f in clusters.tsv cluster_genotypes.vcf ambient_rna.txt; do
    [ -s "$work/out_gex/$f" ] || { echo "souporcell $sample: no $f -- leaving $work for inspection" >&2; exit 1; }
    cp "$work/out_gex/$f" "$out_dir/"
done
rm -rf "$work"
echo "== $(date) souporcell $sample done: $(tail -n +2 "$out_dir/clusters.tsv" | cut -f2 | sort | uniq -c | tr '\n' ' ')"
