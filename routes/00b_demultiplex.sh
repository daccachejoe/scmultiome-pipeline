#!/bin/bash
# Stage 00b (optional): genotype demultiplexing of donor-pooled libraries.
# Runs BEFORE stage 01, from Cell Ranger outputs only, because stage 01's doublets
# step uses the genotype doublets and stage 02's donors step writes patient IDs.
# Only for projects that pool several people per 10x library; single-donor
# projects skip it.
#
# Orchestrator (same pattern as 06_linkpeaks.sh): this route is a lightweight
# job that submits, per library in configs/demultiplexing_paths.csv,
#   souporcell  scripts/demultiplex/souporcell.sh   (skipped when demux_path is set:
#               souporcell was run elsewhere, its clusters.tsv is used as is)
#   rescue      scripts/demultiplex/rescue.sh       (after that library's souporcell)
# and then one
#   match       scripts/demultiplex/06_match_genotypes.py  (after every library)
# which writes output/demultiplex/donor-map-suggested.csv. Copy it to
# configs/donor_map.csv and replace the placeholder IDs with patient IDs.
# Jobs are chained by job ID, not name: done(<name>) never resolves if an older
# job with that name exited.
#
# configs/demultiplexing_paths.csv columns:
#   sampleName      must match configs/samplesheet.csv (its path = the Cell Ranger outs dir)
#   demux_path      blank -> run souporcell here; else an existing souporcell clusters.tsv
#   n_donors        people pooled in the library (souporcell -k); needed when demux_path is blank
#   min_rna_margin  optional; restrict calibration truth to confident souporcell calls (default 0)
#   gex_bam         optional; GEX BAM path, or "none" for an ATAC-only rescue
#                   (default <path>/gex_possorted_bam.bam; missing -> ATAC-only)
# Per-run overrides (env vars, not config keys):
#   DEMUX_STEPS=souporcell,rescue,match   which steps to submit
#   DEMUX_SAMPLES=PSO_L,PSO_NL            only these libraries
#   DEMUX_CORES=16                         cores per souporcell/rescue job
#   DEMUX_FORCE=1                          rerun libraries that already have outputs
# e.g. DEMUX_STEPS=rescue,match DEMUX_SAMPLES=AD_Lesional run/runmultiome demultiplex

set -euo pipefail
SHEET=configs/demultiplexing_paths.csv
STEPS=",${DEMUX_STEPS:-souporcell,rescue,match},"
CORES_CHILD=${DEMUX_CORES:-16}
MEM_PER_CORE=4000
mkdir -p output/demultiplex scripts/outs

if [ ! -f "$SHEET" ] || [ "$(tail -n +2 "$SHEET" | grep -c '[^[:space:]]')" -eq 0 ]; then
    echo "$SHEET has no libraries; nothing to demultiplex" >&2; exit 1
fi

# column value by header name (plain CSV, no quoted commas)
col() {  # <csv line> <column name> <header line>
    awk -F, -v line="$1" -v name="$2" -v header="$3" 'BEGIN{
        n = split(header, h, ","); split(line, v, ",");
        for (i = 1; i <= n; i++) if (h[i] == name) { print v[i]; exit } }'
}
header=$(head -1 "$SHEET" | tr -d '\r')
ss_header=$(head -1 configs/samplesheet.csv | tr -d '\r')

# submit <name> <hours> <cores> <deps (space-separated job IDs, may be empty)> <command...>
# prints the job ID (or "inline")
submit() {
    local name=$1 hours=$2 cores=$3 deps=$4; shift 4
    local log="scripts/outs/demux-${name}-%J.out"
    if [ "$SCHEDULER" == "slurm" ]; then
        local dep=""; [ -n "$deps" ] && dep="--dependency=afterok:$(echo $deps | tr ' ' ':')"
        sbatch --parsable -c "$cores" -p "$slurm_partition" --mem-per-cpu $MEM_PER_CORE -o "$log" \
            -t "$hours:00:00" --export=ALL $dep --wrap="$*"
    elif [ "$SCHEDULER" == "lsf" ]; then
        local dep=()
        [ -n "$deps" ] && dep=(-w "$(for j in $deps; do printf 'done(%s) && ' "$j"; done | sed 's/ && $//')")
        bsub -P "$lsf_project" -q "$lsf_queue" -n "$cores" -R "rusage[mem=$MEM_PER_CORE]" -R span[hosts=1] \
            -o "$log" -W "$hours:00" "${dep[@]}" "$*" | grep -oE '<[0-9]+>' | head -1 | tr -d '<>'
    else
        echo "No job scheduler available; running $name inline" >&2
        bash -c "$*" >&2
        echo inline
    fi
}

all_jobs=""
while IFS= read -r line; do
    line=$(echo "$line" | tr -d '\r')
    [ -z "$(echo "$line" | tr -d '[:space:],')" ] && continue
    sample=$(col "$line" sampleName "$header")
    if [ -n "${DEMUX_SAMPLES:-}" ] && [[ ",$DEMUX_SAMPLES," != *",$sample,"* ]]; then continue; fi
    demux_path=$(col "$line" demux_path "$header")
    n_donors=$(col "$line" n_donors "$header")
    margin=$(col "$line" min_rna_margin "$header"); margin=${margin:-0}
    gex_bam=$(col "$line" gex_bam "$header")
    ss_line=$(awk -F, -v s="$sample" '$1 == s' configs/samplesheet.csv | head -1 | tr -d '\r')
    [ -z "$ss_line" ] && { echo "$sample is in $SHEET but not in configs/samplesheet.csv" >&2; exit 1; }
    outs=$(col "$ss_line" path "$ss_header")

    sp_job=""
    if [ -z "$demux_path" ]; then
        clusters=output/demultiplex/$sample/souporcell/clusters.tsv
        if [[ "$STEPS" == *",souporcell,"* ]]; then
            [ -z "$n_donors" ] && { echo "$sample: n_donors is required to run souporcell ($SHEET)" >&2; exit 1; }
            sp_job=$(submit "souporcell-$sample" 24 "$CORES_CHILD" "" \
                bash scripts/demultiplex/souporcell.sh "$sample" "$outs" "$n_donors" "$CORES_CHILD" ${gex_bam:+"$gex_bam"})
            echo "$sample: souporcell job $sp_job"
            [ "$sp_job" != inline ] && all_jobs="$all_jobs $sp_job"
        fi
    else
        clusters=$demux_path
        echo "$sample: using existing souporcell calls $demux_path"
    fi

    if [[ "$STEPS" == *",rescue,"* ]]; then
        deps=""; [ -n "$sp_job" ] && [ "$sp_job" != inline ] && deps=$sp_job
        rs_job=$(submit "rescue-$sample" 12 "$CORES_CHILD" "$deps" \
            bash scripts/demultiplex/rescue.sh "$sample" "$outs" "$clusters" "$CORES_CHILD" "$margin" ${gex_bam:+"$gex_bam"})
        echo "$sample: rescue job $rs_job${deps:+ (after $deps)}"
        [ "$rs_job" != inline ] && all_jobs="$all_jobs $rs_job"
    fi
done < <(tail -n +2 "$SHEET")

if [[ "$STEPS" == *",match,"* ]]; then
    m_job=$(submit "match" 2 1 "$(echo $all_jobs)" \
        "source $personal_anaconda_path && conda activate ${demux_env_name:-demux-env} && \
         python scripts/demultiplex/06_match_genotypes.py --demux $SHEET --samplesheet configs/samplesheet.csv \
             --out_dir output/demultiplex")
    echo "match job $m_job${all_jobs:+ (after$all_jobs)}"
fi
