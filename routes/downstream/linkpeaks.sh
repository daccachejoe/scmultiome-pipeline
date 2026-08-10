#!/bin/bash
# Downstream branch: link ATAC peaks to nearby genes, in parallel, one job
# per group. Reads stage 05's output (cell types called, peaks re-called
# grouped by cell type) -- independent of the other downstream branches
# (subcluster, SCENIC+), which also branch from stage 05.
#
# Orchestrator pattern: this route runs as a lightweight job itself (queued
# via `run/runmultiome linkpeaks`). It (1) splits the input object by
# --grouping.var into one object per group, (2) submits one LinkPeaks job
# per group, then (3) submits a merge job that depends on every group job
# finishing, to combine the per-group results back into one output.

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

GROUPING_VAR=ct
OUT_PREFIX="${project_prefix}-improved-clust-filtered-full-"
# override to point this stage at a different input, e.g.:
# INPUT_RDS=output/RDS-files/my-variant-05-callpeaks-obj-list.RDS run/runmultiome linkpeaks
INPUT_RDS="${INPUT_RDS:-output/RDS-files/${project_prefix}-grouped-peaks-05-callpeaks-obj-list.RDS}"

echo "Splitting object by ${GROUPING_VAR} for parallel linkpeaks"
scripts/downstream/linkpeaks_split.R "$INPUT_RDS" "$GROUPING_VAR" "$OUT_PREFIX"

GROUP_LIST_FILE="output/RDS-files/${OUT_PREFIX}-linkpeaks-groups.txt"
if [ ! -s "$GROUP_LIST_FILE" ]; then
    echo "Splitting failed: $GROUP_LIST_FILE not found or empty. Exiting."
    exit 1
fi

JOB_IDS=()
while IFS= read -r group; do
    [ -z "$group" ] && continue
    GROUP_RDS="output/RDS-files/${OUT_PREFIX}-linkpeaks-group-${group}-obj.RDS"
    OUT_LOG="scripts/outs/linkpeaks-${group}-%J.out"
    if [ "$SCHEDULER" == "slurm" ]; then
        jid=$(sbatch --parsable -c 6 -p "$slurm_partition" --mem-per-cpu 12000 -o "$OUT_LOG" -t 0-12:00 --export=ALL \
            --wrap="scripts/downstream/linkpeaks_group.R $GROUP_RDS $group $OUT_PREFIX")
        JOB_IDS+=("$jid")
    elif [ "$SCHEDULER" == "lsf" ]; then
        jid=$(bsub -P "$lsf_project" -q "$lsf_queue" -n 6 -R "rusage[mem=12000]" -R span[hosts=1] -o "$OUT_LOG" -W 12:00 \
            "scripts/downstream/linkpeaks_group.R $GROUP_RDS $group $OUT_PREFIX" | grep -oE '[0-9]+' | head -1)
        JOB_IDS+=("$jid")
    else
        echo "No job scheduler available; running group $group inline"
        scripts/downstream/linkpeaks_group.R "$GROUP_RDS" "$group" "$OUT_PREFIX"
    fi
done < "$GROUP_LIST_FILE"

GROUP_CSV=$(paste -sd, "$GROUP_LIST_FILE")
MERGE_LOG="scripts/outs/linkpeaks-merge-%J.out"

if [ "$SCHEDULER" == "slurm" ] && [ ${#JOB_IDS[@]} -gt 0 ]; then
    DEP=$(IFS=:; echo "${JOB_IDS[*]}")
    sbatch -c 4 -p "$slurm_partition" --mem-per-cpu 8000 -o "$MERGE_LOG" -t 0-2:00 \
        --dependency=afterok:$DEP --export=ALL \
        --wrap="scripts/downstream/linkpeaks_merge.R $OUT_PREFIX $GROUP_CSV"
elif [ "$SCHEDULER" == "lsf" ] && [ ${#JOB_IDS[@]} -gt 0 ]; then
    DEP=""
    for jid in "${JOB_IDS[@]}"; do
        DEP="${DEP}done(${jid}) && "
    done
    DEP=${DEP% && }
    bsub -P "$lsf_project" -q "$lsf_queue" -n 4 -R "rusage[mem=8000]" -o "$MERGE_LOG" -W 2:00 -w "$DEP" \
        "scripts/downstream/linkpeaks_merge.R $OUT_PREFIX $GROUP_CSV"
else
    scripts/downstream/linkpeaks_merge.R "$OUT_PREFIX" "$GROUP_CSV"
fi
