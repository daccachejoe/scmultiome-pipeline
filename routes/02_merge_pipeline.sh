#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

# override to point this stage at a different input, e.g.:
# INPUT_RDS=output/RDS-files/my-variant-01-qc-obj-list.RDS run/runmultiome run_merged_pipeline
INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-01-qc-obj-list.RDS}"

# donors (patient IDs from the demultiplex branch + configs/donor_map.csv) runs
# whenever either config has data rows (init creates both with only a header),
# so a pooled project can't produce a merged object without donors
# a single-sample project has nothing to merge or batch-correct: it is
# clustered on its own (cluster stage: WNN + markers, no Harmony) instead
N_SAMPLES=$(tail -n +2 configs/samplesheet.csv | grep -c '[^[:space:]]')
if [ "$N_SAMPLES" -eq 1 ]; then
    LAST_STAGE=cluster
    HARMONY_FLAG=""
else
    LAST_STAGE=merge
    HARMONY_FLAG="--RunHarmony"
fi
STAGES=filter,$LAST_STAGE
for f in configs/demultiplexing_paths.csv configs/donor_map.csv; do
    if [ -f "$f" ] && [ "$(tail -n +2 "$f" | grep -c '[^[:space:]]')" -gt 0 ]; then
        STAGES=filter,donors,$LAST_STAGE
    fi
done

echo "Running merged pipeline: $STAGES"
scripts/seurat_signac_pipeline.R \
        $STAGES \
        configs/samplesheet.csv \
        --qc.sheet configs/qc_df.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        $HARMONY_FLAG \
        -R "$INPUT_RDS"