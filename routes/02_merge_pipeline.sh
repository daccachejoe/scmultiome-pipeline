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

echo Running merged pipeline
scripts/seurat_signac_pipeline.R \
        filter,merge \
        configs/samplesheet.csv \
        --qc.sheet configs/qc_df.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        --RunHarmony \
        -R "$INPUT_RDS"