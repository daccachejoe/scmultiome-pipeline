#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

echo Running merged pipeline
scripts/multiome-processing.R \
        filter,merge \
        configs/samplesheet.csv \
        --qc.sheet configs/qc_df.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        --RunHarmony \
        -R output/RDS-files/$project_prefix-qc-obj-list.RDS