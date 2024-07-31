#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load r/4.1.2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.1.0
else
    echo "No job scheduler available to submit job: $script"
fi

scripts/multiome-processing.R \
        filter,merge \
        configs/samplesheet.csv \
        --qc.sheet configs/qc_df.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        --RunHarmony \
        -R output/RDS-files/$project_prefix-qc-obj-list.RDS