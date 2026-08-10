#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

scripts/multiome-processing.R \
        linkpeaks \
        configs/samplesheet.csv \
        --project_prefix $project_prefix-improved-clust-filtered-full- \
        --grouping.var cell.group \
        -R output/RDS-files/multiome-control-skin-grouped-peaks-callpeaks-obj-list.RDS
