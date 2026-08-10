#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load r/4.1.2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.2.0
else
    echo "No job scheduler available to submit job: $script"
fi

scripts/multiome-processing.R \
        linkpeaks \
        configs/samplesheet.csv \
        --project_prefix $project_prefix-improved-clust-filtered-full- \
        --grouping.var cell.group \
        -R output/RDS-files/multiome-control-skin-grouped-peaks-callpeaks-obj-list.RDS
