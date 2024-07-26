#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load r/4.1.2
    module load macs2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.1.0
    module load macs/2.1.0
else
    echo "No job scheduler available to submit job: $script"
fi

scripts/multiome-processing.R \
    callpeaks \
    data/samplesheet.csv \
    -g ct \
    --project_prefix $project_prefix-grouped-peaks \
    -R output/RDS-files/$project_prefix-annotated-obj-list.RDS \
    -m $my_macs_path 