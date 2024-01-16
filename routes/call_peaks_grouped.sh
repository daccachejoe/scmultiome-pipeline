#!/bin/bash


module load macs2
module load r
scripts/multiome-processing.R \
    callpeaks \
    data/samplesheet.csv \
    -g ct \
    --project-prefix $project_prefix-grouped-peaks \
    -R output/RDS-files/$project_prefix-annotated-obj-list.RDS \
    -m $my_macs_path 