#!/bin/bash
module load r/4.2.2


if [[ $(wc -l < "configs/qc_df.csv") -gt 0 ]]; then
    scripts/multiome-processing.R \
        filter,cluster \
        configs/samplesheet.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        -R output/RDS-files/$project_prefix-qc-obj-list.RDS \
        --qc.sheet configs/qc_df.csv
else
    echo "qc_df.csv does not have more than one line. Please fill it out and try again."
fi


