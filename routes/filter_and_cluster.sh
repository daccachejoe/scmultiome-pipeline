#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

if [[ $(wc -l < "configs/qc_df.csv") -gt 0 ]]; then
    # scripts/multiome-processing.R \
    #     filter,cluster \
    #     configs/samplesheet.csv \
    #     --project_prefix $project_prefix \
    #     -m $my_macs_path \
    #     -R output/RDS-files/$project_prefix-qc-obj-list-improved.RDS \
    #     --qc.sheet configs/qc_df.csv
    scripts/multiome-processing.R \
        subcluster \
        configs/samplesheet.csv \
        --project_prefix $project_prefix-lineage-imrpoved \
        -m $my_macs_path \
        --RunHarmony \
        -R output/RDS-files/$project_prefix-grouped-peaks-callpeaks-obj-list.RDS \
        --qc.sheet configs/qc_df.csv    
else
    echo "qc_df.csv does not have more than one line. Please fill it out and try again."
fi


