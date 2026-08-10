#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

if [[ $(wc -l < "configs/cluster_labels.csv") -gt 1 ]]; then
    echo Using annotating cell types with cluster_labels.csv
    resolution=$(cat configs/resolution_to_use.txt)

    # override to point this stage at a different input, e.g.:
    # INPUT_RDS=output/RDS-files/my-variant-02-merge-obj-list.RDS run/runmultiome label_celltypes
    INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-02-merge-obj-list.RDS}"

    Rscript scripts/04_label_celltypes.R \
        configs/cluster_labels.csv \
        $resolution \
        "$INPUT_RDS" \
        $project_prefix
    echo Celltype labelling is complete. Exiting. 
else
    echo "cluster_labels.csv does not exist or is empty. Exiting."
fi
