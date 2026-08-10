#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load r/4.1.2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.1.0
else
    echo "No job scheduler available to submit job: $script"
fi

if [[ $(wc -l < "configs/cluster_labels.csv") -gt 1 ]]; then
    echo Using annotating cell types with cluster_labels.csv
    resolution=$(cat configs/resolution_to_use.txt)

    Rscript scripts/label_celltypes.R \
        configs/cluster_labels.csv \
        $resolution \
        output/RDS-files/$project_prefix-merged-obj-list.RDS \
        $project_prefix
    echo Celltype labelling is complete. Exiting. 
else
    echo "cluster_labels.csv does not exist or is empty. Exiting."
fi
