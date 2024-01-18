#!/bin/bash

module load r/4.2.2

if [[ $(wc -l < "configs/cluster_labels.csv") -gt 1 ]]; then
    echo Using annotating cell types with cluster_labels.csv
    resolution=$(cat configs/resolution_to_use.txt)

    module load r
    Rscript scripts/label_celltypes.R \
        configs/cluster_labels.csv \
        $resolution \
        output/RDS-files/$project_prefix-cluster-obj-list.h5ad \
        $project_prefix
    echo Celltype labelling is complete. Exiting. 
else
    echo "cluster_labels.csv does not exist or is empty. Exiting."
fi
