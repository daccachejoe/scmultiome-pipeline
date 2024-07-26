#!/bin/bash

# seurat preprocessing chunk
if [ "$SCHEDULER" == "slurm" ]; then
    module load r/4.1.2
    module load macs2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load R/4.1.0
    module load macs/2.1.0
else
    echo "No job scheduler available to submit job: $script"
fi

# pc
# conda activate $conda_env_name

if [ -f "data/raw/souporcell/out_gex/barcode-assignment-df.csv" ]; then
    scripts/multiome-processing.R \
        init,create,callpeaks,qc \
        configs/samplesheet.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path \
        --SoupOrCellDF data/raw/souporcell/out_gex/barcode-assignment-df.csv 
else
    scripts/multiome-processing.R \
        init,create,callpeaks,qc \
        configs/samplesheet.csv \
        --project_prefix $project_prefix \
        -m $my_macs_path
fi



