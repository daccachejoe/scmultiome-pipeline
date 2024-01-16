#!/bin/bash

echo $my_macs_path
echo $project_prefix
# seurat preprocessing chunk
module load r

if [ -f "data/raw/souporcell/out_gex/barcode-assignment-df.csv" ]; then
    scripts/multiome-processing.R \
        init,create,callpeaks,qc \
        data/samplesheet.csv \
        --project-prefix $project_prefix \
        -m $my_macs_path \
        --SoupOrCellDF data/raw/souporcell/out_gex/barcode-assignment-df.csv 
else
    scripts/multiome-processing.R \
        init,create,callpeaks,qc \
        data/samplesheet.csv \
        --project-prefix $project_prefix \
        -m $my_macs_path
fi



