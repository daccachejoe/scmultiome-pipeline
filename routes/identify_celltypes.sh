#!/bin/bash


if [ -s "configs/resolution_to_use.txt" ]; then
    resolution=$(cat configs/resolution_to_use.txt)
    echo Using resolution $resolution to identify cell types

    source $personal_anaconda_path
    conda activate $sceasy_env_name
    
    Rscript scripts/convert_seurat_to_h5ad.R \
        output/RDS-files/$project_prefix-cluster-obj-list.RDS \
        ./output/ucd/$project_prefix-cluster-obj-list.h5ad \
        $sceasy_env_name
    
    conda deactivate 
    conda activate $UCD_env_name

    python scripts/ucd-script.py \
        --resolution $resolution \
        --input_file ./output/ucd/$project_prefix-cluster-obj-list.h5ad \
        --reference_file /gpfs/data/naiklab/SHARED_DATA/Haniffa_Healthy_Only_Updated_Labels_Unicell_Reference.h5ad 
    
    Rscript scripts/plotting-UCD-and-seurat.R /output/ucd/cellmetadata-unbiased.csv $resolution
    echo UCDeconvolve is complete. Exiting. 
else
    echo "resolution_to_use.txt does not exist. Exiting."
fi

