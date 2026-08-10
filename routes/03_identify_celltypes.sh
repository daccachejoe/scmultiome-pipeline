#!/bin/bash


if [ -s "configs/resolution_to_use.txt" ]; then
    resolution=$(cat configs/resolution_to_use.txt)
    echo Using resolution $resolution to identify cell types

    # override to point this stage at a different input, e.g.:
    # INPUT_RDS=output/RDS-files/my-variant-02-merge-obj-list.RDS run/runmultiome identify_celltypes
    INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-02-merge-obj-list.RDS}"

    source $personal_anaconda_path
    conda activate $sceasy_env_name

    Rscript scripts/03_convert_seurat_to_h5ad.R \
        "$INPUT_RDS" \
        ./output/ucd/$project_prefix-cluster-obj-list.h5ad \
        $sceasy_env_name
    
    conda deactivate 
    conda activate $UCD_env_name

    python scripts/03_ucd_deconvolve.py \
        --resolution $resolution \
        --input_file ./output/ucd/$project_prefix-cluster-obj-list.h5ad \
        --reference_file "$celltype_reference_h5ad"
    
    Rscript scripts/03_plot_ucd_results.R output/ucd/cellmetadata-unbiased.csv $resolution
    echo UCDeconvolve is complete. Exiting. 
else
    echo "resolution_to_use.txt does not exist. Exiting."
fi

