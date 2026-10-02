#!/bin/bash

# seurat preprocessing chunk
if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
    module load macs2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

# pc
# conda activate $conda_env_name

# genotype-pooled libraries: run the optional stage 00b (run/runmultiome
# demultiplex) first; doublets reads its calls via configs/demultiplexing_paths.csv
scripts/seurat_signac_pipeline.R \
    init,create,ambient,doublets,callpeaks,qc \
    configs/samplesheet.csv \
    --project_prefix $project_prefix \
    -m $my_macs_path



