#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

# not needed at MSSM
# # install the required packages to make scmultiome executable in the desired directory
# mkdir .lib
# git clone https://bitbucket.org/djhshih/argparser.git .lib/argparser
# mv scripts/install.R .lib/argparser/
# cd .lib/argparser
# Rscript install.R
# R CMD INSTALL .

chmod +x scripts/seurat_signac_pipeline.R
