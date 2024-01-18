#!/bin/bash

export my_macs_path=/gpfs/data/naiklab/jd5457/opt/MACS/bin/macs3
export project_prefix=multiome-control-skin
export personal_anaconda_path=/gpfs/data/ruggleslab/jd5457/miniconda/miniconda/etc/profile.d/conda.sh
export sceasy_env_name=sceasy-joe
export UCD_env_name=ucdenv
export scenicplus_env_name=scenicplus
export create_cistarget_databases_env_name=create_cistarget_databases
export create_cistarget_databases_dir=/gpfs/home/jd5457/.conda/envs/create_cisTarget_databases/

init() {
    echo "Initializing environment"
    bash routes/setup_dirs.sh
    bash routes/install.sh
}

seurat_preprocess() {
    echo "Running seurat preprocessing"
    sbatch -c 32 -p a100_short --mem 300GB -o scripts/outs/preprocessing-%J.out -t 0-23:00:00 --export=ALL routes/seurat_preprocess.sh
}

filter_and_cluster() {
    echo "Running filtering and clustering"
    sbatch -c 32 -p a100_short --mem 300GB -o scripts/outs/filter-and-cluster-%J.out -t 0-23:00:00 --export=ALL  routes/filter_and_cluster.sh
}

identify_celltypes() {
    echo "Helping identify cell types"
    sbatch -c 32 -p a100_short --mem 300GB -o scripts/outs/identify-celltypes-%J.out -t 0-23:00:00 --export=ALL  routes/identify_celltypes.sh
}

label_celltypes() {
    echo "Labelling cell types with user inputs"
    sbatch -c 32 -p a100_short --mem 300GB -o scripts/outs/label-celltypes-%J.out -t 0-23:00:00 --export=ALL  routes/label_celltypes.sh
}

call_peaks_grouped() {
    echo "Calling ATAC-Seq peaks grouped by cell type"
    sbatch -c 32 -p a100_short --mem 300GB -o scripts/outs/callpeaks-grouped-%J.out -t 0-23:00:00 --export=ALL  routes/call_peaks_grouped.sh
}

run_scenicplus() {
    echo "Running SCENIC+"
    sbatch -p gpu8_medium --mem 300GB -o scripts/outs/scenicplus-%J.out -t 0-23:00:00 --export=ALL  routes/run_scenicplus.sh
}

run_multiome() {
    # Check if the argument is provided
    if [ $# -eq 0 ]; then
        echo "Please provide an argument to specify which function to run."
        exit 1
    fi

    # Run the specified function based on the argument
    "$1"
}

# # Call the run_multiome function with the provided argument
# run_multiome "$1"
