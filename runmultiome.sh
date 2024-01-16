#!/bin/bash

export my_macs_path=/gpfs/data/naiklab/jd5457/opt/MACS/bin/macs3
export project_prefix=multiome-control-skin
export personal_anaconda_path=/gpfs/data/ruggleslab/jd5457/miniconda/miniconda/etc/profile.d/conda.sh
export sceasy_env_name=sceasy-joe
export UCD_env_name=ucdenv
export scenicplus_env_name=scenicplus
export create_cistarget_databases_env_name=create_cistarget_databases
export create_cistarget_databases_dir=/gpfs/home/jd5457/.conda/envs/create_cisTarget_databases/


# Check if the argument is provided
if [ $# -eq 0 ]; then
    echo "Please provide an argument to specify which chunk of code to run."
    exit 1
fi

# Define the code chunks
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

# Run the specified chunk based on the argument
eval "$1"

#!/bin/sh
# high memory configuration
#SBATCH --partition gpu8_medium
##SBATCH --partition gpu4_short
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 32
##SBATCH --mem 600GB
#SBATCH --mem 300GB
#SBATCH --time 0-12:00:00
#SBATCH --job-name jup-SCENIC
#SBATCH --output scripts/outs/scenicplus-%J.log
#SBATCH --error scripts/outs/scenicplus-%J.e 