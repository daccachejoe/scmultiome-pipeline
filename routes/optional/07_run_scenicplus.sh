#!/bin/bash

# first we convert the annotated object to h5ad    
source $personal_anaconda_path
conda activate $sceasy_env_name

Rscript scripts/03_convert_seurat_to_h5ad.R \
   output/RDS-files/$project_prefix-grouped-peaks-05-callpeaks-obj-list.RDS \
   data/scenicplus/$project_prefix.h5ad \
   $sceasy_env_name

# then we run a script in R to export the necessary data
Rscript scripts/optional/07a_export_scenicplus_data.R \
    output/RDS-files/$project_prefix-grouped-peaks-05-callpeaks-obj-list.RDS

conda deactivate
conda activate $scenicplus_env_name
# for some reason, sceasy saves the RNA matrix as a sparse matrix in the AnnData.raw slot
# to run SCENIC+, we need to convert this to a AnnData in the .raw slot
# in an interactive is easiest, but run the following lines of code
python scripts/optional/07b_reformat_anndata.py data/scenicplus/$project_prefix.h5ad data/scenicplus/${project_prefix}_new.h5ad
# import anndata as ad
# import scanpy as sc
# obj=sc.read_h5ad("path/to/file")
# obj.raw=ad.AnnData(obj.X)
# new_adata = obj.raw.to_adata()
# new_adata.var_names = obj.var_names #Or whatever your index is called
# new_adata.var_names.name = None
# new_adata.raw = new_adata
# sc.pp.normalize_total(new_adata, target_sum=1e4)
# sc.pp.log1p(new_adata)
# new_adata.write("rna_new.h5ad")

export _JAVA_OPTIONS=-Xmx250g
echo Running SCENIC+ preprocessing pipeline
# make sure the prepcocess config file has compute_topics, indentify_enhancers set to True but not runPycisTarget
config_file=configs/scenicplus-preprocess-config.yml
python scripts/optional/07c_scenicplus_pipeline.py $config_file
echo SCENIC+ preprocessing complete. Moving on to creating cisTarget databases

# deactivate the scenicplus environment and start up the create cisTargetDB env
conda deactivate
module load bedtools
module -q load dataark

echo "Generating fasta's of peaks for cisTarget databases"
REGION_BED="data/raw/macs-peaks/grouped-peaks.bed"
GENOME_FASTA="$genome_fasta"
CHROMSIZES="data/scenicplus/jd_chromsizes"
DATABASE_PREFIX=${project_prefix}_500_bg_padding
SCRIPT_DIR=${create_cistarget_databases_path}

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        data/scenicplus/cisTarget_dbs/pseudobulk_peaks.fa \
        500 \
        yes

# module purge
source $personal_anaconda_path
conda activate $scenicplus_env_name

# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=data/scenicplus/cisTarget_dbs/pseudobulk_peaks.fa

scenicplus_env_dir=$(conda info --base)/envs/${scenicplus_env_name}
# Directory with motifs in Cluster-Buster format. Make sure this data is downloaded
motifs_dir=${scenicplus_env_dir}/aertslab_motif_colleciton/v10nr_clust_public/singletons/
if [ -d "$motifs_dir" ]; then
    echo "Motifs directory exists: ${motifs_dir}"
else
    echo "Motifs directory does not exist: ${motifs_dir}"
    echo "Download the motifs directory using: "
    echo "motif_database_url='https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'"
    echo "wget "${motif_database_url}" --no-check-certificate"
    echo "wget -r https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/ --no-check-certificate "
    exit 1
fi

# File with motif IDs (base name of motif file in ${motifs_dir}).
# motifs_list_filename=data/scenic-plus/cisTarget_dbs/human-motif-names.txt # text file edited from all TFs hg38.txt in R # this did not work
# ls data/scenicplus/cisTarget_dbs/v10nr_clust_public/singletons/*.cb | sed 's/\.cb$//' > data/scenicplus/cisTarget_dbs/prefixes.txt
# motifs_list_filename=data/scenicplus/cisTarget_dbs/prefixes.txt # all .cb files in the v10-clust directory
motifs_list_filename=${scenicplus_env_dir}/motifs.txt
# cisTarget motif database output prefix.
db_prefix=data/scenicplus/$project_prefix
nbr_threads=32
CLUSTER_BUSTER_PATH=${scenicplus_env_dir}/cbust

"${create_cistarget_databases_path}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}" \
    -b 500 \
    -c "${CLUSTER_BUSTER_PATH}"

# check that the necessary output files exist then create a finished.txt file

# # Define the expected output file path
# output_file="${db_prefix}*.feather"

# # Check if the output file exists
# if [ -f "${output_file}" ]; then
#     echo "Output file exists: ${output_file}"
#     touch data/scenic-plus/cisTarget_dbs/outs/finished.txt
# else
#     echo "Output file does not exist: ${output_file}"
#     exit 1
# fi

# NEW SCENIC+ with snakemake pipeline
if [ ! -d "output/scenicplus/scplus_pipeline" ]; then
    mkdir -p output/scenicplus/scplus_pipeline
    scenicplus init_snakemake --out_dir output/scenicplus/scplus_pipeline
    mkdir output/scenicplus/outs
    mkdir output/scenicplus/tmp
    exit 1 
fi

# run through SnakeMake
cd output/scenicplus/scplus_pipeline/Snakemake/
# Detect the number of available CPU cores
nCores=$LSB_DJOB_NUMPROC
echo "Detected number of cores: ${nCores}"

# Run Snakemake with the detected number of cores
snakemake --cores ${nCores} 
