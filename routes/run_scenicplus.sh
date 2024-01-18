#!/bin/bash

# first we convert the annotated object to h5ad    
source $personal_anaconda_path
conda activate $sceasy_env_name

Rscript scripts/convert_seurat_to_h5ad.R \
    output/RDS-files/$project_prefix-grouped-peaks-callpeaks-obj-list.RDS \
    data/scenicplus/$project_prefix.h5ad \
    $sceasy_env_name
    
# then we run a script in R to export the necessary data
Rscript routes/export_scenicplus_data.R \
    output/RDS-files/$project_prefix-grouped-peaks-callpeaks-obj-list.RDS

# then we run the scenicplus preprocessing script
conda deactivate
conda activate $scenicplus_env_name

# make sure the prepcocess config file has compute_topics, indentify_enhancers set to True but not runPycisTarget
config_file=configs/scenicplus-preprocess-config.yaml
python scripts/scenicplus-pipeline.py $config_file

# deactivate the scenicplus environment and start up the create cisTargetDB env
conda deactivate
module load bedtools

# genome.fa = reference file
fasta_filename=/gpfs/data/sequence/cellranger-refdata/refdata-cellranger-arc-GRCh38-2020-A/fasta/genome.fa

# bed file with peaks in scATAC file
bedtools getfasta \
    -fi $fasta_filename \
    -bed data/raw/macs-peaks/grouped-peaks.bed \
    -fo data/scenic-plus/cisTarget_dbs/pseudobulk_peaks.fa
${create_cistarget_databases_dir}create_fasta_with_padded_bg_from_bed.sh $fasta_filename data/raw/macs-peaks/grouped-peaks.bed  500 yes

module purge
module load default-enviroment
source personal_anaconda_path
conda activate $create_cistarget_databases_env_name

# FASTA file with sequences per region IDs / gene IDs.
fasta_filename=data/scenic-plus/cisTarget_dbs/pseudobulk_peaks.fa
# Directory with motifs in Cluster-Buster format. Make sure this data is downloaded
motifs_dir=data/scenicplus/cisTarget_dbs/v10nr_clust_public/singletons
test [ -f $motifs_dir ]; then
    echo "Motifs directory exists: ${motifs_dir}"
else
    echo "Motifs directory does not exist: ${motifs_dir}"
    echo "Download the motifs directory using: "
    echo "motif_database_url='https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'"
    echo "wget "${motif_database_url}" --no-check-certificate"
    exit 1
fi
# File with motif IDs (base name of motif file in ${motifs_dir}).
# motifs_list_filename=data/scenic-plus/cisTarget_dbs/human-motif-names.txt # text file edited from all TFs hg38.txt in R # this did not work
ls data/scenicplus/cisTarget_dbs/v10nr_clust_public/singletons/*.cb | sed 's/\.cb$//' > data/scenicplus/cisTarget_dbs/prefixes.txt
motifs_list_filename=data/scenicplus/cisTarget_dbs/prefixes.txt # all .cb files in the v10-clust directory

# cisTarget motif database output prefix.
db_prefix=$project_prefix
nbr_threads=32

"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"

# check that the necessary output files exist then create a finished.txt file

# Define the expected output file path
output_file="${db_prefix}.feather"

# Check if the output file exists
if test -f "${output_file}"; then
    echo "Output file exists: ${output_file}"
    touch data/scenic-plus/cisTarget_dbs/outs/finished.txt
else
    echo "Output file does not exist: ${output_file}"
    exit 1
fi

# now continue with scenicplus pipeline
conda deactivate
conda activate $scenicplus_env_name

# add db_prefix to the process config file
sed -i "s/db_prefix: ''/db_prefix: '$project_prefix'/g" configs/scenicplus-process-config.yaml

# make sure the downstream process config file has pycisTarget, scenicplus, and load_objects set to true
config_file=data/scenicplus/scenicplus-downstream-config.yaml
python routes/scenicplus-joe-pipeline.py $config_file