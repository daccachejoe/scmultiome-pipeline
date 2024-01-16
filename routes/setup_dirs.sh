# dirsetup for a new multiome project
# . = project_directory

# desired set up
# ./data
# ./data/raw
# .data/scenicplus
# ./data/scenicplus/cisTarget_dbs
# ./scripts
# ./scripts/logs
# ./output
# ./output/RDS-files
# ./output/plots
# ./output/tables
# 
# 
# 
# 
# 
# 
mkdir -p ./data/raw \
         ./data/scenicplus/cisTarget_dbs \
         ./data/raw/macs-peaks \
         ./scripts/logs \
         ./output/RDS-files \
         ./output/plots \
         ./output/tables \
         ./output/ucd \
         ./configs

# init the samplesheet and add required column names
touch configs/samplesheet.csv
echo "sampleName,path" > configs/samplesheet.csv
# init the qc data frame and insert no defaults
touch configs/qc_df.csv
echo "sampleName,cluster.to.remove,vars.to.filter.by,var.filter,filter.direction" > configs/qc_df.csv
# init the resolution to use file, empty
touch configs/resolution_to_use.txt
# init the cluster labels file, with colnames
touch configs/cluster_labels.csv
echo "sampleName,cluster,ct,ct.spec" > configs/cluster_labels.csv
# cp the scenicplus config file, twice, once as the preprocess config and once as the process config
cp scripts/scenicplus-config.yml configs/scenicplus-preprocess-config.yml
cp scripts/scenicplus-config.yml configs/scenicplus-process-config.yml
# change the process conig to set pycisTarget, scenicplus, and load_objects to true
sed -i 's/pycisTarget: false/pycisTarget: true/g' configs/scenicplus-process-config.yml
sed -i 's/scenicplus: false/scenicplus: true/g' configs/scenicplus-process-config.yml
sed -i 's/scenicplus_downstream: false/scenicplus_downstream: true/g' configs/scenicplus-process-config.yml
sed -i 's/load_objects: false/load_objects: true/g' configs/scenicplus-process-config.yml
rm scripts/scenicplus-config.yml

