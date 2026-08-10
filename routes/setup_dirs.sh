#!/bin/bash
# Set up a new project's directory structure and per-project config templates.
# Safe to re-run: existing configs/* files are left untouched.
mkdir -p ./data/raw \
         ./data/scenicplus/cisTarget_dbs \
         ./data/raw/macs-peaks \
         ./scripts/outs \
         ./output/RDS-files \
         ./output/plots \
         ./output/tables \
         ./output/ucd \
         ./configs

# init the samplesheet and add required column names
if [ ! -f configs/samplesheet.csv ]; then
    echo "sampleName,path" > configs/samplesheet.csv
fi
# init the qc data frame and insert no defaults
if [ ! -f configs/qc_df.csv ]; then
    echo "sampleName,cluster.to.remove,vars.to.filter.by,var.filter,filter.direction" > configs/qc_df.csv
fi
# init the resolution to use file, empty
if [ ! -f configs/resolution_to_use.txt ]; then
    touch configs/resolution_to_use.txt
fi
# init the cluster labels file, with colnames
if [ ! -f configs/cluster_labels.csv ]; then
    echo "sampleName,cluster,ct,ct.spec" > configs/cluster_labels.csv
fi
# cp the scenicplus config template, twice, once as the preprocess config and once as the process config
if [ ! -f configs/scenicplus-preprocess-config.yml ]; then
    cp scripts/scenicplus-config-template.yml configs/scenicplus-preprocess-config.yml
    sed -i "s#__SCENICPLUS_TMP_DIR__#${scenicplus_tmp_dir}#g" configs/scenicplus-preprocess-config.yml
fi
if [ ! -f configs/scenicplus-process-config.yml ]; then
    cp scripts/scenicplus-config-template.yml configs/scenicplus-process-config.yml
    sed -i "s#__SCENICPLUS_TMP_DIR__#${scenicplus_tmp_dir}#g" configs/scenicplus-process-config.yml
    # the process config runs the full downstream SCENIC+ steps, unlike preprocess
    sed -i 's/pycisTarget: false/pycisTarget: true/g' configs/scenicplus-process-config.yml
    sed -i 's/scenicplus: false/scenicplus: true/g' configs/scenicplus-process-config.yml
    sed -i 's/scenicplus_downstream: false/scenicplus_downstream: true/g' configs/scenicplus-process-config.yml
    sed -i 's/load_objects: false/load_objects: true/g' configs/scenicplus-process-config.yml
fi

