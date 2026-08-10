#!/bin/bash
# Downstream branch: per-lineage subclustering with Harmony batch
# correction. Reads stage 05's output (cell types called, peaks re-called
# grouped by cell type) -- independent of the other downstream branches
# (linkpeaks, SCENIC+), which also branch from stage 05.

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

if [[ $(wc -l < "configs/qc_df.csv") -gt 0 ]]; then
    # override to point this stage at a different input, e.g.:
    # INPUT_RDS=output/RDS-files/my-variant-05-callpeaks-obj-list.RDS run/runmultiome subcluster
    INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-grouped-peaks-05-callpeaks-obj-list.RDS}"

    scripts/seurat_signac_pipeline.R \
        subcluster \
        configs/samplesheet.csv \
        --project_prefix $project_prefix-lineage-improved \
        -m $my_macs_path \
        --RunHarmony \
        -R "$INPUT_RDS" \
        --qc.sheet configs/qc_df.csv
else
    echo "qc_df.csv does not have more than one line. Please fill it out and try again."
fi


