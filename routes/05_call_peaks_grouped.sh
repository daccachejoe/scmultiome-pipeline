#!/bin/bash

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
    module load macs2
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

# override to point this stage at a different input, e.g.:
# INPUT_RDS=output/RDS-files/my-variant-04-label-celltypes-obj-list.RDS run/runmultiome call_peaks_grouped
INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-04-label-celltypes-obj-list.RDS}"

scripts/seurat_signac_pipeline.R \
    callpeaks \
    configs/samplesheet.csv \
    -g ct \
    --project_prefix $project_prefix-grouped-peaks \
    -R "$INPUT_RDS" \
    -m $my_macs_path