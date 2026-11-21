#!/bin/bash
# Downstream branch: ArchR reproducible peak calling on labelled cell groups.
# Reads stage 05's output by default (cell types labelled) and adds an
# ArchR peak assay next to the existing MACS/Signac one. Independent of the
# other downstream branches. Recipe and rationale: header of
# scripts/downstream/archr/02_call_peaks.R.
#
# Three steps, each its own R process inside this one job, so that no
# Seurat object is in memory while ArchR forks its workers (step 2):
#   1 01_export_cells.R    object -> cells.csv (library, barcode, fragments, peak_group, rep)
#   2 02_call_peaks.R      Arrow files -> pseudo-replicates -> MACS2 -> PeakMatrix
#   3 03_add_peak_assay.R  PeakMatrix -> new ChromatinAssay -> checkpoint RDS
# MACS2 is not module-loaded here: macs2-wrapper.sh loads it only for MACS2
# itself (loading it into the R job breaks Rcpp).

set -e

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler detected; running with the R already on PATH"
fi

# override to point this branch at a different input, e.g.:
# INPUT_RDS=output/RDS-files/my-variant-05-callpeaks-obj-list.RDS run/runmultiome archr_peaks
INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-grouped-peaks-05-callpeaks-obj-list.RDS}"

# Per-run overrides: ARCHR_* env var -> archr_* config key -> default. e.g.
# per cell type x timepoint, with donor-level pseudo-replicates:
#   ARCHR_GROUP_COL=cell_group ARCHR_SPLIT_COL=timepoint ARCHR_REP_COLS=orig.ident,donor \
#     run/runmultiome archr_peaks
GROUP_COL="${ARCHR_GROUP_COL:-${archr_group_col:-ct.spec}}"
SPLIT_COL="${ARCHR_SPLIT_COL:-${archr_split_col:-}}"
REP_COLS="${ARCHR_REP_COLS:-${archr_rep_cols:-orig.ident}}"
MIN_GROUP="${ARCHR_MIN_GROUP:-${archr_min_group:-50}}"
ATAC_ASSAY="${ARCHR_ATAC_ASSAY:-auto}"
ASSAY_NAME="${ARCHR_ASSAY_NAME:-${archr_assay_name:-peaks_archr}}"
# threads: explicit -> config -> the cores this job was given -> 8
THREADS="${ARCHR_THREADS:-${archr_threads:-${LSB_DJOB_NUMPROC:-${SLURM_CPUS_PER_TASK:-8}}}}"

# The output paths encode the input and grouping, so runs with different
# settings never overwrite each other. ARCHR_RUN_TAG sets a shorter name.
sanitize() { echo "$1" | tr -c 'A-Za-z0-9_.-' '_' | sed 's/_*$//'; }
INPUT_TAG="$(basename "$INPUT_RDS" .RDS)"
INPUT_TAG="${INPUT_TAG%-obj-list}"
INPUT_TAG="${INPUT_TAG%-obj}"
AUTO_RUN_TAG="$(sanitize "$INPUT_TAG")_grp-$(sanitize "$GROUP_COL")_split-$(sanitize "${SPLIT_COL:-none}")_rep-$(sanitize "$REP_COLS")_min-$MIN_GROUP"
RUN_TAG="${ARCHR_RUN_TAG:-$AUTO_RUN_TAG}"
OUT_DIR="output/archr/$RUN_TAG"
OUT_RDS="output/RDS-files/$project_prefix-archr-$RUN_TAG-obj-list.RDS"
echo "ArchR run tag: $RUN_TAG"
echo "Input: $INPUT_RDS"
echo "Output: $OUT_DIR, $OUT_RDS"

if [ ! -f "$INPUT_RDS" ]; then
    echo "Input not found: $INPUT_RDS (set INPUT_RDS)"
    exit 1
fi

echo "[archr 1/3] Export cells, peak groups and pseudo-replicates"
if [ -f "$OUT_DIR/cells.csv" ]; then
    echo "  $OUT_DIR/cells.csv exists; reusing it"
else
    Rscript scripts/downstream/archr/01_export_cells.R "$INPUT_RDS" \
        --out_dir "$OUT_DIR" \
        --atac_assay "$ATAC_ASSAY" \
        --group_col "$GROUP_COL" \
        --split_col "$SPLIT_COL" \
        --rep_cols "$REP_COLS" \
        --min_group "$MIN_GROUP"
fi

echo "[archr 2/3] Arrow files -> pseudo-replicates -> MACS2 reproducible peaks -> peak matrix"
if [ -f "$OUT_DIR/peaks.rds" ]; then
    echo "  $OUT_DIR/peaks.rds exists; reusing it"
else
    Rscript scripts/downstream/archr/02_call_peaks.R \
        --out_dir "$OUT_DIR" \
        --threads "$THREADS" \
        --macs2 scripts/downstream/archr/macs2-wrapper.sh
fi

echo "[archr 3/3] Add the ArchR peak matrix to the Seurat object"
Rscript scripts/downstream/archr/03_add_peak_assay.R "$INPUT_RDS" \
    --out_dir "$OUT_DIR" \
    --out_rds "$OUT_RDS" \
    --atac_assay "$ATAC_ASSAY" \
    --assay_name "$ASSAY_NAME"
