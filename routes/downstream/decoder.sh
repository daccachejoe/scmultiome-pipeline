#!/bin/bash
# Downstream branch: constrained multinomial logistic regression decoder.
# Reads stage 05's output (cell types called, peaks re-called grouped by
# cell type) -- independent of the other downstream branches (subcluster,
# SCENIC+), which also branch from stage 05. See README.md's "Decoder
# methodology & caveats" section for the analytical design of this branch.
#
# Steps 1-2 (pseudobulk creation, feature engineering) are R, run under
# this scheduler's R module like the rest of the pipeline. Steps 3-6
# (CV training, peak projection, stability testing, motif enrichment) are
# Python, run under decoder_env_name -- kept in a separate conda env
# because PyTorch pins conflicting dependency versions vs. the other
# Python tools this pipeline uses (see docs/INSTALLATION.md).

if [ "$SCHEDULER" == "slurm" ]; then
    module load "$slurm_r_module"
elif [ "$SCHEDULER" == "lsf" ]; then
    module load "$lsf_r_module"
else
    echo "No job scheduler available to submit job: $script"
fi

# override to point this branch at a different input, e.g.:
# INPUT_RDS=output/RDS-files/my-variant-05-callpeaks-obj-list.RDS run/runmultiome decoder
INPUT_RDS="${INPUT_RDS:-output/RDS-files/$project_prefix-grouped-peaks-05-callpeaks-obj-list.RDS}"

# Schema/env knobs, in precedence order: DECODER_* env var (per-run override,
# survives run/runmultiome's unconditional `source config/pipeline.config`
# since these names aren't config keys) -> decoder_* config/pipeline.config
# value -> the script's own hardcoded default. e.g. to run this branch once
# against an object that labels cell types "cell_type" instead of "ct.spec",
# without editing config/pipeline.config (which every other stage also
# reads):
#   DECODER_CELLTYPE_COL=cell_type run/runmultiome decoder
CELLTYPE_COL="${DECODER_CELLTYPE_COL:-${decoder_celltype_col:-ct.spec}}"
DONOR_COL="${DECODER_DONOR_COL:-${decoder_donor_col:-donor_id}}"
CONDITION_COL="${DECODER_CONDITION_COL:-${decoder_condition_col:-orig.ident}}"
LOCATION_COL="${DECODER_LOCATION_COL:-${decoder_location_col:-location}}"
RNA_ASSAY="${DECODER_RNA_ASSAY:-${decoder_rna_assay:-RNA}}"
ATAC_ASSAY="${DECODER_ATAC_ASSAY:-${decoder_atac_assay:-auto}}"
CONSTRAINT_MATRIX="${DECODER_CONSTRAINT_MATRIX:-configs/decoder/constraint_matrix.csv}"
FOLD_COL="${DECODER_FOLD_COL:-${decoder_fold_col:-donor_id}}"

# Output directory encodes the actual parameter combination this run used,
# so two runs with different schema/fold overrides land in separate,
# non-clobbering trees instead of silently overwriting output/decoder/ --
# what makes this branch modular across objects/schemas rather than
# single-configuration. Set DECODER_RUN_TAG directly for a shorter,
# hand-chosen label instead of the auto-generated one (e.g. for a run you
# expect to reference by name later).
sanitize() { echo "$1" | tr -c 'A-Za-z0-9_.-' '_'; }
CONSTRAINT_MATRIX_TAG="$(basename "$CONSTRAINT_MATRIX" .csv)"
AUTO_RUN_TAG="ct-$(sanitize "$CELLTYPE_COL")_donor-$(sanitize "$DONOR_COL")_cond-$(sanitize "$CONDITION_COL")_loc-$(sanitize "$LOCATION_COL")_fold-$(sanitize "$FOLD_COL")_rna-$(sanitize "$RNA_ASSAY")_atac-$(sanitize "$ATAC_ASSAY")_cmat-$(sanitize "$CONSTRAINT_MATRIX_TAG")"
RUN_TAG="${DECODER_RUN_TAG:-$AUTO_RUN_TAG}"

DECODER_ROOT="output/decoder/$RUN_TAG"
PSEUDOBULK_DIR="$DECODER_ROOT/pseudobulk"
FEATURES_DIR="$DECODER_ROOT/features"
echo "Decoder run tag: $RUN_TAG"
echo "Output root: $DECODER_ROOT"

echo "[decoder 1/6] Pseudobulk creation (cell type x donor)"
Rscript scripts/downstream/decoder/01_pseudobulk_creation.R \
    "$INPUT_RDS" \
    --out_dir "$PSEUDOBULK_DIR" \
    --celltype_col "$CELLTYPE_COL" \
    --donor_col "$DONOR_COL" \
    --condition_col "$CONDITION_COL" \
    --location_col "$LOCATION_COL" \
    --rna_assay "$RNA_ASSAY" \
    --atac_assay "$ATAC_ASSAY"

echo "[decoder 2/6] Feature engineering (distal-peak LSI + RNA PCA)"
Rscript scripts/downstream/decoder/02_feature_engineering.R \
    --pseudobulk_dir "$PSEUDOBULK_DIR" \
    --out_dir "$FEATURES_DIR" \
    --genome_build "${genome:-hg38}"

source "$personal_anaconda_path"
conda activate "$decoder_env_name"

echo "[decoder 3/6] CV grid search + decoder fit (ATAC)"
python scripts/downstream/decoder/03_train_cv.py \
    --embeddings "$FEATURES_DIR/atac_lsi_embeddings.csv" \
    --sample_metadata "$PSEUDOBULK_DIR/sample_metadata.csv" \
    --constraint_matrix "$CONSTRAINT_MATRIX" \
    --modality atac --fold_col "$FOLD_COL" \
    --out_dir "$DECODER_ROOT/decoder_atac"

echo "[decoder 3/6] CV grid search + decoder fit (RNA)"
python scripts/downstream/decoder/03_train_cv.py \
    --embeddings "$FEATURES_DIR/rna_pca_embeddings.csv" \
    --sample_metadata "$PSEUDOBULK_DIR/sample_metadata.csv" \
    --constraint_matrix "$CONSTRAINT_MATRIX" \
    --modality rna --fold_col "$FOLD_COL" \
    --out_dir "$DECODER_ROOT/decoder_rna"

echo "[decoder 4/6] Project ATAC component weights to peak space"
python scripts/downstream/decoder/04_project_to_peak_space.py \
    --loadings "$FEATURES_DIR/atac_lsi_loadings.csv" \
    --cv_fold_weights "$DECODER_ROOT/decoder_atac/cv_fold_W_components.csv" \
    --final_weights "$DECODER_ROOT/decoder_atac/final_decoder_W.csv" \
    --out_dir "$DECODER_ROOT/peak_projections"

echo "[decoder 5/6] Identify consistently positive (stable) peaks per direction"
python scripts/downstream/decoder/05_identify_stable_peaks.py \
    --peak_weights_per_fold "$DECODER_ROOT/peak_projections/peak_weights_per_fold.csv" \
    --fdr_threshold 0.05 --correction_scope per_component \
    --out_dir "$DECODER_ROOT/stable_peaks"

echo "[decoder 6/6] Motif enrichment per regulatory direction (HOMER)"
if [ "$SCHEDULER" == "slurm" ] || [ "$SCHEDULER" == "lsf" ]; then
    module load homer 2>/dev/null || echo "  No 'homer' module found -- relying on decoder_homer_script/\$PATH."
fi
python scripts/downstream/decoder/06_run_motif_enrichment.py \
    --stable_peaks_dir "$DECODER_ROOT/stable_peaks" \
    --peaks_used "$FEATURES_DIR/atac_peaks_used.csv" \
    --genome "${genome:-hg38}" \
    --homer_script "${decoder_homer_script:-findMotifsGenome.pl}" \
    --out_dir "$DECODER_ROOT/motif_enrichment"

conda deactivate
