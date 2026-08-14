#!/usr/bin/env Rscript
# =============================================================================
# Decoder branch, step 1: pseudobulk creation.
#
# Aggregates raw counts per (cell type x donor) for both assays using
# Seurat::AggregateExpression(), and computes per-cell-type ATAC peak
# prevalence (fraction of cells with count > 0) as feature-selection
# metadata for step 2's distal-peak filter.
#
# Reads the object produced by stage 05 (call_peaks_grouped) -- peaks are
# already called grouped by cell type, and cluster->cell-type labels are
# already applied (stage 04). "sample" throughout = one (cell_type,
# donor_id) combination; only samples with >= --min_cells cells are kept.
#
# Column/assay defaults come from config/pipeline.config (decoder_*
# keys) via load_pipeline_config(), same mechanism every other stage in
# this pipeline uses, and can still be overridden per-invocation with the
# corresponding CLI flag.
#
# Outputs (under --out_dir):
#   rna_counts.csv, atac_counts.csv       -- genes/peaks x samples raw sums
#   sample_metadata.csv                    -- one row per sample
#   rna_gene_metadata.csv, atac_peak_metadata.csv
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(Seurat, quietly = TRUE)
  library(Signac, quietly = TRUE)
  library(Matrix, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(stringr, quietly = TRUE)
})

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/seurat_io.R")

p <- arg_parser("Decoder step 1: pseudobulk creation (cell type x donor)")
p <- add_argument(p, "seurat_rds", help = "Path to stage 05 (or later) RDS object", type = "character")
p <- add_argument(p, "--out_dir", help = "Output directory", type = "character",
                   default = "output/decoder/pseudobulk")
p <- add_argument(p, "--rna_assay", help = "RNA assay name", type = "character",
                   default = Sys.getenv("decoder_rna_assay", unset = "RNA"))
p <- add_argument(p, "--atac_assay", help = "ATAC assay name ('auto' = ATAC if present, else peaks)",
                   type = "character", default = Sys.getenv("decoder_atac_assay", unset = "auto"))
p <- add_argument(p, "--celltype_col", help = "Metadata column with cell type labels", type = "character",
                   default = Sys.getenv("decoder_celltype_col", unset = "ct.spec"))
p <- add_argument(p, "--donor_col", help = "Metadata column with donor ID", type = "character",
                   default = Sys.getenv("decoder_donor_col", unset = "donor_id"))
p <- add_argument(p, "--condition_col", help = "Metadata column with experimental condition", type = "character",
                   default = Sys.getenv("decoder_condition_col", unset = "orig.ident"))
p <- add_argument(p, "--location_col", help = "Metadata column with lesional/non-lesional location", type = "character",
                   default = Sys.getenv("decoder_location_col", unset = "location"))
p <- add_argument(p, "--min_cells", help = "Minimum cells required to retain a pseudobulk sample",
                   type = "numeric", default = 10)
p <- add_argument(p, "--peak_prevalence_min_pct", help = "Prevalence threshold reported per cell type (%%); informational only, actual filtering happens in step 2",
                   type = "numeric", default = 5)
p <- add_argument(p, "--skip_peak_prevalence", help = "Skip peak prevalence computation", flag = TRUE)
argv <- parse_args(p)

resolve_atac_assay <- function(obj, requested) {
  if (requested != "auto") return(requested)
  if ("ATAC" %in% Assays(obj)) "ATAC" else "peaks"
}

compute_peak_prevalence <- function(counts_mat, cell_type_vec) {
  cts <- sort(unique(cell_type_vec))
  frac_mat <- matrix(NA_real_, nrow = nrow(counts_mat), ncol = length(cts),
                      dimnames = list(rownames(counts_mat), cts))
  for (ct in cts) {
    idx <- cell_type_vec == ct
    frac_mat[, ct] <- Matrix::rowMeans(counts_mat[, idx, drop = FALSE] > 0)
  }
  frac_mat
}

run_pseudobulk <- function(argv) {
  dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Output directory: ", argv$out_dir)

  message("\n[1/6] Loading Seurat object...")
  obj <- first_seurat(readRDS(argv$seurat_rds))
  atac_assay <- resolve_atac_assay(obj, argv$atac_assay)
  message("  Cells total: ", ncol(obj), " | ATAC assay resolved to: ", atac_assay)

  required_cols <- c(argv$celltype_col, argv$donor_col, argv$condition_col, argv$location_col)
  missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "),
         ". Available: ", paste(colnames(obj@meta.data), collapse = ", "))
  }
  for (assay in c(argv$rna_assay, atac_assay)) {
    if (!assay %in% Assays(obj)) {
      stop("Assay '", assay, "' not found. Available: ", paste(Assays(obj), collapse = ", "))
    }
  }

  message("\n[2/6] Building sample index (cell_type x donor)...")
  meta <- obj@meta.data %>%
    transmute(
      cell_type = .data[[argv$celltype_col]],
      donor_id  = .data[[argv$donor_col]],
      condition = .data[[argv$condition_col]],
      location  = .data[[argv$location_col]]
    ) %>%
    mutate(
      pb_group = paste(
        str_replace_all(as.character(cell_type), "[^A-Za-z0-9]", "_"),
        str_replace_all(as.character(donor_id),  "[^A-Za-z0-9]", "_"),
        sep = "__"
      )
    )

  cell_counts <- meta %>%
    count(pb_group, cell_type, donor_id, condition, location, name = "n_cells")
  retained <- cell_counts %>% filter(n_cells >= argv$min_cells)
  dropped  <- cell_counts %>% filter(n_cells < argv$min_cells)

  if (nrow(dropped) > 0) {
    message("  Dropping ", nrow(dropped), " group(s) with < ", argv$min_cells, " cells:")
    dropped %>% arrange(n_cells) %>%
      mutate(msg = paste0("    ", pb_group, " (", n_cells, " cells)")) %>%
      pull(msg) %>% message()
  }
  message("  Retaining ", nrow(retained), " pseudobulk samples from ",
          n_distinct(retained$cell_type), " cell types x ",
          n_distinct(retained$donor_id), " donors.")

  keep_ids <- retained$pb_group
  keep_idx <- meta$pb_group %in% keep_ids
  obj <- obj[, keep_idx]
  obj$.pb_group <- meta$pb_group[keep_idx]

  message("\n[3/6] Aggregating RNA + ATAC counts via Seurat::AggregateExpression()...")
  agg <- AggregateExpression(object = obj, assays = c(argv$rna_assay, atac_assay),
                              group.by = ".pb_group", return.seurat = FALSE)
  rna_pb  <- as.matrix(agg[[argv$rna_assay]])
  atac_pb <- as.matrix(agg[[atac_assay]])
  message("  Pseudobulk RNA : ", nrow(rna_pb), " genes x ", ncol(rna_pb), " samples")
  message("  Pseudobulk ATAC: ", nrow(atac_pb), " peaks x ", ncol(atac_pb), " samples")

  if (!setequal(colnames(rna_pb), keep_ids) || !setequal(colnames(atac_pb), keep_ids)) {
    warning("AggregateExpression() column names do not exactly match the sanitized pb_group keys. ",
            "Inspect colnames(rna_pb)/colnames(atac_pb) vs. keep_ids before trusting sample_metadata.csv alignment.")
  }

  nz_genes <- rowSums(rna_pb) > 0
  if (any(!nz_genes)) rna_pb <- rna_pb[nz_genes, , drop = FALSE]
  nz_peaks <- rowSums(atac_pb) > 0
  if (any(!nz_peaks)) atac_pb <- atac_pb[nz_peaks, , drop = FALSE]

  write.csv(rna_pb, file.path(argv$out_dir, "rna_counts.csv"))
  write.csv(atac_pb, file.path(argv$out_dir, "atac_counts.csv"))

  message("\n[4/6] Writing gene metadata...")
  gene_meta <- data.frame(gene = rownames(rna_pb), row.names = rownames(rna_pb))
  gene_feature_info <- tryCatch(obj[[argv$rna_assay]]@meta.features[rownames(rna_pb), , drop = FALSE],
                                 error = function(e) NULL)
  if (!is.null(gene_feature_info) && ncol(gene_feature_info) > 0) gene_meta <- cbind(gene_meta, gene_feature_info)
  write.csv(gene_meta, file.path(argv$out_dir, "rna_gene_metadata.csv"))

  message("\n[5/6] Writing peak metadata", if (!argv$skip_peak_prevalence) " + prevalence stats" else "", "...")
  peak_coords <- str_match(rownames(atac_pb), "^(.+)[:\\-_](\\d+)[:\\-_](\\d+)$")
  peak_meta <- data.frame(peak = rownames(atac_pb), chr = peak_coords[, 2],
                           start = as.integer(peak_coords[, 3]), end = as.integer(peak_coords[, 4]),
                           row.names = rownames(atac_pb), stringsAsFactors = FALSE)
  signac_peak_info <- tryCatch(obj[[atac_assay]]@meta.features[rownames(atac_pb), , drop = FALSE],
                                error = function(e) NULL)
  if (!is.null(signac_peak_info) && ncol(signac_peak_info) > 0) peak_meta <- cbind(peak_meta, signac_peak_info)

  if (!argv$skip_peak_prevalence) {
    atac_counts_sc <- tryCatch(LayerData(obj, assay = atac_assay, layer = "counts"),
                                error = function(e) GetAssayData(obj, assay = atac_assay, slot = "counts"))
    atac_counts_sc <- atac_counts_sc[rownames(atac_pb), , drop = FALSE]
    frac_mat <- compute_peak_prevalence(atac_counts_sc, obj@meta.data[[argv$celltype_col]])
    colnames(frac_mat) <- paste0("pct_cells_", colnames(frac_mat))
    frac_mat <- frac_mat * 100
    passes_union <- apply(frac_mat, 1, max, na.rm = TRUE) >= argv$peak_prevalence_min_pct
    peak_meta <- cbind(peak_meta, as.data.frame(frac_mat), prevalence_pass_union = passes_union)
    message("  Peaks passing >", argv$peak_prevalence_min_pct, "% prevalence in >=1 cell type: ",
            sum(passes_union), " / ", nrow(peak_meta))
  }
  write.csv(peak_meta, file.path(argv$out_dir, "atac_peak_metadata.csv"))

  message("\n[6/6] Writing sample metadata...")
  sample_meta <- retained %>%
    filter(pb_group %in% colnames(rna_pb)) %>%
    arrange(match(pb_group, colnames(rna_pb))) %>%
    transmute(sample_id = pb_group, cell_type, donor_id, condition, location, n_cells)
  stopifnot(identical(sample_meta$sample_id, colnames(rna_pb)))
  stopifnot(identical(sample_meta$sample_id, colnames(atac_pb)))
  write.csv(sample_meta, file.path(argv$out_dir, "sample_metadata.csv"), row.names = FALSE)

  message("\nDone. ", nrow(sample_meta), " samples, ", nrow(rna_pb), " genes, ", nrow(atac_pb), " peaks.")
  message("Output in: ", normalizePath(argv$out_dir))
}

run_pseudobulk(argv)
