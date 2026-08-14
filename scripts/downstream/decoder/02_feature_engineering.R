#!/usr/bin/env Rscript
# =============================================================================
# Decoder branch, step 2: distal-peak LSI (ATAC) + PCA (RNA).
#
# Consumes step 1's pseudobulk outputs. ATAC: filters to distal peaks
# (>= --min_distance_to_tss from any TSS) AND peaks passing the prevalence
# union filter, then Signac::RunTFIDF + RunSVD (LSI component 1 dropped
# by default -- sequencing-depth-correlated). RNA: DESeq2 median-of-ratios
# size factors + log1p, then PCA (centered, not scaled) on pseudobulk
# profiles.
#
# Outputs (under --out_dir):
#   atac_lsi_embeddings.csv, atac_lsi_loadings.csv, atac_peaks_used.csv
#   rna_pca_embeddings.csv, rna_pca_loadings.csv
#   feature_engineering_config.json
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(Signac, quietly = TRUE)
  library(Seurat, quietly = TRUE)
  library(Matrix, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(jsonlite, quietly = TRUE)
  library(DESeq2, quietly = TRUE)
})

p <- arg_parser("Decoder step 2: distal-peak LSI (ATAC) + PCA (RNA)")
p <- add_argument(p, "--pseudobulk_dir", help = "Directory containing step 1 outputs", type = "character",
                   default = "output/decoder/pseudobulk")
p <- add_argument(p, "--out_dir", help = "Output directory for this step", type = "character",
                   default = "output/decoder/features")
p <- add_argument(p, "--tss_distance_col", help = "Column in atac_peak_metadata.csv giving distance (bp) to nearest TSS",
                   type = "character", default = "distance")
p <- add_argument(p, "--min_distance_to_tss", help = "Minimum distance (bp) from TSS to count as distal",
                   type = "numeric", default = 2000)
p <- add_argument(p, "--peak_prevalence_min_pct", help = "Prevalence threshold (%%), union across cell types",
                   type = "numeric", default = 5)
p <- add_argument(p, "--skip_prevalence_filter", help = "Skip the prevalence filter (distal filter still applies)", flag = TRUE)
p <- add_argument(p, "--n_lsi_components", help = "Number of LSI components to compute (before dropping dim 1)",
                   type = "numeric", default = 30)
p <- add_argument(p, "--keep_lsi_dim1", help = "Retain LSI component 1 instead of dropping it", flag = TRUE)
p <- add_argument(p, "--n_pca_components", help = "Number of RNA PCA components", type = "numeric", default = 30)
p <- add_argument(p, "--genome_build", help = "Genome build label, recorded for provenance", type = "character",
                   default = Sys.getenv("genome", unset = "hg38"))
argv <- parse_args(p)

drop_dim1 <- !argv$keep_lsi_dim1
dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

read_mat <- function(path) as.matrix(read.csv(path, row.names = 1, check.names = FALSE))

message("[1/5] Loading pseudobulk inputs from ", argv$pseudobulk_dir, "...")
atac_counts <- read_mat(file.path(argv$pseudobulk_dir, "atac_counts.csv"))
rna_counts  <- read_mat(file.path(argv$pseudobulk_dir, "rna_counts.csv"))
peak_meta   <- read.csv(file.path(argv$pseudobulk_dir, "atac_peak_metadata.csv"), row.names = 1, check.names = FALSE)
sample_meta <- read.csv(file.path(argv$pseudobulk_dir, "sample_metadata.csv"), check.names = FALSE)
peak_meta   <- peak_meta[rownames(atac_counts), , drop = FALSE]
stopifnot(identical(colnames(atac_counts), sample_meta$sample_id))
stopifnot(identical(colnames(rna_counts), sample_meta$sample_id))

message("[2/5] Filtering ATAC peaks (distal + prevalence)...")
if (!argv$tss_distance_col %in% colnames(peak_meta)) {
  stop("Column '", argv$tss_distance_col, "' not found in atac_peak_metadata.csv. ",
       "Available: ", paste(colnames(peak_meta), collapse = ", "))
}
dist_vec <- peak_meta[[argv$tss_distance_col]]
distal_pass <- !is.na(dist_vec) & abs(dist_vec) >= argv$min_distance_to_tss
message("  Distal (>= ", argv$min_distance_to_tss, " bp from TSS): ", sum(distal_pass), " / ", length(distal_pass))

if (argv$skip_prevalence_filter) {
  prevalence_pass <- rep(TRUE, nrow(peak_meta))
} else {
  pct_cols <- grep("^pct_cells_", colnames(peak_meta), value = TRUE)
  if (length(pct_cols) == 0) {
    warning("No pct_cells_* columns found; proceeding WITHOUT a prevalence filter.")
    prevalence_pass <- rep(TRUE, nrow(peak_meta))
  } else {
    prevalence_pass <- apply(peak_meta[, pct_cols, drop = FALSE], 1, max, na.rm = TRUE) >= argv$peak_prevalence_min_pct
  }
  message("  Prevalence pass (>= ", argv$peak_prevalence_min_pct, "% in >=1 cell type): ",
          sum(prevalence_pass), " / ", length(prevalence_pass))
}

keep_peaks <- distal_pass & prevalence_pass
message("  Peaks retained after both filters: ", sum(keep_peaks), " / ", length(keep_peaks))
if (sum(keep_peaks) < argv$n_lsi_components) {
  stop("Only ", sum(keep_peaks), " peaks survived filtering, fewer than --n_lsi_components (",
       argv$n_lsi_components, "). Loosen filters or reduce --n_lsi_components.")
}

peaks_used_log <- peak_meta %>%
  mutate(peak = rownames(peak_meta), distal_pass = distal_pass, prevalence_pass = prevalence_pass,
         used_in_lsi = keep_peaks) %>%
  select(peak, !!argv$tss_distance_col, distal_pass, prevalence_pass, used_in_lsi)
write.csv(peaks_used_log, file.path(argv$out_dir, "atac_peaks_used.csv"), row.names = FALSE)

atac_distal <- Matrix::Matrix(atac_counts[keep_peaks, , drop = FALSE], sparse = TRUE)

message("[3/5] Running TF-IDF + SVD (LSI) on distal peaks...")
tfidf_mat <- RunTFIDF(atac_distal, verbose = FALSE)
lsi <- RunSVD(tfidf_mat,
              n = min(argv$n_lsi_components, ncol(atac_distal) - 1, nrow(atac_distal) - 1),
              reduction.key = "LSI_", verbose = FALSE)
lsi_embeddings <- Embeddings(lsi)
lsi_loadings   <- Loadings(lsi)
if (drop_dim1) {
  lsi_embeddings <- lsi_embeddings[, -1, drop = FALSE]
  lsi_loadings   <- lsi_loadings[, -1, drop = FALSE]
}
lsi_embeddings <- lsi_embeddings[sample_meta$sample_id, , drop = FALSE]
message("  Final LSI feature set: ", ncol(lsi_embeddings), " components (dim1 dropped: ", drop_dim1, ")")
write.csv(lsi_embeddings, file.path(argv$out_dir, "atac_lsi_embeddings.csv"))
write.csv(lsi_loadings, file.path(argv$out_dir, "atac_lsi_loadings.csv"))

message("[4/5] Normalizing RNA (DESeq2 median-of-ratios) + PCA...")
size_factors <- estimateSizeFactorsForMatrix(rna_counts)
rna_norm <- sweep(rna_counts, 2, size_factors, FUN = "/")
rna_log  <- log1p(rna_norm)
n_pcs <- min(argv$n_pca_components, ncol(rna_log) - 1, nrow(rna_log) - 1)
pca <- prcomp(t(rna_log), center = TRUE, scale. = FALSE, rank. = n_pcs)
pca_embeddings <- pca$x[sample_meta$sample_id, , drop = FALSE]
pca_loadings   <- pca$rotation
message("  Final PCA feature set: ", ncol(pca_embeddings), " components")
write.csv(pca_embeddings, file.path(argv$out_dir, "rna_pca_embeddings.csv"))
write.csv(pca_loadings, file.path(argv$out_dir, "rna_pca_loadings.csv"))

message("[5/5] Writing provenance record...")
config <- list(
  step = "feature_engineering", pseudobulk_dir = argv$pseudobulk_dir,
  tss_distance_col = argv$tss_distance_col, min_distance_to_tss = argv$min_distance_to_tss,
  peak_prevalence_min_pct = argv$peak_prevalence_min_pct,
  prevalence_filter_applied = !argv$skip_prevalence_filter,
  n_peaks_total = length(keep_peaks), n_peaks_distal = sum(distal_pass),
  n_peaks_used_in_lsi = sum(keep_peaks), n_lsi_components_computed = argv$n_lsi_components,
  lsi_dim1_dropped = drop_dim1, n_lsi_components_final = ncol(lsi_embeddings),
  n_pca_components = ncol(pca_embeddings), rna_normalization = "DESeq2_median_of_ratios_log1p",
  genome_build = argv$genome_build
)
write_json(config, file.path(argv$out_dir, "feature_engineering_config.json"), auto_unbox = TRUE, pretty = TRUE)

message("\nDone. Output in: ", normalizePath(argv$out_dir))
