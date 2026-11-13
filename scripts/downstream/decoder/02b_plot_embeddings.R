#!/usr/bin/env Rscript
# =============================================================================
# Decoder branch, step 2b: check whether the first N LSI/PCA components are
# driven by batch rather than biology.
#
# Runs on step 2's RAW (unstandardized) embeddings -- the same files
# 03_train_cv.py consumes -- against step 1's sample_metadata.csv (which
# now carries a --batch_col value per pseudobulk sample; see
# 01_pseudobulk_creation.R). For each modality (ATAC LSI / RNA PCA) and
# each of the first --n_components components:
#
#   1. Pairwise scatter plots across every pair of those components,
#      colored three ways (batch / condition / cell_type) side by side --
#      a component cleanly separating by batch but not by condition or
#      cell type is the visual signature of a batch-driven artifact.
#   2. A quantitative complement to the same question: variance_explained.csv
#      reports R-squared from a one-way ANOVA (component ~ grouping_var)
#      for each (modality, component, grouping_var) combination. A
#      component with high R^2 against batch and low R^2 against
#      condition/cell_type is the concrete, numeric version of "this
#      looks batch-driven" -- more decisive than eyeballing a scatter plot.
#
# This is diagnostic only -- it doesn't filter or otherwise change any
# decoder input; it exists so a batch effect gets caught by inspection
# before it silently shapes the fitted decoder weights.
#
# USAGE:
#   Rscript 02b_plot_embeddings.R \
#     --pseudobulk_dir output/decoder/pseudobulk \
#     --features_dir output/decoder/features \
#     --n_components 4 \
#     --out_dir output/decoder/embedding_qc
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(patchwork, quietly = TRUE)
})

p <- arg_parser("Decoder step 2b: batch-vs-biology check on the first N LSI/PCA components")
p <- add_argument(p, "--pseudobulk_dir", help = "Directory containing step 1 outputs", type = "character",
                   default = "output/decoder/pseudobulk")
p <- add_argument(p, "--features_dir", help = "Directory containing step 2 outputs", type = "character",
                   default = "output/decoder/features")
p <- add_argument(p, "--n_components", help = "Number of leading components to check (per modality)",
                   type = "numeric", default = 4)
p <- add_argument(p, "--out_dir", help = "Output directory", type = "character",
                   default = "output/decoder/embedding_qc")
argv <- parse_args(p)

dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

sample_meta <- read.csv(file.path(argv$pseudobulk_dir, "sample_metadata.csv"), stringsAsFactors = FALSE)
if (!"batch" %in% colnames(sample_meta)) {
  stop("sample_metadata.csv has no 'batch' column -- re-run 01_pseudobulk_creation.R (it now writes ",
       "one via --batch_col) before this step.")
}
grouping_vars <- c("batch", "condition", "cell_type")

modalities <- list(
  atac = file.path(argv$features_dir, "atac_lsi_embeddings.csv"),
  rna  = file.path(argv$features_dir, "rna_pca_embeddings.csv")
)

variance_rows <- list()

for (modality in names(modalities)) {
  embeddings_path <- modalities[[modality]]
  if (!file.exists(embeddings_path)) {
    message("Skipping ", modality, ": ", embeddings_path, " not found.")
    next
  }
  message("[", modality, "] Loading ", embeddings_path, "...")
  embeddings <- read.csv(embeddings_path, row.names = 1, check.names = FALSE)

  n_comp <- min(argv$n_components, ncol(embeddings))
  if (n_comp < 2) {
    message("  Fewer than 2 components available (", ncol(embeddings), ") -- skipping.")
    next
  }
  comp_names <- colnames(embeddings)[seq_len(n_comp)]
  message("  Checking components: ", paste(comp_names, collapse = ", "))

  df <- embeddings[, comp_names, drop = FALSE]
  df$sample_id <- rownames(df)
  df <- df %>% left_join(sample_meta %>% select(sample_id, all_of(grouping_vars)), by = "sample_id")

  missing_meta <- sum(is.na(df$batch))
  if (missing_meta > 0) {
    warning("  ", missing_meta, " sample(s) in ", embeddings_path, " have no matching row in ",
            "sample_metadata.csv -- dropped from this check.")
    df <- df %>% filter(!is.na(batch))
  }

  # -- Quantitative: R^2 of component ~ grouping_var, one-way ANOVA via lm() --
  #
  # adj_r_squared, not raw r_squared, is the primary flagging statistic:
  # grouping variables with more levels (batch typically has far more
  # levels than condition/cell_type) get a mechanical R^2 boost purely
  # from having more free parameters to fit noise with, independent of any
  # real effect -- comparing raw R^2 across factors of very different
  # cardinality (e.g. 10-level batch vs 5-level condition) systematically
  # favors whichever has more levels. Adjusted R^2 penalizes for degrees
  # of freedom used and is the fair comparison here.
  #
  # cooks_distance/top_influence_sample/r_squared_excl_top_influence exist
  # to catch a DIFFERENT failure mode a global R^2 threshold can't see:
  # one extreme outlier sample dragging an otherwise-unremarkable
  # regression to a high R^2, vs. a real, general batch structure. A large
  # drop between r_squared and r_squared_excl_top_influence means "this is
  # one point's leverage," not "this component tracks batch broadly."
  for (comp in comp_names) {
    for (gv in grouping_vars) {
      if (n_distinct(df[[gv]]) < 2) next  # can't explain variance with a constant grouping variable
      fit <- lm(df[[comp]] ~ factor(df[[gv]]))
      s <- summary(fit)

      cooks_d <- cooks.distance(fit)
      top_idx <- which.max(cooks_d)
      top_sample <- if ("sample_id" %in% colnames(df)) df$sample_id[top_idx] else as.character(top_idx)

      # Refit excluding only that one most-influential point (only meaningful
      # if doing so doesn't collapse a group to zero members).
      r_squared_excl <- NA_real_
      keep <- setdiff(seq_len(nrow(df)), top_idx)
      if (n_distinct(df[[gv]][keep]) >= 2) {
        fit_excl <- lm(df[[comp]][keep] ~ factor(df[[gv]][keep]))
        r_squared_excl <- summary(fit_excl)$r.squared
      }

      variance_rows[[length(variance_rows) + 1]] <- data.frame(
        modality = modality, component = comp, grouping_var = gv,
        r_squared = s$r.squared, adj_r_squared = s$adj.r.squared, n_groups = n_distinct(df[[gv]]),
        max_cooks_distance = max(cooks_d), top_influence_sample = top_sample,
        r_squared_excl_top_influence = r_squared_excl,
        r_squared_drop_from_top_point = s$r.squared - r_squared_excl
      )
    }
  }

  # -- Visual: every pairwise combination of the first n_comp components,
  # one PDF page per pair, three colorings (batch/condition/cell_type)
  # side by side via patchwork so they're directly comparable. --
  pairs <- combn(comp_names, 2, simplify = FALSE)
  pdf_path <- file.path(argv$out_dir, paste0(modality, "_embedding_batch_check.pdf"))
  pdf(pdf_path, width = 15, height = 5)
  for (pr in pairs) {
    panels <- lapply(grouping_vars, function(gv) {
      ggplot(df, aes(x = .data[[pr[1]]], y = .data[[pr[2]]], color = .data[[gv]])) +
        geom_point(size = 1.5, alpha = 0.8) +
        labs(title = paste0("Colored by ", gv), x = pr[1], y = pr[2]) +
        theme_bw() +
        theme(legend.position = "bottom", legend.title = element_blank())
    })
    print(Reduce(`|`, panels) + plot_annotation(title = paste0(toupper(modality), ": ", pr[1], " vs ", pr[2])))
  }
  dev.off()
  message("  Written: ", pdf_path, " (", length(pairs), " page(s))")
}

variance_df <- bind_rows(variance_rows) %>% arrange(modality, component, desc(adj_r_squared))
write.csv(variance_df, file.path(argv$out_dir, "variance_explained.csv"), row.names = FALSE)
message("\nWritten: ", file.path(argv$out_dir, "variance_explained.csv"))

# Two INDEPENDENT flags, deliberately not nested inside each other:
#   - BATCH-DRIVEN (general): adjusted batch R^2 clears both an absolute
#     bar and the biological grouping variables' own adjusted R^2 -- fair
#     across batch's typically-higher level count vs condition/cell_type
#     (see the comment above the R^2 computation loop).
#   - SINGLE-SAMPLE OUTLIER: dropping the single most influential point
#     (by Cook's distance) collapses the batch R^2 substantially. This is
#     checked unconditionally, NOT gated behind the >0.5 threshold above --
#     a component whose aggregate batch R^2 never clears that bar can
#     still be visibly distorted by one outlier sample (this is exactly
#     what a first version of this check missed for LSI_2: its overall
#     batch R^2 was unremarkable, but a single cntrl.1 outlier point was
#     still clearly visible in the scatter plot).
message("\n-- Batch-vs-biology summary (adjusted R^2 per component) ------------------")
for (modality in unique(variance_df$modality)) {
  sub <- variance_df %>% filter(modality == !!modality)
  for (comp in unique(sub$component)) {
    comp_sub <- sub %>% filter(component == comp)
    batch_row <- comp_sub %>% filter(grouping_var == "batch")
    bio_r2 <- comp_sub %>% filter(grouping_var != "batch") %>% pull(adj_r_squared) %>% max(na.rm = TRUE)
    batch_adj_r2 <- if (nrow(batch_row) == 1) batch_row$adj_r_squared else NA_real_

    flags <- character(0)
    if (!is.na(batch_adj_r2) && batch_adj_r2 > 0.5 && batch_adj_r2 > bio_r2) {
      flags <- c(flags, "LIKELY BATCH-DRIVEN (general)")
    }
    if (nrow(batch_row) == 1 && !is.na(batch_row$r_squared_excl_top_influence) &&
        batch_row$r_squared_drop_from_top_point > 0.3) {
      flags <- c(flags, sprintf("SINGLE-SAMPLE OUTLIER (dropping '%s' takes batch R^2 %.3f -> %.3f)",
                                 batch_row$top_influence_sample, batch_row$r_squared,
                                 batch_row$r_squared_excl_top_influence))
    }
    flag_str <- if (length(flags) > 0) paste0(" <-- ", paste(flags, collapse = "; ")) else ""
    message(sprintf("  [%s] %s: batch adj.R^2=%.3f, max(condition/cell_type) adj.R^2=%.3f%s",
                     modality, comp, ifelse(is.na(batch_adj_r2), NA, batch_adj_r2), bio_r2, flag_str))
  }
}
message("----------------------------------------------------------------------------")
message("Done.")
