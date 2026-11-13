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
# donor_id) combination -- or (cell_type, donor_id, timepoint) for
# longitudinal Lesional data, see --timepoint_col -- only samples with
# >= --min_cells cells are kept.
#
# Column/assay defaults come from config/pipeline.config (decoder_*
# keys) via load_pipeline_config(), same mechanism every other stage in
# this pipeline uses, and can still be overridden per-invocation with the
# corresponding CLI flag.
#
# LONGITUDINAL DATA (--timepoint_col, default "timepoint"): Lesional cells
# (--location_col == --lesional_value) with a valid --timepoint_values
# entry get pseudobulked per (cell_type, donor, timepoint) instead of just
# (cell_type, donor), and their recorded 'condition' becomes
# '<condition>_<timepoint>' (e.g. AD_L -> AD_L_T1) to match a longitudinal
# constraint matrix's expanded row names (configs/decoder/
# constraint_matrix_longitudinal.csv). HC/NonLesional cells, and any
# Lesional cell with no valid timepoint value, are completely unaffected --
# fully backward compatible with a static (pre-longitudinal) analysis if
# --timepoint_col isn't present in the object at all.
#
# Outputs (under --out_dir):
#   rna_counts.csv, atac_counts.csv       -- genes/peaks x samples raw sums
#   sample_metadata.csv                    -- one row per sample (includes
#                                              a 'timepoint' column, NA for
#                                              non-longitudinal samples)
#   rna_gene_metadata.csv, atac_peak_metadata.csv
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(Seurat, quietly = TRUE)
  library(Signac, quietly = TRUE)
  library(Matrix, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(stringr, quietly = TRUE)
  library(GenomicRanges, quietly = TRUE)
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
p <- add_argument(p, "--batch_col", help = paste(
                    "Metadata column identifying each cell's raw sample/sequencing batch --",
                    "recorded in sample_metadata.csv so 02b_plot_embeddings.R can check whether",
                    "the LSI/PCA components are driven by batch rather than biology."),
                   type = "character", default = "orig.ident")
p <- add_argument(p, "--timepoint_col", help = paste(
                    "Metadata column with longitudinal timepoint (values in --timepoint_values),",
                    "meaningful only for --location_col == --lesional_value cells -- HC and",
                    "NonLesional samples stay timepoint-less regardless of this column's content.",
                    "For --location_col == --lesional_value cells, pseudobulk aggregates per",
                    "(cell_type, donor, timepoint) instead of just (cell_type, donor), and",
                    "'condition' becomes '<condition>_<timepoint>' (e.g. AD_L -> AD_L_T1) so it",
                    "matches a longitudinal constraint matrix's expanded row names. A Lesional cell",
                    "with no valid timepoint value keeps its plain (non-expanded) condition -- lets",
                    "static (pre-longitudinal) Lesional data coexist with new longitudinal data in",
                    "the same object; if that's NOT what you want (e.g. every Lesional cell should",
                    "have a timepoint and doesn't), validate_constraint_matrix() below will catch",
                    "the resulting stray plain-condition value against a longitudinal constraint",
                    "matrix that no longer has a row for it. If this column doesn't exist in the",
                    "object at all, no expansion happens anywhere -- fully backward compatible with",
                    "a pre-longitudinal (static) analysis."),
                   type = "character", default = "timepoint")
p <- add_argument(p, "--lesional_value", help = "--location_col value that marks a cell as Lesional (subject to --timepoint_col expansion)",
                   type = "character", default = "Lesional")
p <- add_argument(p, "--timepoint_values", help = "Comma-separated valid --timepoint_col values, in temporal order",
                   type = "character", default = "T1,T2,T3")
p <- add_argument(p, "--celltype_filter", help = paste(
                    "Comma-separated --celltype_col values to restrict this run to (e.g. 'FIB'). ",
                    "Applied before pseudobulk aggregation, so this is a real subset of the",
                    "analysis, not a post-hoc filter of the output. Empty (default) = no filtering."),
                   type = "character", default = "")
p <- add_argument(p, "--condition_filter", help = paste(
                    "Comma-separated --condition_col values to restrict this run to (e.g.",
                    "'AD_L,AD_NL' for a within-disease Lesional-vs-NonLesional comparison,",
                    "decoupled from the shared multi-condition model's beta_loc/beta_int",
                    "pooling). Pair with a --constraint_matrix whose rows are exactly this set.",
                    "Empty (default) = no filtering."),
                   type = "character", default = "")
p <- add_argument(p, "--min_cells", help = "Minimum cells required to retain a pseudobulk sample",
                   type = "numeric", default = 10)
p <- add_argument(p, "--peak_prevalence_min_pct", help = "Prevalence threshold reported per cell type (%%); informational only, actual filtering happens in step 2",
                   type = "numeric", default = 5)
p <- add_argument(p, "--skip_peak_prevalence", help = "Skip peak prevalence computation", flag = TRUE)
p <- add_argument(p, "--skip_tss_distance", help = paste(
                    "Skip computing distance-to-nearest-TSS (the 'distance' column step 2's",
                    "distal-peak filter reads). Only needed if the object's ATAC assay has no",
                    "Signac Annotation() to compute it from -- see compute_tss_distance()."),
                   flag = TRUE)
p <- add_argument(p, "--constraint_matrix", help = paste(
                    "Path to the decoder's constraint matrix CSV (conditions x components).",
                    "If given, --condition_col's values are cross-checked against its row names",
                    "right here, before AggregateExpression() (the expensive step below) runs --",
                    "instead of surfacing as a Python ValueError in step 3, after steps 1-2 have",
                    "already spent minutes of compute on data that was always going to fail",
                    "downstream. Pass an empty string to skip this check (e.g. for standalone runs",
                    "with no constraint matrix yet)."),
                   type = "character", default = "configs/decoder/constraint_matrix.csv")
argv <- parse_args(p)

# Fails loudly, with the actual mismatched condition names, if
# --constraint_matrix is given and any condition value found in the data
# is absent from its row names -- see the flag's help text above for why
# this check lives here rather than being left to fail downstream in
# Python. Not a hard error if the file is simply missing/unset, since this
# script is also runnable standalone without a constraint matrix in hand
# yet; 03_train_cv.py will still fail loudly on a real mismatch later if
# this check is skipped.
validate_constraint_matrix <- function(constraint_matrix_path, condition_values) {
  if (is.na(constraint_matrix_path) || !nzchar(constraint_matrix_path)) {
    message("  (--constraint_matrix not given -- skipping condition/constraint-matrix cross-check)")
    return(invisible(NULL))
  }
  if (!file.exists(constraint_matrix_path)) {
    warning("--constraint_matrix '", constraint_matrix_path, "' does not exist -- skipping the ",
            "condition/constraint-matrix cross-check here. (Not a hard error in this script since ",
            "it's runnable standalone; a real mismatch will still fail loudly in step 3 later.)")
    return(invisible(NULL))
  }

  cmat_conditions <- rownames(read.csv(constraint_matrix_path, row.names = 1, check.names = FALSE))
  data_conditions <- sort(unique(as.character(condition_values)))
  missing_from_cmat <- setdiff(data_conditions, cmat_conditions)
  unused_in_data <- setdiff(cmat_conditions, data_conditions)

  if (length(missing_from_cmat) > 0) {
    stop(
      "\n\nConstraint matrix mismatch: condition value(s) found via --condition_col are NOT ",
      "present as row names in the constraint matrix. Fix '", constraint_matrix_path, "'s row ",
      "names (or --condition_col, if it's pointing at the wrong metadata column) and re-run.\n\n",
      "  Condition value(s) in the data but MISSING from the constraint matrix:\n",
      "    ", paste(missing_from_cmat, collapse = ", "), "\n\n",
      "  Row names actually in '", constraint_matrix_path, "':\n",
      "    ", paste(cmat_conditions, collapse = ", "), "\n"
    )
  }

  if (length(unused_in_data) > 0) {
    # Not an error: a constraint-matrix row with zero matching samples is
    # exactly the expected state for a not-yet-demultiplexed reference
    # condition (e.g. HC) -- see decoder_model.py's compute_sample_weights().
    message("  NOTE: constraint matrix row(s) with no matching samples in this data yet (fine if ",
            "not demultiplexed, e.g. a pending HC arm): ", paste(unused_in_data, collapse = ", "))
  }

  message("  Constraint matrix cross-check passed: all ", length(data_conditions),
          " condition value(s) found in the data have a matching row in the constraint matrix.")
}

# "auto" (the default) mirrors functions.R's own
# ifelse("ATAC" %in% names(seu@assays), "ATAC", "peaks") logic, since which
# name a given object uses depends on which Seurat/Signac version created
# it -- older objects in this pipeline used "ATAC", current ones use "peaks".
resolve_atac_assay <- function(obj, requested) {
  if (requested != "auto") return(requested)
  if ("ATAC" %in% Assays(obj)) "ATAC" else "peaks"
}

# Step 2's distal-peak filter needs distance-to-nearest-TSS (bp), which is
# NOT something Signac/MACS peak calling adds to a ChromatinAssay's
# meta.features by default (it only has calling stats like count/percentile
# there) -- so this computes it here, once, while the live object (and its
# embedded gene annotation) is still loaded, rather than requiring the
# caller to have pre-supplied it.
#
# Uses the SAME annotation the object was built with: functions.R's
# CreateMultiomeSeurat() attaches a Signac Annotation() (exon-level GRanges
# from GetGRangesFromEnsDb(), 'chr'-prefixed) to the ChromatinAssay at
# object-creation time (stage 01) -- see functions.R's CreateChromatinAssay
# call. Reusing that instead of reloading EnsDb here keeps this consistent
# with whatever genome build/annotation the rest of the pipeline actually
# used for this object, rather than risking a second, possibly different
# annotation load.
#
# TSS = start of gene body for "+"/"*" strand genes, end of gene body for
# "-" strand genes. Annotation() is typically exon-level (multiple ranges
# per gene), so this collapses it to one TSS point per gene with a plain
# group-by rather than Signac::CollapseToLongestTranscript() -- that
# function turned out not to be exported from this cluster's installed
# Signac version (Error: 'CollapseToLongestTranscript' is not an exported
# object from 'namespace:Signac'), and reaching into Signac's internals
# via ::: would just trade one version-fragility for another. This uses
# only base GenomicRanges/dplyr, which every stage of this pipeline
# already depends on.
compute_tss_distance <- function(obj, atac_assay, peak_meta) {
  ann <- tryCatch(Signac::Annotation(obj[[atac_assay]]), error = function(e) NULL)
  if (is.null(ann) || length(ann) == 0) {
    warning("No Signac Annotation() found on the '", atac_assay, "' assay -- cannot compute ",
            "distance-to-TSS. Step 2's distal filter will need --tss_distance_col pointing at a ",
            "column you supply some other way, or --skip_prevalence_filter is not a substitute here.")
    return(rep(NA_real_, nrow(peak_meta)))
  }

  # row.names = NULL: as.data.frame(GRanges) defaults to using names(ann)
  # (here, exon IDs) as row names, and those aren't guaranteed unique --
  # the same exon ID can legitimately appear under more than one
  # transcript model in Ensembl's annotation, which is exactly what
  # tripped the "duplicate row.names" error here. We don't need row names
  # at all (gene_col/seqnames/strand below are all regular columns), so
  # just skip that entirely instead of trying to de-duplicate them.
  ann_df <- as.data.frame(ann, row.names = NULL)
  gene_col <- intersect(c("gene_id", "gene_name"), colnames(ann_df))
  if (length(gene_col) == 0) {
    warning("Annotation() has no gene_id/gene_name column to group exons into genes by -- ",
            "cannot compute distance-to-TSS. Available columns: ",
            paste(colnames(ann_df), collapse = ", "))
    return(rep(NA_real_, nrow(peak_meta)))
  }
  gene_col <- gene_col[1]

  tss_df <- ann_df %>%
    dplyr::group_by(.data[[gene_col]], seqnames, strand) %>%
    dplyr::summarise(tss_pos = if (strand[1] == "-") max(end) else min(start), .groups = "drop")

  tss_gr <- GenomicRanges::GRanges(
    seqnames = tss_df$seqnames,
    ranges = IRanges::IRanges(start = tss_df$tss_pos, width = 1),
    strand = tss_df$strand
  )
  peaks_gr <- GenomicRanges::GRanges(
    seqnames = peak_meta$chr,
    ranges = IRanges::IRanges(start = peak_meta$start, end = peak_meta$end)
  )

  # ignore.strand = TRUE: "distal from any TSS" should consider TSS points
  # on BOTH strands as candidates for every peak (ATAC peaks have no
  # strand of their own). peaks_gr above is built with no strand=
  # argument, so it already defaults to "*", which GRanges' matching rules
  # treat as compatible with both "+" and "-" -- meaning this already
  # behaved this way even without the explicit flag. Passing it anyway
  # makes that intentional rather than an accident of peaks_gr happening
  # to be unstranded, in case tss_gr or peaks_gr's construction ever
  # changes to carry real strand information later.
  #
  # distanceToNearest() only returns a row for peaks that DO have a nearest
  # feature on a shared seqlevel; a peak on a contig absent from the gene
  # annotation is simply missing from `nearest`, which is why dist_vec
  # starts as all-NA rather than being fully overwritten -- those peaks
  # correctly fail the (!is.na(dist_vec) & ...) distal check in step 2
  # instead of silently defaulting to 0 (which would wrongly mark them distal).
  nearest <- GenomicRanges::distanceToNearest(peaks_gr, tss_gr, ignore.strand = TRUE)
  dist_vec <- rep(NA_real_, nrow(peak_meta))
  dist_vec[S4Vectors::queryHits(nearest)] <- S4Vectors::mcols(nearest)$distance
  dist_vec
}

# Fraction of cells with count > 0, per cell type, for every peak. This is
# feature-selection metadata only (consumed by step 2's prevalence filter,
# and by step 6's motif-enrichment background) -- it does not filter
# anything here. counts_mat stays a sparse dgCMatrix the whole way through
# (Matrix::rowMeans on a `> 0` comparison never densifies it), which matters
# because this runs on the full, unfiltered peak x cell matrix (hundreds of
# thousands of peaks x tens of thousands of cells).
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

  required_cols <- c(argv$celltype_col, argv$donor_col, argv$condition_col, argv$location_col, argv$batch_col)
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

  # --celltype_filter / --condition_filter subset the object BEFORE anything
  # downstream (aggregation, LSI/PCA, the decoder fit itself) ever sees the
  # excluded cells -- a real analysis subset, not a post-hoc filter of
  # already-computed output. Applied here, before validate_constraint_matrix()
  # below, so a --condition_filter run (e.g. just AD_L/AD_NL) is checked
  # against a --constraint_matrix scoped to that same subset, not the full
  # 5-condition one.
  if (nzchar(argv$celltype_filter)) {
    keep_celltypes <- trimws(strsplit(argv$celltype_filter, ",")[[1]])
    unknown_celltypes <- setdiff(keep_celltypes, unique(as.character(obj@meta.data[[argv$celltype_col]])))
    if (length(unknown_celltypes) > 0) {
      stop("--celltype_filter value(s) not found in --celltype_col ('", argv$celltype_col, "'): ",
           paste(unknown_celltypes, collapse = ", "), ". Available: ",
           paste(sort(unique(as.character(obj@meta.data[[argv$celltype_col]]))), collapse = ", "))
    }
    obj <- obj[, obj@meta.data[[argv$celltype_col]] %in% keep_celltypes]
    message("  --celltype_filter '", argv$celltype_filter, "': ", ncol(obj), " cells remain")
  }
  if (nzchar(argv$condition_filter)) {
    keep_conditions <- trimws(strsplit(argv$condition_filter, ",")[[1]])
    unknown_conditions <- setdiff(keep_conditions, unique(as.character(obj@meta.data[[argv$condition_col]])))
    if (length(unknown_conditions) > 0) {
      stop("--condition_filter value(s) not found in --condition_col ('", argv$condition_col, "'): ",
           paste(unknown_conditions, collapse = ", "), ". Available: ",
           paste(sort(unique(as.character(obj@meta.data[[argv$condition_col]]))), collapse = ", "))
    }
    obj <- obj[, obj@meta.data[[argv$condition_col]] %in% keep_conditions]
    message("  --condition_filter '", argv$condition_filter, "': ", ncol(obj), " cells remain")
  }
  if (ncol(obj) == 0) {
    stop("--celltype_filter / --condition_filter left zero cells -- check the filter values ",
         "actually co-occur in this object.")
  }

  message("\n[2/6] Building sample index (cell_type x donor)...")
  # pb_group is the composite key AggregateExpression() will group cells by,
  # and it becomes the pseudobulk sample_id everywhere downstream. It must
  # contain NO underscores anywhere -- not in the sanitized cell_type/donor
  # pieces, and not in the separator joining them. Seurat's own identity-
  # class handling silently rewrites "_" -> "-" wherever this string is used
  # as a group/identity name (its own printed message: "Names of identity
  # class contain underscores, replacing with dashes"), so any underscore
  # we put in here comes back mangled in colnames(rna_pb)/colnames(atac_pb)
  # -- which is exactly why the sample_meta$sample_id == colnames(...)
  # check further down exists, and why it needs to actually pass. Using "-"
  # for BOTH the per-component sanitization AND the separator (never "_")
  # means Seurat has nothing left to rewrite, so the names round-trip
  # unchanged instead of us having to reverse-engineer its substitution.
  san <- function(x) str_replace_all(as.character(x), "[^A-Za-z0-9]", "-")

  timepoint_col_present <- argv$timepoint_col %in% colnames(obj@meta.data)
  timepoint_values <- trimws(strsplit(argv$timepoint_values, ",")[[1]])

  meta <- obj@meta.data %>%
    transmute(
      cell_type = .data[[argv$celltype_col]],
      donor_id  = .data[[argv$donor_col]],
      condition_raw = .data[[argv$condition_col]],
      location  = .data[[argv$location_col]],
      batch     = .data[[argv$batch_col]],
      timepoint_raw = if (timepoint_col_present) as.character(.data[[argv$timepoint_col]]) else NA_character_
    ) %>%
    mutate(
      is_lesional = location == argv$lesional_value,
      # A Lesional cell only gets timepoint-expanded if its value is one of
      # --timepoint_values; anything else (NA, "", an unrecognized label)
      # leaves it as a plain, unexpanded condition -- see --timepoint_col's
      # help text for why that's intentional rather than an error here.
      has_valid_timepoint = is_lesional & !is.na(timepoint_raw) & timepoint_raw %in% timepoint_values,
      timepoint = if_else(has_valid_timepoint, timepoint_raw, NA_character_),
      condition = if_else(has_valid_timepoint, paste(condition_raw, timepoint_raw, sep = "_"), as.character(condition_raw)),
      pb_group = if_else(
        has_valid_timepoint,
        paste(san(cell_type), san(donor_id), san(timepoint_raw), sep = "--"),
        paste(san(cell_type), san(donor_id), sep = "--")
      )
    )

  if (timepoint_col_present) {
    n_lesional <- sum(meta$is_lesional)
    n_lesional_no_tp <- sum(meta$is_lesional & !meta$has_valid_timepoint)
    if (n_lesional_no_tp > 0) {
      message("  NOTE: ", n_lesional_no_tp, " / ", n_lesional, " Lesional cells have no valid ",
              "--timepoint_col value (expected one of: ", paste(timepoint_values, collapse = ", "),
              ") -- kept as a plain (non-timepoint-expanded) condition. If every Lesional cell ",
              "should have a timepoint in this run, this is a real data gap to investigate.")
    }
    unexpected_tp <- setdiff(unique(meta$timepoint_raw[meta$is_lesional]), c(timepoint_values, NA))
    if (length(unexpected_tp) > 0) {
      warning("  Lesional cell(s) have --timepoint_col value(s) not in --timepoint_values: ",
              paste(unexpected_tp, collapse = ", "), " -- treated as missing (plain condition), ",
              "not an error, but likely a typo in the data or in --timepoint_values.")
    }
  }

  message("  Cross-checking condition values against the constraint matrix...")
  validate_constraint_matrix(argv$constraint_matrix, meta$condition)

  cell_counts <- meta %>%
    count(pb_group, cell_type, donor_id, condition, location, timepoint, name = "n_cells")

  # batch is looked up separately (majority value per pb_group) rather than
  # folded into the count() grouping above: unlike condition/location, which
  # are safe to assume constant per (cell_type, donor_id), a demuxed donor_id
  # could in principle collide across two different raw pooled runs (souporcell
  # cluster IDs are only unique within one pool) -- grouping by batch too
  # would silently split such a pb_group into two count() rows and break the
  # "one row per pb_group" invariant everything below this point relies on.
  # Recording the majority batch (and warning on any real inconsistency) is
  # tolerant of that instead of crashing on it -- batch is diagnostic-only
  # (02b_plot_embeddings.R), never load-bearing for the decoder fit itself.
  batch_lookup <- meta %>%
    count(pb_group, batch, name = "n") %>%
    group_by(pb_group) %>%
    summarise(batch = batch[which.max(n)], n_distinct_batch = n_distinct(batch), .groups = "drop")
  inconsistent_batch <- batch_lookup %>% filter(n_distinct_batch > 1)
  if (nrow(inconsistent_batch) > 0) {
    warning(nrow(inconsistent_batch), " pseudobulk group(s) span more than one --batch_col value ",
            "-- recording each group's majority batch. This means a donor_id value is shared ",
            "across raw samples/batches; only affects the diagnostic embedding plots, not the ",
            "decoder fit itself. Groups: ", paste(inconsistent_batch$pb_group, collapse = ", "))
  }
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

  # Subset the object down to cells belonging to a retained (>= min_cells)
  # group BEFORE aggregating, rather than aggregating everything and
  # dropping small-group columns afterward -- keeps AggregateExpression()
  # from wasting time summing counts for groups we're about to discard, and
  # keeps compute_peak_prevalence() below scoped to the same cell set as
  # everything else in this file.
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

  # Fail loudly HERE, with the actual mismatched names, rather than letting
  # a silent mismatch surface three steps later as an opaque
  # stopifnot(identical(...)) failure with no diagnostic content (that's
  # what happened before this check existed: the sample_meta re-ordering
  # at the bottom of this function assumes colnames(rna_pb)/colnames(atac_pb)
  # are exactly the keep_ids set, and if AggregateExpression() renamed
  # anything, that assumption breaks two write.csv() calls after we could
  # have caught it).
  if (!setequal(colnames(rna_pb), keep_ids) || !setequal(colnames(atac_pb), keep_ids)) {
    only_in_rna  <- setdiff(colnames(rna_pb), keep_ids)
    only_in_atac <- setdiff(colnames(atac_pb), keep_ids)
    only_expected <- setdiff(keep_ids, union(colnames(rna_pb), colnames(atac_pb)))
    stop(
      "AggregateExpression() column names do not match the sanitized pb_group keys. ",
      "This almost always means pb_group still contains a character Seurat's identity-name ",
      "handling rewrites (e.g. '_' -> '-'); see the comment above pb_group's construction.\n",
      "  Only in rna_pb colnames : ", paste(head(only_in_rna, 5), collapse = ", "), "\n",
      "  Only in atac_pb colnames: ", paste(head(only_in_atac, 5), collapse = ", "), "\n",
      "  Expected (from pb_group) but missing from both: ", paste(head(only_expected, 5), collapse = ", ")
    )
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
  # Peak feature names look like "chr1-1000-1500", "chr1:1000-1500", or
  # "chr1_1000_1500" depending on how Signac's peak set was built; this
  # regex accepts any of ":"/"-"/"_" as the coordinate separators so it
  # doesn't silently mis-parse peaks from a differently-formatted object.
  peak_coords <- str_match(rownames(atac_pb), "^(.+)[:\\-_](\\d+)[:\\-_](\\d+)$")
  peak_meta <- data.frame(peak = rownames(atac_pb), chr = peak_coords[, 2],
                           start = as.integer(peak_coords[, 3]), end = as.integer(peak_coords[, 4]),
                           row.names = rownames(atac_pb), stringsAsFactors = FALSE)
  signac_peak_info <- tryCatch(obj[[atac_assay]]@meta.features[rownames(atac_pb), , drop = FALSE],
                                error = function(e) NULL)
  if (!is.null(signac_peak_info) && ncol(signac_peak_info) > 0) peak_meta <- cbind(peak_meta, signac_peak_info)

  if (!argv$skip_tss_distance) {
    peak_meta$distance <- compute_tss_distance(obj, atac_assay, peak_meta)
    n_dist_na <- sum(is.na(peak_meta$distance))
    message("  Distance-to-TSS computed for ", nrow(peak_meta) - n_dist_na, " / ", nrow(peak_meta),
            " peaks", if (n_dist_na > 0) paste0(" (", n_dist_na, " NA -- no nearest feature found)") else "")
  }

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
    left_join(batch_lookup %>% select(pb_group, batch), by = "pb_group") %>%
    arrange(match(pb_group, colnames(rna_pb))) %>%
    transmute(sample_id = pb_group, cell_type, donor_id, condition, location, batch, timepoint, n_cells)
  stopifnot(identical(sample_meta$sample_id, colnames(rna_pb)))
  stopifnot(identical(sample_meta$sample_id, colnames(atac_pb)))
  write.csv(sample_meta, file.path(argv$out_dir, "sample_metadata.csv"), row.names = FALSE)

  message("\nDone. ", nrow(sample_meta), " samples, ", nrow(rna_pb), " genes, ", nrow(atac_pb), " peaks.")
  message("Output in: ", normalizePath(argv$out_dir))
}

run_pseudobulk(argv)
