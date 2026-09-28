# Stage: filter
# Removes doublets called in stage 01 (doublet.call == "doublet", see
# scripts/stages/doublets.R; skipped with remove_doublets=false in
# config/pipeline.config), then cells per configs/qc_df.csv (whole clusters
# to drop, plus threshold-based filters on arbitrary metadata columns), and
# saves the filtered per-sample objects as the stage 02 filter checkpoint.
# Sourced by scripts/seurat_signac_pipeline.R.
#
# qc_df.csv format: one or more rows per sampleName. Within a row,
# cluster.to.remove / vars.to.filter.by / var.filter / filter.direction may
# each hold several ";"-separated values; values from all of a sample's rows
# are pooled. filter.direction says which cells to KEEP:
#   nFeature_ATAC,500,greater  -> keep cells with nFeature_ATAC > 500
#   percent.mt,25,less         -> keep cells with percent.mt < 25
# A cell is kept only if it passes every threshold (and has no NA in the
# metrics being tested). Samples with no qc_df row skip these filters, with a
# message, rather than crashing (doublet removal still applies).
#
# Adaptive peak cutoff: var.filter = "auto" on an nFeature_peaks /
# nFeature_ATAC row with direction "greater" sets that sample's cutoff from
# its own ATAC depth:
#   cutoff = max(floor, min(max, round(fraction x median nFeature)))
# (config keys adaptive_peak_cutoff_max / _fraction / _floor, default 200 /
# 0.25 / 50). The median is over the sample's non-doublet cells, before any
# other qc_df filter. Why: a fixed 200 is far stricter for shallow ATAC
# libraries. On the PSO longitudinal v0 library (median 232 peaks) it removed
# 43% of cells whose TSS enrichment, nucleosome signal, RNA content and
# souporcell donor assignment matched the kept cells. The rule holds every
# library to the same relative bar: deep libraries keep 200, only shallower
# ones are lowered. 0.25 is where 200 already sat for a typical library
# (median ratio 0.24 across prelim-long-data); the floor of 50 still drops
# cells with almost no ATAC (v0 cells under 50 peaks had a median of 31 ATAC
# counts). Cutoffs used are written to
# output/tables/{prefix}-02-adaptive-peak-cutoffs.csv.

message("Running Filtering Pipeline")
remove.doublets <- !tolower(Sys.getenv("remove_doublets", unset = "true")) %in% c("false", "f", "0", "no")
qc.df <- read.csv(file = argv$qc.sheet, colClasses = "character", na.strings = c("NA", ""))

ReadFraction <- function(key, default) {
    val <- suppressWarnings(as.numeric(Sys.getenv(key, unset = default)))
    if (is.na(val) || val < 0) stop(key, " in config/pipeline.config must be a non-negative number, got: ", Sys.getenv(key))
    val
}
adaptive.max <- ReadFraction("adaptive_peak_cutoff_max", "200")
adaptive.fraction <- ReadFraction("adaptive_peak_cutoff_fraction", "0.25")
adaptive.floor <- ReadFraction("adaptive_peak_cutoff_floor", "50")
adaptive.vars <- c("nFeature_peaks", "nFeature_ATAC")
adaptive.log <- list()

# pool a qc_df column across a sample's rows, splitting ";"-joined values and
# dropping NA/empty entries
PoolQCField <- function(rows, field) {
    vals <- unlist(strsplit(rows[[field]][!is.na(rows[[field]])], split = ";"))
    vals <- trimws(vals)
    vals[nzchar(vals) & vals != "NA"]
}

obj.list <- lapply(obj.list, function(seu) {
    sample.name <- seu@project.name
    if (remove.doublets) {
        if ("doublet.call" %in% colnames(seu@meta.data)) {
            n.doublets <- sum(seu$doublet.call == "doublet", na.rm = TRUE)
            message("Filtering: ", sample.name, " -- removing ", n.doublets, " doublets called in stage 01")
            if (n.doublets > 0) seu <- subset(seu, cells = colnames(seu)[seu$doublet.call != "doublet"])
        } else {
            message("Filtering: ", sample.name, " -- no doublet.call column (stage 01 ran without the doublets stage)")
        }
    }
    md <- seu@meta.data
    rows <- qc.df[qc.df$sampleName == sample.name, , drop = FALSE]
    if (nrow(rows) == 0) {
        message("Filtering: ", sample.name, " -- no row in ", argv$qc.sheet, ", skipping threshold and cluster filters (", ncol(seu), " cells)")
        return(seu)
    }
    message("Filtering: ", sample.name)

    # whole clusters to remove
    clus.to.remove <- PoolQCField(rows, "cluster.to.remove")
    cells.in.clusters.to.remove <- rownames(md)[as.character(md$seurat_clusters) %in% clus.to.remove]

    # threshold filters: vars/values/directions pair up by position
    vars.to.filter.by <- PoolQCField(rows, "vars.to.filter.by")
    var.filter <- PoolQCField(rows, "var.filter")
    filter.direction <- PoolQCField(rows, "filter.direction")
    if (length(unique(c(length(vars.to.filter.by), length(var.filter), length(filter.direction)))) != 1) {
        stop("qc_df for ", sample.name, ": vars.to.filter.by (", length(vars.to.filter.by), "), var.filter (",
             length(var.filter), ") and filter.direction (", length(filter.direction),
             ") must have the same number of values")
    }
    missing.vars <- setdiff(vars.to.filter.by, colnames(md))
    if (length(missing.vars) > 0) {
        stop("qc_df for ", sample.name, ": metadata column(s) not found: ", paste(missing.vars, collapse = ", "),
             ". Available: ", paste(colnames(md), collapse = ", "))
    }
    bad.directions <- setdiff(filter.direction, c("greater", "less"))
    if (length(bad.directions) > 0) {
        stop("qc_df for ", sample.name, ": filter.direction must be 'greater' or 'less', got: ",
             paste(bad.directions, collapse = ", "))
    }

    # adaptive peak cutoffs, from this sample's non-doublet cells before any
    # other qc_df filter (see header)
    is.adaptive <- tolower(var.filter) == "auto"
    if (any(is.adaptive)) {
        bad <- is.adaptive & !(vars.to.filter.by %in% adaptive.vars & filter.direction == "greater")
        if (any(bad)) {
            stop("qc_df for ", sample.name, ": var.filter 'auto' is only for ", paste(adaptive.vars, collapse = "/"),
                 " with filter.direction 'greater', got: ",
                 paste(vars.to.filter.by[bad], filter.direction[bad], collapse = ", "))
        }
        singlet <- if ("doublet.call" %in% colnames(md)) md$doublet.call != "doublet" else rep(TRUE, nrow(md))
        singlet[is.na(singlet)] <- TRUE
        for (i in which(is.adaptive)) {
            med <- median(md[[vars.to.filter.by[i]]][singlet], na.rm = TRUE)
            cutoff <- max(adaptive.floor, min(adaptive.max, round(adaptive.fraction * med)))
            message("  adaptive ", vars.to.filter.by[i], " cutoff: median ", round(med), " (non-doublet cells) -> ", cutoff)
            adaptive.log[[length(adaptive.log) + 1]] <<- data.frame(
                sample = sample.name, variable = vars.to.filter.by[i], median.non.doublet = med, cutoff = cutoff,
                max = adaptive.max, fraction = adaptive.fraction, floor = adaptive.floor)
            var.filter[i] <- as.character(cutoff)
        }
    }
    thresholds <- suppressWarnings(as.numeric(var.filter))
    if (any(is.na(thresholds))) {
        stop("qc_df for ", sample.name, ": non-numeric var.filter value(s): ",
             paste(var.filter[is.na(thresholds)], collapse = ", "))
    }

    keep <- rep(TRUE, nrow(md))
    for (i in seq_along(vars.to.filter.by)) {
        passes <- if (filter.direction[i] == "greater") {
            md[[vars.to.filter.by[i]]] > thresholds[i]
        } else {
            md[[vars.to.filter.by[i]]] < thresholds[i]
        }
        passes[is.na(passes)] <- FALSE
        message("  keep ", vars.to.filter.by[i], " ", ifelse(filter.direction[i] == "greater", ">", "<"), " ",
                thresholds[i], ": ", sum(!passes), " cells fail")
        keep <- keep & passes
    }
    cells.to.filter <- rownames(md)[!keep]

    cells.to.remove <- union(cells.in.clusters.to.remove, cells.to.filter)
    message("  removing ", length(cells.to.remove), " of ", nrow(md), " cells (",
            length(cells.in.clusters.to.remove), " in clusters ", paste(clus.to.remove, collapse = ";"),
            ", ", length(cells.to.filter), " failing thresholds)")
    if (length(cells.to.remove) == 0) {
        return(seu)
    }
    subset(seu, cells = setdiff(colnames(seu), cells.to.remove))
})
if (length(adaptive.log) > 0) {
    write.csv(dplyr::bind_rows(adaptive.log), row.names = FALSE,
              file = paste0("output/tables/", argv$project_prefix, "-02-adaptive-peak-cutoffs.csv"))
}
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-02-filter-obj-list.RDS"))
