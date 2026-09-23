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

message("Running Filtering Pipeline")
remove.doublets <- !tolower(Sys.getenv("remove_doublets", unset = "true")) %in% c("false", "f", "0", "no")
qc.df <- read.csv(file = argv$qc.sheet, colClasses = "character", na.strings = c("NA", ""))

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
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-02-filter-obj-list.RDS"))
