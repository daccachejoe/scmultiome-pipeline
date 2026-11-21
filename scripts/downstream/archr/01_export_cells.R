#!/usr/bin/env Rscript
# =============================================================================
# ArchR peaks branch, step 1: export the cells to call peaks on.
#
# Reads a Seurat object (stage 05 output by default, or any INPUT_RDS that
# has a cell-group column) and writes one row per cell with everything
# ArchR needs, so step 2 can run in a process with no Seurat object in
# memory (ArchR forks worker processes; forking a job that holds a
# multi-GB Seurat object is what made the original hand-run version hang
# and run out of memory).
#
# Per cell:
#   archr_cell   "<library>#<barcode>", ArchR's cell name
#   library      orig.ident; one Arrow file per library
#   barcode      Cell Ranger barcode, from the object's Fragment objects
#                (not parsed from the cell name, so merge prefixes don't matter)
#   fragments    that library's atac_fragments.tsv.gz
#   group        --group_col value
#   split        --split_col value (e.g. timepoint), or NA
#   rep          pseudo-replicate label = --rep_cols pasted with "|". Use
#                library + donor for donor-pooled libraries, so the two
#                donors in one library count as separate replicates.
#   peak_group   what ArchR calls peaks on: group, or group_split when
#                --split_col is set. Groups with < --min_group cells are
#                pooled into group_all (across splits); groups still below
#                --min_group are left out of peak calling (NA) but still
#                get counts in the final peak matrix.
#
# Inputs:  a Seurat RDS (bare object or a one-object list)
# Outputs: <out_dir>/cells.csv, <out_dir>/peak-groups.csv
# Usage (normally via routes/downstream/archr_peaks.sh):
#   Rscript scripts/downstream/archr/01_export_cells.R <rds> --out_dir output/archr/<tag> \
#     --group_col ct.spec [--split_col timepoint] [--rep_cols orig.ident,donor] [--min_group 50]
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(Seurat, quietly = TRUE)
  library(Signac, quietly = TRUE)
})

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/seurat_io.R")

p <- arg_parser("ArchR peaks step 1: export cells, peak groups and pseudo-replicates")
p <- add_argument(p, "seurat_rds", help = "Seurat RDS (stage 05 or later)", type = "character")
p <- add_argument(p, "--out_dir", help = "Run output directory", type = "character")
p <- add_argument(p, "--atac_assay", help = "ATAC assay holding the fragment files ('auto' = ATAC if present, else peaks)",
                  type = "character", default = "auto")
p <- add_argument(p, "--group_col", help = "Metadata column with the cell groups to call peaks on", type = "character",
                  default = Sys.getenv("archr_group_col", unset = "ct.spec"))
p <- add_argument(p, "--split_col", help = "Optional column to split groups by (e.g. timepoint); '' for none",
                  type = "character", default = Sys.getenv("archr_split_col", unset = ""))
p <- add_argument(p, "--rep_cols", help = "Comma-separated columns that define a pseudo-replicate", type = "character",
                  default = Sys.getenv("archr_rep_cols", unset = "orig.ident"))
p <- add_argument(p, "--min_group", help = "Smallest peak group, in cells", type = "integer",
                  default = as.integer(Sys.getenv("archr_min_group", unset = "50")))
argv <- parse_args(p)

split.col <- if (is.na(argv$split_col)) "" else argv$split_col
rep.cols <- strsplit(argv$rep_cols, ",")[[1]]
dir.create(argv$out_dir, showWarnings = FALSE, recursive = TRUE)

obj.list <- as_seurat_list(readRDS(argv$seurat_rds))
if (length(obj.list) != 1) {
  stop(argv$seurat_rds, " holds ", length(obj.list), " objects; the ArchR branch needs one merged object. ",
       "Point INPUT_RDS at a merged checkpoint (stage 02 or later).")
}
obj <- obj.list[[1]]
rm(obj.list)

assay <- argv$atac_assay
if (assay == "auto") assay <- if ("ATAC" %in% Assays(obj)) "ATAC" else "peaks"
if (!assay %in% Assays(obj)) stop("assay '", assay, "' not in ", argv$seurat_rds, " (assays: ", paste(Assays(obj), collapse = ", "), ")")

meta <- obj@meta.data
need <- unique(c("orig.ident", argv$group_col, if (nzchar(split.col)) split.col, rep.cols))
missing <- setdiff(need, colnames(meta))
if (length(missing)) {
  stop("column(s) ", paste(missing, collapse = ", "), " not in ", argv$seurat_rds,
       "'s metadata. Set ARCHR_GROUP_COL / ARCHR_SPLIT_COL / ARCHR_REP_COLS (or the archr_* config keys).")
}

# cell -> fragment file + barcode, from the object's own Fragment objects
frags <- Fragments(obj[[assay]])
if (!length(frags)) stop("assay '", assay, "' has no fragment files attached")
# subset objects keep a Fragment object for every original library, including ones with no cells left: skip those
fc <- do.call(rbind, lapply(frags, function(f) {
  cells <- GetFragmentData(f, "cells")
  if (!length(cells)) return(NULL)
  data.frame(seurat_cell = names(cells), barcode = unname(cells), fragments = GetFragmentData(f, "path"))
}))
if (is.null(fc)) stop("none of assay '", assay, "''s fragment files has any cells attached")
fc <- fc[fc$seurat_cell %in% rownames(meta), ]
if (anyDuplicated(fc$seurat_cell)) stop("some cells are listed in more than one fragment file")
no.frag <- setdiff(rownames(meta), fc$seurat_cell)
if (length(no.frag)) stop(length(no.frag), " cells have no fragment file (e.g. ", no.frag[1], ")")
bad <- unique(fc$fragments[!file.exists(fc$fragments)])
if (length(bad)) stop("fragment file(s) not found: ", paste(bad, collapse = ", "),
                      ". Fix the paths with Signac::UpdatePath() or re-link the Cell Ranger outputs.")
rm(obj); gc()

out <- data.frame(seurat_cell = rownames(meta), library = as.character(meta$orig.ident))
i <- match(out$seurat_cell, fc$seurat_cell)
out$barcode <- fc$barcode[i]
out$fragments <- fc$fragments[i]
# ArchR makes one Arrow file per library, so library <-> fragment file must be one-to-one
lib.files <- unique(out[, c("library", "fragments")])
if (anyDuplicated(lib.files$library) || anyDuplicated(lib.files$fragments)) {
  stop("orig.ident and fragment files are not one-to-one:\n",
       paste(capture.output(print(lib.files)), collapse = "\n"))
}
out$archr_cell <- paste0(out$library, "#", out$barcode)
if (anyDuplicated(out$archr_cell)) stop("duplicated <library>#<barcode> cell names")

out$group <- as.character(meta[[argv$group_col]])
out$split <- if (nzchar(split.col)) as.character(meta[[split.col]]) else NA_character_
out$rep <- do.call(paste, c(lapply(rep.cols, function(x) as.character(meta[[x]])), sep = "|"))

# peak groups: group (x split); small ones pooled across splits, still-small ones left out of calling
pg <- if (nzchar(split.col)) paste(out$group, out$split, sep = "_") else out$group
pg[is.na(out$group)] <- NA
n <- table(pg)
small <- !is.na(pg) & pg %in% names(n)[n < argv$min_group]
if (nzchar(split.col)) pg[small] <- paste(out$group[small], "all", sep = "_")
n2 <- table(pg)
pg[!is.na(pg) & pg %in% names(n2)[n2 < argv$min_group]] <- NA
# ArchR uses the group name in file paths: keep it to [A-Za-z0-9_-], and make sure that doesn't merge groups
clean <- gsub("[^A-Za-z0-9_-]+", "_", pg)
map <- unique(data.frame(pg = pg, clean = clean))
map <- map[!is.na(map$pg), ]
if (anyDuplicated(map$clean)) stop("peak group names collide after removing special characters: ",
                                   paste(map$pg[map$clean %in% map$clean[duplicated(map$clean)]], collapse = ", "))
out$peak_group <- clean

groups <- as.data.frame(table(peak_group = out$peak_group, rep = out$rep, useNA = "ifany"))
groups <- groups[groups$Freq > 0, ]
print(table(peak_group = out$peak_group, useNA = "ifany"))
message(sum(is.na(out$peak_group)), " cells left out of peak calling (group below ", argv$min_group, " cells or NA)")

write.csv(out, file.path(argv$out_dir, "cells.csv"), row.names = FALSE)
write.csv(groups, file.path(argv$out_dir, "peak-groups.csv"), row.names = FALSE)
message("wrote ", file.path(argv$out_dir, "cells.csv"), " (", nrow(out), " cells, ", nrow(lib.files), " libraries)")
