#!/usr/bin/env Rscript
# =============================================================================
# ArchR peaks branch, step 3: put the ArchR peak counts back into Seurat.
#
# Reads step 2's PeakMatrix (Tn5 insertions per peak, not binarized) and
# adds it to the input object as a new ChromatinAssay (--assay_name,
# default peaks_archr). The existing ATAC assay is left as it is, so MACS
# (Signac) and ArchR peaks can be compared on the same cells. Fragment
# files and gene annotation are copied from that existing assay, so
# CoveragePlot etc. work on the new assay too.
#
# Cells that ArchR didn't count (missing from the Arrow files) are dropped
# from the output object, with a warning naming how many.
#
# Inputs:  the same Seurat RDS step 1 read; <out_dir>/cells.csv; <out_dir>/project/
# Outputs: --out_rds (a one-object list, like the trunk checkpoints)
# Usage (normally via routes/downstream/archr_peaks.sh):
#   Rscript scripts/downstream/archr/03_add_peak_assay.R <rds> --out_dir output/archr/<tag> \
#     --out_rds output/RDS-files/<prefix>-archr-<tag>-obj-list.RDS
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(ArchR, quietly = TRUE)
  library(Seurat, quietly = TRUE)
  library(Signac, quietly = TRUE)
})

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/seurat_io.R")

p <- arg_parser("ArchR peaks step 3: add the ArchR peak matrix to the Seurat object")
p <- add_argument(p, "seurat_rds", help = "Seurat RDS (the one step 1 read)", type = "character")
p <- add_argument(p, "--out_dir", help = "Run output directory (steps 1-2)", type = "character")
p <- add_argument(p, "--out_rds", help = "Output RDS", type = "character")
p <- add_argument(p, "--atac_assay", help = "Existing ATAC assay to copy annotation from ('auto' = ATAC if present, else peaks)",
                  type = "character", default = "auto")
p <- add_argument(p, "--assay_name", help = "Name of the new assay", type = "character",
                  default = Sys.getenv("archr_assay_name", unset = "peaks_archr"))
argv <- parse_args(p)

if (file.exists(argv$out_rds)) stop(argv$out_rds, " already exists. Use a new ARCHR_RUN_TAG.")
proj.dir <- file.path(argv$out_dir, "project")
if (!dir.exists(proj.dir)) stop(proj.dir, " not found: run step 2 (02_call_peaks.R) first")
cells <- read.csv(file.path(argv$out_dir, "cells.csv"), stringsAsFactors = FALSE)

addArchRThreads(threads = 1)
proj <- loadArchRProject(proj.dir, showLogo = FALSE)
pm <- getMatrixFromProject(proj, useMatrix = "PeakMatrix", binarize = FALSE)
pk <- SummarizedExperiment::rowRanges(pm)
pk <- GenomicRanges::GRanges(GenomicRanges::seqnames(pk), IRanges::IRanges(GenomicRanges::start(pk), GenomicRanges::end(pk)))
counts <- SummarizedExperiment::assay(pm, "PeakMatrix")
rownames(counts) <- GRangesToString(pk, sep = c("-", "-"))
colnames(counts) <- cells$seurat_cell[match(colnames(counts), cells$archr_cell)]
if (anyNA(colnames(counts))) stop("ArchR project has cells that are not in ", file.path(argv$out_dir, "cells.csv"))
rm(pm, proj); gc()

obj.list <- as_seurat_list(readRDS(argv$seurat_rds))
if (length(obj.list) != 1) stop(argv$seurat_rds, " holds ", length(obj.list), " objects; expected one merged object")
obj <- obj.list[[1]]
rm(obj.list)
assay <- argv$atac_assay
if (assay == "auto") assay <- if ("ATAC" %in% Assays(obj)) "ATAC" else "peaks"
if (!assay %in% Assays(obj)) stop("assay '", assay, "' not in ", argv$seurat_rds)
if (argv$assay_name %in% Assays(obj)) stop("assay '", argv$assay_name, "' already exists; set ARCHR_ASSAY_NAME")

keep <- intersect(colnames(obj), colnames(counts))
message(length(keep), " of ", ncol(obj), " cells have ArchR peak counts; ", nrow(counts), " peaks")
if (length(keep) < ncol(obj)) {
  warning(ncol(obj) - length(keep), " cells without ArchR peak counts are dropped from the output object")
  obj <- subset(obj, cells = keep)
}
counts <- counts[, colnames(obj)]

cells <- cells[match(colnames(obj), cells$seurat_cell), ]
frags <- lapply(split(cells, cells$fragments), function(d)
  CreateFragmentObject(path = d$fragments[1], cells = setNames(d$barcode, d$seurat_cell), validate.fragments = FALSE))
obj[[argv$assay_name]] <- CreateChromatinAssay(counts = counts, sep = c("-", "-"), fragments = unname(frags),
                                               annotation = Annotation(obj[[assay]]))

dir.create(dirname(argv$out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(obj), argv$out_rds)
message("wrote ", argv$out_rds, " (assay '", argv$assay_name, "')")
