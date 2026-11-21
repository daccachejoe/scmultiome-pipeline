#!/usr/bin/env Rscript
# =============================================================================
# ArchR peaks branch, step 2: ArchR reproducible peak set + peak matrix.
#
# Recipe (from prelim-long-data's hand-run archr-call-peaks.R and
# cross-tissue-build.R, made project-agnostic):
#   - Arrow files from each library's Cell Ranger ARC atac_fragments.tsv.gz,
#     restricted to step 1's cells (validBarcodes). Arrow QC is lenient
#     (minTSS = 0, minFrags = 1) because these cells already passed the
#     pipeline's QC; ArchR's own filters would drop cells we kept on purpose.
#   - Pseudo-replicates from step 1's `rep` (addGroupCoverages sampleLabels),
#     so donor-pooled libraries give one replicate per donor.
#   - MACS2 per peak_group, 501 bp fixed-width peaks (summit +- 250),
#     iterative overlap removal, kept if reproducible in --reproducibility
#     replicates. Cells with peak_group = NA are left out of calling.
#   - The peak set is then added to *all* cells and counted (PeakMatrix).
#
# MACS2 runs through macs2-wrapper.sh (same folder), which loads the MACS2
# module only inside the MACS2 subprocess: loading macs/2.1.0 into this R
# job swaps in an old libstdc++ and breaks Rcpp.
#
# Inputs:  <out_dir>/cells.csv (step 1)
# Outputs: <out_dir>/arrows/<library>.arrow, <out_dir>/project/ (ArchR project),
#          <out_dir>/peaks.rds, peaks.bed (BED6 + peakType, nearestGene,
#          distToTSS, GroupReplicate), peaks-annotated.csv
# Usage (normally via routes/downstream/archr_peaks.sh):
#   Rscript scripts/downstream/archr/02_call_peaks.R --out_dir output/archr/<tag> --threads 8
# =============================================================================

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(ArchR, quietly = TRUE)
})

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/genome.R")
set.seed(1234)

p <- arg_parser("ArchR peaks step 2: Arrow files -> pseudo-replicates -> MACS2 reproducible peaks -> peak matrix")
p <- add_argument(p, "--out_dir", help = "Run output directory (holds step 1's cells.csv)", type = "character")
p <- add_argument(p, "--threads", help = "ArchR threads", type = "integer", default = 8L)
p <- add_argument(p, "--macs2", help = "MACS2 executable or wrapper", type = "character",
                  default = "scripts/downstream/archr/macs2-wrapper.sh")
p <- add_argument(p, "--min_cells", help = "addGroupCoverages minCells", type = "integer", default = 40L)
p <- add_argument(p, "--max_cells", help = "addGroupCoverages maxCells", type = "integer", default = 500L)
p <- add_argument(p, "--reproducibility", help = "replicates a peak must appear in", type = "character", default = "2")
p <- add_argument(p, "--extend_summits", help = "bp each side of the summit (250 = 501 bp peaks)", type = "integer", default = 250L)
p <- add_argument(p, "--peaks_per_cell", help = "addReproduciblePeakSet peaksPerCell", type = "integer", default = 500L)
p <- add_argument(p, "--max_peaks", help = "addReproduciblePeakSet maxPeaks per group", type = "integer", default = 150000L)
argv <- parse_args(p)

out.dir <- normalizePath(argv$out_dir, mustWork = TRUE)
cells.csv <- file.path(out.dir, "cells.csv")
if (!file.exists(cells.csv)) stop(cells.csv, " not found: run step 1 (01_export_cells.R) first")
proj.dir <- file.path(out.dir, "project")
peaks.rds <- file.path(out.dir, "peaks.rds")
if (file.exists(peaks.rds)) stop(peaks.rds, " already exists. Use a new ARCHR_RUN_TAG for a different run.")

macs2 <- normalizePath(argv$macs2, mustWork = TRUE)   # ArchR calls it from other directories
macs.version <- suppressWarnings(system2(macs2, "--version", stdout = TRUE, stderr = TRUE))
if (!length(macs.version) || !grepl("macs", macs.version[1], ignore.case = TRUE)) {
  stop("MACS2 did not run via ", macs2, ": ", paste(macs.version, collapse = " "),
       ". Check archr_macs2_module / archr_macs2_path in config/pipeline.config.")
}
message("macs2: ", macs.version[1], "; threads: ", argv$threads)

g <- load_species_genome()
addArchRThreads(threads = argv$threads)
addArchRGenome(g$archr_genome)

meta <- read.csv(cells.csv, stringsAsFactors = FALSE)
libs <- unique(meta[, c("library", "fragments")])

# ---- Arrow files (restricted to our cells) ----
arrow.dir <- file.path(out.dir, "arrows")
dir.create(arrow.dir, showWarnings = FALSE, recursive = TRUE)
arrows <- file.path(arrow.dir, paste0(libs$library, ".arrow"))
valid <- setNames(lapply(libs$library, function(l) meta$barcode[meta$library == l]), libs$library)
# Arrow files are reused on a rerun (they're the slowest part), but only if they were made from the same cells
valid.csv <- file.path(arrow.dir, "valid-barcodes.csv")
valid.df <- data.frame(library = rep(names(valid), lengths(valid)), barcode = unlist(valid, use.names = FALSE))
if (all(file.exists(arrows))) {
  old <- if (file.exists(valid.csv)) read.csv(valid.csv, stringsAsFactors = FALSE) else NULL
  if (is.null(old) || !identical(sort(paste(old$library, old$barcode)), sort(paste(valid.df$library, valid.df$barcode)))) {
    stop("Arrow files in ", arrow.dir, " were made from different cells. Use a new ARCHR_RUN_TAG.")
  }
  message("reusing Arrow files in ", arrow.dir)
} else {
  createArrowFiles(inputFiles = setNames(libs$fragments, libs$library), sampleNames = libs$library,
                   outputNames = file.path(arrow.dir, libs$library), validBarcodes = valid,
                   QCDir = file.path(out.dir, "QualityControl"),
                   minTSS = 0, minFrags = 1, maxFrags = 1e7, addTileMat = TRUE, addGeneScoreMat = TRUE, force = TRUE)
  write.csv(valid.df, valid.csv, row.names = FALSE)
}

# ---- project + metadata ----
proj <- ArchRProject(ArrowFiles = arrows, outputDirectory = proj.dir, copyArrows = FALSE, showLogo = FALSE)
keep <- intersect(getCellNames(proj), meta$archr_cell)
message(length(getCellNames(proj)), " cells in Arrow files; ", length(keep), " of ", nrow(meta), " cells matched")
# a barcode-format mismatch (e.g. a suffix the fragment file doesn't use) shows up as almost nothing matching
if (length(keep) < 0.9 * nrow(meta)) stop("fewer than 90% of the cells in ", cells.csv, " matched ArchR cell names")
proj <- proj[keep, ]
m <- meta[match(getCellNames(proj), meta$archr_cell), ]
for (col in c("seurat_cell", "group", "split", "rep", "peak_group")) {
  proj <- addCellColData(proj, data = as.character(m[[col]]), name = col, cells = getCellNames(proj), force = TRUE)
}

# ---- pseudo-replicates + reproducible peaks ----
call <- !is.na(m$peak_group)
if (!any(call)) stop("no cells have a peak_group: lower --min_group or check the group column")
proj.call <- proj[getCellNames(proj)[call], ]
proj.call <- addGroupCoverages(proj.call, groupBy = "peak_group", sampleLabels = "rep",
                               minCells = argv$min_cells, maxCells = argv$max_cells,
                               minReplicates = 2, maxReplicates = 5, sampleRatio = 0.8, force = TRUE)
proj.call <- addReproduciblePeakSet(proj.call, groupBy = "peak_group", pathToMacs2 = macs2,
                                    reproducibility = argv$reproducibility, extendSummits = argv$extend_summits,
                                    peaksPerCell = argv$peaks_per_cell, maxPeaks = argv$max_peaks, force = TRUE)
peaks <- getPeakSet(proj.call)
message(length(peaks), " reproducible peaks; by type:")
print(table(peaks$peakType))

proj <- addPeakSet(proj, peakSet = peaks, force = TRUE)   # all cells, incl. those left out of calling
proj <- addPeakMatrix(proj, force = TRUE)
saveArchRProject(proj, outputDirectory = proj.dir, load = FALSE)

saveRDS(peaks, peaks.rds)
bed <- data.frame(chr = as.character(GenomicRanges::seqnames(peaks)), start = GenomicRanges::start(peaks) - 1,
                  end = GenomicRanges::end(peaks), name = paste0("peak", seq_along(peaks)), score = round(peaks$score, 3),
                  strand = ".", peakType = peaks$peakType, nearestGene = peaks$nearestGene, distToTSS = peaks$distToTSS,
                  GroupReplicate = peaks$GroupReplicate)
write.table(bed, file.path(out.dir, "peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(bed, file.path(out.dir, "peaks-annotated.csv"), row.names = FALSE)
message("wrote ", peaks.rds)
