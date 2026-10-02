# Genotype demultiplexing helpers, shared by stage 01 doublets and stage 02 donors.
#
# configs/demultiplexing_paths.csv lists the genotype-pooled libraries
# (sampleName, demux_path, optional n_donors / min_rna_margin / gex_bam). Their
# per-cell calls come from the demultiplex branch (route 00b):
#   output/demultiplex/<sample>/combined_clusters.tsv   souporcell + pooled rescue
# and, when that doesn't exist, from demux_path: a souporcell clusters.tsv made
# outside the pipeline. Both have barcode / status / assignment columns;
# combined_clusters.tsv adds source and souporcell_status.

ReadDemuxSheet <- function(path = "configs/demultiplexing_paths.csv") {
    if (!file.exists(path)) return(NULL)
    sheet <- read.csv(path, stringsAsFactors = FALSE, colClasses = "character", na.strings = c("", "NA"))
    if (nrow(sheet) == 0) return(NULL)
    if (!all(c("sampleName", "demux_path") %in% colnames(sheet))) {
        stop(path, " needs columns sampleName,demux_path (got: ", paste(colnames(sheet), collapse = ", "), ")")
    }
    sheet
}

# per-cell calls for one pooled library, or stop() naming what's missing
ReadDemuxCalls <- function(sample.name, sheet) {
    combined <- file.path("output/demultiplex", sample.name, "combined_clusters.tsv")
    path <- if (file.exists(combined)) combined else sheet$demux_path[sheet$sampleName == sample.name][1]
    if (is.na(path) || !file.exists(path)) {
        stop(sample.name, " is in configs/demultiplexing_paths.csv but has no calls: neither ", combined,
             " (run/runmultiome demultiplex) nor a demux_path file")
    }
    calls <- read.delim(path, colClasses = "character")
    if (!all(c("barcode", "status", "assignment") %in% colnames(calls))) {
        stop(path, " needs barcode/status/assignment columns (souporcell clusters.tsv format)")
    }
    if (!"source" %in% colnames(calls)) calls$source <- ifelse(calls$status == "unassigned", "none", "souporcell")
    if (!"souporcell_status" %in% colnames(calls)) calls$souporcell_status <- calls$status
    attr(calls, "path") <- path
    calls
}
