# Stage: donors
# Writes patient identity into each per-sample object, between filter and
# merge in stage 02, so the merged object carries it. Sourced by
# scripts/seurat_signac_pipeline.R; route 02 adds it to the stage list
# whenever configs/demultiplexing_paths.csv or configs/donor_map.csv has rows,
# so a pooled project can't produce a merged object without donors.
#
# Per sample:
#   genotype-pooled (in configs/demultiplexing_paths.csv): the per-cell calls
#     from scripts/lib/demux.R (the demultiplex branch's combined_clusters.tsv,
#     else demux_path). Singlets get configs/donor_map.csv's patient ID for
#     (sample, cluster); other cells get "doublet" or "unassigned".
#   single-donor: configs/donor_map.csv's row (sample, "*").
# Why per-sample and not after merge: here the cell names are still the Cell
# Ranger barcodes and project.name is the sample, so the join needs no parsing
# of merge's cell-name prefixes. Why stage 02 and not stage 01: donor_map.csv
# is filled in by hand after the demultiplex branch's cross-library matching,
# and editing it should cost a stage 02 rerun, not stage 01's (AMULET etc.).
#
# configs/donor_map.csv: sampleName,cluster,donor (patient IDs consistent
# across libraries that share people; run/runmultiome demultiplex writes a
# suggested one to output/demultiplex/donor-map-suggested.csv).
#
# New metadata columns:
#   donor              patient ID | doublet | unassigned
#   donor.source       souporcell | pooled_rescue | pooled | single-donor library | none
#   souporcell.status  souporcell's own call (NA for single-donor libraries)
# Doublets are normally already gone (filter removes doublet.call ==
# "doublet", which includes genotype doublets), so few "doublet" rows remain.
# Stops if a pooled library's cluster, or a single-donor sample, has no
# donor_map row, or if a cell has no genotype call.
# Outputs: output/tables/{prefix}-02-donor-counts.csv,
#          output/RDS-files/{prefix}-02-donors-obj-list.RDS

message("Running Donor Assignment")
source("scripts/lib/demux.R")
demux.sheet <- ReadDemuxSheet()
donor.map.file <- "configs/donor_map.csv"
if (!file.exists(donor.map.file)) stop(donor.map.file, " not found (sampleName,cluster,donor); see README 'Donor demultiplexing'")
donor.map <- read.csv(donor.map.file, colClasses = "character")
if (!all(c("sampleName", "cluster", "donor") %in% colnames(donor.map))) {
    stop(donor.map.file, " needs columns sampleName,cluster,donor")
}
dup <- duplicated(donor.map[, c("sampleName", "cluster")])
if (any(dup)) stop(donor.map.file, " has duplicate (sampleName, cluster) rows: ",
                   paste(unique(donor.map$sampleName[dup]), collapse = ", "))

obj.list <- lapply(obj.list, function(seu) {
    sample.name <- seu@project.name
    rows <- donor.map[donor.map$sampleName == sample.name, , drop = FALSE]
    pooled <- !is.null(demux.sheet) && sample.name %in% demux.sheet$sampleName
    if (pooled) {
        calls <- ReadDemuxCalls(sample.name, demux.sheet)
        idx <- match(colnames(seu), calls$barcode)
        if (anyNA(idx)) {
            stop(sum(is.na(idx)), " of ", ncol(seu), " ", sample.name, " cells have no row in ", attr(calls, "path"),
                 " (e.g. ", colnames(seu)[is.na(idx)][1], ") -- was it made from this sample's Cell Ranger barcodes?")
        }
        calls <- calls[idx, ]
        clusters <- unique(calls$assignment[calls$status == "singlet"])
        unmapped <- setdiff(clusters, rows$cluster)
        if (length(unmapped) > 0) {
            stop(donor.map.file, " has no row for ", sample.name, " cluster(s) ", paste(unmapped, collapse = ", "),
                 " -- see output/demultiplex/donor-map-suggested.csv")
        }
        donor <- ifelse(calls$status == "singlet", rows$donor[match(calls$assignment, rows$cluster)], calls$status)
        seu$donor <- unname(donor)
        seu$donor.source <- unname(calls$source)
        seu$souporcell.status <- unname(calls$souporcell_status)
    } else {
        if (!"*" %in% rows$cluster) {
            stop(sample.name, " is not genotype-pooled (not in configs/demultiplexing_paths.csv) and has no ",
                 "'", sample.name, ",*,<donor>' row in ", donor.map.file)
        }
        seu$donor <- rep(rows$donor[rows$cluster == "*"], ncol(seu))
        seu$donor.source <- rep("single-donor library", ncol(seu))
        seu$souporcell.status <- rep(NA_character_, ncol(seu))
    }
    counts <- table(seu$donor)
    message("Donors: ", sample.name, " -- ", paste(names(counts), counts, sep = ": ", collapse = ", "))
    seu
})

donor.counts <- dplyr::bind_rows(lapply(obj.list, function(seu) {
    dplyr::count(seu@meta.data, sample = seu@project.name, donor, donor.source, name = "n.cells")
}))
write.csv(donor.counts, file = paste0("output/tables/", argv$project_prefix, "-02-donor-counts.csv"), row.names = FALSE)
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-02-donors-obj-list.RDS"))
