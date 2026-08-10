#!/usr/bin/env Rscript
# Run LinkPeaks on a single group's object (produced by
# scripts/06a_linkpeaks_split.R). Invoked once per group, in
# parallel, by routes/06_linkpeaks.sh.
#
# Usage: scripts/06b_linkpeaks_group.R <group.RDS> <group_name> <out_prefix>

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/genome.R")
source("scripts/lib/seurat_io.R")
species.info <- load_species_genome()
peak.genome <- species.info$genome

source("scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)
input.rds <- args[[1]]
group.name <- args[[2]]
out.prefix <- args[[3]]

obj <- first_seurat(readRDS(input.rds))
obj <- LinkPeaksToGenes(obj,
                         genes = NULL,
                         distance.to.use = 250001,
                         peak.genome = peak.genome)

out.file <- paste0("output/RDS-files/", out.prefix, "-06-linkpeaks-group-", group.name, "-linked-obj.RDS")
saveRDS(obj, file = out.file)
message("Linked group '", group.name, "' written to ", out.file)
