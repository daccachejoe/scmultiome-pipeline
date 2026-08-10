#!/usr/bin/env Rscript
# Run LinkPeaks on a single group's object (produced by
# scripts/downstream/linkpeaks_split.R). Invoked once per group, in
# parallel, by routes/downstream/linkpeaks.sh.
#
# Usage: scripts/downstream/linkpeaks_group.R <group.RDS> <group_name> <out_prefix>

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))

source("scripts/lib/config.R")
load_pipeline_config()
source("scripts/lib/genome.R")
species.info <- load_species_genome()
peak.genome <- species.info$genome

source("scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)
input.rds <- args[[1]]
group.name <- args[[2]]
out.prefix <- args[[3]]

obj <- readRDS(input.rds)
obj <- LinkPeaksToGenes(obj,
                         genes = NULL,
                         distance.to.use = 250001,
                         peak.genome = peak.genome)

out.file <- paste0("output/RDS-files/", out.prefix, "-linkpeaks-group-", group.name, "-linked-obj.RDS")
saveRDS(obj, file = out.file)
message("Linked group '", group.name, "' written to ", out.file)
