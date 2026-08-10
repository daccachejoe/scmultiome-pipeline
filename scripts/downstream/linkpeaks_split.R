#!/usr/bin/env Rscript
# Split a merged/annotated object by a grouping variable into one RDS per
# group, so LinkPeaks can be run in parallel, one job per group (see
# routes/downstream/linkpeaks.sh, which orchestrates split -> per-group -> merge).
#
# Usage: scripts/downstream/linkpeaks_split.R <input.RDS> <grouping.var> <out_prefix>

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
source("scripts/lib/seurat_io.R")

args <- commandArgs(trailingOnly = TRUE)
input.rds <- args[[1]]
grouping.var <- args[[2]]
out.prefix <- args[[3]]

obj.list <- readRDS(input.rds)
obj <- first_seurat(obj.list)

Idents(obj) <- obj[[grouping.var]][, 1]
groups.obj.list <- SplitObject(obj, split.by = grouping.var)

# sanitize group names for use in filenames
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

groups.used <- c()
for (group.name in names(groups.obj.list)) {
  safe.name <- sanitize(group.name)
  out.file <- paste0("output/RDS-files/", out.prefix, "-linkpeaks-group-", safe.name, "-obj.RDS")
  saveRDS(groups.obj.list[[group.name]], file = out.file)
  groups.used <- c(groups.used, safe.name)
  message("Wrote group '", group.name, "' (", ncol(groups.obj.list[[group.name]]), " cells) to ", out.file)
}

groups.file <- paste0("output/RDS-files/", out.prefix, "-linkpeaks-groups.txt")
writeLines(groups.used, groups.file)
message("Wrote group list to ", groups.file)
