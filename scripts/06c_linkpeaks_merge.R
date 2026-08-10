#!/usr/bin/env Rscript
# Merge the per-group linked objects (produced by
# scripts/06b_linkpeaks_group.R) back into a single named list, once
# all group jobs have completed. Invoked by routes/06_linkpeaks.sh
# as a job dependent on every group job finishing.
#
# Usage: scripts/06c_linkpeaks_merge.R <out_prefix> <comma_separated_group_names>

source("scripts/lib/seurat_io.R")

args <- commandArgs(trailingOnly = TRUE)
out.prefix <- args[[1]]
groups <- strsplit(args[[2]], ",")[[1]]
groups <- groups[nzchar(groups)]

linked.list <- lapply(groups, function(group.name) {
  in.file <- paste0("output/RDS-files/", out.prefix, "-06-linkpeaks-group-", group.name, "-linked-obj.RDS")
  if (!file.exists(in.file)) {
    stop("Missing expected group output: ", in.file, " (did the group job for '", group.name, "' fail?)")
  }
  first_seurat(readRDS(in.file))
})
names(linked.list) <- groups

out.file <- paste0("output/RDS-files/", out.prefix, "-06-linkpeaks-obj-list.RDS")
saveRDS(linked.list, file = out.file)
message("Merged ", length(linked.list), " group(s) into ", out.file)
