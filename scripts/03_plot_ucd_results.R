source("scripts/lib/config.R")
load_pipeline_config()
if (nzchar(Sys.getenv("r_libs_personal_path"))) {
  .libPaths(c(.libPaths(), Sys.getenv("r_libs_personal_path")))
}

library(dplyr)
library(ggplot2)
library(SeuratObject) # for NoLegend()
library(argparser, quietly=TRUE)

cli <- commandArgs(trailingOnly = TRUE) 
# args <- strsplit(cli, "=", fixed = TRUE)

infile = cli[[1]]
resolution = cli[[2]]
ucDenv.MD <- read.csv(infile, row.names = 1)
plot <- 
    ucDenv.MD %>%
    # resolution is a column name passed as a string, so look it up with
    # .data[[]] -- aes(x = resolution) would plot the string itself
    ggplot(aes(x = .data[[make.names(resolution)]], 
        y = pred_celltype_ucdbase, 
        color = pred_celltype_ucdbase)) +
    geom_jitter() +
    theme_classic() +
    NoLegend() +
  theme(axis.text = element_text(color = "black"))

pdf("./output/ucd/cluster-to-UCD-unbiased.pdf", height = 8, width = 8)
plot
dev.off()
