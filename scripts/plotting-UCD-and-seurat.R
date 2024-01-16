library(dplyr)
library(ggplot2)
library(argparser, quietly=TRUE)

cli <- commandArgs(trailingOnly = TRUE) 
# args <- strsplit(cli, "=", fixed = TRUE)

infile = cli[[1]]
resolution = cli[[2]]
ucDenv.MD <- read.csv(infile, row.names = 1)
plot <- 
    ucDenv.MD %>%
    ggplot(aes(x = resolution, 
        y = pred_celltype_ucdbase, 
        color = pred_celltype_ucdbase)) +
    geom_jitter() +
    theme_classic() +
    NoLegend() +
  theme(axis.text = element_text(color = "black", size = "black"))

pdf("./output/ucd/cluster-to-UCD-unbiased.pdf", height = 8, width = 8)
plot
dev.off()
