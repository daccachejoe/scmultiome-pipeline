# route: label_celltypes
library(Seurat)
library(dplyr)
library(ggplot2)
library(Signac)

cli <- commandArgs(trailingOnly = TRUE) 
annotations = cli[[1]]
resolution = cli[[2]]
object = cli[[3]]
project_prefix = cli[[4]]

obj <- readRDS(object)

ct.df <- read.csv(annotations)
obj$ct <- ct.df$ct[match(obj@meta.data[[resolution]], ct.df$cluster)]
obj$ct.spec <- ct.df$ct.spec[match(obj@meta.data[[resolution]], ct.df$cluster)]
saveRDS(list(obj), file = paste0("output/RDS-files/",project_prefix,"-annotated-obj-list.RDS"))

# plotting
md <- obj@meta.data
count.table <- 
  md %>%
  group_by(assignment, ct) %>%
  summarise(counts = n()) %>%
  mutate(perc = counts/sum(counts)) %>%
  filter(perc > 0.01)
p1 <- md %>%
  ggplot(aes(x = assignment, fill = ct)) +
  geom_bar(position = "fill") +
  geom_text(data = count.table,
            aes(label = paste0(counts, "\n", round(perc*100, digits = 0), "%"), y = perc),
            position = position_fill(vjust = 0.5)) +
  theme_classic() +
  NoLegend()
p2 <- md %>%
  ggplot(aes(x = assignment, fill = ct)) +
  geom_bar(position = "stack", stat = "count") +
  theme_classic()

pdf("output/plots/proportion-barplot-annotated-object.pdf", height = 8, width = 10)
p1 + p2
dev.off()

pdf("output/plots/dimplots-annotated-object.pdf", height = 8, width = 10)
p1 <- DimPlot(obj, reduction = "pca", group.by = "ct", label = F, label.size = 5, repel = FALSE) + ggtitle("PCA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(obj, reduction = "umap.rna", group.by = "ct", label = F, label.size = 5, repel = FALSE) + ggtitle("RNA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(obj, reduction = "umap.atac", group.by = "ct", label = F, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p4 <- DimPlot(obj, reduction = "wnn.umap", group.by = "ct", label = F, label.size = 5, repel = FALSE) + ggtitle("WNN") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T)

p1 <- DimPlot(obj, reduction = "pca", group.by = "ct.spec", label = F, label.size = 5, repel = FALSE) + ggtitle("PCA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(obj, reduction = "umap.rna", group.by = "ct.spec", label = F, label.size = 5, repel = FALSE) + ggtitle("RNA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p3 <- DimPlot(obj, reduction = "umap.atac", group.by = "ct.spec", label = F, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
p4 <- DimPlot(obj, reduction = "wnn.umap", group.by = "ct.spec", label = F, label.size = 5, repel = FALSE) + ggtitle("WNN") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T)
dev.off()