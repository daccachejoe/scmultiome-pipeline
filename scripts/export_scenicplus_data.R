library(Signac)
library(Seurat)

cli <- commandArgs(trailingOnly = TRUE) 
infile = cli[[1]]

grouped.peaks.obj <- readRDS(infile)
mat <- grouped.peaks.obj[[1]][["peaks"]]@counts
rownames(mat) <- stringr::str_replace(rownames(mat), pattern = "-", replacement = ":")
data.table::fwrite(as.data.frame(mat),
                row.names = T, 
                col.names = T, 
                sep = "\t",
                file=paste0('data/scenicplus/atac.txt'))
data.table::fwrite(grouped.peaks.obj[[1]]@meta.data, 
                file = paste0("data/scenicplus/metadata.txt"),
                sep = "\t", 
                row.names = T, 
                col.names = T)

# making a bed peak file for scenic plus cistarget databases
df <- as.data.frame(grouped.peaks.obj[[1]][["peaks"]]@ranges)
data.table::fwrite(df[,1:3], file = "data/raw/macs-peaks/grouped-peaks.bed", sep="\t", col.names=F)

# exporting the first two pca dimensions for scenicplus
pca_df <- grouped.peaks.obj@reductions[["pca"]]@cell.embeddings
pca_df <- pca_df[,1:2]
data.table::fwrite(pca_df, file = "data/scenicplus/pca.txt", sep="\t", col.names=F)