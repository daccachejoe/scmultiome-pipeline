options(future.globals.maxSize = Inf)
future::plan("multicore", workers = as.numeric(future::availableCores()))
library(Signac)
library(Seurat)
library(readr)

cli <- commandArgs(trailingOnly = TRUE) 
infile = cli[[1]]

obj <- readRDS(infile)

# If type of object is not a list, make it into a list then run all functions as if it were a list
if (is.list(obj)) {
    obj <- obj[[1]]
}

# # Exporting the data for SCENICplus
mat <- obj[["peaks"]]@counts
Matrix::writeMM(mat,file='data/scenicplus/atac.txt')

readr::write_lines(colnames(mat), "data/scenicplus/cell_names.txt")
rownames(mat) <- stringr::str_replace(rownames(mat), pattern = "-", replacement = ":")
readr::write_lines(rownames(mat), "data/scenicplus/region_names.txt")

data.table::fwrite(obj@meta.data, 
                file = paste0("data/scenicplus/metadata.txt"),
                sep = "\t", 
                row.names = T, 
                col.names = T)

# making a bed peak file for scenic plus cistarget databases
df <- as.data.frame(obj[["peaks"]]@ranges)
data.table::fwrite(df[,1:3], file = "data/raw/macs-peaks/grouped-peaks.bed", sep="\t", col.names=F)

# exporting the first two pca dimensions for scenicplus
pca_df <- obj@reductions[["pca"]]@cell.embeddings
pca_df <- pca_df[,1:2]
data.table::fwrite(pca_df, file = "data/scenicplus/pca.txt", sep="\t", col.names=F)