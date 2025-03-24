# Run this in terminal on interactive/batched job
# pc
# conda activate sceasy-joe
# Rscript scripts/convert-seurat-to-h5ad.R arg1
# args1 = object to convert to h5ad format
# args2 = output file name
# args3 = codna_env name

# only export is used right now, better to export meta data in csv file
.libPaths(c(.libPaths(),"/gpfs/data/naiklab/jd5457/R/x86_64-pc-linux-gnu-library/4.2"))
args = commandArgs(trailingOnly=TRUE)

library(sceasy)
library(reticulate)
use_condaenv(args[[3]])
library(Seurat)
library(Signac)

obj <- readRDS(args[[1]])
obj
if(typeof(obj) == "list"){
  obj <- obj[[1]]
}

DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = "RNA", dimreducs = c("pca","wnn.umap"))
sceasy::convertFormat(obj, from="seurat", to="anndata",
                      outFile=args[[2]])
