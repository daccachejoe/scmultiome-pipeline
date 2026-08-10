# Run this in terminal on interactive/batched job
# conda activate <sceasy_env_name>
# Rscript scripts/03_convert_seurat_to_h5ad.R arg1 arg2 arg3
# args1 = object to convert to h5ad format
# args2 = output file name
# args3 = conda env name

# only export is used right now, better to export meta data in csv file
source("scripts/lib/config.R")
load_pipeline_config()
if (nzchar(Sys.getenv("r_libs_personal_path"))) {
  .libPaths(c(.libPaths(), Sys.getenv("r_libs_personal_path")))
}
args = commandArgs(trailingOnly=TRUE)

library(sceasy)
library(reticulate)
use_condaenv(args[[3]])
library(Seurat)
library(Signac)
source("scripts/lib/seurat_io.R")

obj <- first_seurat(readRDS(args[[1]]))

DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = "RNA", dimreducs = c("pca","wnn.umap"))
sceasy::convertFormat(obj, from="seurat", to="anndata",
                      outFile=args[[2]])
