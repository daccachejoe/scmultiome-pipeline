#!/usr/bin/env Rscript
# Pinned R package installation for the Seurat/Signac side of the pipeline
# (conda_env_name in config/pipeline.config). Run this once when setting up
# a new environment, then run renv::init() + renv::snapshot() (see
# docs/INSTALLATION.md) to produce a real, hash-verified renv.lock for your
# environment -- this script alone is not a substitute for that lockfile.

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Seurat v4 (pinned -- Seurat v5 changes several APIs this pipeline relies on)
remotes::install_version("SeuratObject", version = "4.1.4", upgrade = "never")
remotes::install_version("Seurat", version = "4.4.0", upgrade = "never")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Signac",
  "EnsDb.Hsapiens.v86",
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Mmusculus.v79",         # for species=mouse (config/pipeline.config)
  "BSgenome.Mmusculus.UCSC.mm10", # for species=mouse
  "motifmatchr",
  "TFBSTools",
  "JASPAR2020",
  "scDblFinder",                 # stage 01 doublet detection (RNA, ATAC mode, AMULET)
  "SingleCellExperiment"
), update = FALSE, ask = FALSE)

install.packages(c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "stringr",
  "reshape2",
  "readr",
  "clustree",
  "future",
  "future.apply",
  "enrichR",
  "SoupX",  # stage 01 ambient RNA correction
  "hdf5r",  # Read10X_h5 for Cell Ranger raw_feature_bc_matrix.h5
  "MASS" # only needed if scripts/density-scatter-function.R is wired in
))

remotes::install_github("immunogenomics/harmony")

# argparser is vendored at .lib/argparser -- installed by routes/00_install.sh
