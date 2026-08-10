#!/usr/bin/env Rscript

# Entry point / dispatcher for the Seurat+Signac side of the pipeline.
# Parses CLI args, loads config/species/shared objects, then sources the
# requested stage(s) from scripts/stages/. Each stage is a self-contained
# script that reads/writes obj.list (and other shared objects set up below)
# -- see scripts/stages/*.R for the actual stage logic.
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Run single-cell RNA + ATAC Multiomic Analysis from 10X Genomics platform")

# Add required command line arguments
p <- add_argument(p, "pipeline",
                help="Comma delimted combinations of: init, create, callpeaks, qc, filter, cluster, subcluster, merge, linkpeaks, footprint",
                type="character")
p <- add_argument(p, "samplesheet", help="samplesheet in csv format", type="character")

# Add optional command line flags
p <- add_argument(p, "--project_prefix", help="outfile name, no .RDS!", type="character", default="multiome")
p <- add_argument(p, "--grouping.var", help="grouping variable for peak calling algorithm", type = "character", default="NA")
p <- add_argument(p, "--RDS.file.in", help="RDS in-file for the pipeline desired", default="NA")
p <- add_argument(p, "--RunHarmony",flag=TRUE, help="Run Harmony batch correction")
p <- add_argument(p, "--footprint-peaks", help="atac peaks to run footprinting analysis on. txt file", type="character")
p <- add_argument(p, "--my.macs.path", help="path to macs environment. only used if callpaks pipeline is run", type="character")
p <- add_argument(p, "--qc.sheet", help="path to csv file containing clusters to remove per sample", type="character")
p <- add_argument(p, "--SoupOrCellDF", help="path to csv file containing barcodes and their assigned samples", type="character", default="NA")
p <- add_argument(p, "--qc.split", flag=TRUE, help="whether or not to perform qc clustering on objects individually or combined")


# Parse the command line arguments
argv <- parse_args(p)

pipelines.to.run <- unlist(strsplit(argv$pipeline, split = ","))
samplesheet <- argv$samplesheet

library(future, quietly=TRUE)
library(future.apply, quietly=TRUE)

options(future.globals.maxSize = Inf)
options(future.rng.onMisuse = "ignore")
# parrallelize the processing
plan("multicore", workers = as.numeric(future::availableCores()))
plan()

typeof(argv$RunHarmony)

# load pipeline config (config/pipeline.config) for species/genome selection
# and other machinery settings, then set up species-specific annotation
# packages, genome, blacklist, and JASPAR taxon ID.
source("scripts/lib/config.R")
load_pipeline_config()

library(Seurat, quietly=TRUE)
library(Signac, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(ggplot2, quietly=TRUE)
# library(enrichR, quietly=TRUE)

source("scripts/lib/genome.R")
source("scripts/lib/seurat_io.R")
species.info <- load_species_genome()
ensdb.to.use <- species.info$ensdb
genome.to.use <- species.info$genome
blacklist.to.use <- species.info$blacklist
jaspar.taxid <- species.info$jaspar_taxid

# source in wrapper functions
source("scripts/functions.R")

# being processing
# always read in the samplesheet for referemnce
samplesheet <- read.csv(samplesheet)
samples <- samplesheet$sampleName

# load in annotation files for peaks and ranges
annotation <- GetGRangesFromEnsDb(ensdb = ensdb.to.use)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
peak.genome <- genome.to.use
# seqlevels(peak.genome) <- paste0('chr', seqlevels(peak.genome))

# load in data if desired
if(!(argv$RDS.file.in == "NA")){
    message("Loading in: ", argv$RDS.file.in)
    obj.list <- readRDS(argv$RDS.file.in)
    obj.list <- as_seurat_list(obj.list)
}

# Stages run in this fixed order regardless of the order they're listed in
# argv$pipeline -- matches the original monolithic script's behavior.
if("init" %in% pipelines.to.run){
    message("=== Stage: init ===")
    source("scripts/stages/init.R")
}

if("create" %in% pipelines.to.run){
    message("=== Stage: create ===")
    source("scripts/stages/create.R")
}

# Splitting up object(s) by SouporCell called assignment. Not a named
# pipeline stage -- runs whenever --SoupOrCellDF is set, between create
# and callpeaks, same as in the original script.
if(!(argv$SoupOrCellDF == "NA")){
    souporcelldf <- read.csv(argv$SoupOrCellDF)
    obj.list <-
        lapply(obj.list,
        function(obj){
            souporcelldf <- souporcelldf[souporcelldf$barcode %in% colnames(obj), ] # This should be 100% for one object, less if one data frame is for multiple samplesheet instances
            cells.not.in.SoC <- colnames(obj)[!(colnames(obj) %in% souporcelldf$barcode)]
            residual.df <- data.frame(barcode = cells.not.in.SoC, assignment = "NA")
            souporcelldf <- rbind(souporcelldf, residual.df)
            rownames(souporcelldf) <- souporcelldf$barcode
            obj <- AddMetaData(obj, metadata = souporcelldf)
            obj <- subset(obj, assignment == "NA", invert = TRUE)

            # split object into multiple objects based on assignment
            if(argv$qc.split){
                message("Splitting object into ", length(unique(souporcelldf$assignment))-1, " objects")
                mini.obj.list <- SplitObject(obj, split = "assignment")
                mini.obj.list <- lapply(mini.obj.list,
                    function(mini.obj){
                            mini.obj@project.name <- paste0("sample ", mini.obj$assignment[1])
                            mini.obj$orig.ident <- paste0("control.skin.", mini.obj$assignment[1])
                            return(mini.obj)})
                return(mini.obj.list)
            }

            return(obj)
        })
    obj.list <- unlist(obj.list, recursive = FALSE)
}

if("callpeaks" %in% pipelines.to.run){
    message("=== Stage: callpeaks ===")
    source("scripts/stages/callpeaks.R")
}

if("qc" %in% pipelines.to.run){
    message("=== Stage: qc ===")
    source("scripts/stages/qc.R")
}

if("filter" %in% pipelines.to.run){
    message("=== Stage: filter ===")
    source("scripts/stages/filter.R")
}

if("cluster" %in% pipelines.to.run){
    message("=== Stage: cluster ===")
    source("scripts/stages/cluster.R")
}

if("subcluster" %in% pipelines.to.run){
    message("=== Stage: subcluster ===")
    source("scripts/stages/subcluster.R")
}

if("merge" %in% pipelines.to.run){
    message("=== Stage: merge ===")
    source("scripts/stages/merge.R")
}

if("linkpeaks" %in% pipelines.to.run){
    message("=== Stage: linkpeaks ===")
    source("scripts/stages/linkpeaks.R")
}

if("footprint" %in% pipelines.to.run){
    message("=== Stage: footprint ===")
    source("scripts/stages/footprint.R")
}

message("disco complete.")
