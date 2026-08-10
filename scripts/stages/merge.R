# Stage: merge
# Builds a consensus ATAC peak set across obj.list, re-quantifies each
# object against it, merges into one object, and re-runs
# PreprocessAndReduceDims/ConstructWNNGraph on the merged object.
# Sourced by scripts/seurat_signac_pipeline.R.

message("Running Merging Pipeline")
peak.list <- lapply(obj.list, function(seu){
    assay.to.use <- "ATAC"
    if(!(assay.to.use %in% names(seu@assays))){
        assay.to.use <- "peaks"
    }
    peak.granges <- seu[[assay.to.use]]@ranges
    return(peak.granges)
})

# intersecting the peak list
combined.peaks <- reduce(unlist(GRangesList(peak.list)))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")

frag.paths <- lapply(obj.list, function(seu){
    assay.to.use <- "ATAC"
    if(!(assay.to.use %in% names(seu@assays))){
        assay.to.use <- "peaks"
    }
    frag.path <- seu@assays[[assay.to.use]]@fragments[[1]]@path
    return(frag.path)
})

# done post-qc so we know which cells to keep already
frag.list <- lapply(names(frag.paths),function(f.path.name){
    f.path <- frag.paths[[f.path.name]]
    cells.to.keep <- colnames(obj.list[[f.path.name]])
    frag.obj <-
        CreateFragmentObject(path = f.path, cells = cells.to.keep)
    return(frag.obj)
})
names(frag.list) <- names(frag.paths)

# feature-matrix generation for peaks assay
feat.mat.list <- lapply(names(frag.paths),function(f.path.name){
    f.path <- frag.paths[[f.path.name]]
    cells.to.keep <- colnames(obj.list[[f.path.name]])
    frag.obj <- frag.list[[f.path.name]]
    frag.matrix <-
        FeatureMatrix(fragments = frag.obj,
                    features = combined.peaks,
                    cells = colnames(obj.list[[f.path.name]]))
    return(frag.matrix)
})
names(feat.mat.list) <- names(frag.paths)

# re-create the peaks assay using the new counts
obj.list <- lapply(names(frag.paths),function(f.path.name){
    seu <- obj.list[[f.path.name]]
    DefaultAssay(seu) <- "RNA"
    seu <- DietSeurat(seu, assays = c("RNA", "SCT"))
    seu[["peaks"]] <-
        CreateChromatinAssay(
                counts = feat.mat.list[[f.path.name]],
                fragments = frag.list[[f.path.name]],
                annotation = annotation)
    return(seu)
})

# merge the objects now with a combined peaks set
if(length(obj.list) == 2){
    merged.obj <- merge(obj.list[[1]], obj.list[[2]])
} else {
    merged.obj <- merge(obj.list[[1]], obj.list[c(2:length(obj.list))])
}

# create union of variable genes
residual.features.to.use <- Reduce(intersect, lapply(obj.list, function(obj){return(obj[["SCT"]]@var.features)}))

# # create union of variable peaks
# top.feats.list <- lapply(obj.list, function(seu){
#     DefaultAssay(seu) <- "peaks"
#     seu <- RunTFIDF(seu)
#     seu <- FindTopFeatures(seu, min.cutoff = "q50")
#     var.peaks <- seu[["peaks"]]@var.features
#     return(var.peaks)
# })
# residual.peaks.to.use <- Reduce(intersect, top.feats.list)
# message("There are: ", length(residual.peaks.to.use), " residual peaks found.")
residual.peaks.to.use <- NULL

# rerun the clustering on the merged object
if(argv$RunHarmony){
    message("RunHarmony flagged as TRUE")
    my.harmony.vars <- "orig.ident"
} else {
    message("RunHarmony flagged as not TRUE")
    my.harmony.vars <- NULL
}
merged.obj <- PreprocessAndReduceDims(merged.obj,
    harmony = argv$RunHarmony,
    harmony.vars = my.harmony.vars,
    residual.features = residual.features.to.use,
    residual.peaks =  residual.peaks.to.use,
    vars.to.regress = NULL)
merged.obj <- ConstructWNNGraph(merged.obj,
                                harmony = argv$RunHarmony,
                                resolution = seq(0,1,0.1))
merged.obj <- list(merged.obj)
saveRDS(merged.obj, file = paste0("output/RDS-files/", argv$project_prefix,"-02-merge-obj-list.RDS"))
