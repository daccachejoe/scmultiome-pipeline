# Stage: subcluster
# Splits the (single, merged) input object by cell.lineage then by
# orig.ident, filters and re-processes each lineage with Harmony batch
# correction across samples. Sourced by scripts/seurat_signac_pipeline.R.

message("Running sub-clustering Pipeline.")
obj.list <- SplitObject(obj.list[[1]], split.by = "cell.lineage")
obj.list <- lapply(obj.list, function(obj){
    message("Subclustering: ", obj$cell.lineage[1])
    obj@project.name <- obj$cell.lineage[1]
    mini.obj.list <- SplitObject(obj, split.by = "orig.ident")
    for(i in 1:length(mini.obj.list)){
        message("Subclustering: ", mini.obj.list[[i]]$cell.lineage[1], "-", mini.obj.list[[i]]$orig.ident[1])
        mini.obj.list[[i]] <- subset(mini.obj.list[[i]], percent.mt < 30 & nFeature_peaks >= 1000 & nFeature_SCT >= 500)
        mini.obj.list[[i]][["peaks"]]@fragments <- mini.obj.list[[i]][["peaks"]]@fragments[i]
        mini.obj.list[[i]] <- RunTFIDF(mini.obj.list[[i]], assay = "peaks")
        mini.obj.list[[i]] <- FindTopFeatures(mini.obj.list[[i]], min.cutoff = 'q50', assay = "peaks", verbose = FALSE)
    }
    var.peaks <- lapply(
        mini.obj.list, function(mini.obj){
        which.variable <- rownames(mini.obj[["peaks"]]) %in% VariableFeatures(mini.obj)
        return(mini.obj[["peaks"]]@ranges[which.variable, ])
    })
    residual.peaks <- lapply(var.peaks, paste0)
    residual.peaks.to.use <- names(table(unlist(residual.peaks)))[table(unlist(residual.peaks)) == 6]

    # residual.peaks <- GenomicRanges::reduce(unlist(GenomicRanges::GRangesList(var.peaks)))
    # residual.peaks.string <- paste0(residual.peaks)
    residual.peaks.string <- stringr::str_replace(residual.peaks.to.use, pattern = ":", replacement = "-")

    residual.features <- lapply(mini.obj.list,
                                function(mini.obj){
                                # mini.obj <- subset(mini.obj, percent.mt < 30 & nFeature_peaks >= 1000 & nFeature_SCT >= 500)
                                DefaultAssay(mini.obj) <- "RNA"
                                mini.obj <- SCTransform(mini.obj)
                                res.feats <- VariableFeatures(mini.obj)
                                return(res.feats)}
                                )
    residual.features.to.use <- names(table(unlist(residual.features)))[table(unlist(residual.features)) >= 5]

    merged.obj <- PreprocessAndReduceDims(obj,
            harmony = TRUE,
            harmony.vars = c("insitution", "orig.ident"),
            residual.features = residual.features.to.use,
            rna.pcs = 30,
            atac.pcs = 30,
            rna.theta = 10,
            atac.theta = 10,
            residual.peaks =  residual.peaks.string,
            vars.to.regress = c("insitution", "orig.ident"))

    merged.obj <- ConstructWNNGraph(merged.obj,
                                harmony = TRUE,
                                rna.pcs = 30,
                                atac.pcs = 30,
                                resolution = 0.4)

    message("Completed subclustering of ", merged.obj$cell.lineage[1])
    return(merged.obj)
})
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-02b-subcluster-obj-list.RDS"))
