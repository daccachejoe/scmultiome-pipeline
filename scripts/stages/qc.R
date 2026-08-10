# Stage: qc
# Generates QC plots (cross-sample metric violins, density scatter, per-
# object clustering violins), runs SCTransform/PCA/clustering per object.
# Sourced by scripts/seurat_signac_pipeline.R.

pdf(file = paste0("output/plots/", argv$project_prefix, "-qc-plots.pdf"),
    height = 8, width = 12)
list.of.vars <- list("1" = c("nCount_RNA",  "nCount_peaks"),
                 "2" = c("nFeature_RNA","nFeature_peaks"),
                 "3" = "percent.mt",
                 "4" = c("nucleosome_signal" ,"TSS.enrichment"))
# custom ggplot2 (not VlnPlot): compares metadata across obj.list, which
# is still a list of separate per-sample objects at this point (samples
# aren't merged until the "merge" stage), and VlnPlot requires a single
# object. p.list.2 below uses VlnPlot per-object once that's possible.
p.list <-
    lapply(list.of.vars, function(vars.to.plot){
        p <-
            bind_rows(lapply(obj.list, function(seu){return(seu@meta.data)})) %>%
                dplyr::select(orig.ident, all_of(vars.to.plot)) %>%
                reshape2::melt() %>%
                ggplot(aes(x = orig.ident, y = value, fill = orig.ident)) +
                geom_violin() +
                facet_grid(~variable) +
                # scale_fill_manual(values = c("#F7BE9F", "#EA7580")) +
                theme_classic() +
                theme(axis.text = element_text(color = "black"))
                if(!("percent.mt" %in% vars.to.plot | "nucleosome_signal" %in% vars.to.plot)){
                    p <- p +scale_y_log10()
                }
    return(p)
})
print(p.list)

density.scatter.list <-
    lapply(obj.list,
    function(seu){
        p <- FeatureScatter(seu, feature1 = 'nCount_peaks', feature2 = 'TSS.enrichment')
        return(p)
    })
print(density.scatter.list)

# preprocessing
obj.list <- lapply(obj.list,
function(seu){
    DefaultAssay(seu) <- "RNA"
    seu <- SCTransform(seu)
    seu <-
        seu %>%
        RunPCA() %>%
        FindNeighbors() %>%
        FindClusters()
    return(seu)
})

p.list.2 <-
    lapply(obj.list,
        function(seu){
        p <- VlnPlot(seu,
            pt.size = 0,
            features = unlist(list.of.vars),
            group.by = "seurat_clusters",
            stack = T,
            fill.by = "ident",
            log = T) +
            NoLegend() +
            ggtitle(seu@project.name)
        return(p)
        })
print(p.list.2)
dev.off()

saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-01-qc-obj-list.RDS"))
