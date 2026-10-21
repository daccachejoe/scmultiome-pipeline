# Stage: cluster
# Runs PreprocessAndReduceDims + ConstructWNNGraph per object, writes
# per-object cluster resolution plots and marker gene tables.
# Route 02 uses it instead of merge for single-sample projects.
# Sourced by scripts/seurat_signac_pipeline.R.

library(clustree, quietly=TRUE)
message("Running Clustering Pipeline.")
obj.list <- lapply(obj.list, function(obj){
    # lineage objects are named by lineage, single samples by sample
    if("cell.lineage" %in% colnames(obj@meta.data)){
        obj@project.name <- obj$cell.lineage[1]
    } else {
        obj@project.name <- obj$orig.ident[1]
    }
    return(obj)
})

obj.list <-
    lapply(obj.list,
        function(obj){
            obj <- PreprocessAndReduceDims(obj,
                                            harmony = argv$RunHarmony,
                                            harmony.vars = "orig.ident")
            obj <- ConstructWNNGraph(obj,
                                    harmony = argv$RunHarmony,
                                    resolution = seq(0,1,0.1))
            return(obj)
        })
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-02-cluster-obj-list.RDS"))

# plots to help decide resolution to use
lapply(obj.list, function(obj){
    pdf(file = paste0("output/plots/", argv$project_prefix, "-",obj@project.name, "-cluster-plots.pdf"),
        height = 8, width = 12)

    # print(clustree(obj@meta.data, prefix = "wsnn_res."))
    print(DimPlot(obj, reduction = "wnn.umap", group.by = paste0("wsnn_res.", seq(0.1, 1, 0.1)), label = T) & NoLegend())
    p1 <- DimPlot(obj, reduction = "pca", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("PCA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
    p2 <- DimPlot(obj, reduction = "umap.rna", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("RNA")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
    p3 <- DimPlot(obj, reduction = "umap.atac", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
    p4 <- DimPlot(obj, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 5, repel = FALSE) + ggtitle("WNN") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))

    # per-donor panels when stage 02 donors ran (genotype-pooled libraries)
    if("donor" %in% colnames(obj@meta.data)){
        p5 <- DimPlot(obj, reduction = "pca", group.by = "donor", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("PCA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p6 <- DimPlot(obj, reduction = "umap.rna", group.by = "donor", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("RNA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p7 <- DimPlot(obj, reduction = "umap.atac", group.by = "donor", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p8 <- DimPlot(obj, reduction = "wnn.umap", group.by = "donor", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("WNN")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        print(
            ggpubr::ggarrange(
                ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1),
                ggpubr::ggarrange(p5, p6, p7, p8, ncol = 4, nrow = 1),
                ncol = 1, nrow = 2)
        )
    } else {
        print(ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1))
    }
    dev.off()
})

# calculate marker genes at each resolution
lapply(obj.list,
    function(obj){
        obj <- PrepSCTFindMarkers(obj)
        M.list <-
        lapply(seq(0.1, 1, 0.1), function(res){
            Idents(obj) <- paste0("wsnn_res.", res)
            M <- FindAllMarkers(obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
            M <- M %>% mutate(resolution = res) %>% arrange(cluster, p_val_adj)
            write.csv(M, file = paste0("output/tables/", argv$project_prefix, "-", obj@project.name, "-02-cluster-markers-res.", res, ".csv"))
            return(M)
        })
        M <- bind_rows(M.list)
        write.csv(M, file = paste0("output/tables/", argv$project_prefix, "-", obj@project.name, "-02-cluster-markers-bound.csv"))
    })
