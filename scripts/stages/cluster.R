# Stage: cluster
# Runs PreprocessAndReduceDims + ConstructWNNGraph per object, writes
# per-object cluster resolution plots and marker gene tables.
# Sourced by scripts/seurat_signac_pipeline.R.

library(clustree, quietly=TRUE)
message("Running Clustering Pipeline.")
obj.list <- lapply(obj.list, function(obj){
    obj@project.name <- obj$cell.lineage[1]
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
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-02b-cluster-obj-list.RDS"))

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

    if(!(argv$SoupOrCellDF == "NA")){
        p5 <- DimPlot(obj, reduction = "pca", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("PCA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p6 <- DimPlot(obj, reduction = "umap.rna", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("RNA") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p7 <- DimPlot(obj, reduction = "umap.atac", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("ATAC") + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
        p8 <- DimPlot(obj, reduction = "wnn.umap", group.by = "assignment", label = FALSE, label.size = 5, repel = FALSE) + ggtitle("WNN")  + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
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
        M.list <-
        lapply(seq(0.1, 1, 0.1), function(res){
            Idents(obj) <- paste0("wsnn_res.", res)
            obj <- PrepSCTFindMarkers(obj)
            M <- FindAllMarkers(obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
            M <- M %>% mutate(resolution = res)
            write.csv(M, file = paste0("output/tables/cluster-markers-res.", res, ".csv"))
            return(M)
        })
        M <- bind_rows(M.list)
        write.csv(M, file = paste0("output/tables/cluster-markers-bound.csv"))
    })

saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix,"-02b-cluster-obj-list.RDS"))
