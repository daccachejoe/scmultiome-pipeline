# Stage: linkpeaks
# Runs LinkPeaksToGenes on each object in obj.list. For HPC runs, prefer
# the parallel divide-and-conquer path (routes/06_linkpeaks.sh +
# scripts/06a/06b/06c_linkpeaks_*.R) over invoking this stage directly --
# this single-process version is kept for small/interactive use.
# Sourced by scripts/seurat_signac_pipeline.R.

future::plan("sequential")
obj.list <-
    lapply(obj.list,
        function(obj){
                Idents(obj) <- argv$grouping.var
                # DefaultAssay(obj) <- "SCT"
                # obj <- PrepSCTFindMarkers(obj)
                # M <- FindAllMarkers(obj,
                #                     only.pos = TRUE)
                # write.csv(M, file = paste0("output/de-genes-",argv$grouping.var,".csv"), row.names = T)

                # genes.to.link <- M %>% filter(p_val_adj < 0.1) %>% arrange(cluster, desc(avg_log2FC)) %>% pull(gene)
                genes.to.link <- NULL
                obj <- LinkPeaksToGenes(obj,
                                    genes = genes.to.link,
                                    distance.to.use = 250001,
                                    peak.genome = peak.genome)
                gc()
                return(obj)
        })
saveRDS(obj.list, file = paste0("output/RDS-files/",argv$project_prefix, "-06-linkpeaks-obj-list.RDS"))
