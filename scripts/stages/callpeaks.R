# Stage: callpeaks
# Calls ATAC peaks with MACS, either per-object (no --grouping.var) or
# grouped by a metadata column (e.g. cell type). Invoked from stage 01
# (per-sample) and stage 05 (grouped by cell type).
# Sourced by scripts/seurat_signac_pipeline.R.

obj.list <-
    lapply(obj.list,
        function(obj){
            if(argv$grouping.var == "NA"){
                obj <- CallPeaksMACS(obj, my.macs2.path=argv$my.macs.path, my.annotation = annotation, blacklist = blacklist.to.use)
            } else {
                obj <- CallPeaksMACS(obj, grouping.var=argv$grouping.var, my.macs2.path=argv$my.macs.path, my.annotation = annotation, blacklist = blacklist.to.use)
            }
            return(obj)
        })
# callpeaks is invoked both from stage 01 (per-sample, no grouping) and
# stage 05 (grouped by cell type, --grouping.var set)
callpeaks.step <- ifelse(argv$grouping.var == "NA", "01", "05")
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-", callpeaks.step, "-callpeaks-obj-list.RDS"))
