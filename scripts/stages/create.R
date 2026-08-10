# Stage: create
# Builds a Seurat object per sample in the samplesheet, attaches per-sample
# metadata columns and QC metrics (percent.mt, nucleosome signal, TSS
# enrichment). Sourced by scripts/seurat_signac_pipeline.R.

obj.list <-
    lapply(samples,
        function(sample){
            # create objects
            seu <- CreateMultiomeSeurat(data.dir = paste0("data/raw/",sample))

            seu@project.name <- sample
            seu$orig.ident <- sample
            # adding meta data columns to seurat object
            path.col <- grep("path", colnames(samplesheet))
            if(path.col < ncol(samplesheet)){
                for(md.col.to.add in colnames(samplesheet)[c((path.col+1):ncol(samplesheet))]){
                    seu[[paste0(md.col.to.add)]] <-
                        samplesheet[,md.col.to.add][samplesheet$sampleName == sample]
                }
            }
            # assay-specific metrics
            seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
            seu <- NucleosomeSignal(seu, assay = "ATAC")
            seu <- TSSEnrichment(seu, assay = "ATAC")
            return(seu)
        })
names(obj.list) <- samples
saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-01-create-obj-list.RDS"))
