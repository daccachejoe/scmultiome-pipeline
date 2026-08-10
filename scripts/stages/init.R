# Stage: init
# Copies raw CellRanger output into data/raw/<sample>/, combining GEX and
# ATAC fragment files into the layout CreateMultiomeSeurat() expects.
# Sourced by scripts/seurat_signac_pipeline.R when "init" is in argv$pipeline.

message("Begining scMulitome processing")
lapply(samples,
        function(sample){
        data.dir <- samplesheet$path[samplesheet$sampleName == sample]
        new.dir <- paste0("data/raw/", sample)
        dir.create(new.dir)
        CombineDirectories(data.dir, new.dir, sample)
})
