# Stage: footprint
# Footprints TF motif activity at a set of peaks (from --footprint-peaks)
# across every object in obj.list. Sourced by scripts/seurat_signac_pipeline.R.

library(motifmatchr, quietly=TRUE)
library(TFBSTools, quietly=TRUE)
library(JASPAR2020, quietly=TRUE)

peaks.to.footprint <- readr::read_delim(argv[["footprint-peaks"]])
footprint.list <-
    lapply(obj.list,
        function(obj){
            footprint.obj <-
                FootprintMotifs(obj,
                    peaks.to.test = peaks.to.footprint,
                    peak.genome = peak.genome,
                    jaspar.taxid = jaspar.taxid)
        })
saveRDS(footprint.list, file = paste0("output/RDS-files/", argv$project_prefix, "-footprinted-obj-list.RDS"))
