# Stage: filter
# Removes cells per configs/qc_df.csv (whole clusters to drop, plus
# threshold-based filters on arbitrary metadata columns).
# Sourced by scripts/seurat_signac_pipeline.R.
# Note: this stage does not currently saveRDS its result (pre-existing gap,
# not something this reformatting pass changed -- see README "Known
# limitations").

message("Running Filtering Pipeline")
qc.df <- read.csv(file = argv$qc.sheet)

obj.list <- lapply(obj.list, function(seu) {
    md <- seu@meta.data
    message("Filtering: ", seu@project.name)

    # identify clusters to remove entirely
    clus.to.remove <- unique(as.character(qc.df$cluster.to.remove[which(seu@project.name == qc.df$sampleName)]))
    if(!(is.na(clus.to.remove))){
        clus.to.remove <- as.numeric(unlist(strsplit(clus.to.remove, split = ";")))
        cells.in.clusters.to.remove <- rownames(md)[md$seurat_clusters %in% clus.to.remove]
    } else {
        cells.in.clusters.to.remove <- c()
    }


    # identify variables to filter data by
    vars.to.filter.by <- as.character(qc.df$vars.to.filter.by[which(seu@project.name == qc.df$sampleName)])
    vars.to.filter.by <- unlist(strsplit(vars.to.filter.by, split = ";"))
    var.filter <- as.character(qc.df$var.filter[which(seu@project.name == qc.df$sampleName)])
    var.filter <- as.numeric(unlist(strsplit(var.filter, split = ";")))

    # check if each vars.to.filter.by is greater than or less than var.filter in the same position
    filter.direction <- as.character(qc.df$filter.direction[which(seu@project.name == qc.df$sampleName)])
    filter.direction <- as.character(unlist(strsplit(filter.direction, split = ";")))
    filter_list <- mapply(function(var, filter, direction) {
        var.vector <- md[[var]]
        if (direction == "greater") {
            return(var.vector > filter)
        } else {
            return(var.vector < filter)
        }
    }, vars.to.filter.by, var.filter, filter.direction)

    # combine all filters
    filter_vector <- Reduce(`|`, filter_list)
    cells.to.filter <- ifelse(length(filter_vector) > 0, rownames(md)[filter_vector], logical(0))
    if(length(cells.in.clusters.to.remove) == 0) {
        cells.to.remove <- cells.to.filter
    } else if (length(cells.to.filter) == 0) {
        cells.to.remove <- cells.in.clusters.to.remove
    } else {
        cells.to.remove <- c(cells.in.clusters.to.remove, cells.to.filter)
    }

    seu <- subset(seu, cells = cells.to.remove, invert = TRUE)
    return(seu)
})
