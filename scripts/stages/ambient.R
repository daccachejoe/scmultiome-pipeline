# Stage: ambient
# Ambient RNA correction (SoupX) and an ambient DNA (ATAC) estimate, per
# sample (= per 10x capture). Sourced by scripts/seurat_signac_pipeline.R
# right after create and before doublets, so doublet scoring sees the
# corrected RNA counts: ambient transcripts make every droplet look a little
# like a mix of cell types, which is exactly what scDblFinder's simulated
# doublets look like.
#
# Inputs: Cell Ranger's raw_feature_bc_matrix (.h5 or directory), found at
# the samplesheet `path` or in data/raw/<sample>/. Empty droplets in it are
# where the ambient ("soup") profile comes from. Samples without one are
# skipped with a warning and stay uncorrected (listed in the summary table).
# DecontX, which needs no raw matrix, was considered and not added: one
# method across samples keeps the correction comparable.
#
# RNA. SoupX: soup profile from droplets with 0 < nUMI < 100 (SoupX's
# default), contamination fraction rho from autoEstCont() using quick
# per-sample clusters, then adjustCounts(roundToInt = TRUE) so scDblFinder
# and SCTransform get integer counts. rho can be fixed for every sample with
# ambient_rna_rho (use when autoEstCont fails or gives an implausible rho).
# With ambient_rna_correction=true (default) the corrected counts replace
# the RNA assay's counts and the originals go to an "RNA.raw" assay, so
# every later stage uses corrected counts without changes. nCount_RNA /
# nFeature_RNA are recomputed; the originals are kept as *.raw columns.
# percent.mt (from create) stays on the uncorrected counts: it is a
# damaged-cell metric and SoupX's removal would shift it. With
# ambient_rna_correction=false, rho is estimated and reported only.
#
# ATAC (estimate only, counts untouched). There is no established ambient
# correction for ATAC counts: per-cell counts are near-binary, and a
# SoupX-style subtraction on them is unvalidated and would remove real
# accessibility. Instead:
#   ambient.atac.sim  cosine similarity between a cell's binarized peak
#                     profile and the pooled peak profile of the empty
#                     droplets above
#   ambient.atac.z    that similarity as a robust z-score within the cell's
#                     quick cluster -- the raw similarity mostly tracks cell
#                     type (the soup is dominated by abundant types), so
#                     only the within-cluster deviation says "more
#                     ambient-like than its peers" (contaminated, or
#                     empty/low-quality)
# plus, per sample, the fraction of RNA UMIs / ATAC peak counts that fall
# outside the cells (ambient load). Use ambient.atac.z as a qc_df filter
# variable if a sample needs it; nothing is removed automatically.
#
# Outputs:
#   output/tables/{prefix}-01-ambient-summary.csv  per-sample method, rho,
#       % RNA UMIs removed, ambient loads, top soup genes
#   output/plots/{prefix}-01-ambient-plots.pdf     per-sample RNA removal and
#       ATAC ambient z by quick cluster
#   output/RDS-files/{prefix}-01-ambient-obj-list.RDS

suppressPackageStartupMessages(library(SoupX))

# Run this stage sequentially, then restore the engine's plan at the end.
# Under plan("multicore"), forking workers here after create has run in the
# same process corrupts R's state: the job hangs to its run limit with
# "stack imbalance" / "no more error handlers" (reproduced on 5376-AB-2 and
# the PSO longitudinal libraries; the same stages pass on 1 core, with
# R_FUTURE_FORK_ENABLE=false, or when resumed from the create checkpoint).
# Nothing here needs parallelism: quick clustering takes ~15 s per sample.
ambient.engine.plan <- future::plan("sequential")

correct.rna <- !tolower(Sys.getenv("ambient_rna_correction", unset = "true")) %in% c("false", "f", "0", "no")
fixed.rho <- Sys.getenv("ambient_rna_rho", unset = "")
if (nzchar(fixed.rho)) {
    fixed.rho <- as.numeric(fixed.rho)
    if (is.na(fixed.rho) || fixed.rho < 0 || fixed.rho >= 1) {
        stop("ambient_rna_rho in config/pipeline.config must be empty (auto-estimate) or a number in [0, 1), got: ",
             Sys.getenv("ambient_rna_rho"))
    }
} else {
    fixed.rho <- NA
}

# first raw matrix found for a sample, or NULL
FindRawMatrix <- function(sample.name) {
    dirs <- c(samplesheet$path[samplesheet$sampleName == sample.name], file.path("data/raw", sample.name))
    candidates <- c(rbind(file.path(dirs, "raw_feature_bc_matrix.h5"), file.path(dirs, "raw_feature_bc_matrix")))
    found <- candidates[file.exists(candidates)]
    if (length(found) == 0) return(NULL)
    found[1]
}

ReadRawMatrix <- function(path) {
    mat.list <- if (grepl("\\.h5$", path)) Read10X_h5(path) else Read10X(path)
    if (!is.list(mat.list) || !all(c("Gene Expression", "Peaks") %in% names(mat.list))) {
        stop(path, " is not a Cell Ranger ARC raw matrix (needs 'Gene Expression' and 'Peaks' features)")
    }
    # CreateSeuratObject turns "_" into "-" in feature names; match that
    rownames(mat.list[["Gene Expression"]]) <- gsub("_", "-", rownames(mat.list[["Gene Expression"]]))
    mat.list
}

# fine clusters for SoupX's marker search (and for the ATAC z-score); a
# throwaway log-normalized object so seu itself isn't touched
QuickClusters <- function(counts) {
    tmp <- CreateSeuratObject(counts = counts)
    tmp <- NormalizeData(tmp, verbose = FALSE)
    tmp <- FindVariableFeatures(tmp, nfeatures = 2000, verbose = FALSE)
    tmp <- ScaleData(tmp, verbose = FALSE)
    tmp <- RunPCA(tmp, npcs = 30, verbose = FALSE)
    tmp <- FindNeighbors(tmp, dims = 1:30, verbose = FALSE)
    tmp <- FindClusters(tmp, resolution = 1, verbose = FALSE)
    setNames(as.character(Idents(tmp)), colnames(tmp))
}

ambient.summary <- list()

obj.list <- lapply(obj.list, function(seu) {
    sample.name <- seu@project.name
    message("Ambient RNA/DNA: ", sample.name, " (", ncol(seu), " cells)")
    raw.path <- FindRawMatrix(sample.name)
    if (is.null(raw.path)) {
        warning("No raw_feature_bc_matrix(.h5) for ", sample.name, " at its samplesheet path or in data/raw/",
                sample.name, "/ -- ambient correction skipped, RNA counts left uncorrected")
        ambient.summary[[sample.name]] <<- data.frame(sample = sample.name, method = "none (no raw matrix)")
        return(seu)
    }
    message("  raw matrix: ", raw.path)
    raw <- ReadRawMatrix(raw.path)
    raw.gex <- raw[["Gene Expression"]]
    raw.peaks <- raw[["Peaks"]]

    cells <- colnames(seu)
    rna.counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
    missing.cells <- setdiff(cells, colnames(raw.gex))
    missing.genes <- setdiff(rownames(rna.counts), rownames(raw.gex))
    if (length(missing.cells) > 0 || length(missing.genes) > 0) {
        stop(raw.path, " doesn't match ", sample.name, "'s filtered matrix (", length(missing.cells),
             " cell barcodes and ", length(missing.genes), " genes missing) -- check the samplesheet path ",
             "points at the same Cell Ranger run as data/raw/", sample.name, "/")
    }
    raw.gex <- raw.gex[rownames(rna.counts), ]

    set.seed(1234)
    clusters <- QuickClusters(rna.counts)
    message("  ", length(unique(clusters)), " quick clusters")

    # RNA: SoupX
    sc <- SoupChannel(tod = raw.gex, toc = rna.counts, calcSoupProfile = TRUE)
    sc <- setClusters(sc, clusters)
    rho.source <- "fixed (ambient_rna_rho)"
    if (is.na(fixed.rho)) {
        rho.source <- "autoEstCont"
        sc <- tryCatch(autoEstCont(sc, doPlot = FALSE, verbose = FALSE), error = function(e) {
            warning("SoupX autoEstCont failed for ", sample.name, " (", conditionMessage(e), "); RNA left ",
                    "uncorrected -- set ambient_rna_rho in config/pipeline.config to use a fixed fraction")
            NULL
        })
    } else {
        sc <- setContaminationFraction(sc, fixed.rho)
    }
    rho <- if (is.null(sc)) NA_real_ else mean(sc$metaData$rho)
    top.soup <- if (is.null(sc)) NA_character_ else
        paste(head(rownames(sc$soupProfile)[order(sc$soupProfile$est, decreasing = TRUE)], 10), collapse = ";")

    seu$ambient.cluster <- unname(clusters[cells])
    seu$nCount_RNA.raw <- seu$nCount_RNA
    seu$nFeature_RNA.raw <- seu$nFeature_RNA
    seu$ambient.rna.frac.removed <- 0
    corrected <- !is.null(sc) && correct.rna
    if (corrected) {
        adj <- adjustCounts(sc, roundToInt = TRUE)[rownames(rna.counts), cells]
        # CreateAssayObject (v3 assay) works under both Seurat v4 and v5
        seu[["RNA.raw"]] <- CreateAssayObject(counts = rna.counts)
        seu <- SetAssayData(seu, assay = "RNA", slot = "counts", new.data = adj)
        seu$nCount_RNA <- unname(Matrix::colSums(adj))
        seu$nFeature_RNA <- unname(Matrix::colSums(adj > 0))
        seu$ambient.rna.frac.removed <- 1 - seu$nCount_RNA / seu$nCount_RNA.raw
    }
    message("  RNA rho = ", round(rho, 3), " (", rho.source, "), ",
            if (corrected) paste0(round(100 * (1 - sum(seu$nCount_RNA) / sum(seu$nCount_RNA.raw)), 1), "% of UMIs removed")
            else "counts not corrected")

    # ATAC: ambient-likeness from the same empty droplets SoupX used
    raw.umi <- Matrix::colSums(raw.gex)
    empty <- colnames(raw.gex)[raw.umi > 0 & raw.umi < 100 & !colnames(raw.gex) %in% cells]
    ambient.peaks <- Matrix::rowSums(raw.peaks[, empty, drop = FALSE])
    cell.peaks <- raw.peaks[, cells]
    cell.peaks@x[] <- 1
    # binarized cells: ||x|| = sqrt(number of accessible peaks)
    seu$ambient.atac.sim <- unname(as.numeric(Matrix::crossprod(cell.peaks, ambient.peaks)) /
                                   (sqrt(Matrix::colSums(cell.peaks)) * sqrt(sum(ambient.peaks^2))))
    seu$ambient.atac.z <- ave(seu$ambient.atac.sim, seu$ambient.cluster, FUN = function(s) {
        spread <- mad(s, na.rm = TRUE)
        if (!is.finite(spread) || spread == 0) return(rep(NA_real_, length(s)))
        (s - median(s, na.rm = TRUE)) / spread
    })

    ambient.summary[[sample.name]] <<- data.frame(
        sample = sample.name,
        method = if (corrected) "SoupX" else if (is.null(sc)) "none (autoEstCont failed)" else "SoupX estimate only",
        rho = round(rho, 4),
        rho.source = rho.source,
        n.cells = length(cells),
        n.empty.droplets = length(empty),
        pct.rna.umi.removed = round(100 * (1 - sum(seu$nCount_RNA) / sum(seu$nCount_RNA.raw)), 2),
        rna.ambient.load = round(1 - sum(raw.umi[cells]) / sum(raw.umi), 4),
        atac.ambient.load = round(1 - sum(raw.peaks[, cells]) / sum(raw.peaks), 4),
        median.ambient.atac.sim = round(median(seu$ambient.atac.sim, na.rm = TRUE), 4),
        top.soup.genes = top.soup)
    seu
})

write.csv(dplyr::bind_rows(ambient.summary),
          file = paste0("output/tables/", argv$project_prefix, "-01-ambient-summary.csv"), row.names = FALSE)

# quick clusters are unlabeled, per-sample groupings (not cell types), so
# one fixed qualitative palette rather than the lineage color scheme
pdf(file = paste0("output/plots/", argv$project_prefix, "-01-ambient-plots.pdf"), height = 5, width = 12)
for (seu in obj.list) {
    if (!"ambient.cluster" %in% colnames(seu@meta.data)) next
    cluster.levels <- as.character(sort(as.numeric(unique(seu$ambient.cluster))))
    cluster.cols <- setNames(grDevices::hcl.colors(length(cluster.levels), "Dark 3"), cluster.levels)
    seu$ambient.cluster <- factor(seu$ambient.cluster, levels = cluster.levels)
    # stack = TRUE: unstacked VlnPlot errors in Seurat 5.1.0 under ggplot2
    # 4.x (what the Minerva R/4.2.0 module loads) -- "cannot get a slot
    # ("slots") from an object of type NULL"; the stacked form works (as in qc.R)
    print(VlnPlot(seu, features = c("ambient.rna.frac.removed", "ambient.atac.z"), group.by = "ambient.cluster",
                  cols = cluster.cols, pt.size = 0, stack = TRUE, fill.by = "ident") +
          NoLegend() + labs(y = "Quick cluster", title = seu@project.name))
}
dev.off()

saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-01-ambient-obj-list.RDS"))
future::plan(ambient.engine.plan)
