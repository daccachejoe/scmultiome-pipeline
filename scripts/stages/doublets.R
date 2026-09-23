# Stage: doublets
# Scores doublets per sample (= per 10x capture) on both modalities and
# records a consensus call in the metadata. Cells are NOT removed here: the
# filter step of stage 02 drops doublet.call == "doublet" cells, so doublets
# stay visible in stage 01's QC clustering and plots for review.
# Sourced by scripts/seurat_signac_pipeline.R, right after create and before
# any souporcell split, since doublet rates are a property of the capture.
#
# Evidence, per cell:
#   RNA   scDblFinder on RNA counts -> doublet.rna (scDblFinder's own
#         threshold, which scales the expected rate with the capture's cell
#         number, ~1% per 1000 cells).
#   ATAC  scDblFinder in ATAC mode (aggregated peak features) combined with
#         AMULET (>2 fragments at a locus in a diploid cell) by Fisher's
#         method, as the scDblFinder scATAC vignette recommends -- the two
#         catch different doublets (heterotypic vs. homotypic) and neither is
#         best on every dataset. doublet.atac = combined p <
#         doublet_atac_combined_p (config, default 0.05).
#   Genotype  souporcell status == "doublet", for samples listed in
#         configs/demultiplexing_paths.csv (sampleName,demux_path ->
#         clusters.tsv). Only doublets between donors are detectable this way.
#
# doublet.call = "doublet" if souporcell calls it, or if BOTH the RNA and
# ATAC evidence agree. A single modality alone is not enough: requiring
# agreement keeps cells that only look odd in one assay (e.g. large or
# transcriptionally active cells with high RNA complexity), which is where
# single-method callers remove real biology. doublet.evidence records which
# evidence fired, so the kept "rna only"/"atac only" cells can be reviewed.
#
# Souporcell doublets are passed to the RNA scDblFinder run as knownDoublets
# with knownUse = "discard": they're kept out of training rather than used
# as positive examples, because genotype doublets are often homotypic and
# would teach the classifier that singlet-looking profiles are doublets.
#
# Expected doublet rate. scDblFinder thresholds each capture to its expected
# rate (dbr, default 1% per 1000 cells). With souporcell, most doublets are
# already caught by genotype, so leaving dbr at the full rate forces RNA to
# flag another ~10% of an overloaded capture and removes real singlets (on
# healthy-human-multiome cntrl.1, 4 donors, 7.3% of genotype singlets were
# removed vs. ~2.5% expected). For souporcell samples dbr is therefore set to
# the doublets genotype can't see -- same-donor pairs -- i.e. the full rate x
# sum(donor fraction^2), renormalized to the non-genotype-doublet cells.

suppressPackageStartupMessages({
    library(scDblFinder)
    library(SingleCellExperiment)
})

atac.p.threshold <- as.numeric(Sys.getenv("doublet_atac_combined_p", unset = "0.05"))
if (is.na(atac.p.threshold)) {
    stop("doublet_atac_combined_p in config/pipeline.config must be numeric, got: ",
         Sys.getenv("doublet_atac_combined_p"))
}

# souporcell calls, if any samples were genotype-demultiplexed
demux.file <- "configs/demultiplexing_paths.csv"
demux.paths <- if (file.exists(demux.file)) read.csv(demux.file, stringsAsFactors = FALSE) else NULL

# AMULET: exclude mito/sex chromosomes (as recommended) plus the genome's
# blacklist; both UCSC and NCBI names are listed since fragment files vary
amulet.exclude <- suppressWarnings(c(
    GRanges(c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"), IRanges(1L, width = 10^8)),
    granges(blacklist.to.use)))

# Fisher's method on (1 - scDblFinder ATAC score) and the AMULET p-value,
# flooring both at 0.001 so a single extreme value can't dominate (as in the
# vignette); an NA score or AMULET p (not scored, too few fragments)
# contributes no evidence
CombineATACEvidence <- function(atac.score, amulet.p) {
    p <- cbind(ifelse(is.na(atac.score), 1, 1 - atac.score), ifelse(is.na(amulet.p), 1, amulet.p))
    p[p < 0.001] <- 0.001
    pchisq(-2 * rowSums(log(p)), df = 2 * ncol(p), lower.tail = FALSE)
}

# AMULET reads the fragment file one chromosome at a time; spread that over
# the job's cores
bp.param <- BiocParallel::MulticoreParam(max(1, as.numeric(future::availableCores())))

StepTimer <- function() {
    t0 <- Sys.time()
    function(step) {
        message("  ", step, " (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
        t0 <<- Sys.time()
    }
}

obj.list <- lapply(obj.list, function(seu) {
    sample.name <- seu@project.name
    step.done <- StepTimer()
    message("Doublet detection: ", sample.name, " (", ncol(seu), " cells)")
    set.seed(1234)

    # genotype doublets from souporcell
    souporcell.doublet <- rep(FALSE, ncol(seu))
    rna.dbr <- NULL  # scDblFinder's default expected rate unless souporcell ran
    has.souporcell <- !is.null(demux.paths) && sample.name %in% demux.paths$sampleName
    if (has.souporcell) {
        soc <- read.delim(demux.paths$demux_path[demux.paths$sampleName == sample.name][1])
        matched <- mean(colnames(seu) %in% soc$barcode)
        message("  souporcell: ", round(100 * matched, 1), "% of barcodes matched, ",
                sum(soc$status == "doublet"), " genotype doublets called")
        if (matched < 0.5) {
            warning("Fewer than half of ", sample.name, "'s barcodes are in its souporcell clusters.tsv -- ",
                    "check configs/demultiplexing_paths.csv points at the right run")
        }
        souporcell.doublet <- colnames(seu) %in% soc$barcode[soc$status == "doublet"]

        donor.frac <- prop.table(table(soc$assignment[soc$status == "singlet" & soc$barcode %in% colnames(seu)]))
        full.rate <- 0.01 * ncol(seu) / 1000
        same.donor <- sum(donor.frac^2)
        rna.dbr <- full.rate * same.donor / (1 - full.rate * (1 - same.donor))
        message("  ", length(donor.frac), " donors; expected doublet rate ", round(100 * full.rate, 1),
                "% overall, ", round(100 * rna.dbr, 1), "% left after genotype doublets (used as RNA dbr)")
    }

    # scDblFinder takes NULL, not an all-FALSE vector, when there are no known doublets
    known.doublets <- if (any(souporcell.doublet)) souporcell.doublet else NULL

    # RNA
    rna.sce <- SingleCellExperiment(list(counts = GetAssayData(seu, assay = "RNA", slot = "counts")))
    rna.sce <- scDblFinder(rna.sce, knownDoublets = known.doublets, knownUse = "discard", dbr = rna.dbr,
                           verbose = FALSE)
    step.done("RNA scDblFinder done")

    # ATAC: scDblFinder in ATAC mode on the Cell Ranger peak counts, on all
    # cells and without knownDoublets -- scDblFinder 1.12 fails when both
    # knownDoublets and aggregateFeatures = TRUE are given ("index out of
    # bounds: feat1 ..."). Only the score is used (not scDblFinder's ATAC
    # threshold), so dbr doesn't matter here.
    atac.sce <- SingleCellExperiment(list(counts = GetAssayData(seu, assay = "ATAC", slot = "counts")))
    atac.sce <- scDblFinder(atac.sce, aggregateFeatures = TRUE, nfeatures = 25, processing = "normFeatures",
                            verbose = FALSE)
    step.done("ATAC-mode scDblFinder done")

    # ATAC: AMULET on the fragment file
    frag.path <- Fragments(seu[["ATAC"]])[[1]]@path
    amulet.res <- tryCatch(
        amulet(frag.path, barcodes = colnames(seu), regionsToExclude = amulet.exclude,
               BPPARAM = bp.param, verbose = FALSE),
        error = function(e) {
            warning("AMULET failed for ", sample.name, " (", conditionMessage(e),
                    "); using the scDblFinder ATAC score alone for ATAC evidence")
            NULL
        })
    amulet.p <- if (is.null(amulet.res)) rep(NA_real_, ncol(seu)) else amulet.res[colnames(seu), "p.value"]
    step.done("AMULET done")

    # unname(): Seurat's $<- rejects named vectors whose names don't match the cells

    seu$doublet.rna.score <- unname(rna.sce$scDblFinder.score)
    seu$doublet.rna <- unname(rna.sce$scDblFinder.class == "doublet")
    seu$doublet.atac.score <- unname(atac.sce$scDblFinder.score)
    seu$doublet.amulet.p <- unname(amulet.p)
    seu$doublet.atac.combined.p <- CombineATACEvidence(seu$doublet.atac.score, seu$doublet.amulet.p)
    seu$doublet.atac <- seu$doublet.atac.combined.p < atac.p.threshold
    seu$doublet.souporcell <- souporcell.doublet

    seu$doublet.evidence <- dplyr::case_when(
        seu$doublet.souporcell ~ "souporcell",
        seu$doublet.rna & seu$doublet.atac ~ "rna+atac",
        seu$doublet.rna ~ "rna only",
        seu$doublet.atac ~ "atac only",
        TRUE ~ "none")
    seu$doublet.call <- ifelse(seu$doublet.evidence %in% c("souporcell", "rna+atac"), "doublet", "singlet")

    counts <- table(factor(seu$doublet.evidence, levels = c("souporcell", "rna+atac", "rna only", "atac only", "none")))
    message("  ", paste(names(counts), counts, sep = ": ", collapse = ", "),
            " -> ", sum(seu$doublet.call == "doublet"), " doublets (",
            round(100 * mean(seu$doublet.call == "doublet"), 1), "%)")
    seu
})

# per-cell calls and a per-sample summary, for review before stage 02
doublet.cols <- c("doublet.rna.score", "doublet.rna", "doublet.atac.score", "doublet.amulet.p",
                  "doublet.atac.combined.p", "doublet.atac", "doublet.souporcell", "doublet.evidence", "doublet.call")
doublet.md <- dplyr::bind_rows(lapply(obj.list, function(seu) {
    data.frame(sample = seu@project.name, barcode = colnames(seu), seu@meta.data[, doublet.cols])
}))
write.csv(doublet.md, file = paste0("output/tables/", argv$project_prefix, "-01-doublet-calls.csv"), row.names = FALSE)
doublet.summary <- as.data.frame.matrix(table(doublet.md$sample, doublet.md$doublet.evidence))
doublet.summary <- data.frame(sample = rownames(doublet.summary), doublet.summary, check.names = FALSE,
                              pct.removed = round(100 * tapply(doublet.md$doublet.call == "doublet", doublet.md$sample, mean), 2))
write.csv(doublet.summary, file = paste0("output/tables/", argv$project_prefix, "-01-doublet-summary.csv"), row.names = FALSE)

# RNA vs. ATAC evidence per sample; removed cells are the souporcell and
# rna+atac groups, and "rna only"/"atac only" show what is being kept
evidence.cols <- c("souporcell" = "#6A3D9A", "rna+atac" = "#E31A1C", "rna only" = "#FB9A99",
                   "atac only" = "#A6CEE3", "none" = "#BDBDBD")
pdf(file = paste0("output/plots/", argv$project_prefix, "-01-doublet-plots.pdf"), height = 6, width = 8)
for (seu in obj.list) {
    seu$doublet.atac.evidence <- -log10(seu$doublet.atac.combined.p)
    Idents(seu) <- factor(seu$doublet.evidence, levels = names(evidence.cols))
    print(FeatureScatter(seu, feature1 = "doublet.rna.score", feature2 = "doublet.atac.evidence",
                         cols = evidence.cols, pt.size = 0.3, shuffle = TRUE) +
          geom_hline(yintercept = -log10(atac.p.threshold), linetype = "dashed") +
          labs(title = seu@project.name, x = "scDblFinder RNA score", color = "Doublet evidence",
               y = "ATAC evidence, -log10 combined p (scDblFinder ATAC + AMULET)"))
}
dev.off()

saveRDS(obj.list, file = paste0("output/RDS-files/", argv$project_prefix, "-01-doublets-obj-list.RDS"))
