# combine atac_fragments files with the GEX folder
CombineDirectories <- function(data.dir, new.dir, sample){
   dir.create(new.dir)
   file.copy(from = paste0(data.dir, "/atac_fragments.tsv.gz"),
             to = paste0("data/raw/", sample, "/fragments.tsv.gz"))
   
   file.copy(from = paste0(data.dir, "/atac_fragments.tsv.gz.tbi"),
             to = paste0("data/raw/", sample, "/fragments.tsv.gz.tbi"))
   
   gex.files <- list.files(path = paste0(data.dir, "/filtered_feature_bc_matrix"))
   lapply(gex.files, 
          function(file.to.copy){
             file.copy(from = paste0(data.dir, "/filtered_feature_bc_matrix/", file.to.copy),
                       to = paste0("data/raw/", sample, "/", file.to.copy))
          })
}

# load in data and create assays for both Gene Expression and ATAC Peaks
CreateMultiomeSeurat <- function(data.dir, annotation = my.annotation, frag.path=NULL, rename.rows=F){
    require(Seurat)
    require(Signac)
    mat.list <- Read10X(data.dir)
    seu <- CreateSeuratObject(counts = mat.list[["Gene Expression"]])

    if(is.null(frag.path)){
        frag.path <- paste0(data.dir, "/fragments.tsv.gz")
    }

    if(rename.rows){
        rownames(mat.list[["Peaks"]]) <- paste0("chr", rownames(mat.list[["Peaks"]]))
    }

    # create ATAC assay and add it to the object
    seu[["ATAC"]] <- CreateChromatinAssay(
        counts = mat.list$Peaks,
        sep = c(":", "-"),
        fragments = frag.path,
        annotation = my.annotation
    )
    seu <- subset(seu, nFeature_RNA > 250)
    # if(seqlevelsStyle(seu[["ATAC"]]@ranges) != "UCSC"){
    #   seqlevelsStyle(seu[["ATAC"]]@ranges) <- "UCSC"
    #   rownames(seu[["ATAC"]]@counts) <- paste0("chr", rownames(seu[["ATAC"]]@counts))
    #   rownames(seu[["ATAC"]]@data) <- paste0("chr", rownames(seu[["ATAC"]]@data))
    # }
    return(seu)
}

CallMyPeaks <- function(seu, fragpath=NULL,grouping.var=NULL,my.macs2.path=NULL,my.annotation=NULL){
    require(EnsDb.Hsapiens.v86)
    require(Seurat)
    require(Signac)

    DefaultAssay(seu) <- "ATAC"
    # if(is.null(fragpath)){
    #   fragments <- Fragments(seu)
    #   fragpath <- fragments[[1]]@path
    # }
    # call peaks using MACS2
    seqlevelsStyle(blacklist_hg38_unified) <- "NCBI"

    peaks <- CallPeaks(seu, 
      group.by = grouping.var,
      outdir = "./data/raw/macs-peaks",
      macs2.path = my.macs2.path)

    #   peaks <- seu[["ATAC"]]@ranges
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(seu),
      features = peaks,
      cells = colnames(seu)
    )

    # create a new assay using the MACS2 peak set and add it to the Seurat object
    seu[["peaks"]] <- 
    CreateChromatinAssay(
      counts = macs2_counts,
      fragments = Fragments(seu),
      annotation = my.annotation
    )
    return(seu)
}

# performing classical normalization and dimensionality reduction
Preprocess.and.Reduce.Dims <- function(seu, harmony=FALSE, harmony.vars = NULL, vars.to.regress=NULL){
    require(Seurat)
    require(Signac)

    # RNA analysis
    DefaultAssay(seu) <- "RNA"
    seu <- SCTransform(seu, verbose = FALSE, vars.to.regress = vars.to.regress) %>% RunPCA(npc = 30) 

    if(harmony){
        require(harmony)
    message("Running harmony integration on RNA: ", harmony.vars)
    seu <- RunHarmony(seu, group.by.vars = harmony.vars, reduction = "pca", assay.use = "SCT", reduction.save = "rna.harmony", project.dim = FALSE)
    }
    reduction.to.use <- ifelse(harmony, "rna.harmony", "pca")
    seu <- seu %>% RunUMAP(dims = 1:30, reduction = reduction.to.use, reduction.name = "umap.rna", reduction.key = 'rnaUMAP_')

    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    if("peaks" %in% names(seu@assays)){
    DefaultAssay(seu) <- "peaks"
    } else {
    DefaultAssay(seu) <- "ATAC"
    }
    seu <- RunTFIDF(seu)
    seu <- FindTopFeatures(seu, min.cutoff = 'q10')
    seu <- RunSVD(seu, n = 100)

    # if(harmony){
    # message("Running harmony integration on ATAC: ", harmony.vars)
    # # seu <- ScaleData(seu, assay = "ATAC")
    # seu <- RunHarmony(seu, group.by.vars = harmony.vars, reduction = "lsi",assay.use = "ATAC", dims.use = 2:50, reduction.save = "atac.harmony", project.dim = FALSE)
    # }
    # reduction.to.use <- ifelse(harmony, "atac.harmony", "lsi")
    reduction.to.use <- "lsi"
    seu <- RunUMAP(seu, reduction = reduction.to.use, dims = 2:100, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    return(seu)
}

# constructing the WNN grapoh that incorporates both ATAC and RNA reduction
ConstructWNNGraph <- function(seu, harmony = FALSE, resolution = 0.8){
    require(Seurat)
    require(Signac)
  if(harmony){
    # reduction_list <- list("rna.harmony", "atac.harmony")
    reduction_list <- list("rna.harmony", "lsi")
  } else {
    reduction_list <-  list("pca", "lsi")
  }
  seu <- FindMultiModalNeighbors(seu, reduction.list = reduction_list, dims.list = list(1:30, 2:100))
  seu <- RunUMAP(seu, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  seu <- FindClusters(seu, graph.name = "wsnn", algorithm = 3, resolution = resolution, verbose = FALSE)
  return(seu)
}

# links peaks to genes to determine if there is a link between accessibility and expression
LinkMyPeaks <- function(seu, genes = NULL, peak.genome,distance.to.use=1000000){
    require(Seurat)
    require(Signac)
  DefaultAssay(seu) <- "peaks"
  seu <- RegionStats(seu, genome = peak.genome)
  # seu <- LinkPeaks(seu, peak.assay = "peaks", expression.assay = "RNA", distance = 5000, score_cutoff = 0.1)
  # seu <- FindVariableFeatures(seu, assay = "RNA")
  # seu <- LinkPeaks(seu, peak.assay = "peaks", expression.assay = "RNA", genes.use = seu@assays[["RNA"]]@var.features, distance = 10000)
  seu <- 
    LinkPeaks(seu, 
      peak.assay = "peaks", 
      expression.assay = "SCT", 
      genes.use = genes,
      distance = distance.to.use)
  return(seu)
}

# EnrichR function
RunEnrichR <- function(genes.to.test, dbs=NULL, plot=T, title =NULL, top.n = 30){
    require(enrichR)
    require(ggplot2)
  if(is.null(dbs)){
    dbs <- c("GO_Biological_Process_2021" ,
             # "GO_Cellular_Component_2021",  "GO_Molecular_Function_2021", "MSigDB_Hallmark_2020",
             "KEGG_2021_Human",
             "Reactome_2022"
             )
  }
  res <- enrichr(genes = genes.to.test, databases = dbs) 
  Sys.sleep(1) # enrichR gets confused if too quick between different gene lists, must be an API thing on the internet
  res <- res[which(unlist(lapply(res, nrow)) > 0)] %>% dplyr::bind_rows(.id = "db")

  res <-  res %>%
        mutate(nGenes = sapply(strsplit(Overlap, split = "/"), getElement, 1)) %>%
               # db = ifelse(grepl("Biological", db), "BP",
               #             ifelse(grepl("Cellular", db), "CC", "MF"))) %>%
        dplyr::filter(nGenes > 1,
                      Adjusted.P.value < 0.05,
               Combined.Score > 0)
  
  if(nrow(res) == 0){
    warning("nrow == 0, skipping")
    return(NULL)
  }
  
  if(plot){
    p <- 
      res %>%
        group_by(db) %>%
        top_n(n = top.n, wt = -Adjusted.P.value) %>%
        arrange(-Adjusted.P.value) %>%
        mutate(sig = ifelse(Adjusted.P.value < 0.05, "*", "ns"),
               Term = factor(Term, levels = Term)) %>%
        ggplot(aes(x = -log(P.value), y = Term, fill = sig)) +
        geom_bar(stat = "identity") +
      geom_text(aes(label = Overlap)) +
        ggtitle(title) +
        theme_classic() +
        facet_grid(db~., scale = "free") +
        scale_fill_manual(values = c("*" = "#990000","ns" = "#999999")) +
        theme(axis.text = element_text(color = "black"),
              strip.text.y = element_text(color = "black", face = "bold"),
              plot.title = element_text(hjust = 1))
    return(p)
  } 
  return(res)
}


# keep standard chromosomes of a BSgenome object
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

# taken from https://stackoverflow.com/questions/25149520/append-filename-with-date 
date.time.append <- function(str, sep = '-', date.format ="%Y_%m_%d_%H_%M_%S") {
  stopifnot(is.character(str))
  return(paste(str, format(Sys.time(), date.format), sep = sep))  
}

# footprint motifs and specific peaks
FootprintMyPeaks <- function(obj, peaks.to.test, motifs.to.test=NULL, peak.genome=NULL, motif.set=NULL){
  require(motifmatchr)
  require(TFBSTools)
  require(JASPAR2020)

  if(is.null(motif.set)){
    motif.set <- getMatrixSet(
      x = JASPAR2020,
      opts = list(species = 9606, all_versions = TRUE)
    )
    # tf.names <- unlist(lapply(motif.set@listData, function(mini.list){return(mini.list@name)}), use.names = F)
    # tf.in.object <- tf.names[tf.names %in% rownames(obj[["RNA"]])]
    # motif.set[tf.names %in% tf.in.object]
  }

  if(is.null(peak.genome)){
    peak.genome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    seqlevelsStyle(peak.genome) <- "UCSC"
  }

  # # subset to retain only peaks in the standard chromosomes and remove non 'peaks' assays
  # # lose information but footprinting functionality is very memory intensive and this should help it
  # main.chroms <- standardChromosomes(peak.genome)
  # keep.peaks <- as.logical(seqnames(granges(obj)) %in% main.chroms)
  # obj <- obj[keep.peaks, ]

  # add motif information
  DefaultAssay(obj) <- "peaks"
  obj <- AddMotifs(obj, genome = peak.genome, pfm = motif.set)

  if(is.null(motifs.to.test)){
    enriched.motifs <- FindMotifs(
      object = obj,
      features = peaks.to.test
    )
    motifs.to.test <-
      enriched.motifs %>%
      dplyr::filter(p.adjust < 0.001) %>%
      top_n(n = 250, wt = -log(p.adjust)) %>%
      pull(motif) %>% 
      unique()
  }

  obj <- Footprint(
    object = obj,
    motif.name = motifs.to.test,
    genome = peak.genome, 
    in.peaks = TRUE, 
    compute.expected = FALSE,
    upstream = 1000,
    downstream = 1000
  )
  return(obj)
}
