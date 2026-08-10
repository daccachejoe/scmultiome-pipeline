# Species-specific genome/annotation resolution, shared by the main engine
# (scripts/seurat_signac_pipeline.R) and the parallel linkpeaks helpers.
# Reads `species` from config/pipeline.config (via load_pipeline_config()),
# one of "human" (default, hg38) or "mouse" (mm10).
load_species_genome <- function() {
  species <- Sys.getenv("species", unset = "human")

  if (species == "mouse") {
    library(EnsDb.Mmusculus.v79, quietly = TRUE)
    library(BSgenome.Mmusculus.UCSC.mm10, quietly = TRUE)
    list(
      species = species,
      ensdb = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
      genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
      blacklist = Signac::blacklist_mm10,
      jaspar_taxid = 10090
    )
  } else {
    library(EnsDb.Hsapiens.v86, quietly = TRUE)
    library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
    list(
      species = species,
      ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
      genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
      blacklist = Signac::blacklist_hg38_unified,
      jaspar_taxid = 9606
    )
  }
}
