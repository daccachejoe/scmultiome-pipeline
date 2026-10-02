# Installation

This pipeline uses **five separate conda environments** (six with optional demultiplexing) plus your HPC's
module system for R. They're kept separate because SCENIC+, UCDeconvolve,
sceasy, and PyTorch pin conflicting dependency versions.

| Environment (config key)             | Purpose                                   |
| ------------------------------------- | ------------------------------------------ |
| `conda_env_name`                      | Main Seurat/Signac R analysis (trunk stages 00-06, downstream subcluster) |
| `sceasy_env_name`                      | Converts Seurat objects to AnnData (`.h5ad`) |
| `UCD_env_name`                         | UCDeconvolve cell-type calling (stage 03)  |
| `scenicplus_env_name`                  | SCENIC+ regulon inference (downstream branch, extra environments) |
| `decoder_env_name`                     | Constrained decoder (downstream branch): PyTorch, scipy, statsmodels |
| `demux_env_name`                       | Donor demultiplexing (optional stage 00b): pysam, numpy, scipy |

Before starting, copy the config template and fill in your environment names
and paths (see `config/pipeline.config.example` for every key):

```
cp config/pipeline.config.example config/pipeline.config
# edit config/pipeline.config
```

## 1. Main analysis environment (`conda_env_name`)

The main pipeline (trunk stages 00-06, plus the downstream subcluster branch) runs under your HPC's R module,
not a conda R install: SLURM jobs load `slurm_r_module` (default `r/4.1.2`)
and LSF jobs load `lsf_r_module` (default `R/4.2.0`) from `config/pipeline.config`
-- keep these consistent with whatever R module is actually available on
your cluster; **all routes now load the same module for a given scheduler**,
resolving a previous inconsistency where some routes loaded a different LSF
R version than others.

```
conda create -n sc-multiome-pipeline-env
conda activate sc-multiome-pipeline-env
module load <your R module>   # e.g. r/4.1.2 or R/4.2.0
Rscript install_r_packages.R
```

This pins Seurat v4 (`Seurat 4.4.0`, `SeuratObject 4.1.4`) deliberately --
Seurat v5 changes APIs this pipeline relies on. `install_r_packages.R`
documents every package; regenerate a real, hash-verified lockfile for your
environment once it's built:

```r
install.packages("renv")
renv::init()      # first time only
renv::snapshot()  # writes renv.lock
```

The vendored `argparser` package (`.lib/argparser/`) is installed
automatically by `run/runmultiome init` via `routes/00_install.sh`.

## 2. sceasy environment (`sceasy_env_name`)

This is a **separate, self-contained conda R install** (not the HPC module
R above) because `devtools`/`sceasy` need a specific R/library combination:

```
conda create -n sceasy-env r-base=4.2.2
conda activate sceasy-env
conda install -c conda-forge libxml2
conda install -c r r-xml=3.98_1.5
```

In R, temporarily add your HPC's system R library to `.libPaths()` (needed
because `devtools` is otherwise hard to install in a bare conda R):

```r
.libPaths(c(.libPaths(), "/path/to/system/R/library"))  # order matters
devtools::install_github("cellgeni/sceasy")
```

Restart R (so `.libPaths()` reverts to just the conda env), then install
Seurat/Signac in this environment too, since sceasy needs to read the
Seurat objects it's converting:

```r
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "EnsDb.Hsapiens.v86"))
remotes::install_version("SeuratObject", "4.1.4",
    repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0",
    repos = c("https://satijalab.r-universe.dev", getOption("repos")))
setRepositories(ind = 1:3)
install.packages("Signac")
```

## 3. UCDeconvolve environment (`UCD_env_name`)

```
conda create -n ucd-env python=3.10
conda activate ucd-env
uv pip install -e ".[ucd]"
uv lock   # produces a real uv.lock for this environment
```

## 4. SCENIC+ environment (`scenicplus_env_name`)

```
conda create -n scenicplus-env python=3.11
conda activate scenicplus-env
module load cmake

# rust is required by some SCENIC+ dependencies
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# restart your shell session before continuing

git clone https://github.com/aertslab/scenicplus
cd scenicplus && pip install -e . && cd ..

uv pip install -e ".[scenicplus]"
uv lock   # produces a real uv.lock for this environment
```

Then install `create_cisTarget_databases` per the
[SCENIC+ docs](https://scenicplus.readthedocs.io/en/latest/install.html)
(everything is usually already installed except `python-flatbuffers`:
`conda install python-flatbuffers`), plus Cluster-Buster and the motif
collection referenced in `routes/downstream/run_scenicplus.sh`.

Set `create_cistarget_databases_path` and `create_cistarget_databases_dir`
in `config/pipeline.config` to point at your install.

## 5. Decoder environment (`decoder_env_name`)

```
conda create -n decoder-env python=3.11
conda activate decoder-env
uv pip install -e ".[decoder]"
uv lock   # produces a real uv.lock for this environment
```

`scripts/downstream/decoder/06_run_motif_enrichment.py` additionally wraps
HOMER's `findMotifsGenome.pl` (see `decoder_homer_script` in
`config/pipeline.config`) -- installed separately, module-loaded or on
`$PATH`, not part of this conda env.

## 6. Demultiplexing environment (`demux_env_name`, optional stage 00b)

Only for projects that pool several donors per 10x library.

```
conda create -n demux-env -c conda-forge -c bioconda python=3.12.2 pysam=0.22.1 numpy=2.1.0 scipy=1.14.1
```

(or `pip install -e ".[demux]"` into a Python 3.12 env). The versions are
pinned to the ones the ported scripts were checked against byte for byte.

The stage also needs the **souporcell Singularity image** (`souporcell_sif`),
which provides souporcell, samtools/bcftools 1.9 and vartrix 1.1.22, and a
`singularity` module (`singularity_module`). The container only works when
its inputs are copied inside its bind directory, so each library needs work
space for a copy of its GEX BAM (up to ~20 GB) and the reference FASTA: set
`demux_tmp_dir` to scratch.
