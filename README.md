# scMultiome Analysis Pipeline

#### Naik Lab

A single-cell multiome (10x Genomics RNA + ATAC) analysis pipeline: bash
orchestration, R (Seurat/Signac) for the core analysis, and Python for
UCDeconvolve cell-type calling and SCENIC+ regulon inference. See
[docs/INSTALLATION.md](docs/INSTALLATION.md) for environment setup.

## Quick start

```
git clone <this repo> my-project
cd my-project
run/runmultiome init                        # creates config/pipeline.config and configs/*
# edit config/pipeline.config with your HPC/environment values
# edit configs/samplesheet.csv with your samples
run/runmultiome seurat_preprocess            # stage 01, then review output/ before continuing
```

Run `run/runmultiome --help` (or with no argument) at any time for the
full list of stages.

## Pipeline stages

The pipeline has a **trunk** (run in strict order) and **downstream branches**
(run independently, in any order, once the trunk completes). Each stage is
a `run/runmultiome <name>` command that submits a SLURM/LSF job running the
corresponding `routes/*.sh` script.

### Trunk (00 -> 06, strict order)

| # | `run/runmultiome` stage | What it does | Route | Key script(s) |
|---|---|---|---|---|
| 00 | `init` | Set up project directories, bootstrap `config/pipeline.config`, install R deps | `routes/00_setup_dirs.sh`, `routes/00_install.sh` | -- |
| 01 | `seurat_preprocess` | Per-sample: create Seurat objects, QC metrics, MACS peak calling, QC plots | `routes/01_seurat_preprocess.sh` | `scripts/seurat_signac_pipeline.R` |
| 02 | `run_merged_pipeline` | Filter samples per `configs/qc_df.csv`, merge into one object with a consensus peak set + clustering | `routes/02_merge_pipeline.sh` | `scripts/seurat_signac_pipeline.R` |
| 03 | `identify_celltypes` | UCDeconvolve-assisted cell type calling (human then fills in `configs/cluster_labels.csv`) | `routes/03_identify_celltypes.sh` | `scripts/03_convert_seurat_to_h5ad.R`, `scripts/03_ucd_deconvolve.py`, `scripts/03_plot_ucd_results.R` |
| 04 | `label_celltypes` | Apply your cluster -> cell-type labels (`configs/cluster_labels.csv`) | `routes/04_label_celltypes.sh` | `scripts/04_label_celltypes.R` |
| 05 | `call_peaks_grouped` | Re-call ATAC peaks grouped by cell type -- **branch point for the optional stages below** | `routes/05_call_peaks_grouped.sh` | `scripts/seurat_signac_pipeline.R` |
| 06 | `linkpeaks` | Link ATAC peaks to nearby genes, in parallel (one job per group + a dependent merge job) | `routes/06_linkpeaks.sh` | `scripts/06a_linkpeaks_split.R`, `06b_linkpeaks_group.R`, `06c_linkpeaks_merge.R` |

Peak-gene linking (06) is part of the trunk because it's a near-universal
step in multiome analysis, not something most users would skip.

### Downstream branches (optional, independent of each other, run any/all after stage 05)

| `run/runmultiome` stage | What it does | Route | Key script(s) |
|---|---|---|---|
| `filter_and_cluster` | Subcluster within cell lineages (Harmony batch correction) | `routes/downstream/subcluster.sh` | `scripts/seurat_signac_pipeline.R` |
| `run_scenicplus` (extra environments required) | SCENIC+ regulon inference | `routes/downstream/run_scenicplus.sh` | `scripts/downstream/export_scenicplus_data.R`, `reformat_anndata.py`, `scenicplus_pipeline.py` |

These are exploratory or heavyweight extras (a closer look at a specific
lineage, or a separate regulon-inference toolchain), unlike linkpeaks --
so they stay outside the trunk and don't run unless you ask for them.

Output RDS checkpoints follow `{project_prefix}-{step}-{name}-obj[-list].RDS`
in `output/RDS-files/` for trunk stages, including linkpeaks (06)
(self-documenting which stage produced them); downstream branch outputs
drop the step number since they're not part of a fixed sequence.

### Why per-sample objects, merged later?

Samples are created, QC'd, and peak-called **individually** before being
merged (stages 01 -> 02), rather than pooling raw counts up front. This
is standard Seurat/Signac practice, not an oversight: ambient RNA and
doublet rates vary per 10x run, so QC thresholds need to be set per sample;
MACS peak calling on pooled fragments biases toward high-depth samples; and
Harmony batch correction needs per-sample structure to correct against.
Don't "fix" this into an early merge.

## Repository structure

```
.
├── run/
│   └── runmultiome                # single CLI entry point (dispatcher)
├── routes/                        # trunk stages, numbered 00-06 (strict order)
│   ├── 00_setup_dirs.sh
│   ├── 00_install.sh
│   ├── 01_seurat_preprocess.sh
│   ├── 02_merge_pipeline.sh
│   ├── 03_identify_celltypes.sh
│   ├── 04_label_celltypes.sh
│   ├── 05_call_peaks_grouped.sh   # branch point for the optional stages below
│   ├── 06_linkpeaks.sh
│   └── downstream/                # optional branches, run any/all after 05
│       ├── subcluster.sh
│       └── run_scenicplus.sh      # extra environments required
├── scripts/                       # R/Python analysis code
│   ├── seurat_signac_pipeline.R   # engine: arg parsing + config/species setup + dispatch
│   ├── functions.R                # shared Seurat/Signac wrapper functions
│   ├── 03_convert_seurat_to_h5ad.R, 03_ucd_deconvolve.py, 03_plot_ucd_results.R
│   ├── 04_label_celltypes.R
│   ├── 06a_linkpeaks_split.R, 06b_linkpeaks_group.R, 06c_linkpeaks_merge.R
│   ├── stages/                    # one file per pipeline stage, sourced by the engine
│   ├── downstream/                # subcluster helper + SCENIC+ scripts
│   └── lib/                       # shared R helpers (config.R, genome.R, seurat_io.R)
├── config/                        # pipeline MACHINERY config -- edit once
│   ├── pipeline.config.example    # tracked template
│   └── pipeline.config            # your values (git-ignored)
├── configs/                       # per-project RUN INPUTS -- edit per project
│   ├── samplesheet.csv
│   ├── qc_df.csv
│   ├── cluster_labels.csv
│   ├── resolution_to_use.txt
│   └── scenicplus-*-config.yml
├── data/                          # per-project raw data (not tracked in git)
├── output/                        # per-project results (not tracked in git)
└── docs/
    └── INSTALLATION.md            # environment setup
```

`config/` (no trailing 's') vs `configs/` (trailing 's') is a deliberately
confusing pair of names carried over from the original pipeline -- `config/`
holds `pipeline.config`, environment/HPC settings you fill in **once**;
`configs/` holds samplesheet/QC/label files you edit **per project**.

## Configuration

Copy `config/pipeline.config.example` to `config/pipeline.config`
(git-ignored) and fill in your values -- HPC paths, conda env names, and
`species`/`genome` (`human`/`hg38` or `mouse`/`mm10`). `run/runmultiome init`
does this for you automatically if the file doesn't exist yet.

### samplesheet.csv

`configs/samplesheet.csv` is created empty by `init`. The first two columns
**must** be `sampleName` and `path` (case-sensitive); any columns after
`path` are treated as per-sample metadata and attached to each object.

| sampleName | path | cond |
| ---------- | ---- | ---- |
| ctrl.1 | /path/to/cellranger/count-CTRL/outs | ctrl |
| il17.1 | /path/to/cellranger/count-IL-17/outs | il17a |

## Scheduler support

`run/runmultiome` auto-detects SLURM vs LSF (via `$SLURM_JOB_ID` /
`$LSF_ENVDIR`) and submits jobs accordingly; partition/project/queue come
from `config/pipeline.config`.

Each stage submits its job with its own default time/mem/cores. Override
any of them for a single run with `TIME`/`MEM`/`CORES` env vars:

```
CORES=64 MEM=64000 run/runmultiome seurat_preprocess
```

Most trunk/downstream stages also accept an `INPUT_RDS` override to point
at a different input object -- see the comment near the top of each
`routes/*.sh` file for its exact default.

## Version

See [VERSION](VERSION). This is the first structurally reformatted release
(numbered stages, centralized config, species-agnostic genome selection,
Seurat-first plotting, parallel linkpeaks) -- see git log for the full list
of changes from the original single-branch pipeline.
