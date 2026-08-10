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

Run in this order. Each stage is a `run/runmultiome <name>` command that
submits a SLURM/LSF job running the corresponding `routes/*.sh` script.
Stage 07 is optional and branches off after stage 05.

| # | `run/runmultiome` stage | What it does | Route | Key script(s) |
|---|---|---|---|---|
| 00 | `init` | Set up project directories, bootstrap `config/pipeline.config`, install R deps | `routes/00_setup_dirs.sh`, `routes/00_install.sh` | -- |
| 01 | `seurat_preprocess` | Per-sample: create Seurat objects, QC metrics, MACS peak calling, QC plots | `routes/01_seurat_preprocess.sh` | `scripts/seurat_signac_pipeline.R` |
| 02a | `run_merged_pipeline` | Filter samples per `configs/qc_df.csv`, merge into one object with a consensus peak set | `routes/02a_merge_pipeline.sh` | `scripts/seurat_signac_pipeline.R` |
| 02b | `filter_and_cluster` | Subcluster within cell lineages (Harmony batch correction) | `routes/02b_filter_and_subcluster.sh` | `scripts/seurat_signac_pipeline.R` |
| 03 | `identify_celltypes` | UCDeconvolve-assisted cell type calling | `routes/03_identify_celltypes.sh` | `scripts/03_convert_seurat_to_h5ad.R`, `scripts/03_ucd_deconvolve.py`, `scripts/03_plot_ucd_results.R` |
| 04 | `label_celltypes` | Apply your cluster -> cell-type labels (`configs/cluster_labels.csv`) | `routes/04_label_celltypes.sh` | `scripts/04_label_celltypes.R` |
| 05 | `call_peaks_grouped` | Re-call ATAC peaks grouped by cell type | `routes/05_call_peaks_grouped.sh` | `scripts/seurat_signac_pipeline.R` |
| 06 | `linkpeaks` | Link ATAC peaks to nearby genes, in parallel (one job per group + a dependent merge job) | `routes/06_linkpeaks.sh` | `scripts/06a_linkpeaks_split.R`, `scripts/06b_linkpeaks_group.R`, `scripts/06c_linkpeaks_merge.R` |
| 07 (optional) | `run_scenicplus` | SCENIC+ regulon inference | `routes/optional/07_run_scenicplus.sh` | `scripts/optional/07a_export_scenicplus_data.R`, `07b_reformat_anndata.py`, `07c_scenicplus_pipeline.py` |

Output RDS checkpoints follow `{project_prefix}-{step}-{name}-obj[-list].RDS`
in `output/RDS-files/`, so the filename tells you which stage produced it.

### Why per-sample objects, merged later?

Samples are created, QC'd, and peak-called **individually** before being
merged (stages 01 -> 02a/02b), rather than pooling raw counts up front. This
is standard Seurat/Signac practice, not an oversight: ambient RNA and
doublet rates vary per 10x run, so QC thresholds need to be set per sample;
MACS peak calling on pooled fragments biases toward high-depth samples; and
Harmony batch correction needs per-sample structure to correct against.
Don't "fix" this into an early merge.

## Repository structure

```
run/runmultiome              # single CLI entry point (dispatcher)
routes/                      # one shell script per stage, numbered 00-06
routes/optional/             # SCENIC+ branch (stage 07)
scripts/                     # R/Python analysis code, numbered to match routes where 1:1
scripts/optional/            # SCENIC+-specific scripts and config template
scripts/lib/                 # shared R helpers (config.R, genome.R)
config/                      # pipeline MACHINERY config (HPC paths, env names) -- edit once
configs/                     # per-project RUN INPUTS (samplesheet, QC thresholds) -- edit per project
data/, output/                # per-project data and results (not tracked in git)
docs/INSTALLATION.md         # environment setup
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

## Known limitations

These are pre-existing gaps in the pipeline's stage-to-stage file handoffs
that this reformatting pass did not change (see git history for the full
rationale) -- worth knowing about if a stage can't find its expected input:

- `routes/05_call_peaks_grouped.sh` expects an input RDS
  (`-improved-clustering-annotated-filtered.RDS`) that no current route
  actually produces under that name.
- `routes/04_label_celltypes.sh` expects `{prefix}-merged-obj-list.RDS`,
  but stage 02a currently produces `{prefix}-02a-merge-obj-list.RDS`.

If you hit either of these, check `output/RDS-files/` for the actual
filename produced by the previous stage and adjust the route's `-R` flag
(or the input file name) accordingly.

## Version

See [VERSION](VERSION). This is the first structurally reformatted release
(numbered stages, centralized config, species-agnostic genome selection,
Seurat-first plotting, parallel linkpeaks) -- see git log for the full list
of changes from the original single-branch pipeline.
