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
| `subcluster` | Subcluster within cell lineages (Harmony batch correction) | `routes/downstream/subcluster.sh` | `scripts/seurat_signac_pipeline.R` |
| `run_scenicplus` (extra environments required) | SCENIC+ regulon inference | `routes/downstream/run_scenicplus.sh` | `scripts/downstream/export_scenicplus_data.R`, `reformat_anndata.py`, `scenicplus_pipeline.py` |
| `decoder` (extra environment required) | Constrained multinomial logistic regression decoder: pseudobulk -> distal-peak LSI / RNA PCA -> CV-fit decoder -> peak-space projection -> stable-peak identification -> motif enrichment | `routes/downstream/decoder.sh` | `scripts/downstream/decoder/*` |
| `decoder_self_test` | Synthetic end-to-end check of the decoder branch (no real/demultiplexed data needed) | `routes/downstream/decoder_self_test.sh` | `scripts/downstream/decoder/self_test_e2e.py` |

These are exploratory or heavyweight extras (a closer look at a specific
lineage, a separate regulon-inference toolchain, or a separate
condition-decoding model), unlike linkpeaks -- so they stay outside the
trunk and don't run unless you ask for them.

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

### Decoder methodology & caveats

The `decoder` branch trains condition-specific classifiers whose weight
vectors are constrained to be linear combinations of a small set of
biologically interpretable "component" directions (e.g. subclinical
effect, disease identity, shared lesional effect, disease-specific
interaction), via an externally-supplied constraint matrix
(`configs/decoder/constraint_matrix.csv`) -- never hardcoded, so a new
condition set or a longitudinal extension is a matrix swap, not a code
change. A few design decisions worth knowing before trusting its output:

- **Preprocessing is fit per CV fold, not once globally.** Standardizing
  and cell-type-centering features before splitting into CV folds would
  leak held-out-fold statistics into training, biasing the L2
  regularization-strength grid search optimistic.
  `scripts/downstream/decoder/03_train_cv.py` fits
  `preprocessing.fit_standardize_then_center()` fresh on each training
  split (grid search and per-fold refit) and only applies those
  statistics to that split's held-out fold. The one exception is the
  final all-data model (Phase C), which has no held-out set and so no
  leakage concern.
- **CV folds are leave-one-donor-out by default, not leave-one-condition-
  out**, even though this dataset's conditions are sometimes described as
  "experiments." This constraint matrix has components that are nonzero
  for exactly one condition (e.g. the PSO x Lesional interaction term) --
  leaving that condition out of a fold would delete 100% of that
  component's training signal, not run a real held-out test. Donor-level
  folds keep every condition (and therefore every component) represented
  in every training split. Configurable via `decoder_fold_col` /
  `--fold_col` if the donor pool grows enough to revisit this.
- **The stable-peak t-test's p-values are optimistic.** Leave-one-donor-
  out folds share most of their training data, so per-fold peak weights
  are correlated, not independent draws -- inherent to this kind of
  CV-based stability selection, not a bug. `05_identify_stable_peaks.py`
  reports `n_folds_positive`/`frac_folds_positive` alongside the BH-
  corrected q-values as a non-parametric cross-check; don't rely on the
  q-value alone near the `--fdr_threshold` boundary.
- **`beta_sub` (the "vs. healthy control" effect) is provisional until HC
  donors are demultiplexed.** The model still fits with zero reference
  samples -- the reference intercept is fixed at 0 by construction -- but
  the effect is then estimated purely by extrapolation from the other
  conditions, never anchored against an observed healthy baseline.
  `decoder_config.json` records a `reference_caveat` field whenever this
  applies.

Run `run/runmultiome decoder_self_test` any time to check the whole
pseudobulk -> feature engineering -> CV -> peak projection -> stability
chain end-to-end against synthetic data (useful before real
demultiplexed donors exist, or after changing any step's I/O format).

#### Overriding metadata schema per-run

`decoder_celltype_col`/`decoder_donor_col`/`decoder_condition_col`/
`decoder_location_col`/`decoder_rna_assay`/`decoder_atac_assay` in
`config/pipeline.config` set this project's metadata column names once.
Unlike `TIME`/`MEM`/`CORES`/`INPUT_RDS`, those are config keys, so a
same-named env var set before `run/runmultiome decoder` gets silently
overwritten by the config file (`run/runmultiome` sources it
unconditionally). To point the decoder at a *different* object's schema
for a single run without editing the shared config, use the
`DECODER_*`-prefixed variants instead -- these aren't config keys, so
they survive:

```
DECODER_CELLTYPE_COL=cell_type DECODER_DONOR_COL=patient_id \
    run/runmultiome decoder
```

This is what actually makes the decoder branch portable to a Seurat
object from a different project/lab with different metadata column
names, rather than requiring you to edit `config/pipeline.config` (which
every other stage also reads) every time you point it somewhere new.

Output lands under `output/decoder/<run_tag>/`, where `<run_tag>` encodes
the actual combination of column/assay/fold_col/constraint-matrix
parameters used (e.g. `ct-ct.spec_donor-donor_id_..._cmat-constraint_matrix`)
-- two runs with different overrides write to separate directories instead
of one overwriting the other's output, so different schemas/constraint
matrices can be run side by side. Set `DECODER_RUN_TAG=my-label` for a
shorter, hand-chosen directory name instead of the auto-generated one.

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
│       ├── run_scenicplus.sh      # extra environments required
│       ├── decoder.sh             # extra environment required
│       └── decoder_self_test.sh
├── scripts/                       # R/Python analysis code
│   ├── seurat_signac_pipeline.R   # engine: arg parsing + config/species setup + dispatch
│   ├── functions.R                # shared Seurat/Signac wrapper functions
│   ├── 03_convert_seurat_to_h5ad.R, 03_ucd_deconvolve.py, 03_plot_ucd_results.R
│   ├── 04_label_celltypes.R
│   ├── 06a_linkpeaks_split.R, 06b_linkpeaks_group.R, 06c_linkpeaks_merge.R
│   ├── stages/                    # one file per pipeline stage, sourced by the engine
│   ├── downstream/                # subcluster helper, SCENIC+, and decoder scripts
│   │   └── decoder/                # pseudobulk -> features -> CV decoder -> peak projection
│   │                               #   -> stable peaks -> motif enrichment (see README's
│   │                               #   "Decoder methodology & caveats")
│   └── lib/                       # shared R helpers (config.R, genome.R, seurat_io.R)
├── config/                        # pipeline MACHINERY config -- edit once
│   ├── pipeline.config.example    # tracked template
│   └── pipeline.config            # your values (git-ignored)
├── configs/                       # per-project RUN INPUTS -- edit per project
│   ├── samplesheet.csv
│   ├── qc_df.csv
│   ├── cluster_labels.csv
│   ├── resolution_to_use.txt
│   ├── scenicplus-*-config.yml
│   └── decoder/constraint_matrix.csv
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
