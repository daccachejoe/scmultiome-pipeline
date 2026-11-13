# scmultiome-pipeline

Naik Lab single-cell multiome (10x RNA + ATAC) pipeline: bash orchestration, R (Seurat v4 / Signac) for the core analysis, Python for UCDeconvolve, SCENIC+ and the decoder. `README.md` covers the stages, layout and configuration. Read it before changing structure. This file holds the rules for editing the pipeline.

## 1. Where edits happen
- This working copy is the dev repo: branch `dev`, remote `github.com/daccachejoe/scmultiome-pipeline`. Analysis projects are clones of it.
- Pipeline code (`run/`, `routes/`, `scripts/`, `docs/`, `README.md`, `pyproject.toml`, `config/pipeline.config.example`, `install_r_packages.R`) changes **here only**. Then commit, push, and `git pull` in each project.
- If a fix had to be made inside a project clone, copy it back here in the same session.
- The project owns `config/pipeline.config`, `configs/*`, `data/` and `output/`. The tracked `configs/*` here are templates and defaults. Never commit one project's samples, QC thresholds or labels to them.
- Run `git status` here (and in the project, if one is involved) before editing. Uncommitted work in progress is common. Don't overwrite or revert it, and don't mix it into your commit.
- Don't commit, push or merge `dev` → `main` unless asked.

## 2. Architecture to preserve
Call chain: `run/runmultiome <stage>` → sources `config/pipeline.config` with every key exported (`set -a`) → `submit_job` → `bash routes/<route>.sh` on LSF/SLURM → the route loads modules or conda → an R/Python script.

- **Trunk vs. branches.** Stages 00–06 run in strict order. Anything optional or project-shaped goes in `routes/downstream/` + `scripts/downstream/<branch>/` and reads stage 05's output. Only promote something to the trunk if nearly every multiome project needs it (linkpeaks is the precedent). The exception is optional work that must run *before* the trunk because trunk stages consume its output: it gets a lettered pre-trunk stage (`routes/00b_demultiplex.sh` + `scripts/demultiplex/`, feeding stage 01 `doublets` and stage 02 `donors`), not a downstream branch.
- **Per-sample, then merge.** Stage 01 creates, QCs and peak-calls each sample alone. Stage 02 filters and merges. Don't move the merge earlier (the README explains why).
- **The engine is a dispatcher.** `scripts/seurat_signac_pipeline.R` parses args, loads config/genome/annotation, then `source()`s `scripts/stages/<name>.R` in a fixed order, whatever order the comma list uses. Stage files work on shared globals (`obj.list`, `argv`, `annotation`, `peak.genome`, `blacklist.to.use`, `jaspar.taxid`, `samples`, `samplesheet`). A new engine stage needs a file in `scripts/stages/`, a slot in the fixed order, and its name in the `pipeline` arg help text.
- **Standalone scripts** (numbered `scripts/NN_*.R|py`, downstream scripts) must start by calling `source("scripts/lib/config.R"); load_pipeline_config()` when they read config keys. That way they also work in an interactive session.
- **Species.** Genome, EnsDb, blacklist and JASPAR taxon come only from `load_species_genome()` in `scripts/lib/genome.R`. Never hardcode hg38/`EnsDb.Hsapiens`/`9606`/`^MT-` in new code. Add the mouse case to `genome.R` instead. Existing violations (`create.R`'s `^MT-`, `CallPeaksMACS`'s blacklist default, `FootprintMotifs`' GRCh38 fallback) should shrink, not grow.
- **Seurat version mismatch.** `install_r_packages.R` pins Seurat v4, but the Minerva `R/4.2.0` module loads Seurat 5.1.0 / SeuratObject 5.0.2 first, and that's what the jobs actually run. Write code that works on both. `GetAssayData(obj, assay =, slot = "counts")` works in both (it only warns in v5). Avoid v5-only `layer=`/`LayerData()`, and don't add new `@counts`-style slot access.
- Scripts assume the working directory is the project root and use relative paths (`scripts/…`, `output/…`, `configs/…`). Keep it that way. No absolute paths, and no `setwd()` in pipeline code.

## 3. Config, overrides and inputs
- **`config/` vs `configs/`.** `config/pipeline.config` holds environment and HPC settings, filled in once. `configs/` holds per-project run inputs (samplesheet, `qc_df`, cluster labels, resolution, SCENIC+ YAMLs, decoder constraint matrices).
- **New config key:** add it to `config/pipeline.config.example` with a comment and a safe default. `init` appends missing keys to existing project configs. Scripts must still run with a fallback when the key is absent (`${key:-default}` in bash, `Sys.getenv("key", unset=…)` in R).
- **Per-run overrides are env vars that aren't config keys.** Config keys get re-sourced and overwrite any same-named env var, so they can't be overridden this way. Supported overrides:
  - `TIME`/`MEM`/`CORES` for the job
  - `INPUT_RDS` for every route that reads an RDS
  - `DECODER_*` for the decoder branch
  
  New override knobs follow the same pattern: `VAR="${OVERRIDE:-${config_key:-default}}"`, with a commented usage example next to it.
- **New input-reading route:** give it `INPUT_RDS="${INPUT_RDS:-output/RDS-files/<default>}"` and a comment above showing an override example.
- **One project's needs are not a pipeline change.** Handle them with overrides (`MEM=`, `INPUT_RDS=`, `DECODER_*`) or project `configs/`, not by editing route defaults, stage lists or memory. Only change a default when it's wrong for projects in general.
- **Samplesheet contract.** The first two columns are `sampleName,path`. Every column after `path` becomes per-sample metadata. `init` creates the empty `configs/*` headers in `routes/00_setup_dirs.sh`. A new required input file gets its header created there, guarded with `if [ ! -f … ]` so re-running `init` is safe.

## 4. Jobs and scheduling
- Every stage job goes through `submit_job` in `run/runmultiome` (auto-detects SLURM/LSF, invokes `bash <script>` so it doesn't depend on exec bits). Don't put `#BSUB`/`#SBATCH` headers in routes. Don't hardcode `acc_naiklab`/`premium`; use `lsf_project`/`lsf_queue`/`slurm_partition`.
- Routes that submit their own child jobs (06 linkpeaks) must handle all three branches: slurm, lsf, and a none/inline fallback. They must also chain dependencies (`afterok` / `done()`) for the merge/collect step.
- Logs go to `scripts/outs/<name>-%J.out`.
- Module loading in routes uses `"$lsf_r_module"`/`"$slurm_r_module"`. Conda environments: `source "$personal_anaconda_path"; conda activate "$<x>_env_name"`. Python dependencies stay split across the ucd, scenicplus and decoder environments (conflicting pins). New Python deps go in the matching `pyproject.toml` extra and `docs/INSTALLATION.md`. New R deps go in `install_r_packages.R`.
- Don't run R/Python analysis on the login node. Syntax checks are fine (see §8).

## 5. Outputs and naming
- Trunk checkpoints: `output/RDS-files/{project_prefix}-{NN}-{name}-obj[-list].RDS`. Downstream: `{project_prefix}-{name}-…` (no step number). Plots go to `output/plots/`, tables to `output/tables/`, branch outputs to `output/<branch>/`.
- Every stage that changes the object must `saveRDS` a checkpoint.
- RDS files may hold a bare Seurat object or a list. Read them through `as_seurat_list()` / `first_seurat()` (`scripts/lib/seurat_io.R`). Trunk stages save lists.
- Never overwrite an existing RDS or output table from a *different* run configuration. Encode the parameters in the path, as the decoder's `RUN_TAG` directory does, or take a user-supplied tag.
- If you rename an output, update every consumer's `INPUT_RDS` default and the README in the same commit. Grep for the old name first. Downstream routes hardcode upstream filenames (e.g. `-grouped-peaks-05-callpeaks-obj-list.RDS`).
- New files and folders use kebab-case, except numbered step scripts, which follow the existing `NN_snake_case` pattern. To add a step between existing ones, use a letter suffix (`02b_…`, `06a/b/c_…`) rather than renumbering.

## 6. Code style
- **R:** match the surrounding style (`obj.list`, dotted names, `lapply` over objects, `message()` for progress). Shared Seurat/Signac helpers go in `scripts/functions.R` with PascalCase names (`PreprocessAndReduceDims`, `CallPeaksMACS`). Non-Seurat utilities go in `scripts/lib/`. Put CLI parsing in `argparser` (vendored at `.lib/argparser`) for new scripts, not raw `commandArgs` indexing.
- **Python:** argparse with a module docstring holding purpose and USAGE. Shared math goes in an importable sibling module (`preprocessing.py`, `projection.py`) rather than being duplicated. Scripts are run directly, never installed as a package.
- **Comments explain why,** including rejected alternatives and non-obvious traps. Every new script starts with a header block saying what it does, which stage/branch it belongs to, its inputs and outputs, and its usage. Match that density.
- **Plots:**
  - Use the native Seurat/Signac function first (`DimPlot`, `VlnPlot`, `FeatureScatter`, `CoveragePlot`, …) and modify its output with ggplot2 layers.
  - `VlnPlot` only works with `stack = TRUE` (plus `fill.by = "ident"`) under Seurat 5.1.0 + ggplot2 4.x. Unstacked, it fails with "cannot get a slot ("slots") from an object of type NULL" (see the comments in `ambient.R` and `qc.R`).
  - Hand-rolled ggplot2 is only for things Seurat can't do. Leave a comment saying why (see `qc.R`).
  - Never rely on ggplot2's default palette. Set colors explicitly, with one named color vector for cell types (same lineage = same color family) reused across a script's plots.
- **Fail loudly on bad input.** Check required columns, file existence and name alignment (see `projection.py` and `01_pseudobulk_creation.R`'s constraint-matrix check) and `stop()`/`raise` with a message naming the file and the fix. Don't silently drop samples or recycle vectors.
- **No secrets or personal paths in tracked files.** API tokens and user paths go in `config/pipeline.config` (git-ignored) or env vars.

## 7. Decoder branch rules (don't "simplify" these)
- Standardization and per-cell-type centering are fit **inside each CV training fold** (`preprocessing.fit_…`/`apply_…`), never globally before splitting.
- Folds are leave-one-**donor**-out (`--fold_col`), not leave-one-condition-out.
- The constraint matrix is always read from `configs/decoder/*.csv`, never hardcoded. Its row names must equal the pseudobulk `condition` values. The all-zero row is the reference class.
- The stable-peak test is two-sided by default, outputs are split `_up`/`_down`, and `n_folds_positive` is reported next to q-values. Keep all of these.
- Filters (`--celltype_filter`, `--condition_filter`) run before pseudobulking and must be part of `RUN_TAG`.
- If you change a step's input/output format, update `self_test_e2e.py` and rerun `decoder_self_test`.
- Rationale: README "Decoder methodology & caveats"; `../SESSION_EXPORT.md` §1–2 (its file names predate integration).

## 8. Before calling a change done
- Syntax checks (fine on the login node):
  - `bash -n run/runmultiome routes/*.sh routes/downstream/*.sh`
  - `module load R/4.2.0; Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f))' <files>`
  - `python -m py_compile <files>`
- `run/runmultiome --help` still prints, and the new stage appears in `usage()`.
- README stage table, repository tree, and any route/`INPUT_RDS` examples match the code. Update `docs/INSTALLATION.md` for new environments or dependencies.
- Behavior changes: test on a small project (e.g. `inflammatory-skin-disease-data/scmultiome/prelim-long-data/`) with `INPUT_RDS` pointing at an existing checkpoint, submitted as a job. Say plainly if this wasn't done.
- Commits: one logical change each, on `dev`. The subject is an imperative sentence saying what changed and, for fixes, the symptom (see `git log`). Untracked `__pycache__/` and project data never get committed.

## 9. Known issues (fix in place when touching these files; delete entries once fixed)
- The UCDeconvolve API token that used to be hardcoded in `scripts/03_ucd_deconvolve.py` is still in git history. Rotate it (it now comes from `ucd_api_token` in `config/pipeline.config`).
- `scripts/stages/subcluster.R` is project-specific: hardcoded `insitution` (sic) harmony var, `== 6` sample count, `fragments[i]`, fixed QC thresholds. It needs these parameterized before it's general.
- `scripts/stages/cluster.R` reads `cell.lineage`, which only exists after subclustering. No route calls it, and its marker tables overwrite across objects.
- `scripts/stages/merge.R`: `residual.features` uses `intersect`, but the comment says union. Decide which one is intended.
- `routes/06_linkpeaks.sh`:
  - Its output prefix is `{prefix}-improved-clust-filtered-full-`, which gives a double hyphen and doesn't match the checkpoint convention.
  - The child jobs call scripts directly, relying on exec bits, instead of `Rscript`.
- `routes/03_identify_celltypes.sh` loads no R module (it relies on R in the conda environments).
- Plots in `qc.R`, `04_label_celltypes.R` and `03_plot_ucd_results.R` use default ggplot2 colors.
- `scripts/downstream/atac_signal_heatmap/` (untracked) still has `#BSUB` headers, a hardcoded `[1-10]` array, and a personal conda environment (`deeptools-joe`). It needs a route, `runmultiome` stage, config keys and README row before it's integrated.
- In `routes/downstream/decoder.sh`, the step labels (`[decoder 4/6]` etc.) no longer match the file numbers (`02b`, `07`).
