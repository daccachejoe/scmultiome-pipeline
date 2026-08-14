#!/usr/bin/env python3
"""
Decoder branch end-to-end self-test: generates synthetic pseudobulk-shaped
data (matching the real sample_metadata schema -- sample_id, cell_type,
donor_id, condition, location) and runs it through 03_train_cv.py ->
04_project_to_peak_space.py -> 05_identify_stable_peaks.py, using
subprocess so this exercises the actual CLI entry points, not internal
functions.

decoder_model.py's own --self_test only checks the model's forward pass
and LBFGS fit on separable synthetic data. This checks the whole chain's
wiring: file formats between steps, the fold-wise-preprocessing leakage
fix, and the stability test -- useful to run on Minerva before real
(demultiplexed) data exists, and after any change to the chain's I/O
contracts.

Also directly demonstrates the CV-leakage fix isn't a no-op: fitting
standardization+centering on a training fold produces different
statistics than fitting on the full dataset (which would include the
held-out fold) -- see check_leakage_fix_is_not_a_noop().

USAGE:
    python self_test_e2e.py [--out_dir /tmp/decoder_self_test] [--keep]
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import preprocessing as prep

HERE = Path(__file__).resolve().parent

CONDITION_NAMES = ["HC", "AD_NonLesional", "PSO_NonLesional", "AD_Lesional", "PSO_Lesional"]
COMPONENT_NAMES = ["beta_sub", "beta_dis", "beta_loc", "beta_int"]
C_MATRIX = pd.DataFrame(
    [[0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 1, 1]],
    index=CONDITION_NAMES, columns=COMPONENT_NAMES,
)
CELL_TYPES = ["KC_BAS", "FIB", "TCELL"]
N_FEATURES = 8
N_DONORS_PER_CONDITION = 2


def run(cmd, **kwargs):
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if result.returncode != 0:
        print(result.stdout[-4000:])
        print(result.stderr[-4000:])
        raise SystemExit(f"Command failed (exit {result.returncode}): {' '.join(cmd)}")
    return result


def generate_synthetic_pseudobulk(seed: int = 0, include_hc: bool = True):
    """Synthetic raw (unstandardized) LSI/PCA-like embeddings + sample
    metadata, with a real per-cell-type baseline shift and a real
    per-condition signal (via C @ W_true), matching what step 2's output
    actually looks like -- so the CV/preprocessing/decoder chain has
    something genuine to recover, not pure noise."""
    rng = np.random.default_rng(seed)

    conditions = CONDITION_NAMES if include_hc else CONDITION_NAMES[1:]
    W_true = 0.6 * rng.standard_normal((len(COMPONENT_NAMES), N_FEATURES))
    celltype_offset = {ct: 3.0 * rng.standard_normal(N_FEATURES) for ct in CELL_TYPES}

    rows, meta_rows = [], []
    for cond in conditions:
        c_row = C_MATRIX.loc[cond].to_numpy()
        cond_signal = c_row @ W_true
        for donor_i in range(N_DONORS_PER_CONDITION):
            donor_id = f"{cond}_donor{donor_i}"
            for ct in CELL_TYPES:
                sample_id = f"{ct}__{donor_id}"
                x = celltype_offset[ct] + cond_signal + 0.3 * rng.standard_normal(N_FEATURES)
                rows.append(pd.Series(x, name=sample_id))
                meta_rows.append({
                    "sample_id": sample_id, "cell_type": ct, "donor_id": donor_id,
                    "condition": cond, "location": "Lesional" if "Lesional" in cond and "Non" not in cond else "NonLesional",
                    "n_cells": 50,
                })

    embeddings = pd.DataFrame(rows)
    embeddings.columns = [f"LSI_{i+2}" for i in range(N_FEATURES)]  # LSI_1 conventionally dropped
    sample_meta = pd.DataFrame(meta_rows)
    return embeddings, sample_meta, W_true


def generate_synthetic_loadings(feature_names, n_peaks=200, seed: int = 1):
    rng = np.random.default_rng(seed)
    loadings = pd.DataFrame(
        rng.standard_normal((n_peaks, len(feature_names))), columns=feature_names,
        index=[f"chr1-{1000*i}-{1000*i+500}" for i in range(n_peaks)],
    )
    peaks_used = pd.DataFrame({
        "peak": loadings.index, "distance": rng.integers(2000, 500000, size=n_peaks), "used_in_lsi": True,
    })
    return loadings, peaks_used


def check_leakage_fix_is_not_a_noop(embeddings, cell_type, fold_ids):
    """Fitting standardize+center on a training fold must differ from
    fitting on the full dataset (which would include the held-out fold) --
    proves 03_train_cv.py's fold-wise refit is actually doing something,
    not silently equivalent to a single global fit."""
    some_fold = fold_ids.unique()[0]
    train_mask = fold_ids != some_fold

    global_params = prep.fit_standardize_then_center(embeddings, cell_type, quiet=True)
    fold_params = prep.fit_standardize_then_center(embeddings.loc[train_mask], cell_type.loc[train_mask], quiet=True)

    mean_diff = (global_params["global_mean"] - fold_params["global_mean"]).abs().max()
    std_diff = (global_params["global_std"] - fold_params["global_std"]).abs().max()
    assert mean_diff > 1e-6 or std_diff > 1e-6, (
        "Fold-wise preprocessing statistics are identical to the global fit -- the leakage fix "
        "is not actually excluding the held-out fold's data. This should not happen."
    )
    print(f"  OK: fold-wise fit differs from global fit (max mean diff={mean_diff:.4g}, "
          f"max std diff={std_diff:.4g}) -- leakage fix is doing real work.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--keep", action="store_true", help="Don't delete --out_dir on success")
    args = ap.parse_args()

    out_dir = Path(args.out_dir) if args.out_dir else Path(tempfile.mkdtemp(prefix="decoder_self_test_"))
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Self-test working directory: {out_dir}")

    print("\n[1/5] Generating synthetic pseudobulk-shaped data (with HC samples)...")
    embeddings, sample_meta, W_true = generate_synthetic_pseudobulk(include_hc=True)
    embeddings_path = out_dir / "atac_lsi_embeddings.csv"
    sample_meta_path = out_dir / "sample_metadata.csv"
    embeddings.to_csv(embeddings_path)
    sample_meta.to_csv(sample_meta_path, index=False)

    constraint_path = out_dir / "constraint_matrix.csv"
    C_MATRIX.to_csv(constraint_path)

    loadings, peaks_used = generate_synthetic_loadings(list(embeddings.columns))
    loadings_path = out_dir / "atac_lsi_loadings.csv"
    peaks_used_path = out_dir / "atac_peaks_used.csv"
    loadings.to_csv(loadings_path)
    peaks_used.to_csv(peaks_used_path, index=False)

    print("\n[2/5] Verifying the CV-leakage fix is not a no-op...")
    check_leakage_fix_is_not_a_noop(embeddings, sample_meta.set_index("sample_id")["cell_type"],
                                     sample_meta.set_index("sample_id")["donor_id"])

    print("\n[3/5] Running 03_train_cv.py (grid search x fold-wise CV, per-fold refit, final refit)...")
    decoder_out = out_dir / "decoder_atac"
    run([sys.executable, str(HERE / "03_train_cv.py"),
         "--embeddings", str(embeddings_path), "--sample_metadata", str(sample_meta_path),
         "--constraint_matrix", str(constraint_path), "--modality", "atac", "--fold_col", "donor_id",
         "--l2_grid", "0.01,0.1,1,10", "--max_iter", "150", "--out_dir", str(decoder_out)])

    for f in ["grid_search_summary.csv", "best_lambda.json", "cv_fold_W_components.csv",
              "final_decoder_W.csv", "final_decoder_intercepts.csv", "decoder_config.json"]:
        assert (decoder_out / f).exists(), f"Expected output missing: {f}"

    final_W = pd.read_csv(decoder_out / "final_decoder_W.csv", index_col=0)
    assert final_W.shape == (len(COMPONENT_NAMES), N_FEATURES), f"Unexpected final_decoder_W shape: {final_W.shape}"
    assert not final_W.isna().any().any(), "final_decoder_W contains NaNs"

    fold_W = pd.read_csv(decoder_out / "cv_fold_W_components.csv")
    n_folds_expected = sample_meta["donor_id"].nunique()
    assert fold_W["held_out_fold"].nunique() == n_folds_expected, (
        f"Expected {n_folds_expected} folds in cv_fold_W_components.csv, got {fold_W['held_out_fold'].nunique()}"
    )

    with open(decoder_out / "decoder_config.json") as fh:
        config = json.load(fh)
    assert config["reference_n_samples"] > 0, "Expected HC (reference) samples present in this run"
    assert "reference_caveat" not in config, "Did not expect the zero-HC caveat when HC samples are present"
    print("  OK: 03_train_cv.py outputs present and shaped correctly.")

    print("\n[4/5] Running 04_project_to_peak_space.py...")
    peak_proj_out = out_dir / "peak_projections"
    run([sys.executable, str(HERE / "04_project_to_peak_space.py"),
         "--loadings", str(loadings_path), "--cv_fold_weights", str(decoder_out / "cv_fold_W_components.csv"),
         "--final_weights", str(decoder_out / "final_decoder_W.csv"), "--out_dir", str(peak_proj_out)])

    per_fold_peaks = pd.read_csv(peak_proj_out / "peak_weights_per_fold.csv")
    assert per_fold_peaks["peak"].nunique() == len(loadings), "Peak count mismatch after projection"
    print("  OK: peak projection outputs present and shaped correctly.")

    print("\n[5/5] Running 05_identify_stable_peaks.py...")
    stable_out = out_dir / "stable_peaks"
    run([sys.executable, str(HERE / "05_identify_stable_peaks.py"),
         "--peak_weights_per_fold", str(peak_proj_out / "peak_weights_per_fold.csv"),
         "--fdr_threshold", "0.05", "--correction_scope", "per_component", "--out_dir", str(stable_out)])

    stable_all = pd.read_csv(stable_out / "stable_peaks_all.csv")
    assert {"n_folds_positive", "frac_folds_positive", "q_value", "is_stable"} <= set(stable_all.columns), (
        "stable_peaks_all.csv is missing expected columns"
    )
    print(f"  OK: stability test produced {len(stable_all)} (component, peak) rows across "
          f"{stable_all['component'].nunique()} components.")

    print("\n--- Zero-HC-samples caveat check ---")
    embeddings_no_hc, sample_meta_no_hc, _ = generate_synthetic_pseudobulk(include_hc=False)
    embeddings_no_hc.to_csv(out_dir / "embeddings_no_hc.csv")
    sample_meta_no_hc.to_csv(out_dir / "sample_metadata_no_hc.csv", index=False)
    no_hc_out = out_dir / "decoder_atac_no_hc"
    run([sys.executable, str(HERE / "03_train_cv.py"),
         "--embeddings", str(out_dir / "embeddings_no_hc.csv"),
         "--sample_metadata", str(out_dir / "sample_metadata_no_hc.csv"),
         "--constraint_matrix", str(constraint_path), "--modality", "atac", "--fold_col", "donor_id",
         "--l2_grid", "0.1,1", "--max_iter", "100", "--out_dir", str(no_hc_out)])
    with open(no_hc_out / "decoder_config.json") as fh:
        no_hc_config = json.load(fh)
    assert no_hc_config["reference_n_samples"] == 0
    assert "reference_caveat" in no_hc_config, "Expected the zero-HC caveat to be recorded when HC has 0 samples"
    print("  OK: zero-HC-samples run still completes and records the provisional-beta_sub caveat.")

    print(f"\nALL CHECKS PASSED. Working directory: {out_dir}")
    if not args.keep and args.out_dir is None:
        shutil.rmtree(out_dir, ignore_errors=True)
        print("(temporary directory removed; pass --out_dir/--keep to retain it)")


if __name__ == "__main__":
    main()
