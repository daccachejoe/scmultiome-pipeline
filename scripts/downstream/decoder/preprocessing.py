"""
Decoder branch: feature preprocessing -- standardize, then center by cell
type. Per the methodology: "Before model fitting, features were first
standardized, then centered by each cell type."

  1. Standardize (z-score) each component GLOBALLY (mean 0, unit variance).
  2. Center each component BY CELL TYPE (subtract the per-cell-type mean
     of the already-standardized values), removing cell-type baseline
     differences so a single pooled decoder can be fit across cell types.

This order is intentional -- standardizing after centering would let
per-cell-type scale differences leak back in.

fit/apply are split deliberately (rather than one combined
`standardize_then_center`) so callers doing cross-validation can fit on a
training split only and apply those exact statistics to the held-out
split. Fitting standardization+centering on the full dataset before a CV
split leaks held-out information into the training features and biases
regularization-strength selection optimistic -- see 03_train_cv.py, which
calls fit_standardize_then_center() once per training fold rather than
once globally.
"""

import numpy as np
import pandas as pd


def fit_standardize_then_center(embeddings: pd.DataFrame, cell_type: pd.Series,
                                 min_cells_per_celltype: int = 2, quiet: bool = False):
    """
    embeddings: samples x components, index = sample_id
    cell_type : Series indexed by sample_id (aligned to embeddings.index)

    Returns a dict of fitted parameters: global_mean, global_std,
    celltype_means -- everything apply_standardize_then_center() needs.
    """
    global_mean = embeddings.mean(axis=0)
    global_std = embeddings.std(axis=0, ddof=1)

    zero_var = global_std[global_std == 0]
    if len(zero_var) > 0:
        raise ValueError(
            f"Component(s) with zero variance in this split: {list(zero_var.index)}. "
            "Cannot standardize -- check upstream LSI/PCA output, or that this split "
            "has more than one sample."
        )

    standardized = (embeddings - global_mean) / global_std
    celltype_means = standardized.groupby(cell_type).mean()

    if not quiet:
        counts = cell_type.value_counts()
        small = counts[counts < min_cells_per_celltype]
        if len(small) > 0:
            print(f"  WARNING: {len(small)} cell type(s) have fewer than "
                  f"{min_cells_per_celltype} samples in this split -- per-cell-type "
                  f"centering will make their contribution degenerate (exactly 0 if n=1):")
            for ct, n in small.items():
                print(f"    {ct}: {n} sample(s)")

    return {"global_mean": global_mean, "global_std": global_std, "celltype_means": celltype_means}


def apply_standardize_then_center(embeddings: pd.DataFrame, cell_type: pd.Series, params: dict,
                                   quiet: bool = False):
    """
    Applies previously-fit standardization + cell-type centering statistics
    to (possibly different) samples. Used both to transform the split a
    fit came from, and to transform a held-out split using the training
    split's statistics.

    A cell type present here but absent from params["celltype_means"]
    (e.g. a cell type entirely held out of the training fold) is left at
    the globally-standardized value -- centering can't be estimated for
    it from this fit, so it falls back to the global mean/scale only.
    """
    standardized = (embeddings - params["global_mean"]) / params["global_std"]
    centered = standardized.copy()

    known_celltypes = set(params["celltype_means"].index)
    unknown = set(cell_type.unique()) - known_celltypes
    if unknown and not quiet:
        print(f"  WARNING: cell type(s) {sorted(unknown)} not present in the fitted "
              f"centering statistics -- left globally-standardized (not cell-type-centered) "
              f"for these samples.")

    for ct in cell_type.unique():
        if ct not in known_celltypes:
            continue
        mask = cell_type == ct
        centered.loc[mask] = standardized.loc[mask] - params["celltype_means"].loc[ct]

    return centered


def standardize_then_center(embeddings: pd.DataFrame, cell_type: pd.Series,
                             min_cells_per_celltype: int = 2, quiet: bool = False):
    """
    Convenience wrapper: fit and apply on the same data in one call. Only
    appropriate when there is no train/test split to keep separate --
    i.e. the final, all-data (Phase C) decoder fit. Cross-validated fits
    (Phase A/B in 03_train_cv.py) must call fit_/apply_ separately so
    training-fold statistics never see held-out data.
    """
    params = fit_standardize_then_center(embeddings, cell_type, min_cells_per_celltype, quiet)
    final = apply_standardize_then_center(embeddings, cell_type, params, quiet=True)
    return final, params
