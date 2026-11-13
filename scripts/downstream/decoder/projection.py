"""
Decoder branch: shared projection math for mapping a decoder's fitted
component weights (in whatever reduced-dimensionality space the decoder
was trained on -- LSI dims for ATAC, PCA dims for RNA) back to individual
peaks or genes, via the same loadings matrix that created that reduced
space in step 2.

    item_weight[item, component, fold] = sum_c( loadings[item, c] * W[component, c, fold] )

i.e. item_weights = loadings @ W.T -- literally undoing step 2's
dimensionality reduction, using its own loadings matrix as the
reconstruction operator, the same way you'd read PCA loadings to see
which original features drive a component.

Used by both 04_project_to_peak_space.py (ATAC: LSI loadings -> peaks)
and 07_project_to_gene_space.py (RNA: PCA loadings -> genes) -- the math
is identical between modalities, only the loadings matrix and the output
item label ("peak" vs "gene") differ, which is why this lives in one
shared module instead of being duplicated in both scripts.
"""

import pandas as pd


def project_fold_weights(loadings: pd.DataFrame, fold_weights_long: pd.DataFrame, item_label: str = "item"):
    """
    loadings          : items (peaks or genes) x decoder-input-dims (LSI or PCA components)
    fold_weights_long : long format (held_out_fold, component, feature, weight) from
                         03_train_cv.py's cv_fold_W_components.csv -- "feature" here means
                         the decoder's INPUT dimensions (LSI/PCA), not the items being
                         projected to; "component" means the biological regulatory
                         direction (beta_sub, beta_dis, ...).
    item_label        : name for the output item column, e.g. "peak" or "gene"

    Returns long format: held_out_fold, component, <item_label>, weight -- one row per
    (fold, component, item), the input 05_identify_stable_peaks.py's stability test needs.
    """
    cols_in_weights = sorted(fold_weights_long["feature"].unique())
    cols_in_loadings = sorted(loadings.columns)
    if cols_in_weights != cols_in_loadings:
        only_in_weights = set(cols_in_weights) - set(cols_in_loadings)
        only_in_loadings = set(cols_in_loadings) - set(cols_in_weights)
        raise ValueError(
            "Component-space dimension names in cv_fold_weights do not match the loadings "
            "columns exactly. Usually means step 2's settings (e.g. --keep_lsi_dim1) differed "
            "between the loadings export and what fed the decoder, or the two files came from "
            "different runs.\n"
            f"  In weights but not loadings: {sorted(only_in_weights)}\n"
            f"  In loadings but not weights: {sorted(only_in_loadings)}"
        )

    results = []
    # One (fold, component) pair at a time: w is that pair's weight vector in
    # the decoder's input space (length = n LSI/PCA dims), reindexed to
    # loadings' column order so the matrix multiply lines up positionally,
    # not just by having "the same names somewhere".
    for (fold, component), sub in fold_weights_long.groupby(["held_out_fold", "component"]):
        w = sub.set_index("feature")["weight"].reindex(loadings.columns)
        item_weights = loadings.to_numpy() @ w.to_numpy()
        results.append(pd.DataFrame({
            "held_out_fold": fold, "component": component, item_label: loadings.index, "weight": item_weights,
        }))
    return pd.concat(results, axis=0, ignore_index=True)


def project_final_weights(loadings: pd.DataFrame, final_W: pd.DataFrame):
    """Same projection as project_fold_weights(), but for every component at
    once (all rows of final_W multiplied through in one matmul) instead of
    fold-by-fold -- this is the point estimate for reporting, not an input
    to the stability test, so there's no need to keep folds separate here."""
    final_W = final_W[loadings.columns]  # align column order to loadings before the matmul
    item_weights = loadings.to_numpy() @ final_W.to_numpy().T
    return pd.DataFrame(item_weights, index=loadings.index, columns=final_W.index)
