#!/usr/bin/env python3
"""
Decoder branch, step 3: cross-validated regularization grid search + final
decoder fit.

Takes step 2's RAW (unstandardized) LSI/PCA embeddings -- NOT
pre-standardized features -- and fits standardization + per-cell-type
centering fresh inside every training split (see preprocessing.py). Fitting
those statistics once globally before cross-validating leaks held-out-fold
information into the training features and biases the regularization grid
search optimistic; that's why this script owns the fit/apply split instead
of consuming an already-centered features_final.csv like the model's
standalone CLI does.

Fold unit (--fold_col, default donor_id -> leave-one-donor-out): this
constraint matrix has components that are nonzero for only a single
condition (e.g. beta_int for PSO_Lesional). Leave-one-CONDITION-out would,
in the fold that removes that condition, delete 100% of that component's
training signal -- not a real held-out test. Leave-one-DONOR-out keeps at
least one donor from every condition in every training split, so every
regulatory direction stays identifiable in every fold. Externalized as a
CLI flag (and configs/pipeline.config's decoder_fold_col), not hardcoded,
so this can be revisited once more donors exist -- but donor-level is the
correct default given this C matrix's structure, not just what happened to
get built first. See README.md's "Decoder methodology & caveats".

Three phases, all written to --out_dir:

  Phase A -- grid search: for every (lambda, fold) pair, fit standardization
    + centering AND the decoder on all data except that fold, evaluate on
    the held-out fold. Average across folds per lambda -> lowest mean
    held-out loss wins.
      -> grid_search_results.csv, grid_search_summary.csv, best_lambda.json

  Phase B -- per-fold refit at the selected lambda (fold-wise preprocessing
    again): needed because Step 5 needs one W per fold to test "consistently
    positive across folds", not just a final point estimate.
      -> cv_fold_W_components.csv (long format: fold, component, feature, weight)

  Phase C -- final refit on ALL samples (standardization/centering fit once,
    globally -- no CV here, so no leakage concern) at best_lambda.
      -> final_decoder_W.csv, final_decoder_intercepts.csv, decoder_config.json

USAGE:
    python 03_train_cv.py \
        --embeddings output/decoder/features/atac_lsi_embeddings.csv \
        --sample_metadata output/decoder/pseudobulk/sample_metadata.csv \
        --constraint_matrix configs/decoder/constraint_matrix.csv \
        --modality atac --fold_col donor_id \
        --out_dir output/decoder/decoder_atac
"""

import argparse
import json
import os

import pandas as pd
import torch

import decoder_model as dec
import preprocessing as prep


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--embeddings", required=True,
                    help="Path to RAW (unstandardized) atac_lsi_embeddings.csv or "
                         "rna_pca_embeddings.csv from step 2")
    p.add_argument("--sample_metadata", required=True,
                    help="Path to sample_metadata.csv from step 1 (sample_id, cell_type, "
                         "donor_id, condition, location)")
    p.add_argument("--constraint_matrix", required=True,
                    help="Path to constraint matrix C (CSV, conditions x components) -- "
                         "external, never hardcoded")
    p.add_argument("--reference_condition", default=None,
                    help="Reference condition name; auto-detected from C's all-zero row if omitted")
    p.add_argument("--modality", required=True, choices=["rna", "atac"],
                    help="Used for output file naming / provenance only")
    p.add_argument("--fold_col", default=os.environ.get("decoder_fold_col", "donor_id"),
                    help="sample_metadata column defining CV folds (leave-one-value-out)")
    p.add_argument("--min_cells_per_celltype", type=int, default=2,
                    help="Passed through to preprocessing's degenerate-centering warning")
    p.add_argument("--l2_grid", default="0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100",
                    help="Comma-separated L2 lambda grid")
    p.add_argument("--max_iter", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def load_data(args):
    embeddings = pd.read_csv(args.embeddings, index_col=0)
    sample_meta = pd.read_csv(args.sample_metadata).set_index("sample_id")
    sample_meta = sample_meta.loc[embeddings.index]

    C, condition_names, component_names = dec.load_constraint_matrix(args.constraint_matrix)
    reference_idx = dec.find_reference_condition(C, condition_names, args.reference_condition)

    cond_to_idx = {c: i for i, c in enumerate(condition_names)}
    missing_cond = set(sample_meta["condition"].unique()) - set(cond_to_idx)
    if missing_cond:
        raise ValueError(f"Condition(s) in sample_metadata not present in constraint matrix rows: {missing_cond}")

    if args.fold_col not in sample_meta.columns:
        raise ValueError(f"--fold_col '{args.fold_col}' not found in sample_metadata columns: "
                          f"{list(sample_meta.columns)}")
    if "cell_type" not in sample_meta.columns:
        raise ValueError(f"sample_metadata must contain 'cell_type'; found: {list(sample_meta.columns)}")

    y_idx = sample_meta["condition"].map(cond_to_idx)
    fold_ids = sample_meta[args.fold_col]
    cell_type = sample_meta["cell_type"]
    feature_names = list(embeddings.columns)

    return embeddings, y_idx, fold_ids, cell_type, C, condition_names, component_names, reference_idx, feature_names


def preprocess_split(embeddings, cell_type, train_mask, val_mask=None, min_cells_per_celltype=2, quiet=True):
    """Fit standardize+center on the training rows only; apply to train (and
    optionally a disjoint val split) using those training statistics."""
    params = prep.fit_standardize_then_center(
        embeddings.loc[train_mask], cell_type.loc[train_mask],
        min_cells_per_celltype=min_cells_per_celltype, quiet=quiet,
    )
    X_train = prep.apply_standardize_then_center(embeddings.loc[train_mask], cell_type.loc[train_mask],
                                                  params, quiet=True)
    if val_mask is None:
        return X_train, params
    X_val = prep.apply_standardize_then_center(embeddings.loc[val_mask], cell_type.loc[val_mask],
                                                params, quiet=True)
    return X_train, X_val, params


def to_tensor(df):
    return torch.tensor(df.to_numpy(), dtype=torch.float32)


def evaluate_fold(X_val, y_val, C, W, intercepts_full, reference_idx):
    """Unweighted CE + accuracy on a held-out fold. Held-out folds may contain
    only one condition (donors are nested 1:1 within condition in this design),
    so equal-condition-weighting isn't meaningful here -- a plain mean is used."""
    logits = dec.forward(X_val, C, W,
                          intercepts_nonref=torch.cat([
                              intercepts_full[:reference_idx], intercepts_full[reference_idx + 1:]
                          ]),
                          reference_idx=reference_idx)
    log_probs = torch.log_softmax(logits, dim=1)
    ce = -log_probs.gather(1, y_val.unsqueeze(1)).squeeze(1).mean().item()
    pred = logits.argmax(dim=1)
    acc = (pred == y_val).float().mean().item()
    return ce, acc


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    (embeddings, y_idx, fold_ids, cell_type, C, condition_names, component_names,
     reference_idx, feature_names) = load_data(args)
    unique_folds = sorted(fold_ids.unique().tolist())
    l2_grid = [float(x) for x in args.l2_grid.split(",")]

    print(f"[{args.modality.upper()}] {embeddings.shape[0]} samples x {embeddings.shape[1]} raw components")
    print(f"  Conditions ({len(condition_names)}): {condition_names}")
    print(f"  Reference: {condition_names[reference_idx]}")
    print(f"  Fold column: {args.fold_col} -> {len(unique_folds)} folds: {unique_folds}")
    print(f"  L2 grid ({len(l2_grid)} values): {l2_grid}")

    # -- Phase A: grid search (fold-wise preprocessing -- see module docstring) --
    print("\n[Phase A] Grid search over L2 lambda x leave-one-fold-out CV...")
    grid_rows = []
    for lam in l2_grid:
        for fold in unique_folds:
            train_mask = fold_ids != fold
            val_mask = ~train_mask

            X_train_df, X_val_df, _ = preprocess_split(
                embeddings, cell_type, train_mask, val_mask,
                min_cells_per_celltype=args.min_cells_per_celltype,
            )
            X_train, X_val = to_tensor(X_train_df), to_tensor(X_val_df)
            y_train = torch.tensor(y_idx[train_mask].to_numpy(), dtype=torch.long)
            y_val = torch.tensor(y_idx[val_mask].to_numpy(), dtype=torch.long)

            result = dec.fit_decoder(
                X_train, y_train, C, reference_idx, lam,
                max_iter=args.max_iter, seed=args.seed,
                condition_names=condition_names, verbose=False,
            )
            val_ce, val_acc = evaluate_fold(X_val, y_val, C, result["W"], result["intercepts_full"], reference_idx)

            grid_rows.append({
                "lambda": lam, "fold": fold, "train_loss": result["final_loss"],
                "val_ce": val_ce, "val_accuracy": val_acc,
                "n_train": int(train_mask.sum()), "n_val": int(val_mask.sum()),
            })
        print(f"  lambda={lam:g} done ({len(unique_folds)} folds)")

    grid_df = pd.DataFrame(grid_rows)
    grid_df.to_csv(os.path.join(args.out_dir, "grid_search_results.csv"), index=False)

    grid_summary = grid_df.groupby("lambda").agg(
        mean_val_ce=("val_ce", "mean"), std_val_ce=("val_ce", "std"),
        mean_val_accuracy=("val_accuracy", "mean"),
    ).reset_index().sort_values("lambda")
    grid_summary.to_csv(os.path.join(args.out_dir, "grid_search_summary.csv"), index=False)

    best_lambda = float(grid_summary.loc[grid_summary["mean_val_ce"].idxmin(), "lambda"])
    print(f"\n  Best lambda (lowest mean held-out CE): {best_lambda:g}")
    print(grid_summary.to_string(index=False))

    with open(os.path.join(args.out_dir, "best_lambda.json"), "w") as f:
        json.dump({"best_lambda": best_lambda, "selection_criterion": "min mean held-out CE across folds"}, f, indent=2)

    # -- Phase B: per-fold refit at best_lambda (fold-wise preprocessing again) --
    print(f"\n[Phase B] Refitting each leave-one-{args.fold_col}-out fold at lambda={best_lambda:g}...")
    fold_W_rows = []
    for fold in unique_folds:
        train_mask = fold_ids != fold
        X_train_df, _ = preprocess_split(embeddings, cell_type, train_mask,
                                          min_cells_per_celltype=args.min_cells_per_celltype)
        X_train = to_tensor(X_train_df)
        y_train = torch.tensor(y_idx[train_mask].to_numpy(), dtype=torch.long)

        result = dec.fit_decoder(
            X_train, y_train, C, reference_idx, best_lambda,
            max_iter=args.max_iter, seed=args.seed,
            condition_names=condition_names, verbose=False,
        )
        W = result["W"].numpy()
        for ci, comp in enumerate(component_names):
            for fi, feat in enumerate(feature_names):
                fold_W_rows.append({"held_out_fold": fold, "component": comp, "feature": feat, "weight": W[ci, fi]})
        print(f"  fold={fold} done (final_loss={result['final_loss']:.6f})")

    fold_W_df = pd.DataFrame(fold_W_rows)
    fold_W_df.to_csv(os.path.join(args.out_dir, "cv_fold_W_components.csv"), index=False)
    print(f"  Written: cv_fold_W_components.csv "
          f"({len(unique_folds)} folds x {len(component_names)} components x {len(feature_names)} features)")

    # -- Phase C: final refit on all data (global fit, no CV -> no leakage concern) --
    print(f"\n[Phase C] Final refit on all {embeddings.shape[0]} samples at lambda={best_lambda:g}...")
    X_final_df, final_prep_params = prep.standardize_then_center(
        embeddings, cell_type, min_cells_per_celltype=args.min_cells_per_celltype,
    )
    X_final = to_tensor(X_final_df)
    y_final = torch.tensor(y_idx.to_numpy(), dtype=torch.long)

    final_result = dec.fit_decoder(
        X_final, y_final, C, reference_idx, best_lambda,
        max_iter=args.max_iter, seed=args.seed,
        condition_names=condition_names, verbose=True,
    )

    final_W = pd.DataFrame(final_result["W"].numpy(), index=component_names, columns=feature_names)
    final_W.to_csv(os.path.join(args.out_dir, "final_decoder_W.csv"))

    final_intercepts = pd.Series(final_result["intercepts_full"].numpy(), index=condition_names, name="intercept")
    final_intercepts.to_frame().to_csv(os.path.join(args.out_dir, "final_decoder_intercepts.csv"))

    final_prep_params["global_mean"].to_frame("global_mean").to_csv(
        os.path.join(args.out_dir, "final_global_mean.csv"))
    final_prep_params["global_std"].to_frame("global_std").to_csv(
        os.path.join(args.out_dir, "final_global_std.csv"))
    final_prep_params["celltype_means"].to_csv(os.path.join(args.out_dir, "final_celltype_means.csv"))

    n_ref = int((y_idx == reference_idx).sum())
    config = {
        "modality": args.modality, "condition_names": condition_names, "component_names": component_names,
        "feature_names": feature_names, "reference_condition": condition_names[reference_idx],
        "reference_n_samples": n_ref, "fold_col": args.fold_col, "n_folds": len(unique_folds),
        "l2_grid": l2_grid, "best_lambda": best_lambda, "final_train_loss": final_result["final_loss"],
        "max_iter": args.max_iter, "seed": args.seed,
    }
    if n_ref == 0:
        config["reference_caveat"] = (
            "Reference condition has zero samples -- beta_sub (and any component whose "
            "identifiability leans on the reference) is not yet anchored against real "
            "reference-condition data. Treat as provisional. See README.md."
        )
        print("  WARNING: reference condition has 0 samples -- see decoder_config.json's "
              "'reference_caveat' field.")
    with open(os.path.join(args.out_dir, "decoder_config.json"), "w") as f:
        json.dump(config, f, indent=2)

    print(f"\nDone. Output in: {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
