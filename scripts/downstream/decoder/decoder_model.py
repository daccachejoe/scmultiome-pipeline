#!/usr/bin/env python3
"""
Decoder branch: constrained multinomial logistic regression model
(PyTorch/LBFGS). Importable module used by 03_train_cv.py. Also runnable
standalone with --self_test (synthetic data, no real inputs needed) or
--self_test_e2e (exercises the full preprocess -> CV -> peak-projection ->
stability chain against synthetic pseudobulk-shaped data -- see
self_test_e2e.py).

MODEL
-----
Condition-specific decoder weight vectors are NOT fit independently.
Instead:

    condition_weights = C @ W        # (n_conditions, n_features)
    logits            = X @ condition_weights.T + intercepts
                       = X @ (C @ W).T + intercepts

  X : (n_samples, n_features)          standardized+centered LSI/PCA features
  C : (n_conditions, n_components_bio) constraint matrix -- EXTERNAL PARAMETER,
                                        loaded from file, never hardcoded here
  W : (n_components_bio, n_features)   learned component weight vectors
  intercepts : (n_conditions,)         free for non-reference classes,
                                        fixed at 0 for the reference class
                                        (row of C that is all-zero)

LOSS
----
Weighted, optionally masked, softmax cross-entropy + L2 penalty on W only:

    loss = sum_i[ w_i * CE_i(mask_i) ] / sum_i[ w_i ]  +  lambda * sum(W ** 2)

  w_i    : per-sample weight so every condition contributes equally to the
           loss, w_i = 1 / n_samples_in_condition(y_i)
  mask_i : optional per-sample class mask (default: all classes valid)

FITTING
-------
torch.optim.LBFGS, strong-Wolfe line search, full-batch -- this is a
smooth, low-dimensional convex problem (params = n_components_bio *
n_features + n_conditions - 1), not a mini-batch deep-learning problem.

Reference condition: the constraint-matrix row that is entirely zero
(conventionally HC). Auto-detected unless --reference_condition is given.

CAVEAT (HC with zero samples): the model still fits with zero reference
samples -- the reference intercept is fixed at 0 by construction and
beta_sub is estimated from X-dependent variation in the other conditions
-- but until real HC samples exist, beta_sub is never directly anchored
against an observed healthy baseline. Treat it as provisional until then.
See README.md's "Decoder methodology & caveats" section.
"""

import argparse
import json

import numpy as np
import pandas as pd
import torch


# =============================================================================
# Constraint matrix I/O
# =============================================================================

def load_constraint_matrix(path: str):
    """
    Load C from a CSV with conditions as rows (first column = condition name)
    and biological components as columns (header row = component names).

    Returns (C: torch.FloatTensor [n_conditions, n_components], condition_names: list,
             component_names: list)
    """
    df = pd.read_csv(path, index_col=0)
    condition_names = list(df.index)
    component_names = list(df.columns)
    C = torch.tensor(df.to_numpy(), dtype=torch.float32)
    return C, condition_names, component_names


def find_reference_condition(C: torch.Tensor, condition_names: list, explicit: str = None):
    """
    Identify the reference (baseline) condition -- the row of C that is
    all-zero. Raises if ambiguous (0 or >1 all-zero rows) unless an
    explicit condition name is given.
    """
    if explicit is not None:
        if explicit not in condition_names:
            raise ValueError(f"--reference_condition '{explicit}' not found in {condition_names}")
        return condition_names.index(explicit)

    zero_rows = [i for i in range(C.shape[0]) if torch.all(C[i] == 0)]
    if len(zero_rows) == 0:
        raise ValueError(
            "No all-zero row found in the constraint matrix C -- cannot auto-detect the "
            "reference condition. Pass --reference_condition explicitly."
        )
    if len(zero_rows) > 1:
        names = [condition_names[i] for i in zero_rows]
        raise ValueError(
            f"Multiple all-zero rows found in C ({names}) -- ambiguous reference condition. "
            "Pass --reference_condition explicitly."
        )
    return zero_rows[0]


# =============================================================================
# Model
# =============================================================================

def forward(X: torch.Tensor, C: torch.Tensor, W: torch.Tensor,
            intercepts_nonref: torch.Tensor, reference_idx: int):
    """
    X                 : (n_samples, n_features)
    C                 : (n_conditions, n_components_bio)
    W                 : (n_components_bio, n_features)
    intercepts_nonref : (n_conditions - 1,) -- one per non-reference condition,
                         in condition order with the reference index skipped
    reference_idx     : int, index of the reference condition in C's row order

    Returns logits: (n_samples, n_conditions)
    """
    n_cond = C.shape[0]
    condition_weights = C @ W                       # (n_conditions, n_features)
    logits = X @ condition_weights.T                # (n_samples, n_conditions)

    device = X.device
    full_intercepts = torch.zeros(n_cond, dtype=X.dtype, device=device)
    nonref_positions = [i for i in range(n_cond) if i != reference_idx]
    full_intercepts[nonref_positions] = intercepts_nonref

    return logits + full_intercepts


def compute_sample_weights(y_idx: torch.Tensor, n_conditions: int, warn_missing: bool = True,
                            condition_names: list = None, reference_idx: int = None):
    """
    Per-sample weight so every condition contributes equally to the loss:
    w_i = 1 / (# samples with the same condition as sample i).

    A zero-count condition never appears in y_idx by definition, so it's
    harmless here, not an error condition -- notably true for the
    reference (HC) condition while it has no demultiplexed samples yet,
    since its intercept is fixed at 0 by construction and needs no
    training examples to anchor. A non-reference condition with zero
    samples, however, means that condition's contrast is entirely
    unidentifiable from this fit -- flagged as a warning so it isn't
    missed silently.
    """
    counts = torch.bincount(y_idx, minlength=n_conditions).float()
    safe_counts = counts.clone()
    safe_counts[safe_counts == 0] = float("inf")

    if warn_missing:
        missing = [i for i in range(n_conditions) if counts[i] == 0 and i != reference_idx]
        if missing:
            names = [condition_names[i] for i in missing] if condition_names else missing
            print(f"    WARNING: condition(s) with zero training samples in this fold/fit "
                  f"(non-reference, unidentifiable): {names}")

    return 1.0 / safe_counts[y_idx]


def weighted_masked_cross_entropy(logits: torch.Tensor, y_idx: torch.Tensor,
                                   sample_weights: torch.Tensor, mask: torch.Tensor = None):
    """
    logits         : (n_samples, n_conditions)
    y_idx          : (n_samples,) integer class labels
    sample_weights : (n_samples,)
    mask           : optional (n_samples, n_conditions) bool, True = valid
                      target class for that sample.

    Returns scalar weighted mean cross-entropy loss.
    """
    if mask is not None:
        masked_logits = logits.masked_fill(~mask, float("-inf"))
        true_class_masked = ~mask.gather(1, y_idx.unsqueeze(1)).squeeze(1)
        if torch.any(true_class_masked):
            raise ValueError("mask excludes the true class for at least one sample.")
    else:
        masked_logits = logits

    log_probs = torch.log_softmax(masked_logits, dim=1)
    per_sample_ce = -log_probs.gather(1, y_idx.unsqueeze(1)).squeeze(1)

    return (per_sample_ce * sample_weights).sum() / sample_weights.sum()


# =============================================================================
# Fitting
# =============================================================================

def fit_decoder(X: torch.Tensor, y_idx: torch.Tensor, C: torch.Tensor,
                 reference_idx: int, l2_lambda: float,
                 sample_weights: torch.Tensor = None, mask: torch.Tensor = None,
                 W_init: torch.Tensor = None, intercept_init: torch.Tensor = None,
                 max_iter: int = 500, tol: float = 1e-9, seed: int = 0, verbose: bool = False,
                 condition_names: list = None):
    """
    Fit W (n_components_bio x n_features) and non-reference intercepts by
    minimizing weighted masked cross-entropy + L2 on W, via full-batch
    LBFGS with strong-Wolfe line search.
    """
    torch.manual_seed(seed)

    n_components_bio, n_features = C.shape[1], X.shape[1]
    n_cond = C.shape[0]

    if sample_weights is None:
        sample_weights = compute_sample_weights(y_idx, n_cond, condition_names=condition_names,
                                                  reference_idx=reference_idx)

    W = torch.nn.Parameter(
        W_init.clone() if W_init is not None
        else 0.01 * torch.randn(n_components_bio, n_features, dtype=X.dtype)
    )
    intercepts_nonref = torch.nn.Parameter(
        intercept_init.clone() if intercept_init is not None
        else torch.zeros(n_cond - 1, dtype=X.dtype)
    )

    optimizer = torch.optim.LBFGS(
        [W, intercepts_nonref], max_iter=max_iter, tolerance_grad=tol, tolerance_change=tol,
        line_search_fn="strong_wolfe",
    )

    loss_history = []

    def closure():
        optimizer.zero_grad()
        logits = forward(X, C, W, intercepts_nonref, reference_idx)
        ce = weighted_masked_cross_entropy(logits, y_idx, sample_weights, mask)
        reg = l2_lambda * (W ** 2).sum()
        loss = ce + reg
        loss.backward()
        loss_history.append(loss.item())
        return loss

    optimizer.step(closure)
    final_loss = closure().item()

    if verbose:
        print(f"    LBFGS: {len(loss_history)} closure evals, final loss = {final_loss:.6f}")

    full_intercepts = torch.zeros(n_cond, dtype=X.dtype)
    nonref_positions = [i for i in range(n_cond) if i != reference_idx]
    full_intercepts[nonref_positions] = intercepts_nonref.detach()

    return {
        "W": W.detach(), "intercepts_nonref": intercepts_nonref.detach(),
        "intercepts_full": full_intercepts, "final_loss": final_loss,
        "n_closure_evals": len(loss_history), "loss_history": loss_history,
    }


def predict_proba(X: torch.Tensor, C: torch.Tensor, W: torch.Tensor, intercepts_full: torch.Tensor):
    """Return softmax class probabilities (n_samples, n_conditions)."""
    condition_weights = C @ W
    logits = X @ condition_weights.T + intercepts_full
    return torch.softmax(logits, dim=1)


# =============================================================================
# Self-test (synthetic data) -- run after install to sanity-check the
# environment (torch/pandas versions) before touching real data.
# =============================================================================

def _self_test():
    print("Running decoder_model.py self-test on synthetic data...")
    torch.manual_seed(0)
    np.random.seed(0)

    condition_names = ["HC", "AD_NonLesional", "PSO_NonLesional", "AD_Lesional", "PSO_Lesional"]
    component_names = ["beta_sub", "beta_dis", "beta_loc", "beta_int"]
    C = torch.tensor([
        [0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 1, 1],
    ], dtype=torch.float32)

    n_features = 10
    n_per_condition = 20
    W_true = 0.5 * torch.randn(len(component_names), n_features)

    X_list, y_list = [], []
    for i, cond in enumerate(condition_names):
        cw = C[i] @ W_true
        # multiplier 2.0 chosen so this is a reliable smoke test, not a
        # coin flip: at 0.8 (as originally drafted) HC-vs-other separation
        # was too weak (cw=0 for HC by construction) and training accuracy
        # landed around 0.55-0.60 on repeated runs, well under the 0.8 bar.
        Xi = torch.randn(n_per_condition, n_features) + 2.0 * cw.unsqueeze(0)
        X_list.append(Xi)
        y_list.append(torch.full((n_per_condition,), i, dtype=torch.long))

    X = torch.cat(X_list, dim=0)
    y = torch.cat(y_list, dim=0)

    reference_idx = find_reference_condition(C, condition_names)
    assert condition_names[reference_idx] == "HC"

    result = fit_decoder(X, y, C, reference_idx, l2_lambda=0.01, max_iter=200, verbose=True)

    probs = predict_proba(X, C, result["W"], result["intercepts_full"])
    pred = probs.argmax(dim=1)
    acc = (pred == y).float().mean().item()
    print(f"  Training accuracy on separable synthetic data: {acc:.3f} (expect > 0.8)")
    assert acc > 0.8, "Self-test failed: decoder could not fit separable synthetic data well."

    corr = torch.corrcoef(torch.stack([result["W"].flatten(), W_true.flatten()]))[0, 1].item()
    print(f"  Correlation between recovered W and true W: {corr:.3f} (expect > 0.5)")
    assert corr > 0.5, "Self-test failed: recovered W does not correlate with the generating W."

    print("Self-test PASSED.")


# =============================================================================
# CLI (single fit, for debugging -- full CV grid search lives in 03_train_cv.py)
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--self_test", action="store_true",
                    help="Run the synthetic-data self-test and exit (no real data needed)")
    p.add_argument("--features", help="Path to a features_final.csv (already preprocessed)")
    p.add_argument("--sample_metadata", help="Path to sample_metadata.csv from step 1")
    p.add_argument("--constraint_matrix", help="Path to constraint matrix C (CSV, conditions x components)")
    p.add_argument("--reference_condition", default=None,
                    help="Name of the reference condition; auto-detected from C's all-zero row if omitted")
    p.add_argument("--l2_lambda", type=float, default=1.0, help="L2 regularization strength on W")
    p.add_argument("--max_iter", type=int, default=500)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out_json", default=None, help="Optional path to dump fit summary as JSON")
    return p.parse_args()


def main():
    args = parse_args()

    if args.self_test:
        _self_test()
        return

    if not (args.features and args.sample_metadata and args.constraint_matrix):
        raise SystemExit(
            "Provide --features, --sample_metadata, and --constraint_matrix for a real fit, "
            "or pass --self_test to run the synthetic sanity check."
        )

    features = pd.read_csv(args.features, index_col=0)
    sample_meta = pd.read_csv(args.sample_metadata).set_index("sample_id")
    sample_meta = sample_meta.loc[features.index]

    C, condition_names, component_names = load_constraint_matrix(args.constraint_matrix)
    reference_idx = find_reference_condition(C, condition_names, args.reference_condition)

    cond_to_idx = {c: i for i, c in enumerate(condition_names)}
    missing = set(sample_meta["condition"]) - set(cond_to_idx)
    if missing:
        raise ValueError(f"Condition(s) in sample_metadata not present in constraint matrix rows: {missing}")

    y_idx = torch.tensor(sample_meta["condition"].map(cond_to_idx).to_numpy(), dtype=torch.long)
    X = torch.tensor(features.to_numpy(), dtype=torch.float32)

    print(f"Fitting decoder: {X.shape[0]} samples x {X.shape[1]} features, "
          f"{len(condition_names)} conditions, {len(component_names)} biological components")
    print(f"  Reference condition: {condition_names[reference_idx]}")
    print(f"  L2 lambda: {args.l2_lambda}")

    result = fit_decoder(X, y_idx, C, reference_idx, args.l2_lambda,
                          max_iter=args.max_iter, seed=args.seed, verbose=True)
    print(f"Final loss: {result['final_loss']:.6f}")

    if args.out_json:
        summary = {
            "final_loss": result["final_loss"], "n_closure_evals": result["n_closure_evals"],
            "condition_names": condition_names, "component_names": component_names,
            "reference_condition": condition_names[reference_idx], "l2_lambda": args.l2_lambda,
        }
        with open(args.out_json, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Written: {args.out_json}")


if __name__ == "__main__":
    main()
