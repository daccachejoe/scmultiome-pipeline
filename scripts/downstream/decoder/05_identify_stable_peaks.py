#!/usr/bin/env python3
"""
Decoder branch, step 5: identify peaks with consistently positive weights
across CV folds.

For each regulatory direction and each peak, runs a one-sided one-sample
t-test of H0: mean(weight across folds) <= 0 vs H1: mean(weight) > 0, on the
per-fold peak weights from step 4, then applies Benjamini-Hochberg FDR
correction (--correction_scope per_component by default: each regulatory
direction is its own hypothesis family; pass --correction_scope global to
pool all direction x peak p-values into one correction instead).

CAVEAT -- fold non-independence: with leave-one-donor-out folds, any two
folds share ~87.5% of their training data, so per-fold peak weights are
correlated, not iid draws -- the t-test's p-values are anti-conservative
(more "significant" than a genuinely independent-samples test would give).
This is inherent to the published method's CV-based stability-selection
design, not a bug here, so the t-test is kept for fidelity to the
methodology. To avoid depending solely on the optimistic p-values, this
script also reports n_folds_positive / frac_folds_positive -- the raw
count/fraction of folds where a peak's weight was positive, a simpler
statistic that doesn't assume independence and is worth cross-checking
q-value-based calls against, especially near the --fdr_threshold boundary.

A one-sided test in the POSITIVE direction only identifies peaks whose
accessibility increases along a given regulatory direction. Rerun with
--alternative less for the negative/repressed side, or --alternative
two-sided for a non-directional stability test.

USAGE:
    python 05_identify_stable_peaks.py \
        --peak_weights_per_fold output/decoder/peak_projections/peak_weights_per_fold.csv \
        --fdr_threshold 0.05 --correction_scope per_component \
        --out_dir output/decoder/stable_peaks
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peak_weights_per_fold", required=True,
                    help="Path to peak_weights_per_fold.csv from step 4")
    p.add_argument("--fdr_threshold", type=float, default=0.05,
                    help="BH-adjusted p-value (q-value) cutoff for calling a peak stable")
    p.add_argument("--correction_scope", choices=["per_component", "global"], default="per_component")
    p.add_argument("--alternative", choices=["greater", "less", "two-sided"], default="greater")
    p.add_argument("--min_folds", type=int, default=3,
                    help="Minimum folds with data required to test a peak")
    p.add_argument("--peak_coord_regex", default=r"^(?P<chr>.+)[:\-_](?P<start>\d+)[:\-_](?P<end>\d+)$")
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def one_sided_stability_test(peak_weights_per_fold: pd.DataFrame, alternative: str, min_folds: int):
    """
    Returns a DataFrame indexed by (component, peak) with columns:
    mean_weight, std_weight, n_folds, n_folds_positive, frac_folds_positive,
    t_stat, p_value.
    """
    rows = []
    for (component, peak), sub in peak_weights_per_fold.groupby(["component", "peak"]):
        w = sub["weight"].to_numpy()
        n = len(w)
        n_pos = int((w > 0).sum())
        row = {
            "component": component, "peak": peak, "mean_weight": np.mean(w),
            "std_weight": np.std(w, ddof=1) if n > 1 else np.nan, "n_folds": n,
            "n_folds_positive": n_pos, "frac_folds_positive": n_pos / n,
        }
        if n < min_folds:
            row.update({"t_stat": np.nan, "p_value": np.nan})
            rows.append(row)
            continue
        t_stat, p_two_sided = stats.ttest_1samp(w, popmean=0.0)
        if alternative == "greater":
            p_val = p_two_sided / 2 if t_stat > 0 else 1 - p_two_sided / 2
        elif alternative == "less":
            p_val = p_two_sided / 2 if t_stat < 0 else 1 - p_two_sided / 2
        else:
            p_val = p_two_sided
        row.update({"t_stat": t_stat, "p_value": p_val})
        rows.append(row)
    return pd.DataFrame(rows)


def apply_bh_correction(df: pd.DataFrame, scope: str, fdr_threshold: float):
    df = df.copy()
    df["q_value"] = np.nan
    testable = df["p_value"].notna()

    if scope == "global":
        if testable.sum() > 0:
            _, qvals, _, _ = multipletests(df.loc[testable, "p_value"], method="fdr_bh")
            df.loc[testable, "q_value"] = qvals
    else:
        for component in df["component"].unique():
            idx_to_test = df.index[(df["component"] == component) & testable]
            if len(idx_to_test) > 0:
                _, qvals, _, _ = multipletests(df.loc[idx_to_test, "p_value"], method="fdr_bh")
                df.loc[idx_to_test, "q_value"] = qvals

    df["is_stable"] = df["q_value"] < fdr_threshold
    return df


def peak_to_bed(peak_ids, regex):
    pattern = re.compile(regex)
    rows = []
    unparsed = 0
    for pid in peak_ids:
        m = pattern.match(str(pid))
        if not m:
            unparsed += 1
            continue
        gd = m.groupdict()
        rows.append({"chrom": gd["chr"], "start": int(gd["start"]), "end": int(gd["end"]), "peak": pid})
    if unparsed > 0:
        print(f"  WARNING: {unparsed} peak ID(s) could not be parsed with --peak_coord_regex and were dropped.")
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    peak_weights = pd.read_csv(args.peak_weights_per_fold)
    n_components = peak_weights["component"].nunique()
    n_peaks = peak_weights["peak"].nunique()
    n_folds = peak_weights["held_out_fold"].nunique()
    print(f"Input: {n_folds} folds x {n_components} components x {n_peaks} peaks")
    if n_folds > 1:
        print("NOTE: leave-one-fold-out folds overlap heavily in training data, so per-fold "
              "peak weights are correlated -- the t-test below is anti-conservative. "
              "n_folds_positive/frac_folds_positive is reported as a non-parametric cross-check.")

    print(f"\nRunning one-sided (alternative='{args.alternative}') t-tests...")
    test_results = one_sided_stability_test(peak_weights, args.alternative, args.min_folds)

    n_untestable = test_results["p_value"].isna().sum()
    if n_untestable > 0:
        print(f"  {n_untestable} (component, peak) pair(s) had < --min_folds ({args.min_folds}) "
              f"folds and were left untested (p_value = NaN, is_stable = False).")

    print(f"\nApplying BH correction (scope={args.correction_scope}, threshold={args.fdr_threshold})...")
    test_results = apply_bh_correction(test_results, args.correction_scope, args.fdr_threshold)

    all_path = os.path.join(args.out_dir, "stable_peaks_all.csv")
    test_results.to_csv(all_path, index=False)
    print(f"Written: {all_path}")

    for component, sub in test_results.groupby("component"):
        stable = sub[sub["is_stable"]].sort_values("q_value")
        print(f"  {component}: {len(stable)} / {len(sub)} peaks stable at q < {args.fdr_threshold}")

        comp_path = os.path.join(args.out_dir, f"stable_peaks_{component}.csv")
        stable.to_csv(comp_path, index=False)

        bed_df = peak_to_bed(stable["peak"], args.peak_coord_regex)
        bed_path = os.path.join(args.out_dir, f"stable_peaks_{component}.bed")
        if len(bed_df) > 0:
            bed_df[["chrom", "start", "end", "peak"]].to_csv(bed_path, sep="\t", header=False, index=False)
        else:
            open(bed_path, "w").close()

    print("\nDone.")


if __name__ == "__main__":
    main()
