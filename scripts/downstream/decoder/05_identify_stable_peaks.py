#!/usr/bin/env python3
"""
Decoder branch, step 5: identify peaks with a consistently non-zero weight
across CV folds, in either direction.

For each regulatory direction and each peak, runs a t-test on the per-fold
peak weights from step 4 (two-sided by default -- H0: mean weight == 0,
H1: mean weight != 0), then applies Benjamini-Hochberg FDR correction
(--correction_scope per_component by default: each regulatory direction is
its own hypothesis family; pass --correction_scope global to pool all
direction x peak p-values into one correction instead).

Two-sided is the default because a one-sided test only ever looks at one
direction -- the original "greater" default here silently meant every
peak that *decreases* along a direction was untestable by construction (a
decreasing peak's mean is negative, so it can never clear a "greater than
zero" bar no matter how consistent the decrease is). Peaks that open AND
peaks that close are both real biology; --alternative greater/less remain
available if you deliberately want a one-directional test.

Output is split by sign regardless of --alternative (stable_peaks_
<component>_up.csv/.bed for mean_weight > 0, _down for mean_weight < 0),
not just one mixed file per component. This matters beyond bookkeeping:
step 6 (motif enrichment) treats each stable_peaks_*.bed as one HOMER
foreground set, and "which motifs are enriched in a set that mixes
opening and closing peaks" is a muddled question with a muddled answer --
splitting keeps "what's enriched in this direction's opening peaks" and
"...closing peaks" as two clean, separately answerable questions. Under
--alternative greater/less this split is a no-op (the unused direction's
file is just empty), so the file-naming convention stays consistent
regardless of which --alternative a given run used.

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
    p.add_argument("--alternative", choices=["greater", "less", "two-sided"], default="two-sided")
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
        # scipy's ttest_1samp only reports the two-sided p-value, so the
        # one-sided ("greater"/"less") p-value has to be derived from it by
        # hand: halve it when the sample mean already points the direction
        # being tested (t_stat has the expected sign), or take the
        # complement (1 - two_sided/2) when it points the other way -- a
        # peak whose mean weight is negative can never have a small
        # "greater" p-value, no matter how tight its distribution is,
        # which is exactly what the `else` branch of each case enforces.
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
    """
    scope="per_component" (default) runs a separate BH correction for each
    regulatory direction's peaks -- appropriate when "5% FDR among
    beta_dis-associated peaks" and "5% FDR among beta_loc-associated peaks"
    are treated as two independent scientific questions, which is the
    normal framing for this kind of multi-direction decoder. scope="global"
    pools every (component, peak) p-value into one correction instead --
    more conservative, and appropriate only if you want a single "5% FDR
    across all reported peaks, from whichever direction" guarantee.
    """
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

    print(f"\nRunning t-tests (alternative='{args.alternative}')...")
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
        n_up = int((stable["mean_weight"] > 0).sum())
        n_down = int((stable["mean_weight"] < 0).sum())
        print(f"  {component}: {len(stable)} / {len(sub)} peaks stable at q < {args.fdr_threshold} "
              f"({n_up} up, {n_down} down)")

        # Split by sign regardless of --alternative -- see the module
        # docstring for why a mixed up+down foreground set would corrupt
        # step 6's motif enrichment. Under a one-sided --alternative this
        # is a no-op (the unused direction's files are simply empty), so
        # the naming convention is the same either way.
        for direction, direction_df in [("up", stable[stable["mean_weight"] > 0]),
                                         ("down", stable[stable["mean_weight"] < 0])]:
            comp_path = os.path.join(args.out_dir, f"stable_peaks_{component}_{direction}.csv")
            direction_df.to_csv(comp_path, index=False)

            bed_df = peak_to_bed(direction_df["peak"], args.peak_coord_regex)
            bed_path = os.path.join(args.out_dir, f"stable_peaks_{component}_{direction}.bed")
            if len(bed_df) > 0:
                bed_df[["chrom", "start", "end", "peak"]].to_csv(bed_path, sep="\t", header=False, index=False)
            else:
                open(bed_path, "w").close()

    print("\nDone.")


if __name__ == "__main__":
    main()
