#!/usr/bin/env python3
"""Demultiplex step 4: calibrate ATAC allele evidence on held-out singlets, then call every cell.

Branch: demultiplex (route 00b). Per library.
Model: per site, one alt fraction per donor (from the soft genotypes of step 2)
plus a doublet class, the mean likelihood over all donor pairs (for 2 donors,
the single 50/50 mix). Every model is mixed with `ambient` of the pool average
(ambient DNA/RNA). Priors: `doublet_prior` doublet, the rest split evenly over
donors. A cell is called when its posterior >= min_posterior.

Held-out design: fold-B souporcell singlets (split in step 1) are called with
fold-A genotypes and scored against souporcell, so accuracy isn't circular.
ATAC posteriors were overconfident (v1 calls rated >= 0.999 were 98.3% correct
at weight 1), so a weight on ATAC evidence is fit by held-out log loss on
ATAC-only posteriors (GEX overlaps the truth labels), then applied:
combined = w * ATAC + GEX, on the all-singlet genotypes.
Known weakness: log loss is dominated by a few confident disagreements that are
often RNA's errors, pushing w down (cntrl.2 fitted 0.05 though w = 1 was
calibrated). min_rna_margin > 0 restricts truth to confident souporcell calls
(best minus second-best cluster log-likelihood >= margin); cntrl.2 used 20.

Per-donor held-out precision/recall of the combined calls is written for the
thin-donor guard in step 5: a donor with too few cells for a good pooled
genotype attracts another donor's cells (cntrl.1: 33 of 34 held-out donor-2
cells were called donor 0).

Inputs:  <run_dir>/{barcodes.tsv, singlet_folds.tsv, soft_{foldA,all}.af.tsv,
         {fragcount|vartrix}_soft_{set}_atac/, vartrix_soft_{set}_gex/ (absent = ATAC-only)},
         souporcell clusters.tsv
Outputs (tag = <atac_counts>[_rnamargin<m>]) in <run_dir>:
         calls_soft_all_calibrated_<tag>.tsv   final per-cell calls (all-singlet genotypes)
         calls_soft_foldA_calibrated_<tag>.tsv per-cell calls from fold-A genotypes
         heldout_calibrated_<tag>.tsv          held-out summary at the fitted weight and w = 1
         heldout_per_donor_<tag>.tsv           held-out precision/recall per donor (combined)
USAGE:   04_calibrate.py --run_dir DIR --clusters souporcell/clusters.tsv \\
             [--atac_counts fragments] [--min_rna_margin 0] [--min_posterior 0.95] \\
             [--ambient 0.10] [--doublet_prior 0.05]
"""
import argparse
import csv
import itertools
import os

import numpy as np
from scipy.io import mmread
from scipy.special import logsumexp

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run_dir", required=True)
p.add_argument("--clusters", required=True, help="souporcell clusters.tsv (truth for held-out singlets)")
p.add_argument("--atac_counts", choices=["reads", "fragments"], default="fragments")
p.add_argument("--min_rna_margin", default="0", help="kept as given in output names, e.g. 20 -> _rnamargin20")
p.add_argument("--min_posterior", type=float, default=0.95)
p.add_argument("--ambient", type=float, default=0.10)
p.add_argument("--doublet_prior", type=float, default=0.05)
args = p.parse_args()

ERR, AMBIENT, MIN_POST, DOUBLET_PRIOR = 0.01, args.ambient, args.min_posterior, args.doublet_prior
WEIGHTS = np.round(np.arange(0.05, 1.501, 0.05), 2)

run_dir = args.run_dir
counts = args.atac_counts
min_margin = float(args.min_rna_margin)
tag = counts + (f"_rnamargin{args.min_rna_margin}" if min_margin > 0 else "")
atac_prefix = {"reads": "vartrix_soft", "fragments": "fragcount_soft"}[counts]
barcodes = [l.strip() for l in open(os.path.join(run_dir, "barcodes.tsv"))]
sample = os.path.basename(os.path.normpath(run_dir))

# donors = the af table's f_<donor> columns (d0, d1, ... from souporcell cluster 0, 1, ...)
af_header = open(os.path.join(run_dir, "soft_foldA.af.tsv")).readline().rstrip("\n").split("\t")
DONORS = [c[2:] for c in af_header if c.startswith("f_")]
LABELS = [d[1:] for d in DONORS]  # souporcell cluster label, e.g. "0"
NAMES = np.array(LABELS + ["doublet"])
PAIRS = list(itertools.combinations(range(len(DONORS)), 2))
PRIOR = np.log(np.array([(1 - DOUBLET_PRIOR) / len(DONORS)] * len(DONORS) + [DOUBLET_PRIOR]))


def models(set_name):
    """Per-site alt fractions: one per donor, then one per donor pair (50/50 doublet)."""
    header, *rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(run_dir, f"soft_{set_name}.af.tsv"))]
    assert [c[2:] for c in header if c.startswith("f_")] == DONORS, header
    cols = [header.index(f"f_{d}") for d in DONORS]
    p = [np.clip(np.array([float(r[c]) for r in rows]), ERR, 1 - ERR) for c in cols]
    pbar = sum(p) / len(p)
    mixes = p + [(p[i] + p[j]) / 2 for i, j in PAIRS]
    return [np.clip((1 - AMBIENT) * m + AMBIENT * pbar, 1e-4, 1 - 1e-4) for m in mixes]


def collapse_doublets(ll):
    """Donor columns, then one doublet column: mean likelihood over the donor-pair models."""
    k = len(DONORS)
    return np.column_stack([ll[:, :k], logsumexp(ll[:, k:], axis=1) - np.log(len(PAIRS))])


def load(set_name):
    atac_dir = f"{atac_prefix}_{set_name}_atac"
    gex_dir = f"vartrix_soft_{set_name}_gex"
    out = {}
    for key, d in (("atac", atac_dir), ("gex", gex_dir)):
        ms = models(set_name)
        if key == "gex" and not os.path.isdir(os.path.join(run_dir, d)):
            # ATAC-only library (no GEX BAM): GEX contributes no evidence
            out[key] = (np.zeros((len(barcodes), len(DONORS) + 1)), np.zeros(len(barcodes)))
            continue
        ref = mmread(os.path.join(run_dir, d, "ref.mtx")).tocsc()
        alt = mmread(os.path.join(run_dir, d, "alt.mtx")).tocsc()
        assert ref.shape == (len(ms[0]), len(barcodes)), (d, ref.shape)
        ll = np.vstack([ref.T @ np.log(1 - m) + alt.T @ np.log(m) for m in ms]).T
        out[key] = (collapse_doublets(ll), np.asarray((ref + alt).sum(axis=0)).ravel())
    return out


def posterior(ll):
    lp = ll + PRIOR
    lp = lp - lp.max(axis=1, keepdims=True)
    p = np.exp(lp)
    return p / p.sum(axis=1, keepdims=True)


def calls(post, n):
    c = np.where(post.max(axis=1) >= MIN_POST, NAMES[post.argmax(axis=1)], "unassigned")
    return np.where(n == 0, "unassigned", c)


fold = {r["barcode"]: r["fold"] for r in csv.DictReader(open(os.path.join(run_dir, "singlet_folds.tsv")), delimiter="\t")}
soup_rows = {r["barcode"]: r for r in csv.DictReader(open(args.clusters), delimiter="\t") if r["status"] == "singlet"}
soup = {b: r["assignment"] for b, r in soup_rows.items()}


def rna_margin(r):
    ll = sorted((float(v) for k, v in r.items() if k.startswith("cluster") and k[7:].isdigit()), reverse=True)
    return ll[0] - ll[1] if len(ll) > 1 else float("inf")


idx = np.array([i for i, b in enumerate(barcodes) if fold.get(b) == "B" and rna_margin(soup_rows[barcodes[i]]) >= min_margin])
truth = np.array([LABELS.index(soup[barcodes[i]]) for i in idx])  # column index of the true donor
truth_label = np.array(LABELS)[truth] if len(idx) else np.array([], dtype=str)

fA = load("foldA")
atac_ll, atac_n = fA["atac"]
has = atac_n[idx] > 0
print(f"{sample} [{tag}]: {len(idx)} fold-B singlets (RNA margin >= {min_margin:g}), {has.sum()} with ATAC evidence")
best = None
for w in WEIGHTS:
    post = posterior(w * atac_ll[idx][has])
    nll = -np.mean(np.log(np.clip(post[np.arange(has.sum()), truth[has]], 1e-12, 1)))
    if best is None or nll < best[1]:
        best = (w, nll)
w = best[0]
print(f"fitted ATAC weight {w} (held-out log loss {best[1]:.4f}; weight 1.0 gives "
      f"{-np.mean(np.log(np.clip(posterior(atac_ll[idx][has])[np.arange(has.sum()), truth[has]], 1e-12, 1))):.4f})")


def report(label, post, n, rows):
    c = calls(post, n)[rows]
    t = truth_label
    called = np.isin(c, LABELS)
    acc = (c[called] == t[called]).mean() if called.any() else float("nan")
    print(f"  {label:28} called {100 * called.mean():5.1f}%  accuracy {100 * acc:6.2f}%")
    pm = post.max(axis=1)[rows]
    for lo, hi in [(0.95, 0.99), (0.99, 0.999), (0.999, 1.01)]:
        s = called & (pm >= lo) & (pm < hi)
        if s.sum():
            print(f"      post {lo}-{min(hi, 1)}: n={s.sum():5}  correct {100 * (c[s] == t[s]).mean():5.1f}%")


gex_ll, gex_n = fA["gex"]
for label, ww in (("ATAC only, weight 1.0", 1.0), (f"ATAC only, weight {w}", w)):
    report(label, posterior(ww * atac_ll), atac_n, idx)
for label, ww in (("ATAC + GEX, weight 1.0", 1.0), (f"ATAC + GEX, weight {w}", w)):
    report(label, posterior(ww * atac_ll + gex_ll), atac_n + gex_n, idx)

with open(os.path.join(run_dir, f"heldout_calibrated_{tag}.tsv"), "w") as out:
    out.write("reads_used\tatac_weight\tn_cells\tn_called\tn_correct\n")
    for label, ll_, n_, ww in (("atac", w * atac_ll, atac_n, w), ("combined", w * atac_ll + gex_ll, atac_n + gex_n, w),
                               ("atac", atac_ll, atac_n, 1.0), ("combined", atac_ll + gex_ll, atac_n + gex_n, 1.0)):
        c = calls(posterior(ll_), n_)[idx]
        called = np.isin(c, LABELS)
        out.write(f"{label}\t{ww}\t{len(idx)}\t{called.sum()}\t{(c[called] == truth_label[called]).sum()}\n")

# per-donor held-out precision/recall of the combined calls at the fitted weight,
# for the thin-donor guard (05_combine_calls.py)
c_heldout = calls(posterior(w * atac_ll + gex_ll), atac_n + gex_n)[idx]
with open(os.path.join(run_dir, f"heldout_per_donor_{tag}.tsv"), "w") as out:
    out.write("donor\tn_truth\tn_called_as\tn_correct\tprecision\trecall\n")
    for lab in LABELS:
        n_truth = int((truth_label == lab).sum())
        n_called = int((c_heldout == lab).sum())
        n_correct = int(((c_heldout == lab) & (truth_label == lab)).sum())
        prec = n_correct / n_called if n_called else float("nan")
        rec = n_correct / n_truth if n_truth else float("nan")
        out.write(f"{lab}\t{n_truth}\t{n_called}\t{n_correct}\t{prec:.4f}\t{rec:.4f}\n")
        print(f"  donor {lab}: held-out precision {100 * prec:.1f}% ({n_correct}/{n_called}), "
              f"recall {100 * rec:.1f}% ({n_correct}/{n_truth})")

# per-cell calls from the fold-A genotypes (held-out cells are those with fold B in
# singlet_folds.tsv), at weight 1 and the fitted weight
with open(os.path.join(run_dir, f"calls_soft_foldA_calibrated_{tag}.tsv"), "w") as out:
    out.write("barcode\tatac_reads\tatac_call_w1\tatac_post_w1\tatac_call\tatac_post\tcombined_call\tcombined_post\n")
    pw1, pw = posterior(atac_ll), posterior(w * atac_ll)
    pc = posterior(w * atac_ll + gex_ll)
    cw1, cw, cc = calls(pw1, atac_n), calls(pw, atac_n), calls(pc, atac_n + gex_n)
    for i, b in enumerate(barcodes):
        out.write(f"{b}\t{atac_n[i]:g}\t{cw1[i]}\t{pw1[i].max():.4f}\t{cw[i]}\t{pw[i].max():.4f}\t{cc[i]}\t{pc[i].max():.4f}\n")

# final calls on all-singlet genotypes with the fitted weight
al = load("all")
post = {"atac": posterior(w * al["atac"][0]), "gex": posterior(al["gex"][0]),
        "combined": posterior(w * al["atac"][0] + al["gex"][0])}
n = {"atac": al["atac"][1], "gex": al["gex"][1], "combined": al["atac"][1] + al["gex"][1]}
out_path = os.path.join(run_dir, f"calls_soft_all_calibrated_{tag}.tsv")
with open(out_path, "w") as out:
    out.write("barcode\tatac_reads\tgex_reads\tatac_weight\t" +
              "\t".join(f"{m}_call\t{m}_post" for m in ("atac", "gex", "combined")) + "\n")
    cs = {m: calls(post[m], n[m]) for m in post}
    for i, b in enumerate(barcodes):
        out.write(f"{b}\t{n['atac'][i]:g}\t{n['gex'][i]:g}\t{w}\t" +
                  "\t".join(f"{cs[m][i]}\t{post[m][i].max():.4f}" for m in ("atac", "gex", "combined")) + "\n")
print(f"wrote {out_path}")
