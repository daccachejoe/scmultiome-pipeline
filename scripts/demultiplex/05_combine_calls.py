#!/usr/bin/env python3
"""Demultiplex step 5: combine souporcell and pooled-genotype calls into one per-cell table.

Branch: demultiplex (route 00b). Per library.
Priority per barcode:
  souporcell singlet                       -> singlet, souporcell's cluster   (source souporcell)
  souporcell doublet                       -> doublet                          (souporcell)
  souporcell unassigned, pooled doublet    -> doublet                          (pooled)
  souporcell unassigned, pooled donor call -> singlet, that cluster          (pooled_rescue)
      unless the donor is untrusted (below)
  otherwise                                -> unassigned                       (none)
souporcell leaves 13-38% of barcodes unassigned; they pass QC and are real
low-RNA nuclei, which the pooled ATAC + GEX calls (step 4) mostly recover.

Thin-donor guard: a donor whose held-out precision (step 4's
heldout_per_donor_<tag>.tsv, combined calls) is below --min_donor_precision
gets no rescue calls. Its pooled genotype is too thin to be recognized, so its
rescue calls may be another donor's cells. This replaces a hand-kept exception:
in healthy-human cntrl.1, 33 of 34 held-out donor-2 cells (187 cells) were
called donor 0, whose precision was 97.7%. Default 0.98: across the 20 donors
of prelim-long-data every donor that was kept had >= 98.8% (lowest PSO v1
donor 0), so 0.98 reproduces the hand-made exception exactly (it also drops
cntrl.1 donor 2 itself, 1/2 correct, which had no rescue calls). The margin is
narrow; check heldout_per_donor_*.tsv when a new library sits near it.

The output is in souporcell clusters.tsv format (barcode, status, assignment),
which stage 01 doublets reads for genotype doublets, plus source and
souporcell_status columns that stage 02's donors stage writes into the object.
Cluster labels are souporcell's (0, 1, ...); configs/donor_map.csv maps them
to patient IDs.

Inputs:  souporcell clusters.tsv; optionally the step-4 calls and per-donor tables
Outputs: <out> (combined_clusters.tsv)
USAGE:   05_combine_calls.py --clusters souporcell/clusters.tsv --out combined_clusters.tsv \\
             [--calls pooled/calls_soft_all_calibrated_fragments.tsv \\
              --per_donor pooled/heldout_per_donor_fragments.tsv --min_donor_precision 0.98]
"""
import argparse
import csv
import math

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--clusters", required=True)
p.add_argument("--out", required=True)
p.add_argument("--calls", help="step-4 calls_soft_all_calibrated_<tag>.tsv (omit for souporcell only)")
p.add_argument("--per_donor", help="step-4 heldout_per_donor_<tag>.tsv (required with --calls)")
p.add_argument("--min_donor_precision", type=float, default=0.98)
args = p.parse_args()
if args.calls and not args.per_donor:
    raise SystemExit("--per_donor is required with --calls (it drives the thin-donor guard)")

soup = list(csv.DictReader(open(args.clusters), delimiter="\t"))
labels = sorted({r["assignment"] for r in soup if r["status"] == "singlet"})
pooled, untrusted = {}, set()
if args.calls:
    pooled = {r["barcode"]: r["combined_call"] for r in csv.DictReader(open(args.calls), delimiter="\t")}
    for r in csv.DictReader(open(args.per_donor), delimiter="\t"):
        prec = float(r["precision"])
        if not math.isnan(prec) and prec < args.min_donor_precision:
            untrusted.add(r["donor"])
            print(f"donor {r['donor']}: held-out precision {100 * prec:.1f}% < {100 * args.min_donor_precision:g}% "
                  f"({r['n_correct']}/{r['n_called_as']}) -- its rescue calls are not used")
    missing = [r["barcode"] for r in soup if r["barcode"] not in pooled]
    if missing:
        raise SystemExit(f"{len(missing)} souporcell barcodes are missing from {args.calls} "
                         f"(e.g. {missing[0]}) -- were both made from the same Cell Ranger barcodes?")

n = {}
with open(args.out, "w") as out:
    out.write("barcode\tstatus\tassignment\tsource\tsouporcell_status\n")
    for r in soup:
        call = pooled.get(r["barcode"], "unassigned")
        if r["status"] == "singlet":
            status, assignment, source = "singlet", r["assignment"], "souporcell"
        elif r["status"] == "doublet":
            status, assignment, source = "doublet", "doublet", "souporcell"
        elif call == "doublet":
            status, assignment, source = "doublet", "doublet", "pooled"
        elif call in labels and call not in untrusted:
            status, assignment, source = "singlet", call, "pooled_rescue"
        else:
            status, assignment, source = "unassigned", "unassigned", "none"
        out.write(f"{r['barcode']}\t{status}\t{assignment}\t{source}\t{r['status']}\n")
        n[(status, source)] = n.get((status, source), 0) + 1
print(f"wrote {args.out}: " + ", ".join(f"{s}/{src} {k}" for (s, src), k in sorted(n.items())))
