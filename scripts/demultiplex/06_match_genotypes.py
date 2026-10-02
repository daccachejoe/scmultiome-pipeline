#!/usr/bin/env python3
"""Demultiplex step 6: match souporcell clusters across libraries and suggest patient IDs.

Branch: demultiplex (route 00b). Once, after every library's souporcell run.
Libraries that pool the same people (lesional/non-lesional pairs, longitudinal
timepoints) get independent souporcell cluster labels, so patient identity has
to be matched by genotype. For every pair of libraries, each cluster pair's
alt-allele fractions (AO / (AO + RO), sites with >= MIN_DEPTH reads) at variants
in both cluster_genotypes.vcf files are correlated. Same person r ~0.88-0.92;
different people ~0.3-0.56 (prelim-long-data).

Clusters are grouped into people by linking every cluster pair with r >= min_r
that is also each cluster's best match in the other library (reciprocal best
hits), then taking connected components. Each group gets a placeholder ID
(P01, P02, ...), written as a suggested configs/donor_map.csv. Replace the
placeholders with real patient IDs. A library whose clusters don't match anything
gets its own IDs. Single-donor libraries in the samplesheet get a "*" row.

Inputs:  configs/demultiplexing_paths.csv, configs/samplesheet.csv,
         output/demultiplex/<sample>/souporcell/cluster_genotypes.vcf
         (or next to a user-supplied demux_path)
Outputs: <out_dir>/genotype-comparison.tsv  all cluster pairs, with correlations
         <out_dir>/donor-map-suggested.csv   sampleName,cluster,donor,n_singlets,best_r
USAGE:   06_match_genotypes.py --demux configs/demultiplexing_paths.csv \\
             --samplesheet configs/samplesheet.csv --out_dir output/demultiplex [--min_r 0.7]
"""
import argparse
import csv
import itertools
import math
import os

MIN_DEPTH = 5  # reads per cluster at a site to use its allele fraction

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--demux", required=True)
p.add_argument("--samplesheet", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--min_r", type=float, default=0.7)
args = p.parse_args()


def souporcell_dir(row):
    """Where a library's souporcell outputs are: the pipeline's own run, else next to demux_path."""
    own = os.path.join(args.out_dir, row["sampleName"], "souporcell")
    if os.path.isfile(os.path.join(own, "clusters.tsv")):
        return own
    if row.get("demux_path"):
        return os.path.dirname(row["demux_path"])
    raise SystemExit(f"{row['sampleName']}: no souporcell run in {own} and no demux_path")


def read_genotypes(path):
    sites, clusters = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                clusters = fields[9:]
                continue
            fmt = fields[8].split(":")
            ao_i, ro_i = fmt.index("AO"), fmt.index("RO")
            afs = []
            for sample_field in fields[9:]:
                vals = sample_field.split(":")
                try:
                    ao, ro = int(vals[ao_i]), int(vals[ro_i])
                except (ValueError, IndexError):
                    afs.append(None)
                    continue
                afs.append(ao / (ao + ro) if ao + ro >= MIN_DEPTH else None)
            sites[(fields[0], fields[1], fields[3], fields[4])] = afs
    return clusters, sites


def singlet_counts(path):
    counts = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] == "singlet":
            counts[r["assignment"]] = counts.get(r["assignment"], 0) + 1
    return counts


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")


demux = list(csv.DictReader(open(args.demux)))
libs, data, counts, missing_vcf = [], {}, {}, []
for row in demux:
    d = souporcell_dir(row)
    vcf = os.path.join(d, "cluster_genotypes.vcf")
    counts[row["sampleName"]] = singlet_counts(os.path.join(d, "clusters.tsv"))
    if not os.path.isfile(vcf):
        missing_vcf.append(row["sampleName"])
        continue
    libs.append(row["sampleName"])
    data[row["sampleName"]] = read_genotypes(vcf)
if missing_vcf:
    print(f"no cluster_genotypes.vcf (can't be matched, get their own IDs): {', '.join(missing_vcf)}")

rows = []
for a, b in itertools.combinations(libs, 2):
    ca, sa = data[a]
    cb, sb = data[b]
    shared = sa.keys() & sb.keys()
    print(f"{a} vs {b}: {len(shared)} shared variants")
    for i, x in enumerate(ca):
        for j, y in enumerate(cb):
            pairs = [(sa[s][i], sb[s][j]) for s in shared if sa[s][i] is not None and sb[s][j] is not None]
            rows.append((a, x, b, y, len(pairs), pearson([q[0] for q in pairs], [q[1] for q in pairs])))

os.makedirs(args.out_dir, exist_ok=True)
with open(os.path.join(args.out_dir, "genotype-comparison.tsv"), "w") as out:
    out.write("run_a\tcluster_a\trun_b\tcluster_b\tn_sites\tcorrelation\n")
    for row in rows:
        out.write("\t".join(map(str, row[:5])) + f"\t{row[5]:.3f}\n")

# reciprocal best hits with r >= min_r, then connected components
r_of = {}
for a, x, b, y, _, r in rows:
    if not math.isnan(r):
        r_of[((a, x), (b, y))] = r_of[((b, y), (a, x))] = r
best = {}  # (lib, cluster, other lib) -> (other cluster, r)
for ((a, x), (b, y)), r in r_of.items():
    k = (a, x, b)
    if k not in best or r > best[k][1]:
        best[k] = (y, r)
parent = {}


def find(u):
    parent.setdefault(u, u)
    while parent[u] != u:
        parent[u] = parent[parent[u]]
        u = parent[u]
    return u


nodes = [(s, c) for s in counts for c in sorted(counts[s])]
best_r = {n: float("nan") for n in nodes}
for (a, x, b), (y, r) in best.items():
    best_r[(a, x)] = max(r, best_r[(a, x)]) if not math.isnan(best_r[(a, x)]) else r
    if r >= args.min_r and best.get((b, y, a), (None,))[0] == x:
        parent[find((a, x))] = find((b, y))
groups = {}
for n in nodes:
    groups.setdefault(find(n), []).append(n)
ids = {}
for k, members in enumerate(sorted(groups.values(), key=lambda m: sorted(m)[0]), start=1):
    for n in members:
        ids[n] = f"P{k:02d}"

demux_samples = {r["sampleName"] for r in demux}
single = [r["sampleName"] for r in csv.DictReader(open(args.samplesheet)) if r["sampleName"] not in demux_samples]
out_path = os.path.join(args.out_dir, "donor-map-suggested.csv")
with open(out_path, "w") as out:
    out.write("sampleName,cluster,donor,n_singlets,best_r\n")
    for s, c in nodes:
        out.write(f"{s},{c},{ids[(s, c)]},{counts[s][c]},{best_r[(s, c)]:.3f}\n")
    for s in single:
        out.write(f"{s},*,{s},,\n")
for members in sorted(groups.values(), key=lambda m: sorted(m)[0]):
    print(f"{ids[members[0]]}: " + ", ".join(f"{s} cluster {c}" for s, c in sorted(members)))
print(f"wrote {out_path} -- replace the placeholder IDs with patient IDs and save as configs/donor_map.csv")
