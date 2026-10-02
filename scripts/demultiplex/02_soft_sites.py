#!/usr/bin/env python3
"""Demultiplex step 2: pick donor-informative SNPs from pooled allele counts.

Branch: demultiplex (route 00b). Per library, run on the bcftools joint calls
over the per-donor pools from step 1 (fold-A pools, and all singlets).
"Soft" genotypes: instead of hard genotype calls, each donor's alt-allele
fraction is estimated from its pooled reads, shrunk toward 0.5 at low depth,
f = (alt + 1) / (ref + alt + 2). A SNP is kept when at least two pools have
depth >= min_depth and the shrunk fractions of those pools differ by
>= min_diff. Hard genotype filters were tried and rejected: on PSO v0 they kept
~2k SNPs and called 50% of held-out singlets, against ~30k SNPs and 86% at
99.6% accuracy for soft. Any number of donors (the VCF sample columns d0, d1, ...).

Inputs:  <calls.vcf.gz> (bcftools call -m -v with FORMAT/AD)
Outputs: <prefix>.vcf (sites for vartrix / step 3), <prefix>.af.tsv (per-donor
         ref/alt counts, then f_<donor>)
USAGE:   02_soft_sites.py calls/all.vcf.gz soft_all [--min_depth 2] [--min_diff 0.25]
"""
import argparse
import gzip

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("vcf")
p.add_argument("prefix")
p.add_argument("--min_depth", type=int, default=2)
p.add_argument("--min_diff", type=float, default=0.25)
p.add_argument("--min_qual", type=float, default=20)
args = p.parse_args()

kept = total = 0
with gzip.open(args.vcf, "rt") as fh, open(args.prefix + ".vcf", "w") as vcf, open(args.prefix + ".af.tsv", "w") as af:
    for line in fh:
        if line.startswith("##"):
            vcf.write(line)
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = f[9:]
            vcf.write("\t".join(f[:8]) + "\n")
            af.write("\t".join(["chrom", "pos"] + [f"{d}_{x}" for d in samples for x in ("ref", "alt")]
                               + [f"f_{d}" for d in samples]) + "\n")
            continue
        total += 1
        if "INDEL" in f[7] or len(f[3]) != 1 or len(f[4]) != 1 or float(f[5]) < args.min_qual:
            continue
        fmt = f[8].split(":")
        ad = {}
        for name, val in zip(samples, f[9:]):
            d = dict(zip(fmt, val.split(":")))
            r, a = (int(x) for x in d["AD"].split(",")[:2])
            ad[name] = (r, a)
        covered = [d for d in samples if sum(ad[d]) >= args.min_depth]
        if len(covered) < 2:
            continue
        fr = {k: (a + 1) / (r + a + 2) for k, (r, a) in ad.items()}
        if max(fr[d] for d in covered) - min(fr[d] for d in covered) < args.min_diff:
            continue
        kept += 1
        vcf.write("\t".join(f[:8]) + "\n")
        af.write("\t".join([f[0], f[1]] + [str(x) for d in samples for x in ad[d]]
                           + [f"{fr[d]:.4f}" for d in samples]) + "\n")

print(f"{args.vcf}: {total} records, {kept} soft-informative SNPs across {len(samples)} donor pools "
      f"(depth >= {args.min_depth} in >= 2 pools, max diff >= {args.min_diff})")
