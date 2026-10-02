#!/usr/bin/env python3
"""Demultiplex step 1: split a library's BAMs into per-donor, per-fold pools.

Branch: demultiplex (route 00b, run before stage 01). Per library.
Uses souporcell's GEX singlet calls to route each read (by its CB tag) to its
donor. Singlets are split into folds A/B by barcode hash (demux_common.fold_of)
so genotypes built from fold A can be tested on held-out fold B in step 4.
Each output read gets a read group whose SM is its donor, which is how
bcftools mpileup pools reads across files in rescue.sh. BAMs are read in place
(the ATAC BAMs are 24-56 GB), with one worker per (modality, chromosome).

Inputs:  souporcell clusters.tsv (barcode/status/assignment), ATAC and/or GEX BAM
Outputs: <out_dir>/singlet_folds.tsv (barcode, fold, donor d<cluster>)
         <out_dir>/split/<fold>_<donor>_<modality>.<chrom>.bam
USAGE:
  01_split_by_donor.py --clusters clusters.tsv --out_dir DIR --species human \\
      --threads 16 --bam atac=/path/atac_possorted_bam.bam [--bam gex=/path/gex_possorted_bam.bam]
"""
import argparse
import csv
import multiprocessing as mp
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from demux_common import chromosomes, fold_of  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--clusters", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--species", required=True)
p.add_argument("--threads", type=int, default=1)
p.add_argument("--bam", action="append", required=True, help="modality=path, e.g. atac=/path.bam (repeatable)")
args = p.parse_args()

BAMS = dict(b.split("=", 1) for b in args.bam)
CHROMS = chromosomes(args.species)
for modality, path in BAMS.items():
    if not os.path.exists(path):
        raise SystemExit(f"{modality} BAM not found: {path}")

groups = {}
for r in csv.DictReader(open(args.clusters), delimiter="\t"):
    if r["status"] == "singlet":
        groups[r["barcode"]] = (fold_of(r["barcode"]), "d" + r["assignment"])
if not groups:
    raise SystemExit(f"{args.clusters} has no singlets -- check the souporcell run")

DONORS = sorted({d for _, d in groups.values()})  # d0, d1, ... one per souporcell cluster
os.makedirs(os.path.join(args.out_dir, "split"), exist_ok=True)
with open(os.path.join(args.out_dir, "singlet_folds.tsv"), "w") as out:
    out.write("barcode\tfold\tdonor\n")
    for b, (fold, donor) in groups.items():
        out.write(f"{b}\t{fold}\t{donor}\n")


def split(task):
    modality, chrom = task
    bam_in = pysam.AlignmentFile(BAMS[modality])
    header = bam_in.header.to_dict()
    header.pop("PG", None)
    outs = {}
    for fold in "AB":
        for donor in DONORS:
            rg = f"{fold}_{donor}_{modality}"
            header["RG"] = [{"ID": rg, "SM": donor, "LB": modality}]
            outs[(fold, donor)] = pysam.AlignmentFile(os.path.join(args.out_dir, "split", f"{rg}.{chrom}.bam"),
                                                      "wb", header=header)
    n = 0
    for read in bam_in.fetch(chrom):
        if not read.has_tag("CB"):
            continue
        group = groups.get(read.get_tag("CB"))
        if group is None:
            continue
        read.set_tag("RG", f"{group[0]}_{group[1]}_{modality}")
        outs[group].write(read)
        n += 1
    for out in outs.values():
        out.close()
    return modality, chrom, n


tasks = [(m, c) for m in BAMS for c in CHROMS]
with mp.Pool(args.threads) as pool:
    for modality, chrom, n in pool.imap_unordered(split, tasks):
        print(f"{modality} {chrom}: {n} reads kept", flush=True)
print("donors:", " ".join(DONORS))
