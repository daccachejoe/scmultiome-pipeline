#!/usr/bin/env python3
"""Demultiplex step 3: count ATAC ref/alt alleles per cell at the soft sites, once per fragment.

Branch: demultiplex (route 00b). Per library.
vartrix counts reads, so a SNP inside the overlap of a fragment's two mates is
counted twice; with ~50 bp ATAC reads ~60-70% of fragments overlap, which
inflated ATAC evidence (fitted ATAC weight 0.5 with read counts vs 0.65 with
fragment counts on PSO v0). Here each fragment (read name within a cell) gives
one observation per site: agreeing mates count once, disagreeing mates are
dropped. Duplicates, secondary/supplementary and QC-fail reads, MAPQ < 30 and
base quality < 20 are skipped. GEX alleles are still counted by vartrix --umi
(rescue.sh), which already collapses by UMI.

Inputs:  <run_dir>/barcodes.tsv, <run_dir>/<name>.vcf for each site set, ATAC BAM
Outputs: <run_dir>/fragcount_<name>_atac/{ref,alt}.mtx (sites x cells, vartrix layout)
USAGE:   03_count_fragments.py --run_dir DIR --bam atac_possorted_bam.bam --threads 16 \\
             [--sets soft_foldA soft_all]
"""
import argparse
import multiprocessing as mp
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from demux_common import read_barcodes, read_sites  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--run_dir", required=True)
p.add_argument("--bam", required=True)
p.add_argument("--threads", type=int, default=1)
p.add_argument("--sets", nargs="+", default=["soft_foldA", "soft_all"])
args = p.parse_args()

barcodes = read_barcodes(os.path.join(args.run_dir, "barcodes.tsv"))
cell_idx = {b: i for i, b in enumerate(barcodes)}
bam_path = args.bam


def count_chrom(task):
    chrom, sites = task  # sites: list of (site_index, pos, ref, alt)
    bam = pysam.AlignmentFile(bam_path)
    ref_counts, alt_counts = [], []
    stats = {"agree_pairs": 0, "disagree_pairs": 0, "single": 0}
    for si, pos, ref, alt in sites:
        obs = {}
        for col in bam.pileup(chrom, pos - 1, pos, truncate=True, stepper="all",
                              min_base_quality=20, min_mapping_quality=30,
                              ignore_overlaps=False, ignore_orphans=False, max_depth=1000000):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                read = pr.alignment
                if not read.has_tag("CB"):
                    continue
                ci = cell_idx.get(read.get_tag("CB"))
                if ci is None:
                    continue
                obs.setdefault((ci, read.query_name), []).append(read.query_sequence[pr.query_position])
        per_cell = {}
        for (ci, _), bases in obs.items():
            alleles = set(bases)
            if len(bases) > 1:
                stats["agree_pairs" if len(alleles) == 1 else "disagree_pairs"] += 1
            else:
                stats["single"] += 1
            if len(alleles) != 1:
                continue
            base = alleles.pop()
            r, a = per_cell.get(ci, (0, 0))
            if base == ref:
                per_cell[ci] = (r + 1, a)
            elif base == alt:
                per_cell[ci] = (r, a + 1)
        for ci, (r, a) in per_cell.items():
            if r:
                ref_counts.append((si, ci, r))
            if a:
                alt_counts.append((si, ci, a))
    return ref_counts, alt_counts, stats


def write_mtx(path, entries, n_sites):
    with open(path, "w") as out:
        out.write("%%MatrixMarket matrix coordinate integer general\n")
        out.write(f"{n_sites} {len(barcodes)} {len(entries)}\n")
        for si, ci, v in entries:
            out.write(f"{si + 1} {ci + 1} {v}\n")


for set_name in args.sets:
    sites = read_sites(os.path.join(args.run_dir, set_name + ".vcf"))
    by_chrom = {}
    for i, (chrom, pos, ref, alt) in enumerate(sites):
        by_chrom.setdefault(chrom, []).append((i, pos, ref, alt))
    ref_all, alt_all, stats = [], [], {"agree_pairs": 0, "disagree_pairs": 0, "single": 0}
    with mp.Pool(args.threads) as pool:
        for r, a, s in pool.imap_unordered(count_chrom, by_chrom.items()):
            ref_all += r
            alt_all += a
            for k in stats:
                stats[k] += s[k]
    out_dir = os.path.join(args.run_dir, f"fragcount_{set_name}_atac")
    os.makedirs(out_dir, exist_ok=True)
    write_mtx(os.path.join(out_dir, "ref.mtx"), sorted(ref_all), len(sites))
    write_mtx(os.path.join(out_dir, "alt.mtx"), sorted(alt_all), len(sites))
    frags = sum(stats.values())
    print(f"{set_name}: {len(sites)} sites, {frags} fragment observations; "
          f"both mates cover site {100 * (stats['agree_pairs'] + stats['disagree_pairs']) / max(1, frags):.1f}% "
          f"(mates disagree {stats['disagree_pairs']}, dropped)", flush=True)
