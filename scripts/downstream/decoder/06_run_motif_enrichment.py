#!/usr/bin/env python3
"""
Decoder branch, step 6: motif enrichment on stable peaks per regulatory
direction.

Wraps HOMER's findMotifsGenome.pl (config key decoder_homer_script) --
chosen because it's the most commonly available module in genomics HPC
cores and its -bg background-peak-matching is exactly what a "stable
peaks vs. all peaks tested" enrichment design calls for. This is the one
external-tool dependency in the whole decoder branch, isolated to this
file; swap the subprocess call in run_homer() for another tool without
touching steps 1-5.

For each stable_peaks_<component>.bed from step 5, runs HOMER against a
background of ALL peaks that were tested for stability (step 2's
atac_peaks_used.csv, used_in_lsi == True), NOT the whole genome --
answers "which motifs are enriched relative to other peaks we could
plausibly have called stable", not relative to random genomic background.

PREREQUISITES: HOMER with the target genome installed
(perl configureHomer.pl -install hg38), findMotifsGenome.pl on $PATH or
pointed to via --homer_script / decoder_homer_script.

USAGE:
    python 06_run_motif_enrichment.py \
        --stable_peaks_dir output/decoder/stable_peaks \
        --peaks_used output/decoder/features/atac_peaks_used.csv \
        --genome hg38 --out_dir output/decoder/motif_enrichment
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stable_peaks_dir", required=True, help="Directory with stable_peaks_<component>.bed from step 5")
    p.add_argument("--peaks_used", required=True, help="Path to atac_peaks_used.csv from step 2")
    p.add_argument("--peak_coord_regex", default=r"^(?P<chr>.+)[:\-_](?P<start>\d+)[:\-_](?P<end>\d+)$",
                    help="Keep in sync with step 5's --peak_coord_regex")
    p.add_argument("--genome", required=True, help="HOMER genome identifier (e.g. 'hg38') or path to a FASTA")
    p.add_argument("--size", default="given",
                    help="HOMER -size argument: 'given' uses each peak's actual coordinates "
                         "(recommended for variable-width ATAC peaks)")
    p.add_argument("--n_threads", type=int, default=4)
    p.add_argument("--homer_script", default=os.environ.get("decoder_homer_script", "findMotifsGenome.pl"))
    p.add_argument("--out_dir", required=True)
    p.add_argument("--dry_run", action="store_true",
                    help="Print the HOMER commands without executing them")
    return p.parse_args()


def build_background_bed(peaks_used_path: str, regex: str, out_path: str):
    df = pd.read_csv(peaks_used_path)
    if "used_in_lsi" not in df.columns or "peak" not in df.columns:
        raise ValueError(f"{peaks_used_path} must have 'peak' and 'used_in_lsi' columns; "
                          f"found: {list(df.columns)}")
    tested = df.loc[df["used_in_lsi"], "peak"]

    pattern = re.compile(regex)
    rows = []
    unparsed = 0
    for pid in tested:
        m = pattern.match(str(pid))
        if not m:
            unparsed += 1
            continue
        gd = m.groupdict()
        rows.append((gd["chr"], int(gd["start"]), int(gd["end"]), pid))
    if unparsed > 0:
        print(f"  WARNING: {unparsed} background peak ID(s) unparseable, dropped.")

    bg_df = pd.DataFrame(rows, columns=["chrom", "start", "end", "peak"])
    bg_df.to_csv(out_path, sep="\t", header=False, index=False)
    print(f"  Background: {len(bg_df)} peaks -> {out_path}")
    return out_path


def run_homer(foreground_bed, genome, out_dir, background_bed, size, n_threads, homer_script, dry_run):
    cmd = [homer_script, foreground_bed, genome, out_dir, "-bg", background_bed, "-size", str(size), "-p", str(n_threads)]
    print(f"  Command: {' '.join(cmd)}")
    if dry_run:
        return {"command": cmd, "executed": False}

    os.makedirs(out_dir, exist_ok=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    log_path = os.path.join(out_dir, "homer_stdout_stderr.log")
    with open(log_path, "w") as f:
        f.write("STDOUT:\n" + result.stdout + "\n\nSTDERR:\n" + result.stderr)

    if result.returncode != 0:
        print(f"  WARNING: HOMER exited with code {result.returncode} for {foreground_bed}. See {log_path}.")
    return {"command": cmd, "executed": True, "returncode": result.returncode, "log": log_path}


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("Building background BED (all distal, prevalence-passing peaks tested for stability)...")
    background_bed = build_background_bed(args.peaks_used, args.peak_coord_regex,
                                           os.path.join(args.out_dir, "background.bed"))

    bed_files = sorted(glob.glob(os.path.join(args.stable_peaks_dir, "stable_peaks_*.bed")))
    if not bed_files:
        print(f"No stable_peaks_*.bed files found in {args.stable_peaks_dir}. Nothing to run.")
        sys.exit(1)

    manifest = []
    for bed_path in bed_files:
        component = os.path.basename(bed_path).replace("stable_peaks_", "").replace(".bed", "")
        n_lines = sum(1 for _ in open(bed_path))
        if n_lines == 0:
            print(f"\n[{component}] Skipping -- 0 stable peaks (empty BED from step 5).")
            manifest.append({"component": component, "skipped": True, "reason": "0 stable peaks"})
            continue

        print(f"\n[{component}] {n_lines} stable peaks -> running HOMER...")
        comp_out_dir = os.path.join(args.out_dir, component)
        run_info = run_homer(bed_path, args.genome, comp_out_dir, background_bed,
                              args.size, args.n_threads, args.homer_script, args.dry_run)
        run_info["component"] = component
        run_info["n_foreground_peaks"] = n_lines
        manifest.append(run_info)

    manifest_path = os.path.join(args.out_dir, "motif_enrichment_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWritten: {manifest_path}")
    print("Done." if not args.dry_run else "Dry run complete -- no HOMER jobs were executed.")


if __name__ == "__main__":
    main()
