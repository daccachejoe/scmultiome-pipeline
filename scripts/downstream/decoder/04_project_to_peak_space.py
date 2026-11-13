#!/usr/bin/env python3
"""
Decoder branch, step 4: project decoder component weights from LSI space
back to peak space.

ATAC-only step. See 07_project_to_gene_space.py for the RNA-side
equivalent (PCA loadings -> genes, gene_weights = rna_pca_loadings @ W) --
the underlying math (loadings @ W.T) is identical between modalities, so
it lives once in projection.py and both scripts just supply their own
loadings matrix and item label ("peak" vs "gene").

loadings (peaks x LSI components) comes from step 2's
atac_lsi_loadings.csv -- already restricted to distal peaks with LSI dim 1
dropped, so its columns already match W's feature axis one-to-one by name.
projection.py verifies that alignment explicitly rather than assuming it.

USAGE:
    python 04_project_to_peak_space.py \
        --loadings output/decoder/features/atac_lsi_loadings.csv \
        --cv_fold_weights output/decoder/decoder_atac/cv_fold_W_components.csv \
        --final_weights output/decoder/decoder_atac/final_decoder_W.csv \
        --out_dir output/decoder/peak_projections
"""

import argparse
import os

import pandas as pd

from projection import project_fold_weights, project_final_weights


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loadings", required=True, help="Path to atac_lsi_loadings.csv from step 2")
    p.add_argument("--cv_fold_weights", required=True,
                    help="Path to cv_fold_W_components.csv from step 3 "
                         "(long format: held_out_fold, component, feature, weight)")
    p.add_argument("--final_weights", default=None,
                    help="Optional: final_decoder_W.csv from step 3 for a point-estimate projection")
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    loadings = pd.read_csv(args.loadings, index_col=0)
    fold_weights_long = pd.read_csv(args.cv_fold_weights)

    print(f"Loadings: {loadings.shape[0]} peaks x {loadings.shape[1]} LSI components")
    print(f"Fold weights: {fold_weights_long['held_out_fold'].nunique()} folds x "
          f"{fold_weights_long['component'].nunique()} components x "
          f"{fold_weights_long['feature'].nunique()} LSI components")

    per_fold_peak_weights = project_fold_weights(loadings, fold_weights_long, item_label="peak")
    out_path = os.path.join(args.out_dir, "peak_weights_per_fold.csv")
    per_fold_peak_weights.to_csv(out_path, index=False)
    print(f"Written: {out_path} ({len(per_fold_peak_weights)} rows)")

    if args.final_weights:
        final_W = pd.read_csv(args.final_weights, index_col=0)
        final_peak_weights = project_final_weights(loadings, final_W)
        out_path2 = os.path.join(args.out_dir, "peak_weights_final.csv")
        final_peak_weights.to_csv(out_path2)
        print(f"Written: {out_path2} ({final_peak_weights.shape[0]} peaks x "
              f"{final_peak_weights.shape[1]} components)")

    print("Done.")


if __name__ == "__main__":
    main()
