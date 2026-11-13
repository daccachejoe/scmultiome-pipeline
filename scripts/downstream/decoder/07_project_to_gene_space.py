#!/usr/bin/env python3
"""
Decoder branch, step 4 (RNA side): project decoder component weights from
PCA space back to gene space.

RNA counterpart of 04_project_to_peak_space.py -- gene_weights =
rna_pca_loadings @ W, the same loadings @ W.T math (shared in
projection.py), just with genes standing in for peaks and PCA loadings
standing in for LSI loadings. No ATAC-style distal/prevalence filtering
concern here: rna_pca_loadings.csv already covers every gene that went
into step 2's PCA.

This is a standalone diagnostic/interpretation step, not a dependency of
anything else in the branch: unlike the ATAC side (peak projection ->
stability test -> motif enrichment), there's no "gene stability test" or
downstream RNA-only step consuming this output yet. Its purpose is
cross-modality concordance -- for a given regulatory component, do genes
with strong RNA-side signal here overlap with genes near the ATAC-side
peaks that same component flags (via gene_weights_final.csv's ranking
alongside peak_weights_final.csv, or gene_weights_per_fold.csv joined to
peak_weights_per_fold.csv fold-by-fold for a more rigorous version).

USAGE:
    python 07_project_to_gene_space.py \
        --loadings output/decoder/features/rna_pca_loadings.csv \
        --cv_fold_weights output/decoder/decoder_rna/cv_fold_W_components.csv \
        --final_weights output/decoder/decoder_rna/final_decoder_W.csv \
        --out_dir output/decoder/gene_projections
"""

import argparse
import os

import pandas as pd

from projection import project_fold_weights, project_final_weights


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loadings", required=True, help="Path to rna_pca_loadings.csv from step 2")
    p.add_argument("--cv_fold_weights", required=True,
                    help="Path to cv_fold_W_components.csv from step 3's RNA fit "
                         "(long format: held_out_fold, component, feature, weight)")
    p.add_argument("--final_weights", default=None,
                    help="Optional: final_decoder_W.csv from step 3's RNA fit for a point-estimate projection")
    p.add_argument("--out_dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    loadings = pd.read_csv(args.loadings, index_col=0)
    fold_weights_long = pd.read_csv(args.cv_fold_weights)

    print(f"Loadings: {loadings.shape[0]} genes x {loadings.shape[1]} PCA components")
    print(f"Fold weights: {fold_weights_long['held_out_fold'].nunique()} folds x "
          f"{fold_weights_long['component'].nunique()} components x "
          f"{fold_weights_long['feature'].nunique()} PCA components")

    per_fold_gene_weights = project_fold_weights(loadings, fold_weights_long, item_label="gene")
    out_path = os.path.join(args.out_dir, "gene_weights_per_fold.csv")
    per_fold_gene_weights.to_csv(out_path, index=False)
    print(f"Written: {out_path} ({len(per_fold_gene_weights)} rows)")

    if args.final_weights:
        final_W = pd.read_csv(args.final_weights, index_col=0)
        final_gene_weights = project_final_weights(loadings, final_W)
        out_path2 = os.path.join(args.out_dir, "gene_weights_final.csv")
        final_gene_weights.to_csv(out_path2)
        print(f"Written: {out_path2} ({final_gene_weights.shape[0]} genes x "
              f"{final_gene_weights.shape[1]} components)")

    print("Done.")


if __name__ == "__main__":
    main()
