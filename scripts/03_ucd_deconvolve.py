# script for running UCDdeonvolve on a single dataset

# ucdenv script from an h5ad file
import argparse
import os
import sys
import scanpy as sc
import ucdeconvolve as ucd

# parse argument
parser = argparse.ArgumentParser(description = "UCDenvolve on a single sample")
# input yaml file
# named flags, matching how routes/03_identify_celltypes.sh calls this script
parser.add_argument("--resolution", default = "seurat_clusters", help = "metadata column to group cells by for cell type calling")
parser.add_argument("--input_file", required = True, help = "input h5ad file name")
parser.add_argument("--reference_file", default = None, help = "optional reference h5ad; empty or omitted skips the referenced run")
args = parser.parse_args()

# the API token comes from ucd_api_token in config/pipeline.config (exported
# by run/runmultiome) rather than being hardcoded in this tracked file
api_token = os.environ.get("ucd_api_token", "")
if not api_token:
    sys.exit("ucd_api_token is not set -- add it to config/pipeline.config (see pipeline.config.example).")
ucd.api.authenticate(api_token)

def ReadinObjects(args):
    # read in the objects; the reference is optional
    adata = sc.read_h5ad(args.input_file)
    reference = sc.read_h5ad(args.reference_file) if args.reference_file else None
    return(adata, reference)


def run_unbiased(adata):
    # perform base cell type calling 
    ucd.tl.base(adata)

    # visualize results using a heatmap
    plot=ucd.pl.base_clustermap(adata, groupby = args.resolution, n_top_celltypes=100)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap.png")

    # examining raw predictions
    plot=ucd.pl.base_clustermap(adata, groupby = args.resolution, category = 'raw', n_top_celltypes = 75)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-raw.png")

    # assign labels (first pass)
    ucd.utils.assign_top_celltypes(adata, category = "raw", groupby = args.resolution)
    celltypes = ucd.utils.assign_top_celltypes(adata, groupby = args.resolution, category = "raw",  inplace = False)

    # export the annotations
    cellDF=adata.obs
    cellDF.to_csv("./output/ucd/cellmetadata-unbiased.csv")

    # explain the genes driving the cell type prediction
    ucd.tl.explain(adata, celltypes = celltypes, groupby = args.resolution, group_n = 64)
    plot=ucd.pl.explain_clustermap(adata, n_top_genes= 128)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-driving-genes.png")

    for value in celltypes.values():
        plot=ucd.pl.explain_boxplot(adata, key = "ucdexplain", celltypes=value, n_top_genes = 16, ncols = 4, return_fig = True)
        plot.savefig(f"./output/ucd/ucdeconvolve-clustermap-boxplot-{value}-unbiased.png")


def run_referenced(adata,reference):
    ucd.tl.select(adata, reference, reference_key = "cell_types")
    celltypes = ucd.utils.assign_top_celltypes(adata, groupby = args.resolution, category = "raw",  inplace = False)

    # visualize results using a heatmap
    plot=ucd.pl.base_clustermap(adata, groupby = args.resolution, n_top_celltypes=100)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-referenced.png")

    # examining raw predictions
    plot=ucd.pl.base_clustermap(adata, groupby = args.resolution, category = 'raw', n_top_celltypes = 75)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-raw-referenced.png")

    # assign labels (first pass)
    ucd.utils.assign_top_celltypes(adata, category = "raw", groupby = args.resolution)
    celltypes = ucd.utils.assign_top_celltypes(adata, groupby = args.resolution, category = "raw",  inplace = False)

    # explain the genes driving the cell type prediction
    ucd.tl.explain(adata, celltypes = celltypes, groupby = args.resolution, group_n = 64)
    plot=ucd.pl.explain_clustermap(adata, n_top_genes= 128)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-driving-genes-referenced.png")

    for value in celltypes.values():
        plot=ucd.pl.explain_boxplot(adata, key = "ucdexplain", celltypes=value, n_top_genes = 16, ncols = 4, return_fig = True)
        plot.savefig(f"./output/ucd/ucdeconvolve-clustermap-boxplot-{value}-referenced.png")

    cellDF=adata.obs
    cellDF.to_csv("./output/ucd/cellmetadata-referenced.csv")


def __main__():
    adata, reference_data = ReadinObjects(args)
    run_unbiased(adata)

    if reference_data is not None:
        run_referenced(adata, reference_data)


if __name__ == "__main__":
    __main__()

