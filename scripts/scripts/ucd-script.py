# script for running UCDdeonvolve on a single dataset

# ucdenv script from an h5ad file
import scanpy as sc
import ucdeconvolve as ucd

# parse argument
parser = argparse.ArgumentParser(description = "UCDenvolve on a single sample")
# input yaml file
parser.add_argument("resolution",  help = "resolution to group cells by for cell type calling", default = "seurat_clusters")
parser.add_argument("input_file",  help = "input file name")
parser.add_argument("reference_file",  help = "reference file name", default = None)
args = parser.parse_args()

# use api token to log into database
ucd.api.authenticate("454dd14a288146dd86a578f2a7a20ded151f2dafb7610e34da16b07549946033")

def ReadinObjects(args):
    # read in the objects
    adata = sc.read_h5ad(args.input_file)
    hannifa_reference = sc.read_h5ad(args.reference_file)
    return(adata, hannifa_reference)


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
    celltypes = ucd.utils.assign_top_celltypes(adata, groupby = rgs.resolution, category = "raw",  inplace = False)

    # export the annotations
    cellDF=adata.obs
    cellDF.to_csv("data/r-objects/anndata/cellmetadata-unbiased.csv")

    # explain the genes driving the cell type prediction
    ucd.tl.explain(adata, celltypes = celltypes, groupby = args.resolution, group_n = 64)
    plot=ucd.pl.explain_clustermap(adata, n_top_genes= 128)
    plot.savefig("./output/ucd/ucdeconvolve-clustermap-driving-genes.png")

    for value in celltypes.values():
        plot=ucd.pl.explain_boxplot(adata, key = "ucdexplain", celltypes=value, n_top_genes = 16, ncols = 4, return_fig = True)
        plot.savefig(f"./output/ucd/ucdeconvolve-clustermap-boxplot-{value}-.png")


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
    adata, reference_data = ReadinObjects(args.input_file, args.reference_file)
    run_unbiased(adata)
    
    if args.reference_file is not None:
        run_referenced(adata, reference_data)

