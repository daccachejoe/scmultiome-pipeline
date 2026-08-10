
import anndata as ad
import scanpy as sc
import argparse 

parser = argparse.ArgumentParser(description = "Reformatting sceasy h5ad anndata object to correct format")
parser.add_argument("in_file")
parser.add_argument("out_file")

args=parser.parse_args()

# Load in data
obj=sc.read_h5ad(args.in_file)

# reformat object
obj.raw=ad.AnnData(obj.X)
new_adata = obj.raw.to_adata()
new_adata.var_names = obj.var_names #Or whatever your index is called
new_adata.var_names.name = None
new_adata.raw = new_adata
sc.pp.normalize_total(new_adata, target_sum=1e4)
sc.pp.log1p(new_adata)

# write the new object
new_adata.write(args.out_file)