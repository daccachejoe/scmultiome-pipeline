installation steps

1. create new env
```
conda create -n sc-multiome-pipeline-env python=3.
conda activate scmultiome-pipeline-env
```
activate cmake (module loaded here)
```module load cmake```

Install SCENIC+ as defined by the authors (Link here)[https://scenicplus.readthedocs.io/en/latest/install.html]



install cCreateCisTargetDatabase
* everything is already installed except for python-flatbuffers
* ```conda install python-flatbuffers```
Install cluster buster and other subsidiary packages they ask you to in the installation vignette

Installing sceasy for conversions
install r-base using conda (version 4.2.2!!)
```
conda install r-base=4.2.2
conda install -c conda-forge libxml2 # required for BioConductor packages
conda install -c r r-xml=3.98_1.5 # also required?
```




In an alternative R session, find the HPC library for 4.2, NYU's is `/gpfs/share/apps/R/4.2.2/lib64/R/library`, becasue `devtools` is hard to install, we will have to add it to our libpath to install `sceasy` in R.  
In `R`
```
.libPaths(c(.libPaths(), "/gpfs/share/apps/R/4.2.2/lib64/R/library"))
```
The order of the paths is very important as R will search the first directory before the second one!
Now we can install sceasy
```
devtools::install_github("cellgeni/sceasy")
```
Restart the R session to return .libPaths to the original directory and now we will install Seurat (Version 4!!!!)
```
quit("no")
# in command line `R`
# check .libPaths() only contains the conda env library
.libPaths()
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))

BiocManager::install("GenomeInfoDb", "EnsDb.Hsapiens.v86", )

BiocManager::install("EnsDb.Hsapiens.v86")

# overwrite the newer Satija lab packages with the older versions for compatiability reasons
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
```

```
# install Signac
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac") # do not allow the installer to upgrade any other packages
```



in R:
```
R
install.packjages("remotes")
remotes::install


```




1. install rust
```
cd scmultiome-pipeline-env/
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
*Restart your session otherwise the rust installation will not complete*
1. SCENIC+ installation
```
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install -e .
```
1. 




r 4.1.2
Seurat v4 (!!!!)
BiocManager::install("GenomeInfoDb") may be required with changing annotation file to USCS style
