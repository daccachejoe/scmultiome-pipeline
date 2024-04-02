## scMulitome Analysis Pipeline
#### Naik Lab
##### Joseph Daccache (05/02/2023)
This github repo contains the necessary scripts to run the scmultiome wrapper functions.  
In order for the pipeline to work, the user must alllow for R scripts to be executed as typical command line functions.  

### Installation
In the directory of your choosing, clone this repository into a project direcotry name of your choosing. For this example it is `test_dir`
example installation:
```
git clone https://github.com/Naiklab/scmultiome-pipeline.git test_dir
cd test_dir
# initialize the directory set up with the "init" route
chmod +x run/runmultiome
run/runmultiome init
```
Now create your conda enviroment and follow the installation steps script to install the correct versions of the required packages




This should set up your project directory to be able to use `argparser` and have the directories set up in the correct manner
```
[jd5457@bigpurple-ln1 testdir]$ tree
.
├── configs
│   ├── cluster_labels.csv
│   ├── qc_df.csv
│   ├── resolution_to_use.txt
│   ├── samplesheet.csv
│   ├── scenicplus-preprocess-config.yml
│   └── scenicplus-process-config.yml
├── data
│   ├── raw
│   │   └── macs-peaks
│   └── scenicplus
│       └── cisTarget_dbs
├── output
│   ├── plots
│   ├── RDS-files
│   ├── tables
│   └── ucd
├── README.md
├── routes
│   ├── call_peaks_grouped.sh
│   ├── filter_and_cluster.sh
│   ├── identify_celltypes.sh
│   ├── install.sh
│   ├── label_celltypes.sh
│   ├── run_scenicplus.sh
│   ├── setup_dirs.sh
│   ├── seurat_preprocess.sh
│   └── test_filter_and_cluster.sh
├── runmultiome.sh
└── scripts
    ├── convert_seruat_to_h5ad.R
    ├── export_scenicplus_data.R
    ├── functions.R
    ├── label_celltypes.R
    ├── logs
    ├── multiome-processing.R
    ├── plotting-UCD-and-seurat.R
    ├── scenicplus-pipeline.py
    └── ucd-script.py
```
## Ensuring you have the correct packages installed
Seurat And Signac Processing of sc-multiome data (everything except for SCENIC+ analysis)
* Seurat
* Signac
* EnsDb.Hsapiens.v86
* enrichR
* ggplot2
* dplyr
* motifmatchr
* TFBSTools
* JASPAR2020

For SCENIC+ Analysis, three conda environments must be made  
First: *sceasy* is used to convert R objects to h5ad files for use [package link](https://github.com/cellgeni/sceasy)  



### Establishing the required files
samplesheet.csv: the skeleton of your samplesheet is already geenrated for you in `configs/samplesheet.csv`. Let's read it.  
The expected set up of your samplehseet is flexible but the first columns **MUST BE** sampleName and path, and those are written in for you already.
```
[jd5457@bigpurple-ln1 testdir]$ cat configs/samplesheet.csv
sampleName,path
```
Now you can fill in the rows as needed as well as add any additional meta data you would like to include in your single cell experiments. Here is an example with 2 samples and an additional meta data column `cond` added

| sampleName | path                                                                       | cond  |
| ---------- | -------------------------------------------------------------------------- | ----- |
| ctrl.1     | /gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-CTRL/outs  | ctrl  |
| il17.1     | /gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-IL-17/outs | il17a |


### Running the multiome pipeline
Now that our samplesheet is set up, we can run the first step of the pipeline which will be specified by `seurat_preprocess`
```
run/runmultiome seurat_preprocess
```




```
usage: multiome-processing-v0.1.R [--] [--help] [--opts OPTS]
       [--outfilename OUTFILENAME] [--grouping.var GROUPING.VAR]
       [--RDS.file.in RDS.FILE.IN] [--runHarmony RUNHARMONY]
       [--footprint-peaks FOOTPRINT-PEAKS] [--my.macs2.path
       MY.MACS2.PATH] pipeline samplesheet

Run single-cell RNA + ATAC Multiomic Analysis from 10X Genomics
platform

positional arguments:
  pipeline               Comma delimted combinations of: init, create,
                         callpeaks, qc, cluster, merge, linkpeaks
  samplesheet            samplesheet in csv format

flags:
  -h, --help             show this help message and exit

optional arguments:
  -x, --opts             RDS file containing argument values
  -o, --outfilename      outfile name, no .RDS! [default:
                         data/r-objects/multiome-object]
  -g, --grouping.var     grouping variable for peak calling algorithm
                         [default: NA]
  -R, --RDS.file.in      RDS in-file for the pipeline desired [default:
                         NA]
  -r, --runHarmony       whether or not to run Harmony batch correction
                         [default: FALSE]
  -f, --footprint-peaks  atac peaks to run footprinting analysis on.
                         txt file
  -m, --my.macs.path     path to macs environment. only used if
                         callpaks pipeline is run
```
As the output indicates, the scripts excepts 2 mandatory arguments: pipeline(s) and samplesheet.  
The example samplesheet shows what a samplesheet looks like.  
IMPORTANT: the pipeline expects `sampleName` and `path` in a case sensitive manner and all sample-level meta data information to be added to come *after* the `path` column
