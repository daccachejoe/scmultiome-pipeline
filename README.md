## scMulitome Analysis Pipeline
#### Naik Lab
##### Joseph Daccache (05/02/2023)
This github repo contains the necessary scripts to run the scmultiome wrapper functions.  
In order for the pipeline to work, the user must alllow for R scripts to be executed as typical command line functions.  

### Installation
In the directory of your choosing, clone this repository and souece the `install.sh` file  
example installation:
```
git clone path/to/repository
cd respotory name
source install.sh
```

Let's test the installation by running 
```
test.R -h
```

```
test.R 3.14159 
test.R 3.14159 -d 2
```

### Running our scMulitome pipeline
After ensuring argparser is functional, let's test out our pipeline and make sure it is executable and working
```
scripts/multiome-processing.R -h
```
As the output indicates, the scripts excepts 2 mandatory arguments: pipeline(s) and samplesheet.  
The example samplesheet shows what a samplesheet looks like.  
IMPORTANT: the pipeline expects `sampleName` and `path` in a case sensitive manner and all sample-level meta data information to be added to come *after* the `path` column

|sampleName|path|cond|
|---|---|---|
|ctrl.1|/gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-CTRL/outs|ctrl|
|il17.1|/gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-IL-17/outs|il17a|

