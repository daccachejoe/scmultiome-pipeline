## scMulitome Analysis Pipeline
#### Naik Lab
##### Joseph Daccache (05/02/2023)
This github repo contains the necessary scripts to run the scmultiome wrapper functions.  
In order for the pipeline to work, the user must alllow for R scripts to be executed as typical command line functions.  
Follow instructions from: https://bitbucket.org/djhshih/argparser/src/master/  
In short,
* clone their repository
* open R and run some functions
* run `R CMD INSTALL .`
* `chmod +x` any R files you wish to execute

**Running our scMulitome pipeline**
```
# test out the command line interface
scripts/multiome-processing.R -h
```
As the output indicates, the scripts excepts 2 mandatory arguments: pipeline(s) and samplesheet.  
The example samplesheet shows what a decent samplesheet looks like:

|sampleName|path|cond|
|---|---|---|
|ctrl.1|/gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-CTRL/outs|ctrl|
|il17.1|/gpfs/data/sequence/results/naiklab/2023-03-24/cellranger/count-IL-17/outs|il17a|




