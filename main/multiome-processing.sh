#!/bin/bash
#SBATCH --partition a100_short
#SBATCH --nodes 1
#SBATCH --cpus-per-task 32
#SBATCH --mem 320GB
#SBATCH --mail-type=END
#SBATCH --mail-user=joseph.daccache@nyulangone.org
#SBATCH --time 2-23:00:00
#SBATCH --job-name multiome
#SBATCH --output scripts/outs/multiome-%J.out
#SBATCH --error scripts/outs/multiome-%J.e

# module load macs2
module load r
# scripts/multiome-processing-v0.1.R \
#     init,create,callpeaks \
#     data/samplesheet.csv \
#     -o data/r-objects/mulitome-SCENIC \
#     -m /gpfs/data/naiklab/jd5457/opt/MACS/bin/macs3

# scripts/multiome-processing-v0.1.R \
#     qc \
#     data/samplesheet.csv \
#     -o data/r-objects/mulitome-SCENIC \
#     -R data/r-objects/mulitome-SCENIC-callpeaks-obj-list.RDS \
#     -m /gpfs/data/naiklab/jd5457/opt/MACS/bin/macs3

scripts/multiome-processing-v0.1.R \
    filter,merge \
    data/samplesheet.csv \
    -R data/r-objects/mulitome-SCENIC-qc-obj-list.RDS \
    -o data/r-objects/mulitome-SCENIC \
    -q output/clusters.to.remove.csv \
    --RunHarmony
