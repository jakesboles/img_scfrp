#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 03_doubletfinder
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose

module load R/4.4.0
module load gcc/11.2.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/preprocessing/03_doubletfinder.R
