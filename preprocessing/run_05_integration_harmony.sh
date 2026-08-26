#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 05_integration_harmony
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 96G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose
#
# Resources scaled down from the old CCA-oriented job script (164G/48h) --
# Harmony only needs the PCA embedding to integrate, not dense per-sample
# expression access like CCA does, so it's meaningfully lighter/faster.
# Adjust if the actual run needs more.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/preprocessing/05_integration_harmony.R
