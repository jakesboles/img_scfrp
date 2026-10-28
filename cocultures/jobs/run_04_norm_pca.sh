#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 04_norm_pca
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 164G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%j.log
#SBATCH --verbose
#
# Resources copied from monoculture's run_04_norm_pca.sh unchanged -- this
# arm's cohort (combined microglia + neurons, ~146k cells in the
# surviving reference object alone) is comparable in scale or larger, so
# not scaled down.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/cocultures/r_scripts/04_norm_pca.R
