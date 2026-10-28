#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 06_qc2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 96G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%j.log
#SBATCH --verbose
#
# For the batch-runnable back half only, once doublet_clusters has been
# filled in by eye (see r_scripts/06_qc2.R) -- same precedent as
# monoculture's own run_06_qc2.sh.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/cocultures/r_scripts/06_qc2.R
