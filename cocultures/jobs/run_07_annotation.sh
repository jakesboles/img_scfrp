#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_annotation
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%j.log
#SBATCH --verbose
#
# For the batch-runnable back half only, once cluster_labels has been
# filled in by eye (see r_scripts/07_annotation.R) -- same precedent as
# monoculture's run_06_qc2.sh for a script with an eyeballed decision in
# the middle of it.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/cocultures/r_scripts/07_annotation.R
