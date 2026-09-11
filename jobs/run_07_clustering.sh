#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_clustering
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 300G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose
#
# Clusters the whole cohort at two resolutions and runs FindAllMarkers()
# for each. Needs presto installed for FindAllMarkers() to be fast --
# 07_clustering.R warns at the top if it's missing.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/r_scripts/07_clustering.R
