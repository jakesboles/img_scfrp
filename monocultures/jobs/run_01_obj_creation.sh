#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 01_obj_creation
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/r_scripts/01_obj_creation.R
