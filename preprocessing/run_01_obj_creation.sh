#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition short
#SBATCH --job-name 01_obj_creation
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 160GB
#SBATCH --time 4:00:00
#SBATCH --output /projects/p31535/boles/img_scfrp/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.4.0
module load gcc/11.2.0
module load hdf5/1.14.1-2-gcc-12.3.0 

Rscript /projects/p31535/boles/img_scfrp/r_scripts/01_obj_creation.R