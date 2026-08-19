#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name 05_2-3_only
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 164G
#SBATCH --time 48:00:00
#SBATCH --output /projects/p31535/boles/img_scfrp/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.4.0
module load gcc/11.2.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/p31535/boles/img_scfrp/r_scripts/05_integration_harmony_2-3.R
