#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 05_integration_harmony
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 96G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%j.log
#SBATCH --verbose
#
# Resources copied from monoculture's run_05_integration_harmony.sh
# unchanged -- adjust if the actual run needs more (this arm's cohort is
# comparable in scale or larger).

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/cocultures/r_scripts/05_integration_harmony.R
