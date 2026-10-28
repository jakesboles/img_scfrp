#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name cluster_mode_test
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 02:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A.log
#SBATCH --verbose

# Validates slurm.template and cluster-mode job submission before running any
# of the real iMG1-4 jobs. Uses cellranger's built-in tiny test dataset
# (`cellranger testrun`), which exercises the exact same jobmode/sbatch code
# path as a real `multi` run but finishes in minutes and submits ~200 small,
# cheap jobs instead of a multi-day run. See cellranger_scripts/README.md for
# how to run this and what to check.

cd /projects/b1169/boles/img_scfrp/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger testrun \
--id "cluster_mode_test" \
--jobmode "../cellranger_scripts/slurm.template" \
--maxjobs 24 \
--jobinterval 100 \
--localcores 4 \
--localmem 14 \
--disable-ui
