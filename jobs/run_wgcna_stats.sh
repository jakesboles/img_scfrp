#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name wgcna_stats
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 96G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose
#
# Runs after wgcna.R (needs results/wgcna/module_members_consensus.csv).
# UCell module scoring + kNN smoothing over the whole-cohort object,
# genotype x treatment lmer stats per module.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/r_scripts/wgcna_stats.R
