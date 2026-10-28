#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name JSB134-1_cellranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%A.log
#SBATCH --verbose

# Runs in Cell Ranger cluster mode: this job itself only hosts the
# lightweight Martian (mrp) controller process, which submits one sbatch job
# per pipeline stage chunk to b1042/genomics via slurm.template. See
# cellranger_scripts/test_cluster_mode.sh to validate the template before
# running this, and monocultures/cellranger_scripts/README.md for the
# general cluster-mode mechanics (identical here).
#
# RERUN NEEDED (see cellranger_scripts/README.md): an earlier copy-paste
# error in this arm's config files pointed JSB134-1 and JSB134-2 at the
# same underlying FASTQ selection, so the existing cellranger/JSB134-1/
# output on disk does not reflect this pool's current (corrected)
# ../configs/JSB134-1_config.csv. Remove or rename the stale
# cellranger/JSB134-1/ directory before submitting this -- Cell Ranger
# will refuse to run into an existing pipestance directory.

cd /projects/b1169/boles/img_scfrp/cocultures/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger multi \
--id "JSB134-1" \
--csv "../configs/JSB134-1_config.csv" \
--jobmode "../cellranger_scripts/slurm.template" \
--maxjobs 24 \
--jobinterval 100 \
--localcores 4 \
--localmem 14 \
--disable-ui
