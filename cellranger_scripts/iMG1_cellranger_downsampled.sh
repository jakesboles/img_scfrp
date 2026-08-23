#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name iMG1_cellranger_downsampled
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A.log
#SBATCH --verbose

# Runs in Cell Ranger cluster mode: this job itself only hosts the
# lightweight Martian (mrp) controller process, which submits one sbatch job
# per pipeline stage chunk to b1042/genomics via slurm.template. See
# cellranger_scripts/test_cluster_mode.sh to validate the template before
# running this.
#
# Covers only the Gate38_11.22.2024 capture. The Gate40 "redo" capture is a
# separate GEM well/gel bead pool from a different cell aliquot and is run
# independently via iMG1_redo_cellranger.sh -- see that script's comment.

cd /projects/b1169/boles/img_scfrp/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger multi \
--id "iMG1_downsampled" \
--csv "../cellranger_configs/iMG1_config_downsampled.csv" \
--jobmode "../cellranger_scripts/slurm.template" \
--maxjobs 24 \
--jobinterval 100 \
--localcores 4 \
--localmem 14 \
--disable-ui
