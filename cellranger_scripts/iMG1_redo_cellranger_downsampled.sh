#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name iMG1_redo_cellranger_downsampled
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
# Kept separate from iMG1_cellranger.sh: iMG1 (Gate38_11.22.2024) and
# iMG1_redo (Gate40) are two independent GEM wells/gel bead pools captured
# from two different cell aliquots, so their 10x cell barcodes are not
# comparable -- processing them in one `multi` run would risk collapsing
# distinct cells from each pool that happen to share a barcode into one
# called cell.

cd /projects/b1169/boles/img_scfrp/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger multi \
--id "iMG1_redo_downsampled" \
--csv "../cellranger_configs/iMG1_redo_config_downsampled.csv" \
--jobmode "../cellranger_scripts/slurm.template" \
--maxjobs 24 \
--jobinterval 100 \
--localcores 4 \
--localmem 14 \
--disable-ui
