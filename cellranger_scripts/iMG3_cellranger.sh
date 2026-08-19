#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name iMG3_cellranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 350G
#SBATCH --time 96:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A.log
#SBATCH --verbose

cd /projects/b1169/boles/img_scfrp/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger multi \
--id "iMG3" \
--csv "../cellranger_configs/iMG3_config.csv" \
--localcores 16 \
--localmem 340 \
--disable-ui
