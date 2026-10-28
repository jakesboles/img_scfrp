#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name JSB127_cellranger
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
# Not confirmed affected by the JSB134-1/JSB134-2 copy-paste config bug
# (see cellranger_scripts/README.md) -- included here for completeness, not
# because this pool is known to need a rerun.

cd /projects/b1169/boles/img_scfrp/cocultures/cellranger

/projects/p31535/boles/cellranger-10.0.0/cellranger multi \
--id "JSB127" \
--csv "../configs/JSB127_config.csv" \
--jobmode "../cellranger_scripts/slurm.template" \
--maxjobs 24 \
--jobinterval 100 \
--localcores 4 \
--localmem 14 \
--disable-ui
