#!/bin/bash
# Creates required output/log directories (SLURM will not create them for
# #SBATCH --output) and submits all 5 pool jobs at once so they run
# concurrently. iMG1 and iMG1_redo are two independent GEM wells from the
# same batch-1 samples (see iMG1_redo_cellranger.sh) and are run separately.

mkdir -p "/projects/b1169/boles/img_scfrp/logs" "/projects/b1169/boles/img_scfrp/cellranger"

for pool in iMG1 iMG1_redo iMG2 iMG3 iMG4; do
    sbatch "$(dirname "$0")/${pool}_cellranger.sh"
done
