#!/bin/bash
# Creates required output/log directories (SLURM will not create them for
# #SBATCH --output) and submits all 4 pool jobs at once so they run
# concurrently.

mkdir -p "/projects/b1169/boles/img_scfrp/logs" "/projects/b1169/boles/img_scfrp/cellranger"

for pool in iMG1 iMG2 iMG3 iMG4; do
    sbatch "$(dirname "$0")/${pool}_cellranger.sh"
done
