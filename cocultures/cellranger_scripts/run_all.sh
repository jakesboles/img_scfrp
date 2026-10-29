#!/bin/bash
# Creates required output/log directories (SLURM will not create them for
# #SBATCH --output) and submits the pools confirmed affected by the
# JSB134-1/JSB134-2 copy-paste config bug (see README.md) -- NOT all 4
# pools. JSB127 and JSB134-3 are not known to need a rerun; submit
# JSB127_cellranger.sh/JSB134-3_cellranger.sh yourself if that ever
# changes.
#
# IMPORTANT: remove or rename the existing cellranger/JSB134-1/ and
# cellranger/JSB134-2/ directories before running this -- Cell Ranger
# refuses to run into an existing pipestance directory, and those two
# currently hold the stale output from the config bug (see README.md).

mkdir -p "/projects/b1169/boles/img_scfrp/cocultures/logs" "/projects/b1169/boles/img_scfrp/cocultures/cellranger"

for pool in JSB127 JSB134-1 JSB134-2 JSB134-3; do
    sbatch "$(dirname "$0")/${pool}_cellranger.sh"
done
