#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 03_doubletfinder
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/cocultures/logs/%x_%A_%a.log
#SBATCH --verbose
#
# One task per (pool_dir, sample) capture -- 30 total for this arm (10
# genotype-pair conditions x 3 batch replicates each). Run directly (not
# via sbatch) the first time; it computes the array size from
# cellbender_scripts/cell_quantities.csv and resubmits itself as a SLURM
# array, matching the pattern used elsewhere in this repo (seqtk/seqtk.sh,
# cellbender_scripts/run_all_cellbender*.sh, and monoculture's own
# run_03_doubletfinder.sh).

REPO="/projects/b1169/boles/img_scfrp/cocultures"
CSV="$REPO/cellbender_scripts/cell_quantities.csv"

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    if [ ! -f "$CSV" ]; then
        echo "Missing $CSV" >&2
        exit 1
    fi
    n=$(( $(wc -l < "$CSV") - 1 ))
    echo "Submitting ${n} DoubletFinder jobs (1-${n})"
    exec sbatch --array=1-${n} "$0"
fi

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/img_scfrp/cocultures/r_scripts/03_doubletfinder.R
