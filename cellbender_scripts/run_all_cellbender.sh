#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name cellbender
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 30G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A_%a.log
#SBATCH --verbose
#
# Replaces Thomas's CellbenderBatch.sh + CellbenderInd.sh (one sbatch call
# per sample, submitted from a plain bash loop) with a single SLURM array
# job, one task per (pool, sample) row in cell_quantities.csv -- matching
# the seqtk/cellranger pattern used elsewhere in this repo, and giving a
# single throttle point (--array=1-N%k below) for how many A100 jobs run
# concurrently, since the GPU partition is far more contended than CPU
# nodes.
#
# Usage:
#   Rscript cellbender_scripts/make_cellbender_params.R   # writes cell_quantities.csv
#   bash cellbender_scripts/run_all_cellbender.sh  # self-submits the array job

REPO="/projects/b1169/boles/img_scfrp"
CSV="$REPO/cellbender_scripts/cell_quantities.csv"

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    # Not running under SLURM yet: this is the plain `bash run_all_cellbender.sh`
    # invocation. Compute the array size from the CSV (minus its header) and
    # resubmit this same script as an array job.
    if [ ! -f "$CSV" ]; then
        echo "Missing $CSV -- run make_cellbender_params.R first." >&2
        exit 1
    fi
    n=$(( $(wc -l < "$CSV") - 1 ))
    echo "Submitting ${n} CellBender jobs (1-${n})"
    # %8 caps concurrent A100 usage -- raise/lower to match how much of the
    # genomics-gpu partition you can actually use at once.
    exec sbatch --array=1-${n}%8 "$0"
fi

module purge
module load python-miniconda3/4.10.3
module load mamba
source activate "$REPO/envs/cellbender"

# Without this, CellBender's report-generation step fails at the very end
# with UnicodeDecodeError ('ascii' codec can't decode byte ...) any time the
# rendered HTML happens to contain a non-ASCII byte -- Python falls back to
# the locale's default encoding when a file is opened without one specified,
# and this cluster's default locale is plain ASCII/C, not UTF-8. Everything
# before that point (the actual denoised .h5 outputs) still completes either
# way, but this avoids the job ending on a traceback for a report that isn't
# otherwise needed.
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CSV")
IFS=',' read -r POOL_DIR SAMPLE CELLS DROPLETS <<< "$LINE"

# Catches a real bug we hit: cell_quantities.csv briefly had a spurious
# second row per sample (a library-wide total mislabeled as that sample's
# count) that produced huge, wrong --expected-cells/--total-droplets-included
# values and burned GPU time before crashing partway through. Fail fast
# instead of launching cellbender on obviously-bad input.
if ! [[ "$CELLS" =~ ^[0-9]+$ && "$DROPLETS" =~ ^[0-9]+$ && "$CELLS" -gt 0 && "$DROPLETS" -gt "$CELLS" ]]; then
    echo "Bad row for task ${SLURM_ARRAY_TASK_ID}: '${LINE}' (pool_dir=${POOL_DIR} sample=${SAMPLE} cells=${CELLS} droplets=${DROPLETS})" >&2
    exit 1
fi

OUT_DIR="$REPO/cellbender/${POOL_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

# CellBender writes checkpoint files (ckpt.tar.gz) relative to the current
# working directory by default. Running two jobs from the same cwd lets one
# job's checkpoint lock out another's, so every task must cd into its own
# unique output directory before running -- do not remove this.
cd "$OUT_DIR"

echo "Task ${SLURM_ARRAY_TASK_ID}: pool=${POOL_DIR} sample=${SAMPLE} cells=${CELLS} droplets=${DROPLETS}"

cellbender remove-background \
--cuda \
--input "${REPO}/cellranger/${POOL_DIR}/outs/per_sample_outs/${SAMPLE}/sample_raw_feature_bc_matrix.h5" \
--output "${OUT_DIR}/${SAMPLE}" \
--expected-cells "$CELLS" \
--total-droplets-included "$DROPLETS" \
--fpr 0.01 \
--epochs 100
