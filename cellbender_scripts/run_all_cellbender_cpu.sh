#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name cellbender_cpu
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 30G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A_%a_%N.log
#SBATCH --verbose
#
# CPU-only variant of run_all_cellbender.sh, for while the genomics-gpu
# partition's one practically-accessible node (qgpu0102) is oversubscribed
# (see cellbender_scripts/README.md's "GPU troubleshooting" section for the
# full diagnosis). CellBender fully supports CPU-only inference -- omit
# --cuda and pass --cpu-threads instead -- and these datasets (5-10k cells,
# ~15k reads/cell) are well within what CPU inference handles in a
# reasonable time. It's slower per sample than the A100 (likely low
# hours instead of ~15-20 minutes), but runs on the `genomics` CPU
# partition, which has had zero contention issues anywhere else in this
# project (cellranger cluster mode, seqtk) even at high concurrency.
#
# Usage is identical to the GPU version:
#   Rscript cellbender_scripts/make_cellbender_params.R   # writes cell_quantities.csv
#   bash cellbender_scripts/run_all_cellbender_cpu.sh      # self-submits the array job

REPO="/projects/b1169/boles/img_scfrp"
CSV="$REPO/cellbender_scripts/cell_quantities.csv"

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    if [ ! -f "$CSV" ]; then
        echo "Missing $CSV -- run make_cellbender_params.R first." >&2
        exit 1
    fi
    n=$(( $(wc -l < "$CSV") - 1 ))
    echo "Submitting ${n} CellBender (CPU) jobs (1-${n})"
    # %20 -- genomics has handled comparable or higher concurrency
    # elsewhere in this repo without issue; adjust to taste.
    exec sbatch --array=1-${n}%20 "$0"
fi

module purge
module load python-miniconda3/4.10.3
module load mamba
source activate "$REPO/envs/cellbender"

echo "Running on host: $(hostname), ${SLURM_CPUS_PER_TASK} CPUs allocated"

# Same UTF-8 fix as the GPU script: without it, CellBender's report step
# fails at the end with UnicodeDecodeError since this cluster's default
# locale is ASCII, not UTF-8. Cosmetic only -- the real outputs are
# already written by that point -- but avoids ending on a traceback.
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CSV")
IFS=',' read -r POOL_DIR SAMPLE CELLS DROPLETS <<< "$LINE"

if ! [[ "$CELLS" =~ ^[0-9]+$ && "$DROPLETS" =~ ^[0-9]+$ && "$CELLS" -gt 0 && "$DROPLETS" -gt "$CELLS" ]]; then
    echo "Bad row for task ${SLURM_ARRAY_TASK_ID}: '${LINE}' (pool_dir=${POOL_DIR} sample=${SAMPLE} cells=${CELLS} droplets=${DROPLETS})" >&2
    exit 1
fi

OUT_DIR="$REPO/cellbender/${POOL_DIR}/${SAMPLE}"
mkdir -p "$OUT_DIR"

# Same checkpoint-locking reason as the GPU script: cd into a unique
# directory before running so concurrent tasks can't collide on
# ckpt.tar.gz.
cd "$OUT_DIR"

echo "Task ${SLURM_ARRAY_TASK_ID}: pool=${POOL_DIR} sample=${SAMPLE} cells=${CELLS} droplets=${DROPLETS}"

cellbender remove-background \
--cpu-threads "${SLURM_CPUS_PER_TASK}" \
--input "${REPO}/cellranger/${POOL_DIR}/outs/per_sample_outs/${SAMPLE}/sample_raw_feature_bc_matrix.h5" \
--output "${OUT_DIR}/${SAMPLE}" \
--expected-cells "$CELLS" \
--total-droplets-included "$DROPLETS" \
--fpr 0.01 \
--epochs 100
