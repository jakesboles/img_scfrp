#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name seqtk_downsample
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 8G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/img_scfrp/logs/%x_%A_%a.log
#SBATCH --verbose
#
# --array is intentionally not set here: the number of tasks depends on how
# many FASTQ files downsample_manifest.txt lists, which varies as source
# data changes. Run seqtk/build_manifest.sh first, then submit with:
#   sbatch --array=1-N seqtk/seqtk.sh
# where N is the line count build_manifest.sh reports (add a throttle if
# you don't want all files downsampling at once, e.g. --array=1-N%20).

module load seqtk

REPO="/projects/b1169/boles/img_scfrp"
MANIFEST="$REPO/seqtk/downsample_manifest.txt"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
SRC=$(echo "$LINE" | cut -d',' -f1)
FRACTION=$(echo "$LINE" | cut -d',' -f2)
DEST=$(echo "$LINE" | cut -d',' -f3)

mkdir -p "$(dirname "$DEST")"

echo "Task ${SLURM_ARRAY_TASK_ID}: ${SRC} -> ${DEST} (fraction ${FRACTION})"

# Each FASTQ file (R1 and R2 are separate array tasks) is downsampled
# independently with the same fixed seed. seqtk's per-read sampling
# decision depends only on the read's position/content and the seed, so R1
# and R2 of the same lane -- sampled in separate tasks -- still pick the
# same reads and stay paired. The same seed being reused across different
# lanes/pools is fine too: pairing only depends on a lane's R1 and R2 seeing
# the same seed, not on every file in the manifest using a different one.
seqtk sample -s100 "$SRC" "$FRACTION" | gzip > "$DEST"

echo "Done."
