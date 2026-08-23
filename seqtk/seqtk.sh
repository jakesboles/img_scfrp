#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name stk_img
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80GB
#SBATCH --time 8:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output /gpfs/projects/b1169/boles/img_scfrp/logs/%x_%j.log
#SBATCH --verbose
#SBATCH --array=1-5

module load seqtk

cd /gpfs/projects/b1169/boles/img_scfrp/

# 2. Get the specific file and fraction for THIS task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" seqtk/downsample_targets.txt)
FILE=$(echo $LINE | cut -d' ' -f1)
FRACTION=$(echo $LINE | cut -d' ' -f2)

echo "Task ${SLURM_ARRAY_TASK_ID}: Downsampling ${FILE} with fraction ${FRACTION}"

# 3. Run seqtk
# Note: We use the same seed -s100 across all tasks. 
# Because the R1 and R2 for a specific lane are processed with the same seed, 
# the same reads will be picked and pairing will be preserved.
seqtk sample -s100 "$FILE" "$FRACTION" | gzip > "downsampled_fastq/${FILE}"

echo "Done."