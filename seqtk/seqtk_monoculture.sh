#!/bin/bash
#SBATCH --account p31535
#SBATCH --partition normal
#SBATCH --job-name stk_img
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 80GB
#SBATCH --time 8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output /gpfs/projects/b1169/thomas/TDPC_PPA_Spatial/GSEA/logs/%x_oe%j.log
#SBATCH --verbose
#SBATCH --array=1-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.watson1@northwestern.edu

module load seqtk

cd /gpfs/projects/b1042/Gate_Lab/boles/img_scfrp/fastq

# 2. Setup directories
mkdir -p downsampled
mkdir -p logs

# 2. Get the specific file and fraction for THIS task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" batch_list_remaining.txt)
FILE=$(echo $LINE | cut -d' ' -f1)
FRACTION=$(echo $LINE | cut -d' ' -f2)

echo "Task ${SLURM_ARRAY_TASK_ID}: Downsampling ${FILE} with fraction ${FRACTION}"

# 3. Run seqtk
# Note: We use the same seed -s100 across all tasks. 
# Because the R1 and R2 for a specific lane are processed with the same seed, 
# the same reads will be picked and pairing will be preserved.
seqtk sample -s100 "$FILE" "$FRACTION" | gzip > "downsampled/${FILE}"

echo "Done."