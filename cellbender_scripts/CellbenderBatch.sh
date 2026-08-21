#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name CBDR.Batch
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 20GB
#SBATCH --time 1:00:00
#SBATCH --output /gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/logs/%x_oe%j.log
#SBATCH --verbose
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.watson1@northwestern.edu

# State date and name
date
echo 'Thomas Watson'

# Load python
module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

samples="/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/CellQuantities_012126.csv"

sample_arr=()
while IFS=',' read -ra row; do
sample_arr+=("${row[0]}")
done < $samples

sample_arr=("${sample_arr[@]:1}")

# For each sample...
counter=0
echo $"counter: $counter"
# for sample in "${sample_arr[0]}"
for sample in "${sample_arr[@]}"
do
# Update counter and print status
let counter++
  echo "Processing $counter: $sample now"

# Get output sample name, slide ID, and slide slot from metadata
cells=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $4}' "$samples")
echo $Cells
droplets=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $5}' "$samples")
echo $Droplets
pool=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $2}' "$samples")
pool="${pool%\"}"
pool="${pool#\"}"
echo $Pool
sample=$(awk -v sam="$sample" -F, '{ if ($1 == sam) print $3}' "$samples")
sample="${sample%\"}"
sample="${sample#\"}"
echo $sample

# Run spaceranger
sbatch /gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/CellbenderInd.sh \
$sample $cells $droplets $pool

echo "Submitted $counter: $sample"

# Sleep to give scheduler a break
sleep 3
done