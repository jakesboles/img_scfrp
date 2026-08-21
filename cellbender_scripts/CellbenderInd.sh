#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --job-name CBDR.Ind
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 30GB
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/logs/%x_oe%j.log
#SBATCH --verbose
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.watson1@northwestern.edu

# $1 = sample
# $2 = cells
# $3 = droplets
# $4 = pool

# State date and name
date
echo 'Thomas Watson'

# Load python
module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

# Load environment
source activate "/projects/b1169/thomas/CellbenderEnv2/env/Cellbender2"

mkdir "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/$4"
mkdir "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/$4/$1"
mkdir "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/$4/$1/$1"

# DO THIS, or else............................................................................................. checkpoint files will lock other checkpoints out. One job per directory
cd "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/$4/$1/$1"

# Run script
cellbender remove-background \
--cuda \
--input "/projects/b1042/Gate_Lab/thomas/NeilSC/Coculture/$4/outs/per_sample_outs/$1/sample_raw_feature_bc_matrix.h5" \
--output "/gpfs/projects/p31535/thomas/NeilSC/Coculture_Cellbender/$4/$1/$1" \
--expected-cells "$2" \
--total-droplets-included "$3" \
--fpr 0.01 \
--epochs 100