#!/bin/bash

# Load miniconda
module purge
module load python-miniconda3/4.10.3
module load mamba

cd /projects/b1169/thomas/CellbenderEnv2

export PYTHONNOUSERSITE="literallyanyletters"

srun -p genomics-gpu -A b1042 --mem=50G --gres=gpu:a100:1 -N 1 -n 4 -t 2:00:00 --pty bash -l

mamba create --prefix ./env/Cellbender2 python=3.7.12

# Activate environment
mamba activate ./env/Cellbender2

# Install ipykernel
mamba install ipykernel

pip install cellbender

python -m ipykernel install --user --name Cellbender2 --display-name "Cellbender2"

# Generate yml
conda env export --prefix /projects/b1169/thomas/CellbenderEnv2 > Cellbender2.yml