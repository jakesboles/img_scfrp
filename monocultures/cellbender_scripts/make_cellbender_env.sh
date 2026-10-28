#!/bin/bash

# Creates a conda env for running CellBender, under this project's own
# directory rather than a colleague's. Run this interactively on a GPU node
# (pip's cellbender install can pick a CPU-only PyTorch build if no GPU is
# visible at install time, so this keeps that part of Thomas's original
# script) -- but only briefly: env creation itself is a few minutes of
# single-threaded pip/mamba work, not 2 hours of 4-core/50G compute, so the
# allocation below is shrunk accordingly. If you already know a CPU-only
# build resolves fine for pip install cellbender, you can drop the srun
# line entirely and run the rest on a login node instead.

module purge
module load python-miniconda3/4.10.3
module load mamba

ENV_DIR="/projects/b1169/boles/img_scfrp/envs/cellbender"
mkdir -p "$(dirname "$ENV_DIR")"

export PYTHONNOUSERSITE="literallyanyletters"

srun -p genomics-gpu -A b1042 --mem=16G --gres=gpu:a100:1 -N 1 -n 2 -t 0:30:00 --pty bash -l

mamba create --prefix "$ENV_DIR" python=3.7.12

# Activate environment
mamba activate "$ENV_DIR"

# Install ipykernel
mamba install ipykernel

pip install cellbender

python -m ipykernel install --user --name cellbender --display-name "CellBender"

# Generate yml -- NOTE: --prefix here must be the environment path itself
# ($ENV_DIR), not its parent directory, and the redirect target must be an
# absolute path (not a bare filename, whose destination depends on
# whatever directory you happened to be in). Thomas's original script did
# both wrong -- it pointed --prefix at the parent project folder instead of
# env/Cellbender2, which silently produced a yml with an empty
# dependencies list (see CellbenderEnv2/Cellbender2.yml's header comment).
YML_OUT="/projects/b1169/boles/img_scfrp/cellbender_scripts/cellbender_env.yml"
conda env export --prefix "$ENV_DIR" > "$YML_OUT"

echo "Wrote $YML_OUT -- verify it has a real dependencies: list (not just"
echo "name/channels/prefix) before trusting it."
