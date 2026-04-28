#!/bin/bash
# Convert one sweep dataset CSV → PyG using the acorn conda environment.
# Called by run_multiplicity_sweep_job.sh via env -i to ensure a clean
# environment with no ACTS/LCG variables bleeding in.
# Argument: absolute path to data_N directory

set -e

DATA_DIR="$1"

echo "--- Step 3: CSV → PyG (acorn env) ---"
echo "DataDir: ${DATA_DIR}"

source /data/alice/bkuipers/miniconda3/etc/profile.d/conda.sh
conda activate acorn

cd /data/alice/bkuipers/low_pt_gnn_pipeline

python latent_stage_\(1\)/run_sweep_convert.py \
    --data-dir "${DATA_DIR}" \
    --n-events 500

echo "--- Step 3 done ---"
