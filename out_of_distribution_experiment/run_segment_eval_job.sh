#!/bin/bash
# Run segment eval (track building + evaluation) for one data directory.
# Called by segment_eval.sub with the data dir name as argument.

set -e

DATA_DIR="/data/alice/bkuipers/low_pt_gnn_pipeline/${1}"   # e.g. data_50
EXPERIMENT_DIR="/data/alice/bkuipers/out_of_distribution_experiment"
OOD_CONFIG="${EXPERIMENT_DIR}/ood_experiment.yaml"

echo "=========================================="
echo "Segment Eval: $(basename ${DATA_DIR})"
echo "=========================================="
echo "Started: $(date)"
echo "Node:    $(hostname -f)"
echo "DataDir: ${DATA_DIR}"
echo ""

source /data/alice/bkuipers/miniconda3/etc/profile.d/conda.sh
conda activate acorn

python "${EXPERIMENT_DIR}/run_sweep_segment_eval.py" \
    --data-dir "${DATA_DIR}" \
    --ood-config "${OOD_CONFIG}"

echo ""
echo "=========================================="
echo "${DATA_DIR_NAME} done: $(date)"
echo "=========================================="
