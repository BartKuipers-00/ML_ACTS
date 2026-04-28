#!/bin/bash
# Evaluate one multiplicity-sweep dataset through the full GNN pipeline:
#   build latent graphs → GNN inference → segment matching + evaluation
# Argument: HTCondor $(Process) index (0-12), maps to particles_per_vertex 10-130

set -e

PROCESS=${1:-0}
PPV=$(( (PROCESS + 1) * 10 ))
DATA_DIR="/data/alice/bkuipers/low_pt_gnn_pipeline/data_${PPV}"
EXPERIMENT_DIR="/data/alice/bkuipers/out_of_distribution_experiment"
PIPELINE_DIR="/data/alice/bkuipers/low_pt_gnn_pipeline"
OOD_CONFIG="${EXPERIMENT_DIR}/ood_experiment.yaml"

echo "=========================================="
echo "Multiplicity Sweep Eval - particles_per_vertex=${PPV}"
echo "=========================================="
echo "Started: $(date)"
echo "Node:    $(hostname -f)"
echo "CPUs:    $(nproc)"
echo "DataDir: ${DATA_DIR}"
echo ""

source /data/alice/bkuipers/miniconda3/etc/profile.d/conda.sh
conda activate acorn

# Step 1: Build latent-space KNN graphs (feature_store → graph_constructed_latent)
echo "--- Step 1: Latent graph construction ---"
python "${EXPERIMENT_DIR}/run_sweep_build_graphs.py" \
    --data-dir "${DATA_DIR}" \
    --ood-config "${OOD_CONFIG}"

echo ""

# Step 2: GNN edge-score inference (graph_constructed_latent → gnn_stage)
echo "--- Step 2: GNN inference ---"
python "${EXPERIMENT_DIR}/run_sweep_infer_gnn.py" \
    --data-dir "${DATA_DIR}" \
    --ood-config "${OOD_CONFIG}"

echo ""

# Step 3: Track building + evaluation (gnn_stage → track_building + track_evaluation)
# Branch on use_segment_gnn in ood_experiment.yaml
USE_SEGMENT_GNN=$(python -c "import yaml; c=yaml.safe_load(open('${OOD_CONFIG}')); print(c.get('use_segment_gnn', False))")

if [ "${USE_SEGMENT_GNN}" = "True" ]; then
    echo "--- Step 3: GNN segment matching + evaluation ---"
    python "${PIPELINE_DIR}/track_building_stage_(3)/GNN_segmentmatcher_walkthrough.py" \
        testset --data-dir "${DATA_DIR}" --ood-config "${OOD_CONFIG}"
else
    echo "--- Step 3: Helix segment matching + evaluation ---"
    python "${PIPELINE_DIR}/track_building_stage_(3)/segment_matching_track_builder.py" \
        testset --data-dir "${DATA_DIR}" --ood-config "${OOD_CONFIG}"
fi

echo ""
echo "=========================================="
echo "PPV=${PPV} eval done: $(date)"
echo "=========================================="
