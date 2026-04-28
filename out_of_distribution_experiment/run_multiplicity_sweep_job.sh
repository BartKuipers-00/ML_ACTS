#!/bin/bash
# Run one multiplicity-sweep job: simulate → clean → convert to PyG
# Argument: HTCondor $(Process) index (0-12), maps to particles_per_vertex 10-130

set -e

PROCESS=${1:-0}
PPV=$(( (PROCESS + 1) * 10 ))
DATA_DIR="/data/alice/bkuipers/low_pt_gnn_pipeline/data_${PPV}"

echo "=========================================="
echo "Multiplicity Sweep - particles_per_vertex=${PPV}"
echo "=========================================="
echo "Started: $(date)"
echo "Node:    $(hostname -f)"
echo "CPUs:    $(nproc)"
echo "DataDir: ${DATA_DIR}"
echo ""

# Steps 1+2: ACTS environment (LCG + ACTS Python bindings)
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
source /data/alice/bkuipers/acts/build/python/setup.sh

cd /data/alice/bkuipers/low_pt_gnn_pipeline

# Step 1: ACTS simulation (500 events, ppv=${PPV}, output → data_${PPV}/)
echo "--- Step 1: ACTS simulation ---"
python simulation_\(0\)/run_simulation_sweep.py --particles-per-vertex ${PPV}

echo ""

# Step 2: Clean loops and attribute segments
echo "--- Step 2: Cleaning CSV ---"
python simulation_\(0\)/clean_loops_and_attribute_segments.py \
    --data-dir "${DATA_DIR}"

echo ""

# Step 3: CSV → PyG — run in a clean environment (no ACTS/LCG vars)
echo "--- Step 3: CSV → PyG ---"
env -i HOME="${HOME}" \
    bash /data/alice/bkuipers/jobs/run_multiplicity_sweep_convert.sh "${DATA_DIR}"

echo ""
echo "=========================================="
echo "PPV=${PPV} done: $(date)"
echo "=========================================="
