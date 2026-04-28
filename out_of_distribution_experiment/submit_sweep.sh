#!/bin/bash
# Submit simulation or eval sweep jobs.
# Reads max_particles_per_vertex from ood_experiment.yaml to set queue size.
#
# Usage:
#   ./submit_sweep.sh sim     # simulate + clean + convert
#   ./submit_sweep.sh eval    # build latent graphs + GNN inference + track building (old path)
#   ./submit_sweep.sh seg     # segment eval only (track building + evaluation) for all data_* dirs

set -e

EXPERIMENT_DIR="$(cd "$(dirname "$0")" && pwd)"
OOD_CONFIG="${EXPERIMENT_DIR}/ood_experiment.yaml"

N_JOBS=$(python3 -c "
import yaml, math
c = yaml.safe_load(open('${OOD_CONFIG}'))
print(math.ceil(c['max_particles_per_vertex'] / 10))
")

echo "max_particles_per_vertex → ${N_JOBS} jobs (ppv = 10, 20, ..., $((N_JOBS * 10)))"

case "${1}" in
    sim)
        echo "Submitting simulation sweep (${N_JOBS} jobs)..."
        condor_submit "${EXPERIMENT_DIR}/multiplicity_sweep.sub" n_jobs="${N_JOBS}"
        ;;
    eval)
        echo "Submitting eval sweep (${N_JOBS} jobs)..."
        condor_submit "${EXPERIMENT_DIR}/multiplicity_sweep_eval.sub" n_jobs="${N_JOBS}"
        ;;
    seg)
        PIPELINE_ROOT="$(dirname "${EXPERIMENT_DIR}")/low_pt_gnn_pipeline"
        DATA_DIRS_FILE="${EXPERIMENT_DIR}/data_dirs.txt"
        ls -d "${PIPELINE_ROOT}/data_"[0-9]*/ 2>/dev/null \
            | sort -t_ -k2 -n \
            | xargs -I{} basename {} > "${DATA_DIRS_FILE}"
        N_SEG=$(wc -l < "${DATA_DIRS_FILE}")
        if [ "${N_SEG}" -eq 0 ]; then
            echo "No data_* directories found in ${EXPERIMENT_DIR}"
            exit 1
        fi
        echo "Found ${N_SEG} data directories:"
        cat "${DATA_DIRS_FILE}"
        echo ""
        echo "Submitting segment eval (${N_SEG} jobs)..."
        condor_submit "${EXPERIMENT_DIR}/segment_eval.sub"
        ;;
    *)
        echo "Usage: $0 sim|eval|seg"
        exit 1
        ;;
esac
