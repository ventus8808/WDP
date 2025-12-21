#!/bin/bash
# Single task submit script for Delta_bayesian_ridgeline.R (Test Mode)
# Runs a single scenario: C00_C97, Lag 5, Overall Model, K=3, Cluster=1
# Usage: sbatch Code/brms/submit_Delta_bayesian_ridgeline_test.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_delta_ridge_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24
#SBATCH --time=01:00:00
#SBATCH --output=delta_ridge_test_%j.out
#SBATCH --error=delta_ridge_test_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Locate project root ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
  if [ -f "config.yaml" ]; then PROJECT_ROOT="$(pwd -P)"; fi
fi

if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "Cannot determine project root (config.yaml not found)"; exit 1
fi
cd "$PROJECT_ROOT"
log INFO "Project Root: $PROJECT_ROOT"

# --- Activate conda environment ---
ENV_NAME="${ENV_NAME:-brms}"
set +u
if [ -z "${CONDA_DEFAULT_ENV-}" ] || [ "${CONDA_DEFAULT_ENV}" != "$ENV_NAME" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "Conda initialization script not found"; exit 1
  fi
  conda activate "$ENV_NAME" || { log ERROR "Failed to activate conda env: $ENV_NAME"; exit 1; }
fi
set -u

RUNNER="Code/brms/Delta_bayesian_ridgeline.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "R script not found: $RUNNER"; exit 1
fi

# --- Configuration for Test ---
CANCER="C00_C97"
LAG=5
MODEL="overall"
K=3
CLUSTER_ID=1

log INFO "==================================================================="
log INFO "Test Run: $CANCER | Lag $LAG | Model $MODEL | K=$K | Cluster=$CLUSTER_ID"
log INFO "==================================================================="

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export TBB_CXX_TYPE=gcc

# Run with --test flag (reduces iterations to 800/400)
Rscript "$RUNNER" \
  --cancer "$CANCER" \
  --lag "$LAG" \
  --model "$MODEL" \
  --k "$K" \
  --cluster-id "$CLUSTER_ID" \
  --chains 4 \
  --iter 800 \
  --warmup 400 \
  --test

log INFO "✅ Test run complete. Check Result/Ridgeline/ for output."
