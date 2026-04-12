#!/bin/bash
# Slurm array launcher for sensitivity covariate analysis
# 10 tasks: 5 model groups × 2 cancer types, each task runs all 3 lags internally
# Task layout (SLURM_ARRAY_TASK_ID):
#   0-4  → NDD    (G20_G30_G12.2_F01_F03) × [Baseline_Full, Insurance, Doctor, Forest, Monitoring]
#   5-9  → Cancer (C00_C97)               × [Baseline_Full, Insurance, Doctor, Forest, Monitoring]
# A merge job runs automatically after all 10 tasks complete.
#
# Usage:
#   bash Code/brms/submit_cmdstan_sensativity_covar.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_sens_covar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=0-12:00:00
#SBATCH --output=logs/sens_covar_%A_%a.out
#SBATCH --error=logs/sens_covar_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] $2"; }

# ---- Task mapping ----
CANCER_TYPES=("G20_G30_G12.2_F01_F03" "C00_C97")
MODEL_GROUPS=("Baseline_Full" "Insurance" "Doctor" "Forest" "Monitoring")
N_MODELS=${#MODEL_GROUPS[@]}   # 5

# ---- Locate project root ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if   [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then PROJECT_ROOT="$SLURM_SUBMIT_DIR"
elif [ -f "config.yaml" ]; then PROJECT_ROOT="$(pwd -P)"
fi
if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "Cannot determine project root (config.yaml not found)"; exit 1
fi
cd "$PROJECT_ROOT"
mkdir -p logs
log INFO "Project root: $PROJECT_ROOT"

# ---- Activate conda ----
ENV_NAME="${ENV_NAME:-brms}"
set +u
if [ -z "${CONDA_DEFAULT_ENV-}" ] || [ "${CONDA_DEFAULT_ENV}" != "$ENV_NAME" ]; then
  if   [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ];   then source "/opt/anaconda3/etc/profile.d/conda.sh"
  else log ERROR "conda init script not found"; exit 1
  fi
  conda activate "$ENV_NAME" || { log ERROR "Failed to activate conda env: $ENV_NAME"; exit 1; }
fi
set -u

module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/cmdstan_Sensativity_Covar.R"
MERGER="Code/brms/merge_sensativity_covar.R"
if [ ! -f "$RUNNER" ]; then log ERROR "R script not found: $RUNNER"; exit 1; fi

# ---- Controller mode: submit array + merge dependency ----
if [ -z "${SLURM_ARRAY_TASK_ID-}" ] && [ -z "${MERGE_MODE-}" ]; then
  TOTAL_TASKS=$(( ${#CANCER_TYPES[@]} * N_MODELS ))   # 10
  log INFO "Submitting array of ${TOTAL_TASKS} tasks (${#CANCER_TYPES[@]} cancers × ${N_MODELS} model groups)"

  ARRAY_JOB=$(sbatch --array=0-$((TOTAL_TASKS-1)) \
    --export=ALL,ENV_NAME="$ENV_NAME" \
    --parsable \
    "$0")
  log INFO "Array job ID: $ARRAY_JOB"

  # Submit merge job that runs after all array tasks succeed
  sbatch --dependency=afterok:"$ARRAY_JOB" \
    --job-name=WDP_sens_merge \
    --partition=kshctest \
    --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=0-00:30:00 \
    --output=logs/sens_merge_%j.out \
    --error=logs/sens_merge_%j.err \
    --export=ALL,ENV_NAME="$ENV_NAME",MERGE_MODE=1 \
    "$0"
  log INFO "Merge job submitted (runs after array completes)."
  log INFO "Monitor with: squeue -u \$USER"
  exit 0
fi

# ---- Merge mode ----
if [ -n "${MERGE_MODE-}" ]; then
  log INFO "Running merge step..."
  Rscript "$MERGER" --output-dir "Result/brms_Sensativity"
  log INFO "Merge complete."
  exit 0
fi

# ---- Worker mode ----
TASK_ID=$SLURM_ARRAY_TASK_ID
CANCER_IDX=$(( TASK_ID / N_MODELS ))
MODEL_IDX=$(( TASK_ID % N_MODELS ))

CANCER_TYPE="${CANCER_TYPES[$CANCER_IDX]}"
MODEL_GROUP="${MODEL_GROUPS[$MODEL_IDX]}"

if [ -z "$CANCER_TYPE" ] || [ -z "$MODEL_GROUP" ]; then
  log ERROR "Invalid task ID $TASK_ID"; exit 1
fi

log INFO "Task $TASK_ID → Cancer=$CANCER_TYPE | ModelGroup=$MODEL_GROUP"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Different seed per task for better chain jitter
SEED=$(( 1234 + TASK_ID ))

Rscript "$RUNNER" \
  --cancer-types  "$CANCER_TYPE" \
  --model-group   "$MODEL_GROUP" \
  --output-dir    "Result/brms_Sensativity" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "Done: Task=$TASK_ID Cancer=$CANCER_TYPE ModelGroup=$MODEL_GROUP"
