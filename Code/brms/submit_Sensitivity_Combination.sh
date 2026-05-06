#!/bin/bash
# Slurm array launcher for sensitivity combination analysis.
# 7 tasks — one per overall disease outcome.
# Each task runs all 64 covariate combinations × 3 EQI-2000-2005 scenarios internally.
# Task layout (SLURM_ARRAY_TASK_ID → outcome):
#   0 → CLD     (K70_K76_C22)
#   1 → CRD     (J40_J47_J60_J70_J84_D86_C34)
#   2 → CKD     (N00_N29_C64_C65)
#   3 → CVD     (I00_I99)
#   4 → Suicide (X60_X84_Y87.0)
#   5 → NDD     (G20_G30_G12.2_F01_F03)
#   6 → Cancer  (C00_C97)
#
# Usage:
#   bash Code/brms/submit_Sensitivity_Combination.sh
#SBATCH --partition=kshctest02
#SBATCH --job-name=WDP_sens_comb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=cmdstan_sens_comb_%A_%a.out
#SBATCH --error=cmdstan_sens_comb_%A_%a.err
set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] $2"; }
# ---- Overall outcomes (one task each) ----
OUTCOMES=(
  "I00_I99"
  "C00_C97"
)
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
log INFO "Project root: $PROJECT_ROOT"
# ---- Activate micromamba env ----
ENV_NAME="${ENV_NAME:-brms}"
export MAMBA_EXE="$HOME/micromamba/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate "$ENV_NAME" || { log ERROR "Failed to activate micromamba env: $ENV_NAME"; exit 1; }
RUNNER="Code/brms/brms_Sensitivity_Combination.R"
if [ ! -f "$RUNNER" ]; then log ERROR "R script not found: $RUNNER"; exit 1; fi
# ---- Controller mode: submit array ----
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  N=${#OUTCOMES[@]}   # 7
  log INFO "Submitting ${N}-task array (one task per overall outcome)"
  log INFO "Each task runs 64 covariate combinations × 3 scenarios internally"
  ARRAY_JOB=$(sbatch --array=0-$((N-1)) \
    --export=ALL,ENV_NAME="$ENV_NAME" \
    --parsable \
    "$0")
  log INFO "Array job ID: $ARRAY_JOB"
  log INFO "Monitor with: squeue -u \$USER"
  log INFO "Output: Result/brms_Sensitivity_Combination/{CVD,NDD,CLD,CRD,CKD,Suicide,Cancer}_Sensitivity_Combination.csv"
  exit 0
fi
# ---- Worker mode ----
TASK_ID=$SLURM_ARRAY_TASK_ID
OUTCOME="${OUTCOMES[$TASK_ID]}"
if [ -z "$OUTCOME" ]; then
  log ERROR "Invalid task ID $TASK_ID (max index: $((${#OUTCOMES[@]}-1)))"; exit 1
fi
log INFO "Task $TASK_ID → Outcome=$OUTCOME"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
SEED=$(( 1234 + TASK_ID ))
micromamba run -n "$ENV_NAME" Rscript "$RUNNER" \
  --cancer-types "$OUTCOME" \
  --output-dir   "Result/brms_Sensitivity_Combination" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"
log INFO "Done: Task=$TASK_ID Outcome=$OUTCOME"
