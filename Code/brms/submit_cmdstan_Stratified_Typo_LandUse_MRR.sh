#!/bin/bash
# Slurm array launcher for cmdstan_Stratified_Typo_LandUse_MRR.R
# One task per base ICD code; each task runs:
#   - Typology MRR (Farming/Mining/Manufacturing/Government/Services/Nonspecialized)
#   - LandUse MRR  (Natural/Water_Sensitive/Agricultural/Urban)
#   - Stratified (sex/race) MRR, if the disease has strata in the stratified data
# Output directory: Result/brms_Stratified_Typo_LandUse_MRR/
#
# Usage:
#   bash Code/brms/submit_cmdstan_Stratified_Typo_LandUse_MRR.sh   # submit array

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Strat_Typo_LandUse_MRR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=cmdstan_Strat_Typo_LandUse_MRR_%A_%a.out
#SBATCH --error=cmdstan_Strat_Typo_LandUse_MRR_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Full disease list (one task per base ICD code) ---
CANCER_TYPES=(
  "C00_C97"
  "C18_C21"
  "C22"
  "C25"
  "C34"
  "C50"
  "C56"
  "C61"
  "C64_C65"
  "C82_C85"
  "C91_C95"
  "F01"
  "F03"
  "G10"
  "G12.2"
  "G20"
  "G20_G30_G12.2_F01_F03"
  "G30"
  "G30_F01_F03"
)

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
log INFO "Project root: $PROJECT_ROOT"

# --- Activate conda environment ---
ENV_NAME="${ENV_NAME:-brms}"
set +u
if [ -z "${CONDA_DEFAULT_ENV-}" ] || [ "${CONDA_DEFAULT_ENV}" != "$ENV_NAME" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "conda init script not found"; exit 1
  fi
  conda activate "$ENV_NAME" || { log ERROR "Failed to activate conda env: $ENV_NAME"; exit 1; }
fi
set -u

module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"

export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/cmdstan_Stratified_Typo_LandUse_MRR.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "R script not found: $RUNNER"; exit 1
fi

# --- Controller mode: submit array ---
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  N=${#CANCER_TYPES[@]}
  log INFO "Submitting array: 0-$((N-1)) ($N diseases)"
  for i in "${!CANCER_TYPES[@]}"; do
    log INFO "  [$i] ${CANCER_TYPES[$i]}"
  done
  CANCER_LIST=$(IFS=':'; echo "${CANCER_TYPES[*]}")
  sbatch --array=0-$((N-1)) \
    --export=ALL,CANCER_LIST="$CANCER_LIST",ENV_NAME="$ENV_NAME" \
    "$0"
  log INFO "Submitted. Use squeue to monitor progress."
  exit 0
fi

# --- Worker mode ---
log INFO "Starting task $SLURM_ARRAY_TASK_ID"

IFS=':' read -ra CANCERS <<< "$CANCER_LIST"
CANCER_TYPE="${CANCERS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$CANCER_TYPE" ]; then
  log ERROR "Cannot read cancer type for task $SLURM_ARRAY_TASK_ID"; exit 1
fi
log INFO "Cancer type: $CANCER_TYPE"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SEED=$((1234 + SLURM_ARRAY_TASK_ID))

Rscript "$RUNNER" \
  --cancer-type "$CANCER_TYPE" \
  --output-dir "Result/brms_Stratified_Typo_LandUse_MRR" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "Done: $CANCER_TYPE"
log INFO "Output: Result/brms_Stratified_Typo_LandUse_MRR/"
