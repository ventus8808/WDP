#!/bin/bash
# Slurm array launcher for cmdstan_Stratified_All_MRR.R
# 10 tasks: NDD only, one per stratification section
# bash Code/brms/submit_cmdstan_Stratified_All_MRR.sh
# Task list (index: disease section → output file):
#  0: NDD   Male          → G20_G30_G12.2_F01_F03_Male_MRR.csv
#  1: NDD   Female        → G20_G30_G12.2_F01_F03_Female_MRR.csv
#  2: NDD   White         → G20_G30_G12.2_F01_F03_White_MRR.csv
#  3: NDD   Black         → G20_G30_G12.2_F01_F03_Black_MRR.csv
#  4: NDD   Asian         → G20_G30_G12.2_F01_F03_Asian_MRR.csv
#  5: NDD   Indian        → G20_G30_G12.2_F01_F03_Indian_MRR.csv
#  6: NDD   Typology      → G20_G30_G12.2_F01_F03_Typology_MRR.csv
#  7: NDD   RUCC          → G20_G30_G12.2_F01_F03_RUCC_MRR.csv
#  8: NDD   Koppen        → G20_G30_G12.2_F01_F03_Koppen_MRR.csv
#  9: NDD   CensusRegion  → G20_G30_G12.2_F01_F03_CensusRegion_MRR.csv

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Strat_MRR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-12:00:00
#SBATCH --output=cmdstan_Strat_MRR_%A_%a.out
#SBATCH --error=cmdstan_Strat_MRR_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# Task list: "cancer_type:section"  (20 entries, index 0-19)
TASKS=(
  "G20_G30_G12.2_F01_F03:Male"
  "G20_G30_G12.2_F01_F03:Female"
  "G20_G30_G12.2_F01_F03:White"
  "G20_G30_G12.2_F01_F03:Black"
  "G20_G30_G12.2_F01_F03:Asian"
  "G20_G30_G12.2_F01_F03:Indian"
  "G20_G30_G12.2_F01_F03:Typology"
  "G20_G30_G12.2_F01_F03:RUCC"
  "G20_G30_G12.2_F01_F03:Koppen"
  "G20_G30_G12.2_F01_F03:CensusRegion"
)

# --- Locate project root ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
elif [ -f "config.yaml" ]; then
  PROJECT_ROOT="$(pwd -P)"
fi
if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "Cannot determine project root (config.yaml not found)"; exit 1
fi
cd "$PROJECT_ROOT"
log INFO "Project root: $PROJECT_ROOT"

# --- Conda ---
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

RUNNER="Code/brms/cmdstan_Stratified_All_MRR.R"
if [ ! -f "$RUNNER" ]; then log ERROR "R script not found: $RUNNER"; exit 1; fi

# --- Controller: submit array ---
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  N=${#TASKS[@]}
  log INFO "Submitting array: 0-$((N-1)) ($N tasks)"
  for i in "${!TASKS[@]}"; do
    IFS=':' read -r ct sec <<< "${TASKS[$i]}"
    log INFO "  [$i] $ct  section=$sec"
  done
  TASK_LIST=$(IFS='|'; echo "${TASKS[*]}")
  sbatch --array=0-$((N-1)) \
    --export=ALL,TASK_LIST="$TASK_LIST",ENV_NAME="$ENV_NAME" \
    "$0"
  log INFO "Submitted. Use squeue to monitor."
  exit 0
fi

# --- Worker ---
IFS='|' read -ra ALL_TASKS <<< "$TASK_LIST"
TASK="${ALL_TASKS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$TASK" ]; then
  log ERROR "Cannot read task for index $SLURM_ARRAY_TASK_ID"; exit 1
fi
IFS=':' read -r CANCER_TYPE SECTION <<< "$TASK"
log INFO "Task $SLURM_ARRAY_TASK_ID → cancer=$CANCER_TYPE  section=$SECTION"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SEED=$((1234 + SLURM_ARRAY_TASK_ID))

Rscript "$RUNNER" \
  --cancer-type "$CANCER_TYPE" \
  --section     "$SECTION" \
  --output-dir  "Result/Stratified_MRR" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "Done: $CANCER_TYPE / $SECTION → Result/Stratified_MRR/"
