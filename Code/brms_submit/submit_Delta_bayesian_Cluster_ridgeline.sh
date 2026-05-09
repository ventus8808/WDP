#!/bin/bash
# Slurm submission script for Delta_bayesian_Cluster_ridgeline.R
# Specific run for C00_C97 (All Cancers) with k=3 (National + Clusters 0,1,2)
#
# Usage:
#   sbatch Code/brms/submit_Delta_bayesian_Cluster_ridgeline.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_ridgeline_C00_C97
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=ridgeline_C00_C97_%j.out
#SBATCH --error=ridgeline_C00_C97_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Locate project root (expects config.yaml there) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
  # Fallback to cwd if it contains config.yaml
  if [ -f "config.yaml" ]; then PROJECT_ROOT="$(pwd -P)"; fi
fi

if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "无法确定项目根目录 (找不到 config.yaml)"; exit 1
fi
cd "$PROJECT_ROOT"
log INFO "项目根目录: $PROJECT_ROOT"

# --- Activate conda environment (default: brms) ---
ENV_NAME="${ENV_NAME:-brms}"
set +u
if [ -z "${CONDA_DEFAULT_ENV-}" ] || [ "${CONDA_DEFAULT_ENV}" != "$ENV_NAME" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "找不到conda初始化脚本"; exit 1
  fi
  conda activate "$ENV_NAME" || { log ERROR "激活conda环境失败: $ENV_NAME"; exit 1; }
fi
set -u

RUNNER="Code/brms/Delta_bayesian_Cluster_ridgeline.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Configuration
CANCER_TYPE="C00_C97"
K_VALUE="3"

log INFO "==================================================================="
log INFO "任务: Ridgeline Data Generation"
log INFO "疾病: $CANCER_TYPE"
log INFO "K值: $K_VALUE (National + Clusters 0,1,2)"
log INFO "CPU: ${SLURM_CPUS_PER_TASK:-16}"
log INFO "==================================================================="

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
# Fix CmdStan compilation issue on non-standard compilers
export TBB_CXX_TYPE=gcc

# Run analysis
# Note: The R script automatically runs National analysis and then iterates through clusters for the given k
Rscript "$RUNNER" \
  --cancer-types "$CANCER_TYPE" \
  --k "$K_VALUE" \
  --chains 4 \
  --iter 2000 \
  --warmup 1000 \
  --adapt-delta 0.95 \
  --max-treedepth 12 \
  --seed 12345

log INFO "✅ 任务完成: $CANCER_TYPE (k=$K_VALUE)"
