#!/bin/bash
# Test script for Delta_bayesian_Cluster analysis - flexible version
# Usage: bash Code/brms/test_Delta_bayesian_Cluster.sh [cancer_types] [k_values]
# Example: bash Code/brms/test_Delta_bayesian_Cluster.sh "C00_C97" "3"

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Delta_Cluster_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=test_Delta_bayesian_Cluster_%j.out
#SBATCH --error=test_Delta_bayesian_Cluster_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Parse command line arguments ---
CANCER_TYPES="${1:-C00_C97}"
K_VALUES="${2:-3}"

log INFO "测试参数: 癌症类型=$CANCER_TYPES, k值=$K_VALUES"

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
  log ERROR "无法确定项目根目录 (找不到 config.yaml)"; exit 1
fi
cd "$PROJECT_ROOT"
log INFO "项目根目录: $PROJECT_ROOT"

# --- Activate conda environment ---
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

# Load devtoolset for newer g++ on CentOS
module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"

# Set environment variables for CmdStan
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/Delta_bayesian_Cluster.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

log INFO "开始测试: Delta_bayesian_Cluster $CANCER_TYPES with k=$K_VALUES"

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use test mode with reduced iterations
SEED=1234

# Run the Delta_bayesian_Cluster analysis
Rscript "$RUNNER" \
  --cancer-types "$CANCER_TYPES" \
  --k "$K_VALUES" \
  --chains 4 --iter 800 --warmup 300 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --test \
  --seed "$SEED"

log INFO "✅ 测试完成: Delta_bayesian_Cluster $CANCER_TYPES with k=$K_VALUES"