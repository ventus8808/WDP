#!/bin/bash
# Slurm launcher for the cmdstan main interval-censored mixed model pipeline
# Single-cancer runner restricted to C00_C97, keeping other settings identical to submit_cmdstan_main.sh.
#
# Usage:
#   sbatch Code/brms/submit_cmdstan_main_C00_C97_single.sh
#   # Optional overrides:
#   # ENV_NAME=brms CANCER_TYPE=C00_C97 SEED=20242 sbatch Code/brms/submit_cmdstan_main_C00_C97_single.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_cmdstan_C00_C97
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_main_C00_C97_%j.out
#SBATCH --error=cmdstan_main_C00_C97_%j.err

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

# --- Activate conda environment (default: brms; override via ENV_NAME) ---
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
    log ERROR "找不到conda初始化脚本"; exit 1
  fi
  conda activate "$ENV_NAME" || { log ERROR "激活conda环境失败: $ENV_NAME"; exit 1; }
fi
set -u

# Load devtoolset for newer g++ on CentOS
module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"

# Set environment variables for CmdStan
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/cmdstan_main.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# --- Single-cancer configuration (C00_C97 by default) ---
CANCER_TYPE="${CANCER_TYPE:-C00_C97}"
if [ -z "$CANCER_TYPE" ]; then
  log ERROR "CANCER_TYPE 为空"; exit 1
fi
log INFO "仅运行癌症类型: $CANCER_TYPE"

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use a fixed/override-able seed
SEED="${SEED:-1234}"

# Run the main interval-censored pipeline for this cancer type
# This runs Overall + RUCC1-4 stratifications internally
log INFO "启动 R 脚本: $RUNNER"
Rscript "$RUNNER" \
  --cancer-types "$CANCER_TYPE" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "✅ 完成: Cancer=$CANCER_TYPE"
