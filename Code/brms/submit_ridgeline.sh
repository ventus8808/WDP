#!/bin/bash
# Slurm array launcher for ridgeline posterior extraction (C00_C97 only)
# Each task runs one lag scenario (5, 10, or 15 years) with both Overall and Multi-domain models
# Usage:
#   sbatch --array=1-3 Code/brms/submit_ridgeline.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_ridgeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=ridgeline_%A_%a.out
#SBATCH --error=ridgeline_%A_%a.err

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

RUNNER="Code/brms/cmdstan_main_ridgeline.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# --- Map array task ID to scenario ---
TASK_ID=${SLURM_ARRAY_TASK_ID}
case $TASK_ID in
  1) SCENARIO=1; DESC="Lag5" ;;
  2) SCENARIO=2; DESC="Lag10" ;;
  3) SCENARIO=3; DESC="Lag15" ;;
  *)
    log ERROR "无效的 SLURM_ARRAY_TASK_ID: $TASK_ID (应为 1-3)"; exit 1
    ;;
esac

log INFO "开始处理任务 $TASK_ID: $DESC"

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use a different seed per task for better chain jitter
SEED=$((1234 + TASK_ID))

CANCER="C00_C97"
CHAINS=4
ITER=2000
WARMUP=1000
ADAPT_DELTA=0.95
MAX_TREEDEPTH=12

# --- Run Overall model ---
log INFO "运行 Overall EQI 模型 (Scenario=$SCENARIO)"
Rscript "$RUNNER" \
  --scenario "$SCENARIO" \
  --model "overall" \
  --cancer "$CANCER" \
  --chains "$CHAINS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH" \
  --seed "$SEED"

if [ $? -ne 0 ]; then
  log ERROR "Overall 模型失败 (Scenario=$SCENARIO)"; exit 1
fi
log INFO "✓ Overall 模型完成"

# --- Run Multi-domain model ---
log INFO "运行 Multi-domain 模型 (Scenario=$SCENARIO)"
Rscript "$RUNNER" \
  --scenario "$SCENARIO" \
  --model "multi" \
  --cancer "$CANCER" \
  --chains "$CHAINS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH" \
  --seed "$SEED"

if [ $? -ne 0 ]; then
  log ERROR "Multi-domain 模型失败 (Scenario=$SCENARIO)"; exit 1
fi
log INFO "✓ Multi-domain 模型完成"

log INFO "✅ 完成任务 $TASK_ID: $DESC (Overall + MultiDomain)"
