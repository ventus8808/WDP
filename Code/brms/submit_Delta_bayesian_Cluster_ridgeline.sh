#!/bin/bash
# Slurm array launcher for Delta cluster ridgeline RDS generation
# Runs 4 tasks for C00_C97 with k=3:
#   1: national
#   2: cluster 0
#   3: cluster 1
#   4: cluster 2
# Each task fits EQI change and single-domain change models and saves RDS for ridgeline plots.
#
# Usage:
#   sbatch --array=1-4 Code/brms/submit_Delta_bayesian_Cluster_ridgeline.sh
#
# Optional environment overrides:
#   ENV_NAME=brms     # conda env name
#   CHAINS=4 ITER=2000 WARMUP=1000 ADAPT_DELTA=0.95 MAX_TREEDEPTH=12
#   OUTPUT_DIR="Result/Delta_Ridgeline"

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_delta_ridge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=delta_ridge_%A_%a.out
#SBATCH --error=delta_ridge_%A_%a.err

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

# Load devtoolset for newer g++ on CentOS (if available)
module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"

# Set environment variables for CmdStan
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/Delta_bayesian_Cluster_ridgeline.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# --- Map array task ID to cluster selection ---
TASK_ID=${SLURM_ARRAY_TASK_ID}
case $TASK_ID in
  1) CLUSTER="national"; DESC="National" ;;
  2) CLUSTER="0";        DESC="Cluster 0" ;;
  3) CLUSTER="1";        DESC="Cluster 1" ;;
  4) CLUSTER="2";        DESC="Cluster 2" ;;
  *)
    log ERROR "无效的 SLURM_ARRAY_TASK_ID: $TASK_ID (应为 1-4)"; exit 1
    ;;
esac

# Parameters
CANCER="C00_C97"
K=3
CHAINS="${CHAINS:-4}"
ITER="${ITER:-2000}"
WARMUP="${WARMUP:-1000}"
ADAPT_DELTA="${ADAPT_DELTA:-0.95}"
MAX_TREEDEPTH="${MAX_TREEDEPTH:-12}"
OUTPUT_DIR="${OUTPUT_DIR:-Result/Delta_Ridgeline}"

log INFO "开始任务 $TASK_ID: $DESC | Cancer=$CANCER Lag=$LAG k=$K Cluster=$CLUSTER"

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use a different seed per task for better chain jitter
SEED=$((1234 + TASK_ID))

# Run
Rscript "$RUNNER" \
  --cancer "$CANCER" \
  --k "$K" \
  --cluster "$CLUSTER" \
  --output-dir "$OUTPUT_DIR" \
  --chains "$CHAINS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH" \
  --seed "$SEED"

RC=$?
if [ $RC -ne 0 ]; then
  log ERROR "任务失败: $DESC (退出码 $RC)"; exit $RC
fi

log INFO "✅ 完成任务: $DESC | 输出目录: $OUTPUT_DIR"
