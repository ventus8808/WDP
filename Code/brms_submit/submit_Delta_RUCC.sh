#!/bin/bash
# Slurm array launcher for the RUCC-stratified delta mixed model pipeline.
# One task per outcome; each task runs all lag × RUCC × domain models inside the R runner.
# Usage:
#   bash Code/brms_submit/submit_Delta_RUCC.sh

#SBATCH --partition=wzhctest
#SBATCH --job-name=WDP_Delta_RUCC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=0-12:00:00
#SBATCH --output=Delta_RUCC_%A_%a.out
#SBATCH --error=Delta_RUCC_%A_%a.err

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

export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/Delta_RUCC.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: submit array from outcome_overall.list, then exit.
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  OUTCOME_LIST_FILE="outcome_overall.list"
  if [ ! -f "$OUTCOME_LIST_FILE" ]; then
    log ERROR "outcome_overall.list not found in project root"; exit 1
  fi
  N=$(wc -l < "$OUTCOME_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "outcome.list is empty"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个outcomes)"
  sbatch --array=0-$((N-1)) \
    --export=ALL,OUTCOME_FILE="$PROJECT_ROOT/$OUTCOME_LIST_FILE",ENV_NAME="$ENV_NAME" \
         "$0"
  log INFO "提交完成。使用 squeue 查看进度。"
  exit 0
fi

# Worker mode: run the actual job
log INFO "开始处理任务 $SLURM_ARRAY_TASK_ID"
OUTCOME=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$OUTCOME_FILE")
if [ -z "$OUTCOME" ]; then log ERROR "无法读取任务 $SLURM_ARRAY_TASK_ID 的outcome"; exit 1; fi
log INFO "处理outcome: $OUTCOME"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

Rscript "$RUNNER" \
  --outcomes "$OUTCOME"

log INFO "✅ 完成: Outcome=$OUTCOME"
