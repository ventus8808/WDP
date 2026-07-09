#!/bin/bash
# Slurm array launcher for ANAL sensitivity-combination analysis.
# One task per overall outcome from outcome_overall.list; each task runs:
#   4 ANAL exposure metrics x 42 covariate combinations
# internally for pooled 10-year lag.
#
# Usage:
#   bash Code/brms_ANAL_submit/submit_ANAL_Sensitivity_Combination.sh
#
# Optional overrides:
#   ENV_NAME=brms
#   ANAL_EXPOSURES=popw_mean_rad,mean_rad,sol,lit_area_km2
#   ANAL_LAG=10
#   ANAL_MAX_COVARIATES=3
#   ANAL_OUTPUT_DIR=Result/brms_ANAL_Sensitivity_Combination

#SBATCH --partition=wzhctest
#SBATCH --job-name=ANAL_Sensitivity_Combination
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=2-00:00:00
#SBATCH --output=ANAL_Sensitivity_Combination_%A_%a.out
#SBATCH --error=ANAL_Sensitivity_Combination_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Locate project root (expects config.yaml there) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if [ -f "${SCRIPT_DIR}/../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
elif [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
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

# Set environment variables for CmdStan
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms_ANAL/ANAL_brms_Sensitivity_Combination.R"
DATA_PATH="Data/Processed/df_ANAL.csv"
OUTCOME_LIST_FILE="${OUTCOME_LIST_FILE:-outcome_overall.list}"
ANAL_EXPOSURES="${ANAL_EXPOSURES:-popw_mean_rad,mean_rad,sol,lit_area_km2}"
ANAL_LAG="${ANAL_LAG:-10}"
ANAL_MAX_COVARIATES="${ANAL_MAX_COVARIATES:-3}"
ANAL_OUTPUT_DIR="${ANAL_OUTPUT_DIR:-Result/brms_ANAL_Sensitivity_Combination}"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi
if [ ! -f "$DATA_PATH" ]; then
  log ERROR "找不到输入数据: $DATA_PATH"; exit 1
fi

# Controller mode: if not running as an array worker, submit one task per outcome.
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  if [ ! -f "$OUTCOME_LIST_FILE" ]; then
    log ERROR "$OUTCOME_LIST_FILE not found in project root"; exit 1
  fi
  N=$(grep -cve '^[[:space:]]*$' "$OUTCOME_LIST_FILE")
  if [ "$N" -le 0 ]; then log ERROR "$OUTCOME_LIST_FILE is empty"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个overall outcomes)"
  log INFO "数据源: $DATA_PATH"
  log INFO "输出目录: $ANAL_OUTPUT_DIR"
  log INFO "暴露指标: $ANAL_EXPOSURES"
  log INFO "Lag: $ANAL_LAG; max covariates: $ANAL_MAX_COVARIATES"
  sbatch --array=0-$((N-1)) \
    --export=ALL,OUTCOME_FILE="$PROJECT_ROOT/$OUTCOME_LIST_FILE",ENV_NAME="$ENV_NAME",ANAL_EXPOSURES="$ANAL_EXPOSURES",ANAL_LAG="$ANAL_LAG",ANAL_MAX_COVARIATES="$ANAL_MAX_COVARIATES",ANAL_OUTPUT_DIR="$ANAL_OUTPUT_DIR" \
    "$0"
  log INFO "提交完成。使用 squeue 查看进度。"
  exit 0
fi

# Worker mode: run the actual job
log INFO "开始处理任务 $SLURM_ARRAY_TASK_ID"
OUTCOME=$(grep -ve '^[[:space:]]*$' "$OUTCOME_FILE" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
if [ -z "$OUTCOME" ]; then log ERROR "无法读取任务 $SLURM_ARRAY_TASK_ID 的outcome"; exit 1; fi
log INFO "处理outcome: $OUTCOME"
log INFO "数据源: $DATA_PATH"
log INFO "输出目录: $ANAL_OUTPUT_DIR"
log INFO "暴露指标: $ANAL_EXPOSURES"
log INFO "Lag: $ANAL_LAG; max covariates: $ANAL_MAX_COVARIATES"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SEED=$((1234 + SLURM_ARRAY_TASK_ID))

Rscript "$RUNNER" \
  --data "$DATA_PATH" \
  --output-dir "$ANAL_OUTPUT_DIR" \
  --outcomes "$OUTCOME" \
  --exposures "$ANAL_EXPOSURES" \
  --lag "$ANAL_LAG" \
  --max-covariates "$ANAL_MAX_COVARIATES" \
  --chains "${SLURM_CPUS_PER_TASK:-6}" \
  --seed "$SEED" \
  --overwrite

log INFO "✅ 完成: Outcome=$OUTCOME"
