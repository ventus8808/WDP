#!/bin/bash
# Slurm array launcher for Delta_bayesian.R - EQI change vs cancer mortality change analysis
# One task per Cancer_Type; each task runs both Lag=5 and Lag=10 analyses
# Usage:
#   bash Code/brms/submit_Delta_bayesian_array.sh         # auto-discovers diseases and submits array
#   # or, advanced: sbatch --array=0-<N-1> Code/brms/submit_Delta_bayesian_array.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_delta_bayes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=delta_bayesian_%A_%a.out
#SBATCH --error=delta_bayesian_%A_%a.err

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

DATA_CSV="Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
RUNNER="Code/brms/Delta_bayesian.R"

if [ ! -f "$DATA_CSV" ]; then
  log ERROR "找不到数据文件: $DATA_CSV"; exit 1
fi
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: if not running as an array worker, discover diseases and submit array
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  CANCER_LIST_FILE="cancer_types_delta.list"
  log INFO "发现所有 Cancer_Type 并生成任务列表: $CANCER_LIST_FILE"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
  dt <- fread(path, select = "Cancer_Type")
  u <- sort(unique(dt$Cancer_Type))
  if (length(u) == 0) stop("No Cancer_Type found in ", path)
  writeLines(u, "cancer_types_delta.list")
  cat(length(u), "cancer types found\n")
RS
  N=$(wc -l < "$CANCER_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "未找到任何Cancer_Type"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个疾病)"
  
  # Export list path and env name to workers
  sbatch --array=0-$((N-1)) \
    --export=ALL,CANCERS_FILE="$PROJECT_ROOT/$CANCER_LIST_FILE",ENV_NAME="$ENV_NAME" \
         "$0"
  log INFO "提交完成。每个任务分析一个癌症类型的Lag=5和Lag=10。"
  log INFO "使用 squeue 查看进度，结果保存至 Result/brms_delta/{ICD_Code}_delta.csv"
  exit 0
fi

# Worker mode (inside Slurm allocation)
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  log ERROR "SLURM_ARRAY_TASK_ID 未设置；请用 bash 直接运行脚本让其自提交"
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID}
CANCERS_FILE=${CANCERS_FILE:-"$PROJECT_ROOT/cancer_types_delta.list"}

if [ ! -f "$CANCERS_FILE" ]; then
  log WARN "未发现 CANCERS_FILE，在线生成列表"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
  dt <- fread(path, select = "Cancer_Type")
  u <- sort(unique(dt$Cancer_Type))
  if (length(u) == 0) stop("No Cancer_Type found in ", path)
  writeLines(u, "cancer_types_delta.list")
RS
  CANCERS_FILE="$PROJECT_ROOT/cancer_types_delta.list"
fi

if ! CANCER=$(sed -n "$((task_id+1))p" "$CANCERS_FILE"); then
  log ERROR "读取疾病列表失败 (index=$task_id)"; exit 1
fi
if [ -z "$CANCER" ]; then
  log WARN "索引 $task_id 超出疾病列表范围，跳过"; exit 0
fi

log INFO "==================================================================="
log INFO "任务ID=$task_id  疾病=$CANCER  CPU=${SLURM_CPUS_PER_TASK:-NA}"
log INFO "==================================================================="

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use different seed per task
SEED=$((1234 + task_id))

# Run Delta analysis for this cancer (both Lag=5 and Lag=10)
# Output: Result/brms_delta/{ICD_Code}_delta.csv
Rscript "$RUNNER" \
  --cancer-types "$CANCER" \
  --chains 4 \
  --iter 2000 \
  --warmup 1000 \
  --adapt-delta 0.95 \
  --max-treedepth 12 \
  --min-n 50 \
  --min-n-rucc 30 \
  --seed "$SEED"

log INFO "✅ 完成: $CANCER (输出: Result/brms_delta/$(echo $CANCER | sed 's/^delta_AAMR_//')_delta.csv)"
