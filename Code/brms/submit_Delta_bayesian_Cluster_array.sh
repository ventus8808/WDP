#!/bin/bash
# Slurm array launcher for Delta_bayesian_Cluster.R - EQI change vs cancer mortality change analysis with cluster stratification
# One task per Cancer_Type x K combination
# Usage:
#   bash Code/brms/submit_Delta_bayesian_Cluster_array.sh         # auto-discovers diseases and k values, submits array
#   # or, advanced: sbatch --array=0-<N-1> Code/brms/submit_Delta_bayesian_Cluster_array.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_delta_cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=delta_bayesian_cluster_%A_%a.out
#SBATCH --error=delta_bayesian_cluster_%A_%a.err

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
RUNNER="Code/brms/Delta_bayesian_Cluster.R"

if [ ! -f "$DATA_CSV" ]; then
  log ERROR "找不到数据文件: $DATA_CSV"; exit 1
fi
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: if not running as an array worker, discover diseases and k values, submit array
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  TASK_LIST_FILE="delta_cluster_tasks.list"
  log INFO "发现所有 Cancer_Type 和 k 值组合并生成任务列表: $TASK_LIST_FILE"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
  dt <- fread(path, select = "Cancer_Type")
  cancers <- sort(unique(dt$Cancer_Type))
  if (length(cancers) == 0) stop("No Cancer_Type found in ", path)
  
  # Define k values (clusters)
  k_values <- c(3, 4)
  
  # Generate all combinations
  tasks <- expand.grid(cancer = cancers, k = k_values, stringsAsFactors = FALSE)
  task_lines <- paste0(tasks$cancer, ",", tasks$k)
  
  writeLines(task_lines, "delta_cluster_tasks.list")
  cat(length(task_lines), "tasks generated (", length(cancers), " cancers x ", length(k_values), " k values)\n")
RS
  N=$(wc -l < "$TASK_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "未找到任何任务"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个任务)"
  
  # Export list path and env name to workers
  sbatch --array=0-$((N-1)) \
    --export=ALL,TASKS_FILE="$PROJECT_ROOT/$TASK_LIST_FILE",ENV_NAME="$ENV_NAME" \
         "$0"
  log INFO "提交完成。每个任务分析一个癌症类型和k值的组合。"
  log INFO "使用 squeue 查看进度，结果保存至 Result/brms_delta_cluster/"
  exit 0
fi

# Worker mode (inside Slurm allocation)
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  log ERROR "SLURM_ARRAY_TASK_ID 未设置；请用 bash 直接运行脚本让其自提交"
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID}
TASKS_FILE=${TASKS_FILE:-"$PROJECT_ROOT/delta_cluster_tasks.list"}

if [ ! -f "$TASKS_FILE" ]; then
  log WARN "未发现 TASKS_FILE，在线生成列表"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
  dt <- fread(path, select = "Cancer_Type")
  cancers <- sort(unique(dt$Cancer_Type))
  if (length(cancers) == 0) stop("No Cancer_Type found in ", path)
  
  k_values <- c(3, 4)
  tasks <- expand.grid(cancer = cancers, k = k_values, stringsAsFactors = FALSE)
  task_lines <- paste0(tasks$cancer, ",", tasks$k)
  writeLines(task_lines, "delta_cluster_tasks.list")
RS
  TASKS_FILE="$PROJECT_ROOT/delta_cluster_tasks.list"
fi

if ! TASK_LINE=$(sed -n "$((task_id+1))p" "$TASKS_FILE"); then
  log ERROR "读取任务列表失败 (index=$task_id)"; exit 1
fi
if [ -z "$TASK_LINE" ]; then
  log WARN "索引 $task_id 超出任务列表范围，跳过"; exit 0
fi

# Parse task line: CANCER,K
CANCER=$(echo "$TASK_LINE" | cut -d',' -f1)
K=$(echo "$TASK_LINE" | cut -d',' -f2)

if [ -z "$CANCER" ] || [ -z "$K" ]; then
  log ERROR "解析任务失败: $TASK_LINE"; exit 1
fi

log INFO "==================================================================="
log INFO "任务ID=$task_id  疾病=$CANCER  k=$K  CPU=${SLURM_CPUS_PER_TASK:-NA}"
log INFO "==================================================================="

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Fix CmdStan compilation issue on non-standard compilers
export TBB_CXX_TYPE=gcc
log INFO "Set TBB_CXX_TYPE=gcc for CmdStan compilation"

# Use different seed per task
SEED=$((1234 + task_id))

# Run Delta Cluster analysis for this cancer and k
# Output: Result/brms_delta_cluster/{ICD_Code}_k{K}_delta.csv
Rscript "$RUNNER" \
  --cancer-types "$CANCER" \
  --k "$K" \
  --chains 4 \
  --iter 2000 \
  --warmup 1000 \
  --adapt-delta 0.95 \
  --max-treedepth 12 \
  --seed "$SEED"

log INFO "✅ 完成: $CANCER (k=$K) (输出: Result/brms_delta_cluster/$(echo $CANCER | sed 's/^delta_AAMR_//')_k${K}_delta.csv)"