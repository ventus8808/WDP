#!/bin/bash
# Slurm array launcher for the cmdstan stratified interval-censored mixed model pipeline
# One task per stratification combination; each task runs all scenarios and layers inside the R runner.
# Usage:
#   bash Code/brms/submit_cmdstan_stratified_array.sh         # auto-discovers stratifications and submits an array
#   # or, advanced: sbatch --array=0-<N-1> Code/brms/submit_cmdstan_stratified_array.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_cmdstan_stratified
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_stratified_%A_%a.out
#SBATCH --error=cmdstan_stratified_%A_%a.err

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

RUNNER="Code/brms/cmdstan_stratified.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: if not running as an array worker (either outside Slurm or a non-array sbatch),
# discover stratifications and submit an array, then exit.
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  TASKS_LIST_FILE="triangulation_stratified_tasks.list"
  log INFO "发现所有 Cancer_Type 组合并生成任务列表: $TASKS_LIST_FILE (位于项目根目录)"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Stratifed.csv"
  if (!file.exists(path)) stop("Unified stratified CSV not found: ", path)
  dt <- fread(path, select=c("Cancer_Type"))
  if (!"Cancer_Type" %in% names(dt)) stop("Cancer_Type column missing")
  u <- unique(dt$Cancer_Type)
  keep <- grepl("^(C00_C97|NDD)_(Male|Female|White|Black|Asian|Indian)$", u)
  u <- u[keep]
  if (length(u) == 0) stop("No Cancer_Type tasks discovered")
  writeLines(u, "triangulation_stratified_tasks.list")
  cat(length(u))
RS
  N=$(wc -l < "$TASKS_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "未找到任何任务"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个任务)"
  sbatch --array=0-$((N-1)) \
    --export=ALL,TASKS_FILE="$PROJECT_ROOT/$TASKS_LIST_FILE",ENV_NAME="$ENV_NAME" \
         "$0"
  log INFO "提交完成。使用 squeue 查看进度。"
  exit 0
fi

# Worker mode (inside Slurm allocation)
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  log ERROR "SLURM_ARRAY_TASK_ID 未设置；请用 bash 直接运行脚本让其自提交，或使用 --array 提交"
  exit 1
fi
task_id=${SLURM_ARRAY_TASK_ID}
TASKS_FILE=${TASKS_FILE:-"$PROJECT_ROOT/triangulation_stratified_tasks.list"}
if [ ! -f "$TASKS_FILE" ]; then
  log WARN "未发现 TASKS_FILE=$TASKS_FILE，回退到在线生成列表 (写入项目根目录)"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  path <- "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Stratifed.csv"
  if (!file.exists(path)) stop("Unified stratified CSV not found: ", path)
  dt <- fread(path, select=c("Cancer_Type"))
  if (!"Cancer_Type" %in% names(dt)) stop("Cancer_Type column missing")
  u <- unique(dt$Cancer_Type)
  keep <- grepl("^(C00_C97|NDD)_(Male|Female|White|Black|Asian|Indian)$", u)
  u <- u[keep]
  if (length(u) == 0) stop("No Cancer_Type tasks discovered")
  writeLines(u, "triangulation_stratified_tasks.list")
RS
  TASKS_FILE="$PROJECT_ROOT/triangulation_stratified_tasks.list"
fi

if ! LINE=$(sed -n "$((task_id+1))p" "$TASKS_FILE"); then
  log ERROR "读取任务列表失败 (index=$task_id)"; exit 1
fi
if [ -z "$LINE" ]; then
  log WARN "索引 $task_id 超出任务列表范围，跳过。"; exit 0
fi

CANCER_TYPE="$LINE"

log INFO "任务ID=$task_id  Cancer_Type=$CANCER_TYPE  CPU=${SLURM_CPUS_PER_TASK:-NA}"

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use a different seed per task for better chain jitter
SEED=$((1234 + task_id))

# Run the stratified interval-censored pipeline for this Cancer_Type
Rscript "$RUNNER" \
  --cancer-type "$CANCER_TYPE" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "✅ 完成: $CANCER_TYPE"
