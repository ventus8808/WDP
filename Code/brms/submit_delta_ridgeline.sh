#!/bin/bash
# Slurm array launcher for delta ridgeline posterior extraction
# Each task runs one combination of (cancer, k, lag, model)
# Usage:
#   sbatch --array=0-N Code/brms/submit_delta_ridgeline.sh
#   # or simply: bash Code/brms/submit_delta_ridgeline.sh (auto-submit mode)

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_delta_ridge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=delta_ridgeline_%A_%a.out
#SBATCH --error=delta_ridgeline_%A_%a.err

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

DATA_CSV="Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"
RUNNER="Code/brms/cmdstan_delta_ridgeline.R"

if [ ! -f "$DATA_CSV" ]; then
  log ERROR "找不到数据文件: $DATA_CSV"; exit 1
fi
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: if not running as an array worker, discover tasks and submit array
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  TASK_LIST_FILE="delta_ridgeline_tasks.list"
  log INFO "生成任务列表: $TASK_LIST_FILE"

  # Generate task combinations: cancer × k × lag × model
  # Focus on C00_C97 with k=3; lag=5,10,15; model=overall,multi
  cat > "$TASK_LIST_FILE" <<'TASKS'
C00_C97 3 5 overall
C00_C97 3 5 multi
C00_C97 3 10 overall
C00_C97 3 10 multi
C00_C97 3 15 overall
C00_C97 3 15 multi
TASKS

  N=$(wc -l < "$TASK_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "未生成任何任务"; exit 1; fi

  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个任务)"
  log INFO "任务配置: Cancer=C00_C97, k=3, lag=5,10,15, model=overall,multi"

  # Export task list path and env name to workers
  sbatch --array=0-$((N-1)) \
    --export=ALL,TASKS_FILE="$PROJECT_ROOT/$TASK_LIST_FILE",ENV_NAME="$ENV_NAME" \
    "$0"

  log INFO "提交完成。每个任务生成一个 ridgeline RDS 文件。"
  log INFO "使用 squeue 查看进度，结果保存至 Result/Ridgeline_Delta/"
  exit 0
fi

# Worker mode (inside Slurm allocation)
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  log ERROR "SLURM_ARRAY_TASK_ID 未设置；请用 bash 直接运行脚本让其自提交"
  exit 1
fi

task_id=${SLURM_ARRAY_TASK_ID}
TASKS_FILE=${TASKS_FILE:-"$PROJECT_ROOT/delta_ridgeline_tasks.list"}

if [ ! -f "$TASKS_FILE" ]; then
  log ERROR "找不到任务列表文件: $TASKS_FILE"; exit 1
fi

# Read task configuration
if ! TASK_LINE=$(sed -n "$((task_id+1))p" "$TASKS_FILE"); then
  log ERROR "读取任务列表失败 (index=$task_id)"; exit 1
fi
if [ -z "$TASK_LINE" ]; then
  log WARN "索引 $task_id 超出任务列表范围，跳过"; exit 0
fi

# Parse task: CANCER K LAG MODEL
read -r CANCER K_VAL LAG_VAL MODEL <<< "$TASK_LINE"

if [ -z "$CANCER" ] || [ -z "$K_VAL" ] || [ -z "$LAG_VAL" ] || [ -z "$MODEL" ]; then
  log ERROR "任务配置解析失败: '$TASK_LINE'"; exit 1
fi

log INFO "==================================================================="
log INFO "任务ID=$task_id  Cancer=$CANCER  k=$K_VAL  Lag=$LAG_VAL  Model=$MODEL"
log INFO "==================================================================="

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use different seed per task for better chain jitter
SEED=$((1234 + task_id))

CHAINS=4
ITER=2000
WARMUP=1000
ADAPT_DELTA=0.95
MAX_TREEDEPTH=12

log INFO "运行 Delta Ridgeline 提取 (Cancer=$CANCER, k=$K_VAL, Lag=$LAG_VAL, Model=$MODEL)"

Rscript "$RUNNER" \
  --cancer "$CANCER" \
  --k "$K_VAL" \
  --lag "$LAG_VAL" \
  --model "$MODEL" \
  --chains "$CHAINS" \
  --iter "$ITER" \
  --warmup "$WARMUP" \
  --adapt-delta "$ADAPT_DELTA" \
  --max-treedepth "$MAX_TREEDEPTH" \
  --seed "$SEED"

if [ $? -ne 0 ]; then
  log ERROR "Delta Ridgeline 失败 (Cancer=$CANCER, k=$K_VAL, Lag=$LAG_VAL, Model=$MODEL)"; exit 1
fi

log INFO "✅ 完成任务 $task_id: $CANCER k=$K_VAL Lag=$LAG_VAL Model=$MODEL"
log INFO "输出文件: Result/Ridgeline_Delta/${CANCER}_k${K_VAL}_Lag${LAG_VAL}_${MODEL^^}.rds"
