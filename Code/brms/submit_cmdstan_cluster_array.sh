#!/bin/bash
# Slurm array launcher for the cmdstan cluster interval-censored mixed model pipeline
# One task per cluster; each task runs all scenarios inside the R runner.
# Usage:
#   bash Code/brms/submit_cmdstan_cluster_array.sh         # auto-discovers clusters and submits an array
#   # or, advanced: sbatch --array=0-<N-1> Code/brms/submit_cmdstan_cluster_array.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_cmdstan_cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_cluster_%A_%a.out
#SBATCH --error=cmdstan_cluster_%A_%a.err

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

RUNNER="Code/brms/cmdstan_cluster.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Controller mode: if not running as an array worker (either outside Slurm or a non-array sbatch),
# discover clusters and submit an array, then exit.
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  CLUSTER_LIST_FILE="clusters.list"
  log INFO "发现所有cluster并生成任务列表: $CLUSTER_LIST_FILE (位于项目根目录)"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  # Load clustered data
  cluster_path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval_Clustered.csv"
  if (file.exists(cluster_path)) {
    dt <- fread(cluster_path, select = "Cluster")
    u <- unique(dt[!is.na(Cluster), .(Cluster)])
    u <- u[order(Cluster)]
  } else {
    stop("Clustered data not found")
  }
  if (nrow(u) == 0) stop("No clusters found")
  # Write cluster IDs
  writeLines(as.character(u$Cluster), "clusters.list")
  cat(nrow(u))
RS
  N=$(wc -l < "$CLUSTER_LIST_FILE" | tr -d ' ')
  if [ "$N" -le 0 ]; then log ERROR "未找到任何cluster"; exit 1; fi
  log INFO "将提交数组任务: 0-$((N-1)) (共 $N 个cluster)"
  # Export list path and env name to workers
  sbatch --array=0-$((N-1)) \
    --export=ALL,CLUSTER_FILE="$PROJECT_ROOT/$CLUSTER_LIST_FILE",ENV_NAME="$ENV_NAME" \
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
CLUSTER_FILE=${CLUSTER_FILE:-"$PROJECT_ROOT/clusters.list"}
if [ ! -f "$CLUSTER_FILE" ]; then
  log WARN "未发现 CLUSTER_FILE=$CLUSTER_FILE，回退到在线生成列表 (写入项目根目录)"
  Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  # Load clustered data
  cluster_path <- "Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval_Clustered.csv"
  if (file.exists(cluster_path)) {
    dt <- fread(cluster_path, select = "Cluster")
    u <- unique(dt[!is.na(Cluster), .(Cluster)])
    u <- u[order(Cluster)]
  } else {
    stop("Clustered data not found")
  }
  if (nrow(u) == 0) stop("No clusters found")
  # Write cluster IDs
  writeLines(as.character(u$Cluster), "clusters.list")
RS
  CLUSTER_FILE="$PROJECT_ROOT/clusters.list"
fi

if ! CLUSTER_ID=$(sed -n "$((task_id+1))p" "$CLUSTER_FILE"); then
  log ERROR "读取cluster列表失败 (index=$task_id)"; exit 1
fi
if [ -z "$CLUSTER_ID" ]; then
  log WARN "索引 $task_id 超出cluster列表范围，跳过。"; exit 0
fi

log INFO "任务ID=$task_id  Cluster=$CLUSTER_ID  CPU=${SLURM_CPUS_PER_TASK:-NA}"

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Use a different seed per task for better chain jitter
SEED=$((1234 + task_id))

# Run the cluster interval-censored pipeline for this cluster
Rscript "$RUNNER" \
  --cluster-ids "$CLUSTER_ID" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "✅ 完成: Cluster=$CLUSTER_ID"