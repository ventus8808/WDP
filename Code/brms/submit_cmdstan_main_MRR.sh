#!/bin/bash
# Slurm array launcher for cmdstan_main_MRR.R
# One task per NDD disease; each task runs 3 EQI0005 scenarios (lags 5/10/15),
# Overall EQI model only, and outputs MRD + MRR + lag test to Result/brms_MRR_lag/.
# Usage:
#   bash Code/brms/submit_cmdstan_main_MRR.sh    # submit array to Slurm
#   # or, advanced: sbatch --array=0-<N-1> Code/brms/submit_cmdstan_main_MRR.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_MRR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_MRR_%A_%a.out
#SBATCH --error=cmdstan_MRR_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Fixed NDD disease list ---
# NDD_DISEASES=(
#  "G20_G30_G12.2_F01_F03"
#  "G30_F01_F03"
#  "G30"
#  "G20"
#  "G12.2"
#  "G10"
#)
NDD_DISEASES=("F01")

# --- Locate project root ---
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

# --- Activate conda environment ---
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

export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/cmdstan_main_MRR.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# --- Controller mode: submit array ---
if [ -z "${SLURM_ARRAY_TASK_ID-}" ]; then
  N=${#NDD_DISEASES[@]}
  log INFO "提交数组任务: 0-$((N-1)) (共 $N 个疾病)"
  for i in "${!NDD_DISEASES[@]}"; do
    log INFO "  [$i] ${NDD_DISEASES[$i]}"
  done
  # Pass the disease list as a colon-separated env var to workers
  DISEASE_LIST=$(IFS=':'; echo "${NDD_DISEASES[*]}")
  sbatch --array=0-$((N-1)) \
    --export=ALL,DISEASE_LIST="$DISEASE_LIST",ENV_NAME="$ENV_NAME" \
    "$0"
  log INFO "提交完成。使用 squeue 查看进度。"
  exit 0
fi

# --- Worker mode ---
log INFO "开始处理任务 $SLURM_ARRAY_TASK_ID"

# Reconstruct array from colon-separated env var
IFS=':' read -ra DISEASES <<< "$DISEASE_LIST"
CANCER_TYPE="${DISEASES[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$CANCER_TYPE" ]; then
  log ERROR "无法读取任务 $SLURM_ARRAY_TASK_ID 的疾病类型"; exit 1
fi
log INFO "处理疾病: $CANCER_TYPE"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

SEED=$((1234 + SLURM_ARRAY_TASK_ID))

Rscript "$RUNNER" \
  --cancer-types "$CANCER_TYPE" \
  --output-dir "Result/brms_MRR_lag" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED"

log INFO "✅ 完成: Cancer=$CANCER_TYPE"
log INFO "输出目录: Result/brms_MRR_lag/"
