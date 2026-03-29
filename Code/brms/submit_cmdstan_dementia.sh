#!/bin/bash
# Single-job Slurm launcher for Dementia (G30_F01_F03) main analysis.
# Runs all 5 EQI×AAMR scenarios with Overall + RUCC1-4 layers.
# Usage:
#   bash Code/brms/submit_cmdstan_dementia.sh        # submit via sbatch
#   bash Code/brms/submit_cmdstan_dementia.sh --run  # run directly (no sbatch)

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_cmdstan_dementia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_dementia_%j.out
#SBATCH --error=cmdstan_dementia_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

CANCER_TYPE="G30_F01_F03"

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

# --- Submit via sbatch if called directly (not inside a Slurm job) ---
if [ -z "${SLURM_JOB_ID-}" ] && [ "${1-}" != "--run" ]; then
  log INFO "提交 Dementia (${CANCER_TYPE}) 任务到 Slurm..."
  sbatch "$0"
  log INFO "提交完成。使用 squeue 查看进度。"
  exit 0
fi

# --- Activate conda environment ---
ENV_NAME="${ENV_NAME:-brms}"
set +u
if [ -z "${CONDA_DEFAULT_ENV-}" ] || [ "${CONDA_DEFAULT_ENV}" != "$ENV_NAME" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

RUNNER="Code/brms/cmdstan_main.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

log INFO "开始处理: $CANCER_TYPE (Overall + RUCC1-4, 5 scenarios)"

Rscript "$RUNNER" \
  --cancer-types "$CANCER_TYPE" \
  --data "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate_Typology_LandUse.csv" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed 1234

log INFO "✅ 完成: Cancer=$CANCER_TYPE"
log INFO "输出文件: Result/brms/${CANCER_TYPE}_main.csv"
