#!/bin/bash
# SLURM submit script for cmdstan_main_Ridgeline_Test.R
# Single job for C00-C97, lag=5, Overall EQI only
# Extracts full MCMC posterior draws for ridgeline visualization

#SBATCH --partition=kshctest
#SBATCH --job-name=Ridge_Test_C00_C97
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=Ridge_Test_C00-C97_%j.out
#SBATCH --error=Ridge_Test_C00-C97_%j.err</parameter>

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

# Ensure output directories exist
LOG_DIR="$PROJECT_ROOT/Result/log"
RIDGE_DIR="$PROJECT_ROOT/Result/Ridgeline"

if [ ! -d "$LOG_DIR" ]; then
  mkdir -p "$LOG_DIR"
  log INFO "创建日志目录: $LOG_DIR"
fi

if [ ! -d "$RIDGE_DIR" ]; then
  mkdir -p "$RIDGE_DIR"
  log INFO "创建输出目录: $RIDGE_DIR"
fi

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

log INFO "Conda环境: $CONDA_DEFAULT_ENV"

# Load devtoolset for newer g++ on CentOS
module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system g++"

# Set environment variables for CmdStan
export TBB_CXX_TYPE=gcc

# Limit threading to allocation
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

log INFO "========================================="
log INFO "Ridgeline Test: Posterior Extraction"
log INFO "========================================="
log INFO "Cancer:     C00_C97"
log INFO "Lag:        5 (2000-2005 EQI → 2006-2010 AAMR)"
log INFO "Model:      Overall EQI only"
log INFO "Iterations: 800 (warmup: 400)"
log INFO "Chains:     4"
log INFO "CPUs:       ${SLURM_CPUS_PER_TASK:-1}"
log INFO "Memory:     32G"
log INFO "========================================="

RUNNER="Code/brms/cmdstan_main_Ridgeline_Test.R"

if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Run the ridgeline test script
log INFO "开始执行 Ridgeline Test..."
Rscript "$RUNNER"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
  log INFO "========================================="
  log INFO "✅ Ridgeline Test 完成"
  log INFO "========================================="
  log INFO "输出文件: Result/Ridgeline/C00_C97_Ridge_Test.rds"
  log INFO ""
  log INFO "后续步骤:"
  log INFO "  1. 加载数据: data <- readRDS('Result/Ridgeline/C00_C97_Ridge_Test.rds')"
  log INFO "  2. 查看摘要: print(data\$summary)"
  log INFO "  3. 绘制山脊图:"
  log INFO "     library(ggridges)"
  log INFO "     ggplot(data\$draws_long, aes(x=effect, y=quintile, fill=quintile)) +"
  log INFO "       geom_density_ridges(alpha=0.7)"
  log INFO "========================================="
else
  log ERROR "Ridgeline Test 失败 (退出码: $EXIT_CODE)"
  exit $EXIT_CODE
fi
