#!/bin/bash
# Single job for testing cmdstan interval-censored mixed model: C00_C97 (non-stratified)
# Run locally:
#   bash Code/brms/test_cmdstan_interval_C00_C97.sh
# Or submit to Slurm:
#   sbatch Code/brms/test_cmdstan_interval_C00_C97.sh

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_cmdstan_interval_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=cmdstan_interval_test_%j.out
#SBATCH --error=cmdstan_interval_test_%j.err

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

# Load devtoolset for newer g++ on CentOS (safe no-op elsewhere)
module load devtoolset-8 2>/dev/null || log WARN "Could not load devtoolset-8, using system compiler"

# Set environment variables for CmdStan/TBB
export TBB_CXX_TYPE=gcc

RUNNER="Code/brms/cmdstan_interval_runner.R"
if [ ! -f "$RUNNER" ]; then
  log ERROR "找不到R脚本: $RUNNER"; exit 1
fi

# Limit threading to allocation to be polite on shared nodes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Fixed seed for testing (can override via SEED)
SEED="${SEED:-1234}"

# Optional fast mode: FAST=1 adds --test to runner (fewer iters)
if [ "${FAST:-0}" = "1" ] || [ "${FAST:-false}" = "true" ]; then
  TEST_FLAG="--test"
  log INFO "启用快速测试模式 (--test)"
else
  TEST_FLAG=""
fi

# Optionally override via CANCER_TYPES env; default try C00_C97
CANCER_TYPES="${CANCER_TYPES:-C00_C97}"

# If the default interval data exists, verify the requested cancer type is present; otherwise, fall back to first available
DATA_PATH="$PROJECT_ROOT/Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval.csv"
if [ -f "$DATA_PATH" ]; then
  export CANCER_TYPES DATA_PATH
  SELECTED=$(Rscript - <<'RS'
  suppressPackageStartupMessages({library(data.table)})
  ct <- Sys.getenv("CANCER_TYPES", unset = "C00_C97")
  fp <- Sys.getenv("DATA_PATH", unset = "")
  if (file.exists(fp)) {
    dt <- tryCatch(fread(fp, select = "Cancer_Type"), error = function(e) NULL)
    if (!is.null(dt) && "Cancer_Type" %in% names(dt)) {
      u <- unique(dt$Cancer_Type)
      if (ct %in% u) {
        cat(ct)
      } else if (length(u)) {
        cat(u[1])
      } else {
        cat(ct)
      }
    } else {
      cat(ct)
    }
  } else {
    cat(ct)
  }
RS
)
  # Trim whitespace/newlines
  SELECTED="$(echo -n "$SELECTED" | tr -d '\n' | xargs)"
else
  SELECTED="$CANCER_TYPES"
fi

log INFO "开始测试: Cancer=${SELECTED}"

# Run the interval-censored pipeline
Rscript "$RUNNER" \
  --cancer-types "$SELECTED" \
  --chains 4 --iter 2000 --warmup 1000 \
  --adapt-delta 0.95 --max-treedepth 12 \
  --seed "$SEED" ${TEST_FLAG}

log INFO "✅ 测试完成: Cancer=${SELECTED}"
