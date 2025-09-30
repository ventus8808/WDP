#!/bin/bash
# Test2: 单模型深度分析 (Single Model Deep Dive)
# 目标：对一个复杂的模型进行完整的、高质量的采样，用于结果验证。
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Test2_SingleFull
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=1-00:00:00  # 申请1天时间以确保完成
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 项目根目录检测
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""

if [ -f "${SCRIPT_DIR}/../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
fi

if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "无法确定项目根目录"
  exit 1
fi

cd "$PROJECT_ROOT" || exit 1
log INFO "项目根目录: $PROJECT_ROOT"

# 日志将输出到当前目录

# 激活conda环境
set +u
if [ "${CONDA_DEFAULT_ENV-}" != "pymc" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "找不到conda"; exit 1
  fi
  conda activate pymc || { log ERROR "激活pymc失败"; exit 1; }
fi

log INFO "Test2: 开始单模型深度分析..."

# 使用生产级参数，对一个指定的复杂模型进行完整采样
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M5_SVI" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "production" \
  --draws 4000 --tune 2000 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-16} --target-accept 0.95 \
  --config-path "config.yaml" --verbose

log INFO "✅ Test2: 单模型深度分析完成！"