#!/bin/bash
# Test1: 快速冒烟测试 (Quick Smoke Test)
# 目标：用最少的资源和时间，验证代码、数据和环境的连通性。
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Test1_Smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=03:00:00
#SBATCH --output=logs/Test1_Smoke_%j.out
#SBATCH --error=logs/Test1_Smoke_%j.err

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

# 确保logs目录存在
mkdir -p logs

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

log INFO "Test1: 开始快速冒烟测试..."

# 使用最简化的参数，目标是快速完成
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M0" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 100 --tune 50 --chains 1 --cores ${SLURM_CPUS_PER_TASK:-4} \
  --config-path "config.yaml" --verbose

log INFO "✅ Test1: 快速冒烟测试完成！"