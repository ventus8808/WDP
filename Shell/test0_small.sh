#!/bin/bash
# 小规模烟囱测试：单化合物+少量模型，短采样，快速验证能跑通
#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_PyMC_Smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/WONDER_Smoke-%j.out
#SBATCH --error=logs/WONDER_Smoke-%j.err

set -eo pipefail

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 切到项目根目录（脚本位于 WDP/Shell 下）
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"
log INFO "项目根目录: $PROJECT_ROOT"

# 激活conda（稳健处理：避免hook未绑定变量；已在pymc则跳过）
set +u || true
export ADDR2LINE="${ADDR2LINE-}"
export CONDA_BACKUP_CXX="${CONDA_BACKUP_CXX-}"
if [ "${CONDA_DEFAULT_ENV-}" != "pymc" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "找不到conda初始化脚本"; exit 1
  fi
  conda activate pymc || { log ERROR "激活pymc失败"; exit 1; }
fi
log INFO "Conda Python: $(which python)"

mkdir -p logs

# 运行最小测试
sbatch_args=(
  --partition=kshctest
  --job-name=WONDER_PyMC_Smoke_Sub
  --nodes=1 --ntasks=1 --cpus-per-task=4
  --mem-per-cpu=2G --time=01:00:00
)

log INFO "提交最小测试：C81-C96 | compound=2 | models=M5_SVI,M6_ENV1"
srun python Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M5_SVI,M6_ENV1" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 200 --tune 100 --chains 2 --cores ${SLURM_CPUS_PER_TASK:-4} --target-accept 0.9 \
  --config-path "config.yaml" --verbose

log INFO "完成"
