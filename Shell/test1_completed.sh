#!/bin/bash
# 完整单模型测试：覆盖全部基础与交互模型于一个化合物，较长采样
#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_PyMC_OneChem_AllModels
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --output=logs/WONDER_OneChem_%j.out
#SBATCH --error=logs/WONDER_OneChem_%j.err

set -euo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 切到项目根目录
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"
log INFO "项目根目录: $PROJECT_ROOT"

# 激活conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/anaconda3/etc/profile.d/conda.sh"
else
  log ERROR "找不到conda初始化脚本"; exit 1
fi
conda activate pymc || { log ERROR "激活pymc失败"; exit 1; }

mkdir -p logs

DISEASE=${DISEASE:-"C81-C96"}
COMPOUND=${COMPOUND:-"2"}
MODELS=${MODELS:-"M0,M1,M2,M3,M5_SVI,M6_ENV1"}

log INFO "开始完整单模型测试：$DISEASE | compound=$COMPOUND | models=$MODELS"
python Code/PYMC/main.py \
  --disease "$DISEASE" \
  --compound "$COMPOUND" \
  --model "$MODELS" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 1000 --tune 1000 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-16} --target-accept 0.9 \
  --config-path "config.yaml" --verbose

log INFO "完成"
