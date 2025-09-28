#!/bin/bash
# 完整化合物批量：一个或多个化合物 × 全部模型配置，生产采样设置
#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_PyMC_AllChem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=98G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/WONDER_AllChem_%j.out
#SBATCH --error=logs/WONDER_AllChem_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 切到项目根目录
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

DISEASE=${DISEASE:-"C81-C96"}
COMPOUNDS=${COMPOUNDS:-"2,9,cat21,cat33"}
MODELS=${MODELS:-"M0,M1,M2,M3,M5_SVI,M6_ENV1"}

log INFO "开始完整化合物批量：$DISEASE | compounds=$COMPOUNDS | models=$MODELS"
python Code/PYMC/main.py \
  --disease "$DISEASE" \
  --compound "$COMPOUNDS" \
  --model "$MODELS" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg,max" \
  --sampling-mode "production" \
  --draws 4000 --tune 2000 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-32} --target-accept 0.95 \
  --config-path "config.yaml" --verbose

log INFO "完成"
