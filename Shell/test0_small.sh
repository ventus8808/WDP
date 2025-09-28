#!/bin/bash
# 小规模烟囱测试：单化合物+少量模型，短采样，快速验证能跑通
#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_PyMC_Smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -eo pipefail

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 切到项目根目录（严格以脚本所在目录的上一级为准，忽略环境同名变量）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
cd "$PROJECT_ROOT" || { log ERROR "无法切换到项目根目录: $PROJECT_ROOT"; exit 1; }
log INFO "项目根目录: $PROJECT_ROOT"

# 将 Slurm 的默认 .out/.err 在作业结束时移动到项目根目录
move_logs_to_root() {
  # 仅在 Slurm 环境下执行
  if [ -n "${SLURM_JOB_ID-}" ] && [ -n "${SLURM_JOB_NAME-}" ]; then
    local submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
    local out_src="${submit_dir}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    local err_src="${submit_dir}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
    local out_dst="${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    local err_dst="${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
    [ -f "$out_src" ] && mv -f "$out_src" "$out_dst" || true
    [ -f "$err_src" ] && mv -f "$err_src" "$err_dst" || true
  fi
}
trap move_logs_to_root EXIT

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

# 基础文件存在性检查
if [ ! -f "Code/PYMC/main.py" ]; then
  log ERROR "未找到 Code/PYMC/main.py，请检查项目根目录是否正确: $PROJECT_ROOT"; exit 1;
fi



# 运行最小测试
sbatch_args=(
  --partition=kshctest
  --job-name=WONDER_PyMC_Smoke_Sub
  --nodes=1 --ntasks=1 --cpus-per-task=4
  --mem-per-cpu=2G --time=01:00:00
)

log INFO "提交最小测试：C81-C96 | compound=2 | models=M5_SVI,M6_ENV1"
python Code/PYMC/main.py \
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
