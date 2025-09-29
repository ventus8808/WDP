#!/bin/bash
# 生产级批量分析：多化合物×全模型配置，高质量采样
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 切到项目根目录（增强路径检测逻辑，处理集群环境下的符号链接）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""

# 方法1: 脚本所在目录的上一级（优先）
if [ -f "${SCRIPT_DIR}/../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
# 方法2: 使用git仓库根目录
elif command -v git >/dev/null 2>&1; then
  GIT_ROOT="$(cd "${SCRIPT_DIR}" && git rev-parse --show-toplevel 2>/dev/null || true)"
  if [ -n "$GIT_ROOT" ] && [ -f "$GIT_ROOT/config.yaml" ]; then
    PROJECT_ROOT="$GIT_ROOT"
  fi
fi

# 方法3: 如果在Slurm环境，使用提交目录
if [ -z "$PROJECT_ROOT" ] && [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
fi

# 方法4: 最后回退到当前工作目录
if [ -z "$PROJECT_ROOT" ] && [ -f "$PWD/config.yaml" ]; then
  PROJECT_ROOT="$PWD"
fi

# 验证项目根目录
if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "无法确定项目根目录。请确保在WDP仓库根目录下提交作业。"
  log ERROR "当前脚本目录: $SCRIPT_DIR"
  log ERROR "SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR-未设置}"
  log ERROR "当前工作目录: $PWD"
  exit 1
fi

cd "$PROJECT_ROOT" || { log ERROR "无法切换到项目根目录: $PROJECT_ROOT"; exit 1; }
log INFO "项目根目录: $PROJECT_ROOT"

# 将 Slurm 的默认 .out/.err 在作业结束时移动到项目根目录
move_logs_to_root() {
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

# 基础文件检查
if [ ! -f "Code/PYMC/main.py" ]; then
  log ERROR "未找到 Code/PYMC/main.py，请检查项目根目录是否正确: $PROJECT_ROOT"; exit 1;
fi



DISEASE=${DISEASE:-"C81-C96"}
COMPOUNDS=${COMPOUNDS:-"2,9,cat21,cat33"}
MODELS=${MODELS:-"M0,M1,M2,M3"}

log INFO "开始生产级批量分析：$DISEASE | compounds=$COMPOUNDS | models=$MODELS"
python Code/PYMC/main.py \
  --disease "$DISEASE" \
  --compound "$COMPOUNDS" \
  --model "$MODELS" \
  --lag "5,10" \
  --measure "Weight,Density" \
  --estimate "avg,max" \
  --sampling-mode "production" \
  --draws 4000 --tune 2000 --chains 2 --cores ${SLURM_CPUS_PER_TASK:-16} --target-accept 0.95 \
  --config-path "config.yaml" --verbose

log INFO "完成"
