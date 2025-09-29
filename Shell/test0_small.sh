#!/bin/bash
# 基础验证测试：双阶段模型验证，短采样，2小时内完成
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Basic_Test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
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

# 设置Python输出不缓冲，确保实时显示进度
export PYTHONUNBUFFERED=1

log INFO "提交最小测试：C81-C96 | compound=2 | models=M5_SVI,M6_ENV1"
log INFO "开始执行Python分析..."

# 使用更保守的采样参数，避免数值问题
# 先测试简单的M0模型，成功后再尝试交互模型
log INFO "开始第一阶段：基础模型M0测试..."
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M0" \
  --lag "5" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 500 --tune 200 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-32} --target-accept 0.85 \
  --config-path "config.yaml" --verbose 2>&1 | tee -a "smoke_test_${SLURM_JOB_ID:-$$}.log"

PHASE1_EXIT_CODE=$?
if [ $PHASE1_EXIT_CODE -eq 0 ]; then
  log INFO "✅ 第一阶段成功，开始第二阶段：社会脆弱性模型测试..."
  python -u Code/PYMC/main.py \
    --disease "C81-C96" \
    --compound "2" \
    --model "M1" \
    --lag "5" \
    --measure "Weight" \
    --estimate "avg" \
    --sampling-mode "test" \
    --draws 500 --tune 200 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-32} --target-accept 0.85 \
    --config-path "config.yaml" --verbose 2>&1 | tee -a "smoke_test_${SLURM_JOB_ID:-$$}.log"
  
  PHASE2_EXIT_CODE=$?
  if [ $PHASE2_EXIT_CODE -ne 0 ]; then
    log ERROR "第二阶段失败，但第一阶段成功。交互模型可能有问题。"
    exit $PHASE2_EXIT_CODE
  fi
else
  log ERROR "第一阶段（基础模型）失败，退出码: $PHASE1_EXIT_CODE"
  exit $PHASE1_EXIT_CODE
fi

# 检查Python命令的退出状态
PYTHON_EXIT_CODE=$?
if [ $PYTHON_EXIT_CODE -ne 0 ]; then
  log ERROR "Python分析失败，退出码: $PYTHON_EXIT_CODE"
  log ERROR "请检查详细日志: smoke_test_${SLURM_JOB_ID:-$$}.log"
  exit $PYTHON_EXIT_CODE
fi

log INFO "完成"
