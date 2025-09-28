#!/bin/bash
# =================================================================
# WONDER PyMC 稳健型任务提交脚本 (v2.2)
#
# 更新日志:
# v2.2: 增加 --estimate 参数控制 min/avg/max
# v2.1: 支持混合ID输入 (e.g., "2,cat21")
# v2.0: 重构为稳健的启动器，将批量逻辑交给main.py
# =================================================================

# --- SLURM 配置 ---
#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_PyMC_Batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=24:00:00
#SBATCH --output=logs/WONDER_PyMC-%x-%j.out
#SBATCH --error=logs/WONDER_PyMC-%x-%j.err

# --- 帮助函数 ---
log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"
}

# --- 可编辑参数 ---
# 批量分析的核心参数 (支持逗号分隔的列表)
DISEASE_CODE=${DISEASE_CODE:-"C81-C96"}
# 使用ID指定具体化合物 (如 "2"), 使用 "cat" 前缀指定类别 (如 "cat21")
COMPOUND=${COMPOUND:-"2,9,cat21,cat33"} 
MODEL_TYPE=${MODEL_TYPE:-"M0,M1,M2,M3,M5_SVI,M6_ENV1"}
LAG_YEARS=${LAG_YEARS:-"5,10"}
MEASURE_TYPE=${MEASURE_TYPE:-"Weight"}
# <--- 新增：控制使用 min, avg, 还是 max 估算值 --->
ESTIMATE_TYPE=${ESTIMATE_TYPE:-"avg,max"} # 支持逗号分隔

# 生产级采样设置
SAMPLING_MODE=${SAMPLING_MODE:-"production"}
DRAWS=${DRAWS:-"4000"}
TUNE=${TUNE:-"2000"}
TARGET_ACCEPT=${TARGET_ACCEPT:-"0.95"}
CHAINS=${CHAINS:-"4"}
CORES=${CORES:-"${SLURM_CPUS_PER_TASK:-16}"}

# --- 命令行参数覆盖 (可选) ---
while [[ $# -gt 0 ]]; do
  key="$1"; shift
  case $key in
    --disease) DISEASE_CODE="$1"; shift ;;
    --compound) COMPOUND="$1"; shift ;;
    --model) MODEL_TYPE="$1"; shift ;;
    --lag) LAG_YEARS="$1"; shift ;;
    --measure) MEASURE_TYPE="$1"; shift ;;
    --estimate) ESTIMATE_TYPE="$1"; shift ;; # <--- 新增
    --draws) DRAWS="$1"; shift ;;
    --tune) TUNE="$1"; shift ;;
    --target-accept) TARGET_ACCEPT="$1"; shift ;;
    --chains) CHAINS="$1"; shift ;;
    --cores) CORES="$1"; shift ;;
  --sampling-mode) SAMPLING_MODE="$1"; shift ;;
    *) log "WARN" "未知参数: $key" ;;
  esac
done

set -eo pipefail

# ==================== 任务开始 ====================
log "INFO" "WONDER PyMC 分析任务启动"
mkdir -p logs

# --- 1. 环境激活 ---
log "INFO" "激活 Conda 环境: pymc"
# 一些集群的激活钩子在 set -u 下会触发未绑定变量（如 ADDR2LINE），临时关闭 -u 更稳妥
set +u || true
# 预定义可能被钩子引用的变量，避免未绑定
export ADDR2LINE="${ADDR2LINE-}"
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/anaconda3/etc/profile.d/conda.sh"
else
  log "ERROR" "无法找到 Conda 初始化脚本。"
  exit 1
fi
conda activate pymc || { log "ERROR" "激活 'pymc' 环境失败。"; exit 1; }
set -u || true
log "INFO" "Conda 环境激活成功. Python路径: $(which python)"

# --- 2. 路径设置 ---
# 解析项目根目录：兼容脚本位于仓库根或 Code/PYMC 目录两种情况
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
if [ -f "${SCRIPT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="${SCRIPT_DIR}"
elif [ -f "${SCRIPT_DIR}/../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
elif [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
else
  PROJECT_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || echo "${SCRIPT_DIR}")"
fi
PYMC_DIR="${PROJECT_ROOT}/Code/PYMC"
CONFIG_PATH="${PROJECT_ROOT}/config.yaml"

cd "${PROJECT_ROOT}"
log "INFO" "项目根目录: ${PROJECT_ROOT}"
log "INFO" "配置文件:   ${CONFIG_PATH}"
log "INFO" "PyMC目录:   ${PYMC_DIR}"

# --- 3. 打印分析摘要 ---
log "INFO" "========== 分析参数摘要 =========="
log "INFO" "疾病 (Disease):    ${DISEASE_CODE}"
log "INFO" "暴露 (Compound):   ${COMPOUND}"
log "INFO" "估算 (Estimate):   ${ESTIMATE_TYPE}" # <--- 新增
log "INFO" "模型 (Model):      ${MODEL_TYPE}"
log "INFO" "滞后 (Lag):        ${LAG_YEARS}"
log "INFO" "测量 (Measure):    ${MEASURE_TYPE}"
log "INFO" "采样模式:          ${SAMPLING_MODE}"
log "INFO" "采样设置: Draws=${DRAWS}, Tune=${TUNE}, Chains=${CHAINS}, Cores=${CORES}"
log "INFO" "======================================"

# --- 4. 构建并执行命令 ---
CMD=(python "${PYMC_DIR}/main.py"
  --disease "${DISEASE_CODE}"
  --compound "${COMPOUND}"
  --model "${MODEL_TYPE}"
  --lag "${LAG_YEARS}"
  --measure "${MEASURE_TYPE}"
  --estimate "${ESTIMATE_TYPE}" # <--- 新增
  --sampling-mode "${SAMPLING_MODE}"
  --config-path "${CONFIG_PATH}"
  --draws "${DRAWS}"
  --tune "${TUNE}"
  --target-accept "${TARGET_ACCEPT}"
  --chains "${CHAINS}"
  --cores "${CORES}"
  --verbose
)

log "INFO" "执行命令: ${CMD[*]}"

if ! "${CMD[@]}"; then
  log "ERROR" "Python 分析脚本执行失败。请检查错误日志。"
  exit 1
fi

log "INFO" "WONDER PyMC 分析任务成功完成"
exit 0