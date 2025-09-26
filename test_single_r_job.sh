#!/bin/bash
# ========================
# WONDER R分析 单任务CPU测试脚本 (v6 - 最终路径修正版)
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_R_Test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

# ========================
# Conda环境设置
# ========================
echo "清理并设置环境模块..."
module purge
module load compiler/dtk/23.10

echo "激活Conda环境..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA
echo "Conda环境 'INLA' 已激活。"
echo "Rscript路径: $(which Rscript)"
echo "[诊断] Rscript环境与INLA包可用性："
Rscript -e "cat('INLA available:', suppressMessages(require(INLA, quietly=TRUE)), '\n')"
Rscript -e "INLA:::inla.version()"

# ========================
# 运行R分析命令
# ========================

# 原始项目根目录
PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
# 包含R脚本和utils文件夹的目录
SCRIPT_DIR="${PROJECT_ROOT}/Code/INLA"

echo "切换工作目录至: ${SCRIPT_DIR}"
cd ${SCRIPT_DIR}

# --- 定义测试参数 ---
DISEASE_CODE="C81-C96"
COMPOUND_ID="1"

echo "开始执行R分析脚本..."
echo "  疾病代码: ${DISEASE_CODE}"
echo "  化合物ID: ${COMPOUND_ID}"
echo "-------------------------------------"

R_SCRIPT_NAME="INLA_Main.R"
# --- 修正配置文件路径 ---
# 使用项目根目录构建配置文件的绝对路径，确保总能找到它
CONFIG_PATH="${PROJECT_ROOT}/Code/INLA/INLA_Config/analysis_config.yaml"

echo "使用配置文件: $(realpath ${CONFIG_PATH})"
# 优先使用节点本地临时目录，避免网络文件系统导致的INLA临时文件问题
export TMPDIR="${SLURM_TMPDIR:-/tmp}/${USER}/wdp_inla_${SLURM_JOB_ID:-manual}"
mkdir -p "$TMPDIR"
echo "使用临时目录: $TMPDIR"

# 若 SLURM_TMPDIR 未定义，则与 TMPDIR 对齐，确保 R 脚本选择相同本地目录
if [ -z "${SLURM_TMPDIR}" ]; then
  export SLURM_TMPDIR="$TMPDIR"
fi

# 诊断：打印与校验本地临时目录可写性与 R 端视角
echo "[诊断] SLURM_TMPDIR=${SLURM_TMPDIR}"
echo "[诊断] TMPDIR=${TMPDIR}"
echo "[诊断] df -h (SLURM_TMPDIR 与当前目录)"
df -h "${SLURM_TMPDIR}" || true
df -h . || true
echo "[诊断] 检查 TMPDIR 写权限"
touch "${TMPDIR}/_wdp_tmp_write_check" && echo "  ✓ TMPDIR 可写" || echo "  ✗ TMPDIR 不可写"
rm -f "${TMPDIR}/_wdp_tmp_write_check" || true

echo "[诊断] R 环境与临时目录视角"
Rscript - <<'RS_EOF'
cat("R version: ", R.version.string, "\n", sep="")
cat("INLA available: ")
ok <- suppressMessages(require(INLA, quietly=TRUE)); cat(ok, "\n")
cat("TMPDIR=", Sys.getenv("TMPDIR"), "\n", sep="")
cat("SLURM_TMPDIR=", Sys.getenv("SLURM_TMPDIR"), "\n", sep="")
cat("tempdir()=", tempdir(), "\n", sep="")
tf <- tempfile(pattern = "wdp_inla_check_")
writeLines("ok", tf)
cat("tempfile exists=", file.exists(tf), " (", tf, ")\n", sep="")
unlink(tf, force=TRUE)
RS_EOF

echo "[运行模式] 最小化组合：Weight / avg / 5y / M0"
Rscript ${R_SCRIPT_NAME} \
  --config ${CONFIG_PATH} \
  --disease-code "${DISEASE_CODE}" \
  --pesticide-category "compound:${COMPOUND_ID}" \
  --measure-type "Weight" \
  --estimate-types "avg" \
  --lag-years "5" \
  --model-types "M0" \
  --verbose

# --- 如需恢复全组合，请改用下方命令并注释掉上方“最小化组合” ---
# Rscript ${R_SCRIPT_NAME} \
#   --config ${CONFIG_PATH} \
#   --disease-code "${DISEASE_CODE}" \
#   --pesticide-category "compound:${COMPOUND_ID}" \
#   --measure-type "Weight,Density" \
#   --estimate-types "avg" \
#   --lag-years "5,10" \
#   --model-types "M0,M1,M2,M3" \
#   --verbose

status=$?
if [ $status -ne 0 ]; then
    echo "！！！R脚本执行失败，请检查错误日志: ${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
else
    echo "--- R脚本执行成功 ---"
    echo "请检查项目根目录下的日志文件: ${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.log"
fi
