#!/bin/bash
# ========================
# WONDER R分析 单任务CPU测试脚本 (v3 - 最终优化版)
# ========================

#SBATCH --partition=kshdtest
#SBATCH --job-name=WONDER_R_Test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8           # 维持8个CPU核心
#SBATCH --mem-per-cpu=2G            # (新) 采纳您的成功配置，共16G内存
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

# --- 已移除DCU申请，因为R/INLA是纯CPU程序 ---
# #SBATCH --gres=dcu:1

# ========================
# Conda环境设置
# ========================
echo "清理并设置环境模块..."
module purge

echo "激活Conda环境..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "Conda环境 'WDP' 已激活。"
echo "Rscript路径: $(which Rscript)"

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

echo "开始执行R分析脚本，测试参数如下:"
echo "  疾病代码: ${DISEASE_CODE}"
echo "  化合物ID: ${COMPOUND_ID}"
echo "-------------------------------------"

R_SCRIPT_NAME="BYM_INLA_Production.R"
CONFIG_PATH="../../config.yaml"

Rscript ${R_SCRIPT_NAME} \
  --config ${CONFIG_PATH} \
  --disease-code "${DISEASE_CODE}" \
  --pesticide-category "compound:${COMPOUND_ID}" \
  --measure-type "Weight,Density" \
  --estimate-types "avg" \
  --lag-years "5,10" \
  --model-types "M0,M1,M2,M3" \
  --verbose

status=$?
if [ $status -ne 0 ]; then
    echo "！！！R脚本执行失败，请检查错误日志: ${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
else
    echo "--- R脚本执行成功 ---"
    echo "请检查项目根目录下的日志文件: ${PROJECT_ROOT}/${SLURM_JOB_NAME}-${SLURM_JOB_ID}.log"
fi
