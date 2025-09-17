#!/bin/bash
# ========================
# WONDER R分析 单任务CPU测试脚本 (含conda环境激活)
# ========================

#SBATCH --partition=kshdtest        # 您的分区
#SBATCH --job-name=WONDER_R_Test    # 清晰的任务名
#SBATCH --nodes=1                   # 单节点任务
#SBATCH --ntasks-per-node=1         # 单个R进程
#SBATCH --cpus-per-task=8           # 为R分配8个CPU核心 (INLA会使用它们)
#SBATCH --mem=32G                   # 32GB内存 (每个核心4G，符合集群策略)
#SBATCH --time=1:00:00              # 测试任务，运行1小时足矣
#SBATCH --output=%x-%j.log          # 日志文件，例如 WONDER_R_Test-12345.log
#SBATCH --error=%x-%j.err           # 错误日志

# --- 关键：我们不需要申请DCU，因为R/INLA用不上 ---
# #SBATCH --gres=dcu:1  <-- 已注释掉此行

# ========================
# Conda环境设置
# ========================
echo "清理并设置环境模块..."
module purge
# 注意：我们不需要加载 dtk 或 rocm 模块，因为它们是DCU专用环境
# 如果您的R环境需要特定编译器，可以在这里加载，否则此步可省略

echo "激活Conda环境..."
# 确保您的conda初始化脚本路径正确
source ~/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "Conda环境 'WDP' 已激活。"

# ========================
# 运行R分析命令
# ========================
echo "当前Rscript路径: $(which Rscript)"
echo "项目工作目录: $(pwd)"

# --- 定义测试参数 ---
DISEASE_CODE="C81-C96"
COMPOUND_ID="1"

echo "开始执行R分析脚本，测试参数如下:"
echo "  疾病代码: ${DISEASE_CODE}"
echo "  化合物ID: ${COMPOUND_ID}"
echo "-------------------------------------"


# 定义项目路径和脚本路径
PROJECT_DIR=$(pwd) # 脚本在WDP目录下运行，所以pwd就是项目目录
SCRIPT_PATH="${PROJECT_DIR}/Code/BYM_INLA_Production.R"
CONFIG_PATH="${PROJECT_DIR}/config.yaml"


Rscript ${SCRIPT_PATH} \
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
    echo "！！！R脚本执行失败，请检查错误日志: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
else
    echo "--- R脚本执行成功 ---"
fi
