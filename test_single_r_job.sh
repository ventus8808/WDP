#!/bin/bash
#SBATCH --partition=kshdtest
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=dcu:1
#SBATCH --mem-per-cpu=2G      # 8核×2G=16G（通常安全值）
#SBATCH --time=1:00:00


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
