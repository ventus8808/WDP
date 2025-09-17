#!/bin/bash
# ========================
# WONDER R分析 深度调试脚本 (v8)
# ========================

#SBATCH --partition=kshdtest
#SBATCH --job-name=WONDER_DEBUG
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=0:30:00 # 调试任务，30分钟足矣
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err
#SBATCH --gres=dcu:1

# --- 调试步骤：开启INLA详细日志 ---
echo "--- 正在进入深度调试模式 ---"
# 切换到包含R脚本的目录
cd /public/home/acf4pijnzl/WDP/Code/INLA

# 备份原始的 model_fitting.R 文件
cp utils/model_fitting.R utils/model_fitting.R.bak

echo "修改 model_fitting.R 以开启INLA详细日志..."
# 1. 强制开启 verbose=TRUE
sed -i "s/verbose = config\$model_fitting\$inla\$verbose/verbose = TRUE/g" utils/model_fitting.R
# 2. 注释掉所有可能隐藏底层C语言错误的 sink() 命令
sed -i "s/sink(sink_connection, type = \"message\")/cat('--- sink() disabled for debugging ---\\n') # sink(sink_connection, type = \"message\")/g" utils/model_fitting.R
sed -i "s/sink(type = \"message\")/cat('--- sink() disabled for debugging ---\\n') # sink(type = \"message\")/g" utils/model_fitting.R
echo "--- 调试模式已开启 ---"

# --- 正常执行环境设置和R脚本 ---
module purge
module load compiler/dtk/23.10
source ~/miniconda3/etc/profile.d/conda.sh
conda activate WDP

Rscript BYM_INLA_Production.R \
  --config config/analysis_config.yaml \
  --disease-code "C81-C96" \
  --pesticide-category "compound:1" \
  --measure-type "Weight" \
  --estimate-types "avg" \
  --lag-years "5" \
  --model-types "M0" \
  --verbose

# --- 调试步骤：恢复原始文件 ---
echo "--- 正在退出深度调试模式 ---"
echo "恢复原始 model_fitting.R 文件..."
mv utils/model_fitting.R.bak utils/model_fitting.R
echo "--- 调试结束 ---"
