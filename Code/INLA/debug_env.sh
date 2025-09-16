#!/bin/bash

# --- Slurm设置 ---
#SBATCH --job-name=Env_Debug
#SBATCH --output=debug_manual_log.out  # 指定一个备用日志名
#SBATCH --error=debug_manual_log.err
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

# --- 关键调试步骤 ---
# 强制将所有后续命令的输出(标准和错误)都重定向到项目根目录的一个文件里
exec &> "${SLURM_SUBMIT_DIR}/debug_output.log"

echo "=== 侦察兵脚本开始 ==="
date
echo ""

echo "--- 1. 初始环境检查 ---"
echo "当前工作目录 (pwd): $(pwd)"
echo "运行节点 (hostname): $(hostname)"
echo "用户 (whoami): $(whoami)"
echo "初始PATH环境变量:"
echo "$PATH"
echo ""

echo "--- 2. 测试 Conda 初始化 ---"
echo "尝试 source conda.sh..."
source /public/home/acf4pijnzl/miniconda3/etc/profile.d/conda.sh
CONDA_INIT_EXIT_CODE=$? # 捕获上一条命令的退出代码
echo "source conda.sh 退出代码: ${CONDA_INIT_EXIT_CODE} (0代表成功)"
echo ""

echo "--- 3. 测试 Conda 环境激活 ---"
echo "尝试 conda activate WDP..."
conda activate WDP
CONDA_ACTIVATE_EXIT_CODE=$?
echo "conda activate WDP 退出代码: ${CONDA_ACTIVATE_EXIT_CODE} (0代表成功)"
echo ""

echo "--- 4. 激活后的环境检查 ---"
echo "激活后的PATH环境变量:"
echo "$PATH"
echo ""

echo "检查 Rscript 路径 (which Rscript)..."
which Rscript
RSCRIPT_WHICH_EXIT_CODE=$?
echo "which Rscript 退出代码: ${RSCRIPT_WHICH_EXIT_CODE} (0代表成功)"
echo ""

echo "=== 侦察兵脚本结束 ==="
