#!/bin/bash

# --- Slurm设置 ---
#SBATCH --job-name=R_Run_Debug
#SBATCH --output=debug_r_run.out
#SBATCH --error=debug_r_run.err
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00

# --- 关键调试步骤 ---
# 强制将所有后续命令的输出重定向到根目录的文件
exec &> "${SLURM_SUBMIT_DIR}/debug_r_run_output.log"

echo "=== 终极侦察兵脚本开始 ==="
date
echo ""

echo "--- 1. 激活环境 (已知成功) ---"
source /public/home/acf4pijnzl/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "Conda环境 'WDP' 已激活"
echo ""

echo "--- 2. 尝试执行R脚本 ---"
echo "当前工作目录 (pwd): $(pwd)"
echo "尝试进入 'Code/INLA/' 目录..."
cd Code/INLA/
CD_EXIT_CODE=$?
echo "cd 命令退出代码: ${CD_EXIT_CODE} (0代表成功)"
echo "进入后，当前工作目录 (pwd): $(pwd)"
echo ""

echo "尝试用 --help 运行R脚本，并捕获其输出..."
Rscript BYM_INLA_Production.R --help
RSCRIPT_EXIT_CODE=$?
echo "Rscript 命令退出代码: ${RSCRIPT_EXIT_CODE} (0代表成功)"
echo ""

echo "=== 侦察兵脚本结束 ==="
