#!/bin/bash
# ========================
# WONDER R分析 渐进式调试脚本
# 从简单到复杂逐步测试INLA模型
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_R_Progressive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA

PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
cd ${PROJECT_ROOT}

echo "=== 渐进式INLA模型测试 ==="

# 步骤1: 测试基础INLA功能
echo "步骤1: 基础INLA测试"
sbatch test_inla_simple.sh
sleep 5

# 步骤2: 测试无时空效应的主分析
echo "步骤2: 主分析脚本 (无时空效应)"
export TMPDIR="${SLURM_TMPDIR:-/tmp}/${USER}/wdp_inla_${SLURM_JOB_ID:-manual}_debug"
mkdir -p "$TMPDIR"

Rscript Code/INLA/INLA_Main.R \
  --config Code/INLA/INLA_Config/analysis_config.yaml \
  --disease-code "C81-C96" \
  --pesticide-category "compound:1" \
  --measure-type "Weight" \
  --estimate-types "avg" \
  --lag-years "5" \
  --model-types "M0" \
  --verbose

debug_status=$?

if [ $debug_status -eq 0 ]; then
    echo "✓ 步骤2成功 - 可以尝试启用时空效应"
else
    echo "❌ 步骤2失败 - 需要进一步调试基础模型"
fi

echo "调试完成，状态码: $debug_status"