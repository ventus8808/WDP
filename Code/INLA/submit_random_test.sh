#!/bin/bash

#================================================================================
# Slurm SBATCH Directives for a Random Compound Full Test (v3 - No Docker)
#================================================================================
# -- Job name
#SBATCH --job-name=WDP_RandomTest
#
# -- Output and error files
#SBATCH --output=slurm_logs/random_test-%j.out
#SBATCH --error=slurm_logs/random_test-%j.err
#
# -- Resource allocation
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

#================================================================================
# Job Execution
#================================================================================
echo "========================================================"
echo "Job Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "========================================================"

# --- 1. 准备参数 ---
mkdir -p slurm_logs

# 正确的映射文件路径
MAPPING_FILE="Data/Processed/Pesticide/mapping.csv" 

if [ ! -f "$MAPPING_FILE" ]; then
    echo "❌ 错误: 化合物映射文件未找到: $MAPPING_FILE"
    exit 1
fi

# 随机选择化合物
RANDOM_LINE=$(tail -n +2 "$MAPPING_FILE" | shuf -n 1)
COMPOUND_ID=$(echo "$RANDOM_LINE" | cut -d',' -f1)
COMPOUND_NAME=$(echo "$RANDOM_LINE" | cut -d',' -f2)

# 设置分析参数
DISEASE_CODE="C81-C96"
MEASURE_TYPES="Weight,Density"
ESTIMATE_TYPES="min,avg,max"
LAG_YEARS="5,10"
MODEL_TYPES="M0,M1,M2,M3"

# 创建结果目录
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SAFE_COMPOUND_NAME=$(echo "$COMPOUND_NAME" | sed 's/[^a-zA-Z0-9]/-/g')
RESULTS_DIR="Result/RandomTest_${SAFE_COMPOUND_NAME}_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "📋 本次随机测试参数:"
echo "--------------------------"
echo "  化合物: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})"
echo "  疾病代码: ${DISEASE_CODE}"
echo "  所有模型: ${MODEL_TYPES}"
echo "  结果将保存在: ${RESULTS_DIR}"
echo "--------------------------"

# --- 2. 激活Conda环境并直接运行R脚本 ---
echo "🚀 激活 Conda 环境 'WDP'..."
source /public/home/acf4pijnzl/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "   ✅ Conda 环境已激活"

echo "🚀 直接运行 R 脚本..."
Rscript Code/INLA/BYM_INLA_Production.R \
  --config Code/INLA/config/analysis_config.yaml \
  --pesticide-category compound:${COMPOUND_ID} \
  --measure-type "${MEASURE_TYPES}" \
  --estimate-types "${ESTIMATE_TYPES}" \
  --lag-years "${LAG_YEARS}" \
  --model-types "${MODEL_TYPES}" \
  --disease-code "${DISEASE_CODE}" \
  --output-file "${RESULTS_DIR}/Results_${SAFE_COMPOUND_NAME}.csv" \
  --verbose

echo "========================================================"
echo "Job Finished: $(date)"
echo "查看结果目录: ${RESULTS_DIR}"
echo "========================================================"