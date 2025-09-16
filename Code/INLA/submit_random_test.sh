#!/bin/bash

#================================================================================
# Slurm SBATCH Directives for a Random Compound Full Test (v4 - Final)
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

# (修正) 获取项目根目录的绝对路径，确保输出路径正确
PROJECT_ROOT=$(pwd)

MAPPING_FILE="Data/Processed/Pesticide/mapping.csv" 
if [ ! -f "$MAPPING_FILE" ]; then
    echo "❌ 错误: 化合物映射文件未找到: $MAPPING_FILE"
    exit 1
fi

RANDOM_LINE=$(tail -n +2 "$MAPPING_FILE" | shuf -n 1)
COMPOUND_ID=$(echo "$RANDOM_LINE" | cut -d',' -f1)
COMPOUND_NAME=$(echo "$RANDOM_LINE" | cut -d',' -f2)

DISEASE_CODE="C81-C96"
MEASURE_TYPES="Weight,Density"
ESTIMATE_TYPES="min,avg,max"
LAG_YEARS="5,10"
MODEL_TYPES="M0,M1,M2,M3"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SAFE_COMPOUND_NAME=$(echo "$COMPOUND_NAME" | sed 's/[^a-zA-Z0-9]/-/g')
# (修正) 使用绝对路径定义结果目录
RESULTS_DIR="${PROJECT_ROOT}/Result/RandomTest_${SAFE_COMPOUND_NAME}_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "📋 本次随机测试参数:"
echo "--------------------------"
echo "  化合物: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})"
echo "  疾病代码: ${DISEASE_CODE}"
echo "  结果将保存在: ${RESULTS_DIR}"
echo "--------------------------"

# --- 2. 激活Conda环境并直接运行R脚本 ---
echo "🚀 激活 Conda 环境 'WDP'..."
source /public/home/acf4pijnzl/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "   ✅ Conda 环境已激活"

# (修正) 进入R脚本所在的目录，以解决 source() 路径问题
echo "🚀 进入脚本目录并运行 R 脚本..."
cd Code/INLA/

# (修正) R脚本的路径和配置文件路径也相应简化
Rscript BYM_INLA_Production.R \
  --config config/analysis_config.yaml \
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