#!/bin/bash
#SBATCH --job-name=wdp_single_compound_test
#SBATCH --output=slurm_logs/single_compound_test-%j.out
#SBATCH --error=slurm_logs/single_compound_test-%j.err
#SBATCH --partition=kshdtest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --gres=dcu:1

# WDP Single Compound Test SLURM Submission Script
# Tests all model combinations for a single compound (default: 2,4-D)
# Expected runtime: ~30-60 minutes
# Expected output: ~24 model results (3 estimates × 2 measures × 4 models)

echo "============================================================"
echo "SLURM Job Information"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start Time: $(date)"
echo ""

# 确保日志目录存在
mkdir -p slurm_logs

# 加载 R 环境
echo "Loading R environment..."
module load R/4.2.2

# 检查R环境
echo "R version:"
R --version | head -n 1

# 设置环境变量
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 显示当前工作目录
echo ""
echo "Working directory: $(pwd)"
echo "Available space:"
df -h $(pwd)

# 检查关键文件是否存在
echo ""
echo "Checking required files..."
KEY_FILES=(
    "Code/INLA/BYM_INLA_Production.R"
    "Code/INLA/config/analysis_config.yaml"
    "Code/INLA/run_single_compound_test.sh"
    "Data/Processed/Pesticide/PNSP.csv"
    "Data/Processed/CDC/C81-C96.csv"
)

for file in "${KEY_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ❌ $file (MISSING!)"
    fi
done

echo ""
echo "============================================================"
echo "Starting Analysis"
echo "============================================================"

# 执行单化合物测试脚本
# 可以通过修改这些参数来测试不同的化合物和条件
bash Code/INLA/run_single_compound_test.sh \
    "2" \           # Compound ID (2 = 2,4-D)
    "2,4-D" \       # Compound Name
    "C81-C96" \     # Disease Code
    "Weight,Density" \     # Measure Types
    "min,avg,max" \        # Estimate Types
    "5" \           # Lag Years
    "M0,M1,M2,M3"   # Model Types

# 获取执行结果
EXIT_CODE=$?

echo ""
echo "============================================================"
echo "Job Summary"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "End Time: $(date)"
echo "Exit Code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "Status: ✅ SUCCESS"
    
    # 显示输出文件信息
    echo ""
    echo "Output files in Result/Filter/:"
    ls -la Result/Filter/Results_C81-C96_2,4-D_*.csv 2>/dev/null | tail -n 3
    
    # 显示最新结果文件的统计信息
    LATEST_FILE=$(ls -t Result/Filter/Results_C81-C96_2,4-D_*.csv 2>/dev/null | head -n 1)
    if [ -f "$LATEST_FILE" ]; then
        echo ""
        echo "Latest results file: $(basename "$LATEST_FILE")"
        echo "Number of results: $(($(wc -l < "$LATEST_FILE") - 1))"
    fi
    
else
    echo "Status: ❌ FAILED"
    echo ""
    echo "Check the error log for details:"
    echo "  slurm_logs/single_compound_test-${SLURM_JOB_ID}.err"
fi

echo ""
echo "Log files:"
echo "  Output: slurm_logs/single_compound_test-${SLURM_JOB_ID}.out"
echo "  Error:  slurm_logs/single_compound_test-${SLURM_JOB_ID}.err"
echo "============================================================"

exit $EXIT_CODE