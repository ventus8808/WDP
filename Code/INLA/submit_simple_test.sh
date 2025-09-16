#!/bin/bash

#================================================================================
# Slurm SBATCH Script for the Minimal INLA BYM2 Test (v2 - Absolute Path)
#================================================================================
# -- Job name
#SBATCH --job-name=INLA_BYM2_Test
#
# -- Output and error files
#SBATCH --output=slurm_logs/simple_test-%j.out
#SBATCH --error=slurm_logs/simple_test-%j.err
#
# -- Resource allocation
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00

#================================================================================
# Job Execution
#================================================================================
echo "========================================================"
echo "Starting Minimal INLA BYM2 Test (v2 - Using Absolute Path)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time: $(date)"
echo "========================================================"

# --- 直接使用我们WDP环境里Rscript的绝对路径来运行R脚本 ---
# --- 这样就完全绕开了 conda activate 在批处理环境中的不确定性 ---

echo "🚀 Using absolute path to run the R test script..."
/public/home/acf4pijnzl/miniconda3/envs/WDP/bin/Rscript Code/INLA/simple_bym2_test.R

echo "========================================================"
echo "Test script finished."
echo "Time: $(date)"
echo "========================================================"
