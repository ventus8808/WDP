#!/bin/bash

#================================================================================
# Slurm SBATCH Script for the "Hello World" R Test
#================================================================================
# -- Job name
#SBATCH --job-name=HelloWorldTest
#
# -- Output and error files
#SBATCH --output=slurm_logs/hello_world-%j.out
#SBATCH --error=slurm_logs/hello_world-%j.err
#
# -- Resource allocation (minimal)
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:02:00 # 2 minutes is plenty

#================================================================================
# Job Execution
#================================================================================
echo "========================================================"
echo "Starting 'Hello World' R Test"
echo "Job ID: $SLURM_JOB_ID"
echo "Time: $(date)"
echo "========================================================"

# 使用我们确认过的、最稳健的绝对路径方法来运行R脚本
echo "🚀 Using absolute path to run the hello_world.R script..."
/public/home/acf4pijnzl/miniconda3/envs/WDP/bin/Rscript Code/INLA/hello_world.R

echo "========================================================"
echo "Test script finished."
echo "Time: $(date)"
echo "========================================================"
