#!/bin/bash
#SBATCH --job-name=hello_world
#SBATCH --output=slurm_logs/hello_world-%j.out
#SBATCH --error=slurm_logs/hello_world-%j.err
#SBATCH --partition=kshdtest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --gres=dcu:1
# 确保日志目录存在
mkdir -p slurm_logs

# 加载 R 环境
module load R/4.2.2

# 执行 R 脚本
Rscript Code/INLA/hello_world.R
