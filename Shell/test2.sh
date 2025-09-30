#!/bin/bash
# Test2: 单模型深度分析 (Single Model Deep Dive)
# 目标：对一个复杂的模型进行完整的、高质量的采样，用于结果验证。
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Test2_SingleFull
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00  # 申请1天时间以确保完成
#SBATCH --output=logs/Test2_SingleFull_%j.out
#SBATCH --error=logs/Test2_SingleFull_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- (激活conda和切换目录的代码，请从您之前的脚本复制) ---
# (此处省略)
# ...

log INFO "Test2: 开始单模型深度分析..."

# 使用生产级参数，对一个指定的复杂模型进行完整采样
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M5_SVI" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "production" \
  --draws 4000 --tune 2000 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-16} --target-accept 0.95 \
  --config-path "config.yaml" --verbose

log INFO "✅ Test2: 单模型深度分析完成！"