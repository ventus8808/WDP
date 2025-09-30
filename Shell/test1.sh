#!/bin/bash
# Test1: 快速冒烟测试 (Quick Smoke Test)
# 目标：用最少的资源和时间，验证代码、数据和环境的连通性。
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Test1_Smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G   # 4*4G=16G
#SBATCH --time=01:00:00    # 1小时足够
#SBATCH --output=logs/Test1_Smoke_%j.out
#SBATCH --error=logs/Test1_Smoke_%j.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- (激活conda和切换目录的代码，请从您之前的脚本复制) ---
# (此处省略)
# ...

log INFO "Test1: 开始快速冒烟测试..."

# 使用最简化的参数，目标是快速完成
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M0" \
  --lag "10" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 100 --tune 50 --chains 1 --cores ${SLURM_CPUS_PER_TASK:-4} \
  --config-path "config.yaml" --verbose

log INFO "✅ Test1: 快速冒烟测试完成！"