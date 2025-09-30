#!/bin/bash
# Test3: 单化合物全模型并行分析 (Full Compound Analysis via Job Array)
# 目标：利用作业数组，为一个化合物并行运行所有模型、测量方式和估算类型的组合。
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Test3_CompoundArray
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00  # 为单个任务申请1天时间

# --- 关键：定义作业数组 ---
# 模型(6) * 测量(2) * 估算(3) = 36个任务
#SBATCH --array=0-35

# --- 将日志按子任务分开存放 ---
#SBATCH --output=logs/Test3_Array_%A_%a.out  # %A=主任务ID, %a=子任务ID
#SBATCH --error=logs/Test3_Array_%A_%a.err

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- (激活conda和切换目录的代码，请从您之前的脚本复制) ---
# (此处省略)
# ...

# --- 参数空间定义 ---
DISEASE="C81-C96"
COMPOUND="2"  # 指定要进行完整分析的化合物
LAG_YEARS="10"
# 根据你的config.yaml定义模型列表
MODELS_LIST=("M0" "M1" "M2" "M3" "M5_SVI" "M6_ENV1")
MEASURES_LIST=("Weight" "Density")
ESTIMATES_LIST=("avg" "min" "max")

# --- 任务映射逻辑 ---
NUM_MODELS=${#MODELS_LIST[@]}
NUM_MEASURES=${#MEASURES_LIST[@]}
NUM_ESTIMATES=${#ESTIMATES_LIST[@]}

# 从 SLURM_ARRAY_TASK_ID 映射到具体的参数索引
task_id=$SLURM_ARRAY_TASK_ID
m_idx=$(( (task_id / (NUM_MEASURES * NUM_ESTIMATES)) % NUM_MODELS ))
me_idx=$(( (task_id / NUM_ESTIMATES) % NUM_MEASURES ))
e_idx=$(( task_id % NUM_ESTIMATES ))

# 根据索引获取当前任务的具体参数
MODEL=${MODELS_LIST[$m_idx]}
MEASURE=${MEASURES_LIST[$me_idx]}
ESTIMATE=${ESTIMATES_LIST[$e_idx]}

log INFO "Test3: 作业数组任务ID: $task_id"
log INFO "--> 参数: 化合物=$COMPOUND, 模型=$MODEL, 测量=$MEASURE, 估算=$ESTIMATE, 滞后=$LAG_YEARS年"

# --- 执行分析 ---
# 使用与Test2相同的生产级采样参数
python -u Code/PYMC/main.py \
  --disease "$DISEASE" \
  --compound "$COMPOUND" \
  --model "$MODEL" \
  --lag "$LAG_YEARS" \
  --measure "$MEASURE" \
  --estimate "$ESTIMATE" \
  --sampling-mode "production" \
  --draws 4000 --tune 2000 --chains 4 --cores ${SLURM_CPUS_PER_TASK:-16} --target-accept 0.95 \
  --config-path "config.yaml" --verbose

log INFO "✅ Test3: 作业数组任务 $task_id 完成！"