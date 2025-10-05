#!/bin/bash
# BRMS Complete Analysis - Simple Version
# 简单版本：C00-C97完整EQI分析
#SBATCH --partition=kshctest
#SBATCH --job-name=BRMS_Complete
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=6:00:00
#SBATCH --array=0-71
#SBATCH --output=logs/BRMS_%A_%a.out
#SBATCH --error=logs/BRMS_%A_%a.err

set -e

# 进入项目目录
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p logs

echo "任务 $SLURM_ARRAY_TASK_ID 开始: $(date)"

# 激活环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pymc

# 参数数组
DOMAINS=(total air water land built sociodemographic)
RUCC_TYPES=(All Urban Rural)
PERIODS=(2006-2010 2011-2015)  
LAGS=(5 10)

# 计算当前任务参数 (6×3×2×2=72)
task=$SLURM_ARRAY_TASK_ID
domain=${DOMAINS[$((task/12))]}
rucc=${RUCC_TYPES[$((task%12/4))]}
period=${PERIODS[$((task%4/2))]}
lag=${LAGS[$((task%2))]}

echo "参数: 域=$domain, RUCC=$rucc, 时间=$period, 滞后=${lag}年"

# 设置RUCC筛选器
case $rucc in
  All) rucc_filter="null";;
  Urban) rucc_filter="[1,2,3]";;
  Rural) rucc_filter="[4,5,6,7,8,9]";;
esac

# 设置EQI期间
case $period in
  2006-2010) eqi_period="0005";;
  2011-2015) eqi_period="0610";;
esac

scenario_name="AllCancer_${domain}_${rucc}_${period//-/}_Lag${lag}"

echo "场景名称: $scenario_name"

# 数据准备 (只在任务0执行)
if [ $task -eq 0 ]; then
  echo "准备数据..."
  Rscript Code/brms/01_prepare_data.R
fi

# 等待数据准备完成
if [ $task -ne 0 ] && [ $task -lt 10 ]; then
  sleep $((task * 5))
fi

# Python脚本执行分析
python3 << EOF
import yaml
import subprocess
import sys

# 读取配置
with open('config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# RUCC筛选器处理
rucc_filter = "$rucc_filter"
if rucc_filter == "null":
    rucc_value = None
else:
    rucc_value = eval(rucc_filter)

# 创建场景
scenario = {
    'name': '$scenario_name',
    'active': True,
    'cancer_type': 'C00_C97',
    'eqi_period': '$eqi_period',
    'time_period': '$period',
    'lag_years': $lag,
    'rucc_filter': rucc_value,
    'domain': '$domain',
    'formula_key': 'total_eqi_quintile',
    'family': 'gaussian',
    'use_midpoint': False,
    'output_suffix': ''
}

# 更新配置
config['brms_analysis']['scenarios'] = [scenario]

# 保存临时配置
temp_config = '/tmp/brms_$task.yaml'
with open(temp_config, 'w') as f:
    yaml.dump(config, f)

print(f"开始分析: $scenario_name")

try:
    # 运行R分析
    result = subprocess.run([
        'Rscript', 'Code/brms/02_run_brms_model.R',
        '--scenario', '$scenario_name'
    ], check=True)
    print("✅ 分析成功")
    
except Exception as e:
    print(f"❌ 分析失败: {e}")
    sys.exit(1)
finally:
    import os
    if os.path.exists(temp_config):
        os.remove(temp_config)
EOF

echo "任务 $task ($scenario_name) 完成: $(date)"

# 最后一个任务做汇总
if [ $task -eq 71 ]; then
  echo "等待其他任务完成..."
  sleep 120
  echo "执行结果汇总..."
  Rscript Code/brms/04_process_results.R || echo "汇总失败但继续"
  echo "🎉 全部完成!"
fi