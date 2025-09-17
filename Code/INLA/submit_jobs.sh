#!/bin/bash

#================================================================================
# WONDER 研究计划 - 批量任务提交脚本
#
# 使用方法:
# 1. 在下面的 DISEASES 和 COMPOUNDS 数组中填入您想分析的条目。
# 2. 保存文件。
# 3. 在终端运行 ./submit_jobs.sh
#================================================================================

echo "开始批量提交 WONDER 分析任务..."

# --- 在这里配置您要分析的疾病和化合物 ---

# 疾病代码列表 (参考 WDP/Data/Processed/CDC/ 目录下的文件名)
DISEASES=(
  "C81-C96"
  "C25"
  "C50"
  "C34"
)

# 化合物ID列表 (参考 WDP/Data/Processed/Pesticide/mapping.csv 文件中的 'compound_id' 列)
COMPOUNDS=(
  1
  5
  10
  23
  45
)
# ----------------------------------------------------

# 循环遍历所有疾病和化合物的组合
for disease in "${DISEASES[@]}"; do
  for compound in "${COMPOUNDS[@]}"; do
    # 为每个任务动态生成一个唯一的、信息丰富的任务名
    JOB_NAME="WDP-${disease}-C${compound}"

    echo "正在提交任务: ${JOB_NAME}"

    # 使用 sbatch 命令提交工作脚本
    # --job-name 参数会覆盖工作脚本中默认的任务名
    # 参数 ${disease} 和 ${compound} 会被传递给 run_single_job.sh
    sbatch --job-name="${JOB_NAME}" run_single_job.sh "${disease}" "${compound}"

    # 短暂暂停1秒，避免瞬间提交大量任务给 Slurm 控制器造成压力
    sleep 1
  done
done

echo "-------------------------------------------------"
echo "所有任务已提交完毕。"
echo "使用 'squeue -u $USER' 命令查看您的任务队列。"
echo "日志文件将保存在 'WDP/slurm_logs' 目录下。"
echo "分析结果将保存在 'WDP/Result' 目录下。"
