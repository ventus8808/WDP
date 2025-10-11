#!/bin/bash
# Quick Reference Commands for Delta Bayesian Analysis
# ======================================================================

# --------------------------------------------------------------------
# 1. 本地测试单个癌症类型
# --------------------------------------------------------------------
bash Code/brms/test_Delta_bayesian_local.sh C00_C97

# 测试其他癌症类型
bash Code/brms/test_Delta_bayesian_local.sh C18_C21
bash Code/brms/test_Delta_bayesian_local.sh C50


# --------------------------------------------------------------------
# 2. 集群提交全部13个癌症类型
# --------------------------------------------------------------------
# 方法A: 自动发现并提交（推荐）
bash Code/brms/submit_Delta_bayesian_array.sh

# 方法B: 手动指定数组任务数
sbatch --array=0-12 Code/brms/submit_Delta_bayesian_array.sh


# --------------------------------------------------------------------
# 3. 监控任务状态
# --------------------------------------------------------------------
# 查看队列中的任务
squeue -u $USER

# 查看特定任务的实时输出
tail -f delta_bayesian_<JOB_ID>_<TASK_ID>.out

# 查看所有任务状态
sacct -j <JOB_ID> --format=JobID,State,ExitCode,Elapsed

# 统计完成情况
sacct -j <JOB_ID> --format=State | grep -c COMPLETED


# --------------------------------------------------------------------
# 4. 检查输出结果
# --------------------------------------------------------------------
# 列出所有生成的结果文件
ls -lh Result/brms_delta/*_delta.csv

# 统计每个文件的行数（应该约122行）
wc -l Result/brms_delta/*_delta.csv

# 查看某个结果文件的前几行
head -20 Result/brms_delta/C00_C97_delta.csv

# 检查所有文件的收敛性（Rhat_max列）
for f in Result/brms_delta/*_delta.csv; do
  echo "File: $(basename $f)"
  awk -F',' 'NR>1 {print $9}' "$f" | sort -rn | head -5
  echo ""
done


# --------------------------------------------------------------------
# 5. 故障排查
# --------------------------------------------------------------------
# 如果遇到TBB_CXX_TYPE编译错误
export TBB_CXX_TYPE=gcc
bash Code/brms/submit_Delta_bayesian_array.sh

# 如果遇到内存不足
sbatch --mem=64G Code/brms/submit_Delta_bayesian_array.sh

# 如果时间超限
sbatch --time=2-00:00:00 Code/brms/submit_Delta_bayesian_array.sh

# 查看错误日志
cat delta_bayesian_<JOB_ID>_<TASK_ID>.err

# 重新提交失败的任务（假设任务3失败）
sbatch --array=3 Code/brms/submit_Delta_bayesian_array.sh


# --------------------------------------------------------------------
# 6. 高级用法：仅运行特定lag
# --------------------------------------------------------------------
# 仅运行Lag=5分析
Rscript Code/brms/Delta_bayesian.R \
  --cancer-types=C00_C97 \
  --lag=5 \
  --chains=4 --iter=2000 --warmup=1000

# 运行多个癌症类型
Rscript Code/brms/Delta_bayesian.R \
  --cancer-types="C00_C97,C18_C21,C50" \
  --chains=4 --iter=2000 --warmup=1000


# --------------------------------------------------------------------
# 7. 结果汇总和分析
# --------------------------------------------------------------------
# 合并所有结果到单个文件（可选）
head -1 Result/brms_delta/C00_C97_delta.csv > Result/brms_delta/ALL_delta_combined.csv
for f in Result/brms_delta/*_delta.csv; do
  if [ "$f" != "Result/brms_delta/ALL_delta_combined.csv" ]; then
    tail -n +2 "$f" >> Result/brms_delta/ALL_delta_combined.csv
  fi
done
echo "Combined $(wc -l < Result/brms_delta/ALL_delta_combined.csv) rows"

# 提取显著结果（包含*的行）
grep '\*' Result/brms_delta/C00_C97_delta.csv

# 按Model类型统计
for model in EQI Air Water Land Built Social; do
  echo "$model models:"
  grep "\"$model\"" Result/brms_delta/C00_C97_delta.csv | wc -l
done
