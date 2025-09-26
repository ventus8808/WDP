#!/bin/bash
# 检查SLURM分区和提交任务脚本

echo "=== 检查SLURM分区 ==="
sinfo

echo ""
echo "=== 当前用户作业状态 ==="
squeue -u $USER

echo ""
echo "=== 如果kshctest分区不可用，可以尝试以下命令 ==="
echo "1. 查看默认分区: scontrol show config | grep DefaultPartition"
echo "2. 查看可用分区: sinfo -o \"%P %a %l %F\""
echo "3. 临时提交命令（无分区指定）:"
echo "   sbatch --job-name=WONDER_R_Test --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --mem-per-cpu=3G --time=1:00:00 test_single_r_job.sh"

echo ""
echo "=== 尝试提交任务 ==="
if sbatch test_single_r_job.sh; then
    echo "✓ 任务提交成功"
else
    echo "❌ 任务提交失败，请检查分区设置"
    echo "尝试查看系统默认分区："
    scontrol show config | grep -i partition || echo "无法获取默认分区信息"
fi