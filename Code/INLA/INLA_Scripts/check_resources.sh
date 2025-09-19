#!/bin/bash
# ========================
# SLURM Resource Information and Optimization Script
# 帮助确定最优的资源配置
# ========================

echo "🔍 SLURM集群资源信息查询"
echo "=========================="

# 1. 查看可用的分区和资源限制
echo "📊 [1] 可用分区信息："
echo "sinfo -o \"%20P %5a %10l %6D %6t %N\""
sinfo -o "%20P %5a %10l %6D %6t %N" 2>/dev/null || echo "需要在SLURM集群上运行"

echo ""
echo "📊 [2] 分区详细限制："
echo "sinfo -Nel"
sinfo -Nel 2>/dev/null || echo "需要在SLURM集群上运行"

echo ""
echo "📊 [3] 当前队列状态："
echo "squeue -u \$USER"
squeue -u $USER 2>/dev/null || echo "需要在SLURM集群上运行"

echo ""
echo "📊 [4] 节点详细信息："
echo "scontrol show partition kshdtest"
scontrol show partition kshdtest 2>/dev/null || echo "需要在SLURM集群上运行"

echo ""
echo "📊 [5] 内存和CPU信息："
echo "sinfo -o \"%20N %5c %10m %25f %10G\""
sinfo -o "%20N %5c %10m %25f %10G" 2>/dev/null || echo "需要在SLURM集群上运行"

echo ""
echo "🎯 资源配置建议："
echo "=================="

echo "基于INLA分析的特点，推荐配置："
echo ""
echo "💡 [轻量级] - 快速测试："
echo "   --cpus-per-task=4"
echo "   --mem-per-cpu=2G"
echo "   --time=1:00:00"
echo "   总内存: 8GB, 适合单个模型快速验证"
echo ""
echo "💪 [标准配置] - 生产分析："
echo "   --cpus-per-task=8" 
echo "   --mem-per-cpu=3G"
echo "   --time=2:00:00"
echo "   总内存: 24GB, 适合完整分析流程"
echo ""
echo "🚀 [高性能] - 复杂模型："
echo "   --cpus-per-task=16"
echo "   --mem-per-cpu=4G" 
echo "   --time=4:00:00"
echo "   总内存: 64GB, 适合大规模数据或复杂空间模型"
echo ""

echo "🔧 INLA特定优化："
echo "=================="
echo "• INLA主要是单线程计算，多核收益有限"
echo "• 内存需求取决于数据大小和空间复杂度"
echo "• 建议 CPU:Memory 比例为 1:2GB 到 1:4GB"
echo "• 空间模型需要更多内存用于邻接矩阵"

echo ""
echo "📝 使用方法："
echo "============"
echo "1. 在集群上运行此脚本查看实际资源限制"
echo "2. 根据数据规模选择合适的配置"
echo "3. 修改 run_single_compound.sh 中的 SBATCH 参数"
echo "4. 监控作业运行情况调整配置"