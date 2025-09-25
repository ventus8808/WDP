#!/bin/bash

# Slurm集群配置快速检查脚本

echo "========================================="
echo "     Slurm集群配置检查"
echo "========================================="

# 1. 基本集群信息
echo -e "\n🏗️  集群基本信息:"
echo "Slurm版本: $(sinfo --version 2>/dev/null || echo '无法获取')"
echo "当前用户: $USER"
echo "默认账户: $(sacctmgr show user $USER format=account%20 -n 2>/dev/null | head -1 | xargs || echo '无法获取')"

# 2. 分区信息
echo -e "\n📊 可用分区信息:"
printf "%-15s %-8s %-10s %-12s %-15s %-10s\n" "分区名" "状态" "节点数" "时间限制" "默认内存/CPU" "最大内存/CPU"
echo "--------------------------------------------------------------------------------"

sinfo --format="%P %a %D %l" --noheader 2>/dev/null | while read partition state nodes timelimit; do
    partition=${partition%*}  # 移除末尾的*号
    
    # 获取分区的内存配置
    partition_info=$(scontrol show partition $partition 2>/dev/null)
    def_mem=$(echo "$partition_info" | grep -o "DefMemPerCPU=[0-9]*" | cut -d'=' -f2)
    max_mem=$(echo "$partition_info" | grep -o "MaxMemPerCPU=[0-9]*" | cut -d'=' -f2)
    
    # 格式化输出
    def_mem_fmt="${def_mem:-N/A}MB"
    max_mem_fmt="${max_mem:-N/A}MB"
    
    printf "%-15s %-8s %-10s %-12s %-15s %-10s\n" "$partition" "$state" "$nodes" "$timelimit" "$def_mem_fmt" "$max_mem_fmt"
done

# 3. 节点类型统计
echo -e "\n🖥️  节点配置统计:"
echo "节点总数: $(sinfo -h -o "%D" | awk '{sum+=$1} END {print sum}')"

# 按节点类型分组显示
sinfo --format="%n %c %m %f" --noheader 2>/dev/null | \
awk '{
    key = $2 "核心_" int($3/1024) "GB"
    if ($4 != "(null)") key = key "_" $4
    count[key]++
} 
END {
    for (config in count) {
        print config ": " count[config] "个节点"
    }
}' | sort

# 4. 资源限制检查
echo -e "\n⚖️  资源限制:"
if command -v sacctmgr &> /dev/null; then
    echo "用户资源限制:"
    sacctmgr show user $USER format=user,account,maxcpus,maxwall,maxmem -n 2>/dev/null || echo "无法获取用户限制信息"
    
    echo -e "\nQOS限制:"
    sacctmgr show qos format=name%15,maxwall%12,maxcpus%8,maxmem%10 2>/dev/null | head -5 || echo "无法获取QOS信息"
fi

# 5. 当前队列状态
echo -e "\n📋 当前队列状态:"
running_jobs=$(squeue -h -t RUNNING | wc -l)
pending_jobs=$(squeue -h -t PENDING | wc -l)
echo "运行中作业: $running_jobs"
echo "等待中作业: $pending_jobs"

if [ $pending_jobs -gt 0 ]; then
    echo -e "\n等待作业的主要原因:"
    squeue -t PENDING -o "%.10i %.15u %.12P %.20r" -h 2>/dev/null | \
    awk '{print $4}' | sort | uniq -c | sort -nr | head -5
fi

# 6. 推荐的内存配置
echo -e "\n💡 内存配置建议:"

# 计算每个分区的推荐内存配置
sinfo --format="%P" --noheader 2>/dev/null | sort | uniq | while read partition; do
    partition=${partition%*}  # 移除末尾的*号
    
    # 获取该分区节点的典型配置
    typical_mem=$(sinfo -p $partition --format="%m" --noheader 2>/dev/null | head -1)
    typical_cpu=$(sinfo -p $partition --format="%c" --noheader 2>/dev/null | head -1)
    
    if [ ! -z "$typical_mem" ] && [ ! -z "$typical_cpu" ]; then
        # 计算推荐的per-cpu内存 (留20%给系统)
        mem_gb=$((typical_mem/1024))
        recommended_per_cpu=$((mem_gb*800/typical_cpu))  # 80% * 1000MB/GB
        
        echo "分区 $partition:"
        echo "  - 典型配置: ${typical_cpu}核心, ${mem_gb}GB内存"
        echo "  - 推荐 --mem-per-cpu: ${recommended_per_cpu}M 或 $(echo "scale=1; $recommended_per_cpu/1024" | bc)G"
        echo "  - 推荐 --mem (整节点): $((mem_gb*4/5))G"
    fi
done

# 7. 常用命令提醒
echo -e "\n🔧 常用检查命令:"
echo "查看作业状态: squeue -u \$USER"
echo "查看历史作业: sacct -u \$USER --format=JobID,JobName,Partition,AllocCPUS,State,MaxRSS,Elapsed"
echo "监控运行中作业: sstat -j <job_id> --format=JobID,MaxRSS,AveCPU"
echo "取消作业: scancel <job_id>"
echo "查看节点详情: scontrol show node <node_name>"

# 8. 脚本模板推荐
echo -e "\n📝 根据您的集群，推荐的脚本模板:"

# 找到最常用的分区
common_partition=$(sinfo --format="%P %D" --noheader 2>/dev/null | sort -k2 -nr | head -1 | awk '{print $1}' | sed 's/\*$//')

if [ ! -z "$common_partition" ]; then
    echo "
# 推荐的基础模板 (适用于分区: $common_partition)
#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --partition=$common_partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# 您的程序命令
./your_program
"
fi

echo -e "\n========================================="
echo "      配置检查完成"
echo "========================================="