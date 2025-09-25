#!/bin/bash

# 服务器系统信息快速查看脚本
# 适用于HCP集群环境

echo "========================================="
echo "        服务器系统信息概览"
echo "========================================="

# 1. 基本系统信息
echo -e "\n🖥️  基本系统信息:"
echo "主机名: $(hostname)"
echo "系统: $(uname -s)"
echo "内核版本: $(uname -r)"
echo "架构: $(uname -m)"

# 2. 发行版信息
echo -e "\n📋 发行版信息:"
if [ -f /etc/os-release ]; then
    source /etc/os-release
    echo "发行版: $PRETTY_NAME"
    echo "版本ID: $VERSION_ID"
elif [ -f /etc/redhat-release ]; then
    echo "发行版: $(cat /etc/redhat-release)"
elif [ -f /etc/debian_version ]; then
    echo "发行版: Debian $(cat /etc/debian_version)"
fi

# 3. 系统运行时间和负载
echo -e "\n⏱️  系统运行信息:"
uptime_info=$(uptime)
echo "运行时间和负载: $uptime_info"

# 4. CPU信息
echo -e "\n🔧 CPU信息:"
cpu_model=$(grep "model name" /proc/cpuinfo | head -1 | cut -d':' -f2 | sed 's/^ *//')
cpu_cores=$(nproc)
echo "CPU型号: $cpu_model"
echo "CPU核心数: $cpu_cores"

# 5. 内存信息
echo -e "\n💾 内存信息:"
mem_info=$(free -h | grep "Mem:")
echo "内存: $mem_info"

# 6. 磁盘使用情况
echo -e "\n💿 磁盘使用情况:"
df -h | grep -E "(^/dev/|Filesystem)" | grep -v tmpfs

# 7. 网络接口信息
echo -e "\n🌐 网络接口:"
ip addr show | grep -E "(^[0-9]+:|inet )" | grep -v "127.0.0.1"

# 8. HCP相关服务检查
echo -e "\n☁️  HCP相关服务状态:"
services=("consul" "nomad" "vault" "docker" "containerd")
for service in "${services[@]}"; do
    if systemctl is-active --quiet $service 2>/dev/null; then
        status="✅ 运行中"
    elif systemctl list-unit-files --type=service | grep -q "^$service.service"; then
        status="❌ 已停止"
    else
        status="⚪ 未安装"
    fi
    printf "%-12s: %s\n" "$service" "$status"
done

# 9. 容器运行时检查
echo -e "\n🐳 容器信息:"
if command -v docker &> /dev/null; then
    echo "Docker版本: $(docker --version 2>/dev/null || echo '未运行')"
    echo "运行的容器数: $(docker ps -q 2>/dev/null | wc -l)"
fi

if command -v kubectl &> /dev/null; then
    echo "Kubectl版本: $(kubectl version --client --short 2>/dev/null || echo '配置异常')"
fi

# 10. 最近的系统日志错误
echo -e "\n📝 最近的系统错误日志 (最近10条):"
if command -v journalctl &> /dev/null; then
    journalctl -p err -n 10 --no-pager -q 2>/dev/null | head -10 || echo "无权限或无错误日志"
else
    tail -10 /var/log/messages 2>/dev/null | grep -i error || echo "无法访问系统日志"
fi

echo -e "\n========================================="
echo "           信息收集完成"
echo "========================================="

