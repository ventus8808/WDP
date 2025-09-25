# 集群环境分析报告

## 📋 基本信息概览

| 项目 | 详情 |
|------|------|
| **主机名** | login04 |
| **操作系统** | CentOS Linux 7 (Core) |
| **内核版本** | 3.10.0-957.el7.x86_64 |
| **架构** | x86_64 |
| **Slurm版本** | 22.05.8-2.3.2-102-20250620 |
| **当前用户** | acf4pijnzl |

## 🖥️ 硬件配置

### CPU配置
- **型号**: Hygon C86 7185 32-core Processor
- **当前节点核心数**: 2 (登录节点)
- **计算节点核心数**: 32核心/节点

### 内存配置
- **登录节点总内存**: 251G
- **已用内存**: 46G
- **可用内存**: 20G
- **计算节点内存**: 123GB/节点

### 存储信息
| 挂载点 | 文件系统 | 总大小 | 已用 | 可用 | 使用率 |
|--------|----------|--------|------|------|--------|
| / | /dev/sda5 | 327G | 103G | 224G | 32% |
| /boot | /dev/sda1 | 477M | 216M | 232M | 49% |
| /opt | /dev/sda2 | 1.4T | 131G | 1.2T | 11% |
| /data | /dev/sdb | 37T | 349G | 35T | 1% |

## 🌐 网络配置

### 网络接口
- **管理网络**: enp97s0f0 (10.15.200.4/16)
- **InfiniBand**: ib0 (11.2.200.4/8) - 高速互联网络

### 系统状态
- **运行时间**: 126天 2小时15分钟
- **当前负载**: 8.34, 8.63, 8.66 (1分钟, 5分钟, 15分钟平均)
- **当前用户数**: 148个用户在线

## 🚀 Slurm集群配置

### 分区信息
| 分区名 | 状态 | 节点数 | 时间限制 | 默认内存/CPU | 节点配置 |
|--------|------|--------|----------|--------------|----------|
| **kshctest** | up | 258 | 3天 | 3500MB | 32核心/123GB |
| **kshdtest** | up | 59 | 3天 | 3500MB | 32核心/123GB |

### 集群规模
- **总节点数**: 317个计算节点
- **总CPU核心数**: 317 × 32 = 10,144 核心
- **总内存**: 317 × 123GB ≈ 39TB

## 💡 内存分配建议

### 推荐配置策略

#### 1. 使用 per-cpu 内存分配（推荐）
```bash
#SBATCH --mem-per-cpu=3G    # 推荐值，接近默认的3500MB
```

#### 2. 使用整节点内存分配
```bash
#SBATCH --mem=98G           # 为整个节点分配内存（预留25GB给系统）
```

### 常见作业类型配置

#### 🔬 小型测试作业
```bash
#!/bin/bash
#SBATCH --job-name=test_job
#SBATCH --partition=kshctest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G    # 总内存: 4×2G=8G
#SBATCH --time=01:00:00
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err
```

#### 🧮 中等规模计算
```bash
#!/bin/bash
#SBATCH --job-name=compute_job
#SBATCH --partition=kshctest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G    # 总内存: 16×3G=48G
#SBATCH --time=12:00:00
#SBATCH --output=compute_%j.out
#SBATCH --error=compute_%j.err
```

#### 🚄 高性能计算（单节点最大配置）
```bash
#!/bin/bash
#SBATCH --job-name=hpc_job
#SBATCH --partition=kshctest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=98G           # 使用几乎整个节点的内存
#SBATCH --time=1-00:00:00   # 1天
#SBATCH --output=hpc_%j.out
#SBATCH --error=hpc_%j.err
```

#### 🔄 MPI并行作业
```bash
#!/bin/bash
#SBATCH --job-name=mpi_job
#SBATCH --partition=kshctest
#SBATCH --nodes=4           # 4个节点
#SBATCH --ntasks-per-node=32 # 每节点32个MPI进程
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G    # 总内存: 4×32×3G=384G
#SBATCH --time=06:00:00
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err

mpirun ./my_mpi_program
```

## 🛠️ 实用命令速查

### 作业管理
```bash
# 提交作业
sbatch my_script.sh

# 查看队列状态
squeue -u $USER

# 查看作业详情
scontrol show job <job_id>

# 取消作业
scancel <job_id>

# 查看历史作业
sacct -u $USER --format=JobID,JobName,Partition,AllocCPUS,State,MaxRSS,Elapsed

# 监控运行中作业内存使用
sstat -j <job_id> --format=JobID,MaxRSS,AveCPU
```

### 集群信息查询
```bash
# 查看分区信息
sinfo -l

# 查看节点详情
scontrol show node <node_name>

# 查看可用资源
sinfo -N -l
```

## ⚠️ 重要注意事项

### 内存配置原则
1. **不能同时使用** `--mem` 和 `--mem-per-cpu`
2. **默认推荐**: `--mem-per-cpu=3G`（接近系统默认的3500MB）
3. **最大单节点内存**: 约98-100GB（为系统预留20-25GB）
4. **测试先行**: 从小内存开始测试，逐步增加到合适值

### 时间限制
- **默认最大时间**: 3天
- **建议策略**: 根据实际需求设置，避免占用资源过长

### 分区选择
- **kshctest**: 主要计算分区（258节点）
- **kshdtest**: 较小测试分区（59节点）

## 📊 资源使用建议

### 内存效率最大化
- 对于内存密集型应用：优先使用 `--mem` 指定总内存
- 对于CPU密集型应用：使用 `--mem-per-cpu` 保持灵活性
- 总是为系统预留20-25%的内存

### 作业优化
1. **先小规模测试**，确定合适的资源需求
2. **使用 `sacct` 分析**已完成作业的实际资源使用
3. **避免资源浪费**，不要申请过多不必要的资源
4. **合理使用时间限制**，避免长时间占用节点

---

*报告生成时间: 2025年9月25日*  
*集群状态: 正常运行，当前无排队作业*