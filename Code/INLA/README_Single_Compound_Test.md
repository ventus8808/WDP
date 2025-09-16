# WDP 单化合物测试脚本使用说明

## 文件说明

1. **`run_single_compound_test.sh`** - 主要运行脚本（服务器版本，无Docker）
2. **`submit_single_compound_test.sh`** - SLURM批处理提交脚本

## 快速开始

### 1. 基本使用（推荐）
```bash
# 提交到SLURM队列（推荐方式）
sbatch Code/INLA/submit_single_compound_test.sh
```

### 2. 直接运行（用于调试）
```bash
# 使用默认参数（2,4-D化合物）
bash Code/INLA/run_single_compound_test.sh

# 或指定其他化合物
bash Code/INLA/run_single_compound_test.sh "5" "Atrazine" "C81-C96"
```

## 参数说明

运行脚本支持以下参数（按顺序）：

1. **COMPOUND_ID** - 化合物ID（默认: 2）
2. **COMPOUND_NAME** - 化合物名称（默认: 2,4-D）
3. **DISEASE_CODE** - 疾病代码（默认: C81-C96）
4. **MEASURE_TYPES** - 测量类型（默认: Weight,Density）
5. **ESTIMATE_TYPES** - 估计类型（默认: min,avg,max）
6. **LAG_YEARS** - 滞后年数（默认: 5）
7. **MODEL_TYPES** - 模型类型（默认: M0,M1,M2,M3）

## 预期结果

### 输出文件
- 位置：`Result/Filter/`
- 格式：`Results_C81-C96_2,4-D_YYYYMMDD_HHMMSS.csv`
- 内容：包含所有模型组合的分析结果

### 预期结果数量
对于单个化合物的完整测试：
- 3 估计类型 × 2 测量类型 × 4 模型类型 = **24条结果**

### 运行时间
- 预计：30-60分钟
- 内存需求：8GB
- CPU：4核

## 资源配置

SLURM作业配置：
```bash
#SBATCH --cpus-per-task=4    # 4个CPU核心
#SBATCH --mem=8G             # 8GB内存
#SBATCH --time=02:00:00      # 2小时时间限制
#SBATCH --partition=kshdtest # 测试分区
```

## 监控和调试

### 查看作业状态
```bash
# 查看队列中的作业
squeue -u $(whoami)

# 查看作业详情
scontrol show job <JOB_ID>
```

### 查看日志
```bash
# 实时查看输出日志
tail -f slurm_logs/single_compound_test-<JOB_ID>.out

# 查看错误日志
tail -f slurm_logs/single_compound_test-<JOB_ID>.err
```

### 检查结果
```bash
# 查看最新的结果文件
ls -la Result/Filter/Results_C81-C96_2,4-D_*.csv | tail -n 3

# 查看结果统计
LATEST=$(ls -t Result/Filter/Results_C81-C96_2,4-D_*.csv | head -n 1)
echo "结果数量: $(($(wc -l < "$LATEST") - 1))"
head -n 5 "$LATEST"  # 显示前5行
```

## 故障排除

### 常见问题

1. **R包缺失**
   ```bash
   # 检查R包是否可用
   R -e "library(INLA); library(dplyr)"
   ```

2. **数据文件不存在**
   ```bash
   # 检查关键数据文件
   ls -la Data/Processed/Pesticide/PNSP.csv
   ls -la Data/Processed/CDC/C81-C96.csv
   ```

3. **内存不足**
   - 增加SLURM内存限制：`#SBATCH --mem=16G`
   - 或减少并行度：`#SBATCH --cpus-per-task=2`

4. **权限问题**
   ```bash
   # 确保脚本可执行
   chmod +x Code/INLA/run_single_compound_test.sh
   chmod +x Code/INLA/submit_single_compound_test.sh
   ```

### 调试模式

如果需要详细调试信息，可以修改R脚本调用：
```bash
# 在运行脚本中添加 --verbose 参数
Rscript ... --verbose
```

## 下一步

测试成功后，可以：

1. **测试多个化合物** - 修改脚本参数测试不同化合物
2. **批量处理** - 创建数组作业处理多个化合物
3. **扩展到其他疾病** - 修改DISEASE_CODE参数
4. **优化资源** - 根据实际运行情况调整SLURM参数

## 联系信息

如有问题，请检查：
- SLURM日志文件
- R错误输出
- 数据文件完整性