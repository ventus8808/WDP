# Delta Bayesian Analysis - Cluster Deployment Guide

## 概述

Delta_bayesian.R 实现了EQI变化与癌症死亡率变化之间关系的三阶段贝叶斯分析。本文档说明如何在集群上批量运行分析。

## 文件说明

### 核心脚本
- **Delta_bayesian.R** - 主分析脚本，支持命令行参数
- **submit_Delta_bayesian_array.sh** - Slurm数组任务提交脚本（集群使用）
- **test_Delta_bayesian_local.sh** - 本地测试脚本（快速验证）

### 输入数据
- `Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv` - 81,549行，包含Lag=5和Lag=10数据

### 输出结果
- `Result/brms_delta/{ICD_Code}_delta.csv` - 每个癌症类型一个文件
- 例如：`C00_C97_delta.csv`, `C18_C21_delta.csv`

## 使用方法

### 1. 本地测试（推荐先测试）

测试单个癌症类型（使用较少迭代）：

```bash
# 在本地MacBook/Linux上测试
bash Code/brms/test_Delta_bayesian_local.sh C00_C97

# 或指定其他癌症类型
bash Code/brms/test_Delta_bayesian_local.sh C18_C21
```

测试参数：
- Chains: 2
- Iterations: 500 (warmup: 250)
- 运行时间：约5-10分钟

### 2. 集群批量提交

在集群上提交所有13个癌症类型的分析：

```bash
# 方法1: 直接运行脚本（推荐）
# 脚本会自动发现所有癌症类型并提交数组任务
bash Code/brms/submit_Delta_bayesian_array.sh

# 方法2: 手动sbatch提交
sbatch --array=0-12 Code/brms/submit_Delta_bayesian_array.sh
```

### 3. 监控任务进度

```bash
# 查看队列中的任务
squeue -u $USER

# 查看特定任务的输出（实时）
tail -f delta_bayesian_<JOB_ID>_<TASK_ID>.out

# 查看所有任务的状态
sacct -j <JOB_ID> --format=JobID,State,ExitCode,Elapsed
```

## 集群资源配置

### 默认配置（submit_Delta_bayesian_array.sh）
```bash
#SBATCH --partition=kshctest
#SBATCH --cpus-per-task=16        # 16核心（自动使用80% = 12-13核）
#SBATCH --mem=48G                 # 48GB内存
#SBATCH --time=1-00:00:00         # 24小时时限
```

### 采样参数
- Chains: 4
- Iterations: 2000 (warmup: 1000)
- adapt_delta: 0.95
- max_treedepth: 12
- 每个癌症预计运行时间：8-12小时

## 预期输出

### 输出文件结构
每个癌症类型生成一个CSV文件，包含：
- **行数**：约122行（61行/lag × 2个lag值）
- **列数**：9列

### 列说明
1. `ICD_Code` - 癌症类型代码（如C00_C97）
2. `Model` - 模型名称（EQI, Air, RUCC1_EQI, Air_Multi等）
3. `Lag` - 滞后年数（5或10）
4. `MRD_Q_Improved` - Improved组效应（均值和95%置信区间）
5. `MRD_Q_Worsened` - Worsened组效应
6. `Intercept_Baseline_Change` - 基线死亡率变化
7. `Control_delta_Smoking_Rate` - 吸烟率控制变量效应
8. `N_Counties` - 样本县数
9. `Rhat_max` - 最大Rhat值（收敛诊断）

### 模型类型分布（每个癌症×每个lag）
- **Phase 1 (National)**: 11行
  - 1行 EQI
  - 5行 单领域（Air, Water, Land, Built, Social）
  - 5行 Multi模型（Air_Multi, Water_Multi等）
  
- **Phase 2 (Stratified)**: 44行
  - 4个RUCC层级 × 11个模型
  
- **Phase 3 (Within-RUCC)**: 6行
  - 1行 Within_RUCC_EQI
  - 5行 Within_RUCC单领域

**总计**: 61行/lag × 2个lag = 122行/癌症

## 示例输出

```csv
ICD_Code,Model,Lag,MRD_Q_Improved,MRD_Q_Worsened,...
C00_C97,EQI,5,"-0.19(-1.80,1.43)","1.50(-0.32,3.28)",...
C00_C97,Air,5,"-1.15(-2.45,0.35)","0.44(-1.53,2.80)",...
C00_C97,Air_Multi,5,"-0.53(-1.71,0.88)","-0.74(-2.73,1.76)",...
C00_C97,RUCC1_EQI,5,"-0.67(-2.23,0.95)","2.41(0.08,4.48)*",...
C00_C97,Within_RUCC_EQI,5,"0.15(-1.32,1.62)","0.89(-0.68,2.34)",...
C00_C97,EQI,10,"0.45(-0.98,1.87)","2.15(0.62,3.68)**",...
...
```

## 故障排查

### 任务失败常见原因

1. **内存不足**
   ```bash
   # 增加内存限制
   sbatch --mem=64G Code/brms/submit_Delta_bayesian_array.sh
   ```

2. **时间超限**
   ```bash
   # 延长时间限制
   sbatch --time=2-00:00:00 Code/brms/submit_Delta_bayesian_array.sh
   ```

3. **收敛问题（Rhat > 1.1）**
   - 检查输出CSV中的Rhat_max列
   - 如果Rhat > 1.1，考虑增加迭代次数或调整adapt_delta

### 查看错误日志

```bash
# 查看stderr输出
cat delta_bayesian_<JOB_ID>_<TASK_ID>.err

# 查看stdout输出
cat delta_bayesian_<JOB_ID>_<TASK_ID>.out
```

## 高级用法

### 仅运行特定lag值

```bash
# 仅Lag=5
Rscript Code/brms/Delta_bayesian.R --cancer-types=C00_C97 --lag=5

# 仅Lag=10
Rscript Code/brms/Delta_bayesian.R --cancer-types=C00_C97 --lag=10
```

### 运行多个癌症类型

```bash
Rscript Code/brms/Delta_bayesian.R \
  --cancer-types="C00_C97,C18_C21,C50" \
  --chains=4 --iter=2000 --warmup=1000
```

### 自定义输出目录

```bash
Rscript Code/brms/Delta_bayesian.R \
  --cancer-types=C00_C97 \
  --output-dir="Result/brms_delta_test"
```

## 数据说明

### Lag设计
- **Lag=5**: EQI变化(2000-2005→2006-2010) → AAMR变化(2006-2010→2011-2015)
- **Lag=10**: EQI变化(2000-2005→2006-2010) → AAMR变化(2011-2015→2016-2020)

### 变化分类
- **Improved**: EQI quintile上升（环境质量改善）
- **Stable**: EQI quintile不变
- **Worsened**: EQI quintile下降（环境质量恶化）

### Within-RUCC相对变化
Phase 3使用RUCC内相对排名变化，而非全国绝对变化。例如：
- 某县全国EQI排名Stable（Q3→Q3）
- 但在RUCC=1内排名Worsened（Q2→Q4，被同类城市超越）
- Phase 3检验这种"相对剥夺"效应对健康的影响

## 参考资料

- 主脚本文档：查看 `Code/brms/Delta_bayesian.R` 开头的注释
- Interval分析参考：`Code/brms/cmdstan_interval_runner.R`
- 数据清洗脚本：`Code/Clean/EQI_AAMR_Delta.py`

## 预期完成时间

- **单个癌症类型**: 8-12小时（包含Lag=5和Lag=10）
- **全部13个癌症**: 并行运行约12-16小时（取决于集群负载）
- **总输出**: 13个CSV文件，约1,586行数据
