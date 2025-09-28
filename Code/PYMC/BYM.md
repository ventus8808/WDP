# BYM2 Spatiotemporal Model for WDP

This document describes the Bayesian spatiotemporal model implementation using PyMC for the WDP (WONDER Data Pipeline) project.

## Model Overview

The BYM2 (Besag-York-Mollié version 2) model is a Bayesian spatiotemporal model that accounts for:
- Spatial autocorrelation between neighboring counties
- Temporal trends across years  
- Both structured and unstructured spatial effects
- Covariate effects from PCA components

## Mathematical Formulation

For each county $i$ and time $t$, we model the observed deaths $y_{i,t}$ as:

$$y_{i,t} \sim \text{Poisson}(\lambda_{i,t})$$

where $\lambda_{i,t} = E_{i,t} \cdot \exp(\eta_{i,t})$ and $E_{i,t}$ is the expected number of deaths.

The linear predictor is:
$$\eta_{i,t} = \alpha + \beta \cdot \text{pesticide}_{i,t-lag} + \mathbf{X}_{i,t}^T \boldsymbol{\gamma} + u_i + v_i + \phi_t$$

Where:
- $\alpha$: intercept
- $\beta$: lagged pesticide exposure effect
- $\mathbf{X}_{i,t}^T \boldsymbol{\gamma}$: PCA covariate effects
- $u_i$: structured spatial random effect (CAR model)
- $v_i$: unstructured spatial random effect  
- $\phi_t$: temporal random effect (RW1 model)

## PCA Covariates (Based on README_PCA.md)

### Social Vulnerability Index (SVI)
- **SVI_PC1**: First principal component capturing socioeconomic vulnerability
  - Higher values = greater vulnerability
  - Loadings: poverty (+), income (-), education (-)

### Environmental Factors  
- **ENV_PC1**: Climate/precipitation component (33.1% variance)
  - High loadings: humidity, precipitation (seasonal)
- **ENV_PC2**: Temperature/radiation component (28.3% variance)  
  - High loadings: temperature, solar radiation
- **ENV_PC3**: Wind/seasonal component (13.3% variance)
  - High loadings: wind speed, seasonal precipitation patterns

## Model Variants

Based on PCA results, four model variants are supported:

1. **M0**: Base model
   - Pesticide exposure + spatial + temporal effects only
   
2. **M1**: M0 + Social vulnerability  
   - Adds: SVI_PC1
   
3. **M2**: M0 + Environmental factors
   - Adds: ENV_PC1, ENV_PC2
   
4. **M3**: Full model
   - Adds: SVI_PC1, ENV_PC1, ENV_PC2

## File Structure

- `Utils_Data.py`: Data loading and preprocessing
- `Utils_Model.py`: BYM2 model fitting (M0-M3)
- `Utils_Result.py`: Results extraction and formatting  
- `Utils_Others.py`: Utility functions
- `Main.py`: Command-line interface and workflow orchestration


## Usage Example

```bash
# Run analysis for C81-C96 with 2,4-D pesticide
python Main.py --disease C81-C96 --compound 24D --model M3 --lag 5

# Batch analysis for multiple models
python Main.py --disease C81-C96 --compound 24D --model M0,M1,M2,M3 --lag 5,10
```

## 数学模型

### BYM2空间时间模型

对于县域 $i$，年份 $t$：

$$O_{i,t} \sim \text{Poisson}(\lambda_{i,t})$$

$$\log(\lambda_{i,t}) = \log(E_{i,t}) + \eta_{i,t}$$

$$\eta_{i,t} = \alpha + \beta \log(X_{i,t-lag} + c) + \sum_j \gamma_j Z_{j,i,t} + \phi_i + \theta_i + \tau_t$$

其中：
- $O_{i,t}$: 观测死亡数
- $E_{i,t}$: 预期死亡数（偏移量）
- $X_{i,t-lag}$: 滞后农药暴露
- $\phi_i$: 结构化空间效应（CAR先验）
- $\theta_i$: 非结构化空间效应（IID先验）
- $\tau_t$: 时间效应（RW1先验）

### 先验分布

```python
# PyMC先验设置
beta_exposure = pm.Normal('beta_exposure', mu=0, sigma=1)           # 暴露效应
alpha = pm.Normal('intercept', mu=0, sigma=10)                     # 截距
tau_spatial = pm.Gamma('tau_spatial', alpha=1, beta=0.01)          # 空间精度
sigma_temporal = pm.HalfNormal('sigma_temporal', sigma=1)          # 时间标准差
```

## 使用指南

### 环境设置

```bash
# 在WDP根目录下创建PyMC环境
cd /path/to/WDP
chmod +x setup_pymc_env.sh
./setup_pymc_env.sh

# 激活环境
source activate_pymc.sh

# 或手动激活
conda create -n pymc -y -c conda-forge \
  python=3.11 \
  pymc=5.25.1 arviz=0.22.0 pytensor=2.25.1 \
  numpy>=1.24 scipy>=1.10 pandas>=2.0 \
  networkx>=3.0 geopandas>=0.13 \
  matplotlib>=3.6 seaborn>=0.12 \
  tqdm>=4.64 pyyaml>=6.0 statsmodels>=0.14

conda activate pymc
```

### 命令行使用

```bash
# 基本分析
python main.py --disease-code C81-C96 --pesticide-category "2,4-D" --measure-type Weight

# 完整分析
python main.py \
  --disease-code C81-C96 \
  --pesticide-category compound:2 \
  --measure-type Weight,Density \
  --estimate-types min,avg,max \
  --lag-years 5,10 \
  --model-types M0,M1,M2,M3 \
  --verbose

# 批量分析
python main.py --pesticide-category ALL --model-types M0,M3
```

### 在集群上提交任务（SLURM）

```bash
# 单个化合物/疾病提交（参数顺序：疾病 化合物 模型 滞后 [测量] [draws] [tune] [target_accept] [chains] [cores]）
sbatch Code/PYMC/submit_pymc_single.sbatch C81-C96 24D M3 5 Weight 4000 2000 0.9 4 8

# 或在本地/登录节点直接运行：
bash Code/PYMC/run_pymc_single.sh "C81-C96" "24D" "M3" "5" "Weight" 1000 1000 0.9 4 4
```

### 参数说明

| 参数 | 说明 | 默认值 | 示例 |
|------|------|--------|------|
| `--disease-code` | 疾病ICD编码 | C81-C96 | C81-C96, C50, C34 |
| `--pesticide-category` | 农药类别/化合物 | TEST | ALL, TOP5, "Atrazine", compound:2 |
| `--measure-type` | 暴露测量类型 | Weight | Weight, Density |
| `--estimate-types` | 暴露估计方法 | avg | min, avg, max |
| `--lag-years` | 滞后年数 | 5 | 5, 10 |
| `--model-types` | 模型类型 | M0,M1,M2,M3 | M0, M1, M2, M3 |
| `--output-dir` | 输出目录 | Result/PyMC | 自定义路径 |
| `--verbose` | 详细输出 | False | 启用调试信息 |

## 输出格式

### CSV结果文件

与INLA版本兼容的CSV格式：

```csv
Timestamp,Disease,Exposure,Category,Measure,Estimate,Lag,Model,Dose_Response_Type,RR_Per_SD,RR_Per_SD_Lower,RR_Per_SD_Upper,RR_P90_vs_P10,RR_P90_vs_P10_Lower,RR_P90_vs_P10_Upper,P_Value,DIC,WAIC,N_Counties,N_Records,Status_Message
2025-09-26 14:30:25,C81-C96,2,4-D,Herbicides,Weight,avg,5,M3,Linear,1.0234,0.9876,1.0612,1.1234,1.0123,1.2456,0.034*,2345.67,2342.12,3108,98765,SUCCESS
```

### 列说明

- **Timestamp**: 分析时间戳
- **Disease**: 疾病编码
- **Exposure**: 农药化合物名称
- **Category**: 农药类别
- **Measure**: 暴露测量方法
- **Estimate**: 暴露估计类型
- **Lag**: 滞后年数
- **Model**: 模型类型(M0-M3)
- **RR_Per_SD**: 每标准差增加的相对风险
- **RR_P90_vs_P10**: 90th vs 10th百分位数相对风险
- **P_Value**: 显著性检验p值(含显著性标记)
- **DIC/WAIC**: 模型拟合指标
- **N_Counties/N_Records**: 样本量信息

## 性能优化

### 计算效率
- **采样策略**: 使用NUTS采样器，自动调整步长
- **并行化**: 支持多链并行采样
- **内存管理**: 流式数据处理，减少内存占用
- **缓存机制**: 重复分析时复用数据预处理结果

### 推荐设置
```python
# 快速测试
with model:
    trace = pm.sample(draws=1000, tune=500, chains=2, cores=2)

# 生产分析
with model:
    trace = pm.sample(draws=4000, tune=2000, chains=4, cores=4)
```

## 与INLA版本的差异

### 优势
1. **环境简化**: 仅需Python环境，无需R和复杂依赖
2. **调试友好**: Python生态系统的调试工具更丰富
3. **内存效率**: 避免INLA的临时文件系统
4. **扩展性**: 易于添加新的模型组件和功能

### 限制
1. **计算速度**: PyMC采样可能比INLA近似推断稍慢
2. **收敛诊断**: 需要更仔细的收敛监控
3. **先验选择**: 需要手动调整超参数先验

### 迁移指南
- INLA脚本参数完全兼容
- 输出格式保持一致
- 可并行运行进行结果验证

## 故障排除

### 常见问题

1. **收敛警告**
```python
# 增加采样数和调优步数
trace = pm.sample(draws=6000, tune=3000)
```

2. **内存不足**
```python
# 减少批量大小或使用单一测量类型
--measure-type Weight  # 而不是 Weight,Density
```

3. **空间矩阵问题**
```python
# 检查邻接矩阵的连通性
assert nx.is_connected(adjacency_graph)
```

## 开发和扩展

### 添加新模型
```python
# 在bym_model.py中添加
def build_custom_model(data, model_config):
    with pm.Model() as model:
        # 自定义模型组件
        pass
```

### 自定义输出
```python
# 在results.py中扩展
def extract_custom_results(trace, model):
    # 自定义结果提取
    pass
```

### 测试框架
```bash
# 运行测试套件
python test_pymc.py

# 验证与INLA结果的一致性
python validate_results.py --compare-with-inla
```

## 参考文献

1. Salvatier J, Wiecki TV, Fonnesbeck C (2016). Probabilistic programming in Python using PyMC3. *PeerJ Computer Science* 2:e55.

2. Betancourt M (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint* arXiv:1701.02434.

3. Riebler A, Sørbye SH, Simpson D, Rue H (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. *Statistical Methods in Medical Research* 25(4):1145-1165.

---

**联系方式**: WDP Analysis Team  
**最后更新**: 2025-09-26  
**版本**: 1.0  