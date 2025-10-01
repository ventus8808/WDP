# EQI LMM 多重插补分析Pipeline设计

## 1. 数据路径结构

```
Data/
├── df/
│   ├── EQI_LMM_Delete_df.csv           # 完整案例分析数据 (已有)
│   ├── EQI_LMM_MI_df.csv               # 多重插补原始数据
│   └── MI_Datasets/                     # 插补数据集目录
│       ├── MI_dataset_1.csv
│       ├── MI_dataset_2.csv
│       ├── ...
│       └── MI_dataset_m.csv            # m个完整插补数据集
└── MI_Diagnostics/                      # 插补诊断结果
    ├── missing_patterns.csv
    ├── convergence_plots.png
    └── imputation_summary.csv

Result/
└── EQI_LMM_MI/                         # 多重插补分析结果
    ├── Individual_Results/              # 每个插补数据集的结果
    │   ├── MI_1_Results/
    │   ├── MI_2_Results/
    │   └── ...
    ├── Pooled_Results/                  # 池化结果
    │   ├── MI_C00_C97_Pooled_Results.csv
    │   ├── MI_C15_C26_Pooled_Results.csv
    │   └── ...
    └── MI_Diagnostics/                  # 插补质量诊断
        ├── between_imputation_variance.csv
        └── fraction_missing_info.csv
```

## 2. 核心模块设计

### 2.1 数据准备模块 (LMM_MI_Data.py)

**功能**:
- 生成包含缺失值的完整数据框
- 分析缺失模式 (MAR/MCAR/MNAR)
- 为插补准备辅助变量

**关键方法**:
```python
class LMMMIDataPreprocessor:
    def create_mi_dataset()           # 生成MI基础数据
    def analyze_missing_patterns()    # 分析缺失模式
    def prepare_auxiliary_variables() # 准备辅助变量
    def validate_mi_assumptions()     # 验证MI假设
```

### 2.2 多重插补模块 (LMM_MI_Imputation.py)

**功能**:
- 使用MICE算法进行多重插补
- 生成m个完整数据集 (默认m=20)
- 插补质量诊断

**关键方法**:
```python
class LMMMIImputation:
    def run_mice_imputation()         # MICE插补算法
    def generate_multiple_datasets()  # 生成多个完整数据集
    def check_convergence()           # 收敛性诊断
    def validate_imputation_quality() # 插补质量检查
```

### 2.3 模型分析模块 (LMM_MI_Model.py)

**功能**:
- 对每个插补数据集运行完整LMM分析
- 收集所有插补结果
- 继承现有LMM分析框架

**关键方法**:
```python
class LMMMIAnalyzer:
    def analyze_single_imputation()   # 分析单个插补数据集
    def run_all_imputations()         # 运行所有插补分析
    def collect_results()             # 收集结果
```

### 2.4 结果池化模块 (LMM_MI_Result.py)

**功能**:
- 使用Rubin规则池化参数估计
- 计算池化标准误和置信区间
- 输出最终结果

**关键方法**:
```python
class LMMMIResultPooler:
    def pool_estimates_rubin()        # Rubin规则池化
    def calculate_pooled_se()         # 池化标准误
    def compute_mi_statistics()       # MI统计量
    def format_final_results()        # 格式化最终结果
```

## 3. 分析流程设计

### 3.1 Phase 1: 数据准备
1. 读取原始数据 (包括有缺失值的记录)
2. 生成 `EQI_LMM_MI_df.csv`
3. 分析缺失模式和机制
4. 准备辅助变量

### 3.2 Phase 2: 多重插补
1. 设置插补参数 (m=20, iterations=10)
2. 运行MICE算法
3. 生成20个完整数据集
4. 插补质量诊断

### 3.3 Phase 3: 模型分析
1. 对每个插补数据集运行LMM分析
2. 4个场景 × 13个癌症类型 × 20个插补数据集
3. 收集所有结果

### 3.4 Phase 4: 结果池化
1. 使用Rubin规则池化估计值
2. 计算池化标准误
3. 生成最终置信区间
4. 输出池化结果

## 4. 配置文件更新

添加到 `config.yaml`:
```yaml
# 多重插补分析配置
multiple_imputation:
  # 基础设置
  n_imputations: 20              # 插补次数
  max_iterations: 10             # MICE最大迭代次数
  random_seed: 42               # 随机种子
  
  # 数据路径
  data_paths:
    mi_source: "Data/df/EQI_LMM_MI_df.csv"
    mi_datasets_dir: "Data/df/MI_Datasets"
    diagnostics_dir: "Data/MI_Diagnostics"
  
  # 结果路径
  result_paths:
    base_dir: "Result/EQI_LMM_MI"
    individual_results: "Result/EQI_LMM_MI/Individual_Results"
    pooled_results: "Result/EQI_LMM_MI/Pooled_Results"
    mi_diagnostics: "Result/EQI_LMM_MI/MI_Diagnostics"
  
  # 插补设置
  imputation_settings:
    method: "mice"                # 插补方法
    predictive_mean_matching: true # 使用PMM
    auxiliary_variables:          # 辅助变量
      - "population_density"
      - "urbanization_level" 
      - "geographic_region"
```

## 5. 实施优先级

### Phase 1 (立即实施):
1. 创建 `LMM_MI_Data.py` - 数据准备
2. 生成 `EQI_LMM_MI_df.csv`
3. 分析缺失模式

### Phase 2 (核心功能):
1. 创建 `LMM_MI_Imputation.py` - MICE插补
2. 创建 `LMM_MI_Model.py` - 模型分析
3. 创建 `LMM_MI_Result.py` - 结果池化

### Phase 3 (集成测试):
1. 创建 `LMM_MI_Main.py` - 主控制脚本
2. 端到端测试
3. 性能优化

## 6. 技术考量

### 6.1 计算复杂度
- 20个插补数据集 × 13个癌症类型 × 4个场景 × 30个模型 = 31,200个模型
- 预计运行时间: 4-6小时
- 建议: 并行化处理和增量保存

### 6.2 内存管理
- 分批处理插补数据集
- 及时释放中间结果
- 使用生成器模式

### 6.3 质量控制
- 收敛性诊断
- 插补合理性检查
- 敏感性分析

## 7. 输出格式

最终池化结果表格格式:
```csv
ICD_Code,EQI_Period,AAMR_Period,Lag,Model,Q1_Est,Q1_SE,Q1_CI_Lower,Q1_CI_Upper,Q1_FMI,Q2_Est,Q2_SE,Q2_CI_Lower,Q2_CI_Upper,Q2_FMI,Q3_Est,Q3_SE,Q3_CI_Lower,Q3_CI_Upper,Q3_FMI,Q4_Est,Q4_SE,Q4_CI_Lower,Q4_CI_Upper,Q4_FMI,Q5_Est,Q5_SE,Q5_CI_Lower,Q5_CI_Upper,Q5_FMI
```

其中:
- Est: 池化估计值
- SE: 池化标准误 
- CI_Lower/Upper: 池化置信区间
- FMI: 缺失信息分数 (Fraction of Missing Information)