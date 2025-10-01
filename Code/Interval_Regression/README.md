# 区间回归分析系统 (Interval Regression Analysis System)

## 概述

本系统实现了基于CDC WONDER癌症死亡率数据的区间回归分析，相比传统的插补方法，直接利用置信区间信息进行统计建模，提供更准确的不确定性量化。

## 核心优势

### 📊 **统计学优势**
- **保留原始不确定性**：直接使用CDC提供的95%置信区间，避免插补假设
- **理论基础更坚实**：区间回归是处理区间数据的标准方法
- **多重比较校正**：支持Benjamini-Hochberg FDR等多重校正方法
- **贝叶斯框架**：使用brms包提供完整的不确定性量化

### 🔧 **技术特点**
- **完整工作流**：从数据处理到结果可视化的端到端解决方案
- **模块化设计**：各功能模块独立，易于维护和扩展
- **自动化分析**：支持批量分析多种癌症类型和场景
- **质量控制**：内置数据验证和结果检查机制

## 系统架构

```
Code/Interval_Regression/
├── interval_data_loader.py      # 数据加载和预处理
├── interval_regression_runner.py   # 分析运行控制器
├── interval_regression_results.py  # 结果处理和可视化
├── demo_and_test.py            # 演示和测试脚本
├── README.md                   # 本文档
├── data/                       # 临时数据文件
├── results/                    # 分析结果
└── figures/                    # 图表输出
```

## 快速开始

### 1. 环境准备

**Python依赖项：**
```bash
pip install pandas numpy matplotlib seaborn pyyaml
```

**R环境和包：**
```r
# 在R中安装
install.packages(c("brms", "dplyr", "readr", "jsonlite", "ggplot2"))
```

### 2. 运行演示

```bash
# 进入项目目录
cd /Users/ventus/Repository/WDP/Code/Interval_Regression

# 运行演示程序（包含依赖项检查）
python demo_and_test.py
```

### 3. 完整分析

```python
from interval_regression_runner import IntervalRegressionAnalyzer

# 创建分析器
analyzer = IntervalRegressionAnalyzer()

# 运行综合分析
results = analyzer.run_comprehensive_analysis()
```

## 详细使用说明

### 数据加载器 (IntervalRegressionDataLoader)

负责加载和预处理区间数据：

```python
from interval_data_loader import IntervalRegressionDataLoader

loader = IntervalRegressionDataLoader()

# 加载原始数据
data = loader.load_raw_data()

# 准备特定癌症类型的分析数据
analysis_data = loader.prepare_analysis_data(
    cancer_types=['C34', 'C50'],  # 肺癌、乳腺癌
    analysis_scenario='EQI0610_AAMR2016_2020'
)

# 导出供R分析
loader.export_for_r_analysis('analysis_data.csv')
```

**主要功能：**
- ✅ 数据质量验证（区间有效性、缺失值等）
- ✅ EQI五分位数创建
- ✅ 变量标准化
- ✅ R兼容格式导出

### 分析运行器 (IntervalRegressionAnalyzer)

协调Python数据处理和R模型拟合：

```python
from interval_regression_runner import IntervalRegressionAnalyzer

analyzer = IntervalRegressionAnalyzer()

# 单个分析场景
result = analyzer.run_analysis_scenario(
    scenario_name='lung_cancer_eqi',
    cancer_types=['C34'],
    analysis_type='total_eqi'
)

# 城乡分层分析
result = analyzer.run_analysis_scenario(
    scenario_name='urban_cancer_analysis',
    cancer_types=['C00_C97'],
    analysis_type='total_eqi',
    rucc_filter=[1, 2, 3]  # 仅城市地区
)
```

**支持的分析类型：**
- `total_eqi`: 总EQI评分分析
- `domain_specific`: EQI领域特异性分析
- `rucc_stratified`: 城乡分层分析

### 结果处理器 (IntervalRegressionResultProcessor)

处理R分析结果并生成可视化：

```python
from interval_regression_results import IntervalRegressionResultProcessor

processor = IntervalRegressionResultProcessor()

# 处理单个场景结果
result = processor.process_scenario_results('lung_cancer_eqi')

# 处理所有场景结果
summary = processor.process_all_scenarios()
```

**输出内容：**
- 📈 EQI效应估计图表
- 📋 结果摘要表格
- 📊 显著性统计
- 🔍 模型诊断信息

## 分析场景示例

### 1. 主要癌症类型分析

```python
scenarios = [
    {
        'name': 'primary_cancers_total_eqi',
        'cancer_types': ['C00_C97', 'C34', 'C50', 'C61'],
        'description': '主要癌症类型 - 总EQI分析'
    }
]
```

### 2. 消化系统癌症分析

```python
scenarios = [
    {
        'name': 'digestive_cancers_total_eqi',
        'cancer_types': ['C15_C26', 'C18_C21', 'C25'],
        'description': '消化系统癌症 - 总EQI分析'
    }
]
```

### 3. 城乡对比分析

```python
# 城市地区分析
urban_result = analyzer.run_analysis_scenario(
    scenario_name='urban_analysis',
    cancer_types=['C00_C97'],
    rucc_filter=[1, 2, 3]
)

# 农村地区分析  
rural_result = analyzer.run_analysis_scenario(
    scenario_name='rural_analysis',
    cancer_types=['C00_C97'],
    rucc_filter=[4, 5, 6, 7, 8, 9]
)
```

## 数据结构说明

### 输入数据格式

系统使用的主要数据文件：`Data/df/EQI_LMM_Interval.csv`

**关键列：**
- `AAMR_lower`, `AAMR_upper`: 年龄调整死亡率置信区间
- `EQI_total_score`: EQI总评分
- `Smoking_Rate`: 吸烟率
- `Cancer_Type`: 癌症类型代码
- `State`: 州代码
- `RUCC`: 城乡连续代码

### 输出结果格式

**分析结果文件：** `results/{scenario_name}_combined_results.csv`

**结果列：**
- `Scenario`: 分析场景名称
- `Cancer_Type`: 癌症类型
- `Parameter`: 模型参数名称
- `Estimate`: 效应估计
- `Lower_CI`, `Upper_CI`: 95%置信区间
- `Rhat`: 收敛诊断统计量

## 模型技术细节

### 区间回归模型

使用brms包实现的贝叶斯多层区间回归：

```r
# 模型公式示例
formula <- bf(
  AAMR_response | cens(cens, AAMR_lower, AAMR_upper) ~ 
    EQI_quintile + Smoking_Rate_std + (1 | State),
  family = gaussian()
)
```

**模型特点：**
- **响应变量**：区间删失的癌症死亡率
- **固定效应**：EQI五分位数、标准化吸烟率
- **随机效应**：州层面的随机截距
- **先验设置**：使用brms默认弱信息先验

### 统计推断

**效应估计：**
- EQI五分位数2-5相对于五分位数1的效应
- 95%可信区间不包含0时认为显著

**模型诊断：**
- Rhat < 1.1 表示收敛良好
- 有效样本量 > 400

## 结果解读指南

### EQI效应图表解读

📈 **效应估计图表特点：**
- **纵轴**：效应估计值（相对于EQI最低五分位数）
- **横轴**：EQI五分位数（2-5）
- **误差线**：95%可信区间
- **颜色**：红色=显著，蓝色=不显著
- **参考线**：零效应线（虚线）

**解读要点：**
- 正值：该EQI水平相对于最低水平增加死亡率
- 负值：该EQI水平相对于最低水平降低死亡率
- 置信区间不跨越0：统计学显著

### 结果表格解读

📋 **摘要表格包含：**
- **效应估计**：点估计值
- **95% CI**：可信区间下限和上限
- **显著性**：基于置信区间是否包含0
- **R_hat**：模型收敛诊断

## 质量控制

### 数据质量检查

✅ **自动验证项目：**
- 区间有效性（下限 ≤ 上限）
- 缺失值检查
- 数值范围合理性
- 样本量充足性

### 模型诊断

✅ **收敛性检查：**
- Rhat统计量 < 1.1
- 有效样本量充足
- 链混合良好

### 结果可靠性

✅ **质量指标：**
- 模型拟合诊断
- 预测区间合理性
- 敏感性分析支持

## 故障排除

### 常见问题及解决方案

**1. R环境问题**
```bash
# 问题：R或brms包未安装
# 解决：
R --version  # 检查R是否安装
# 在R中：install.packages("brms")
```

**2. 数据文件缺失**
```bash
# 问题：找不到 EQI_LMM_Interval.csv
# 解决：确保运行了数据预处理脚本
python Code/Clean/EQI_Interval_Format_Matcher.py
```

**3. 内存不足**
```r
# 问题：R分析内存溢出
# 解决：减少chains数或iter数
params$chains <- 2
params$iter <- 1000
```

**4. 收敛问题**
```r
# 问题：Rhat > 1.1
# 解决：增加迭代次数或调整control参数
control = list(adapt_delta = 0.99, max_treedepth = 15)
```

### 调试技巧

**🔍 启用详细输出：**
```python
# 在演示脚本中查看详细日志
python demo_and_test.py
```

**🔍 检查中间结果：**
```python
# 查看数据加载步骤
loader = IntervalRegressionDataLoader()
data = loader.load_raw_data()
loader.print_data_summary()
```

## 扩展功能

### 添加新的分析类型

1. **在分析运行器中添加新公式：**
```python
# 在 _prepare_r_parameters 方法中
if params['analysis_type'] == "new_analysis":
    # 定义新的分析参数
    pass
```

2. **在R脚本中添加对应的模型公式：**
```r
# 在 interval_regression_analysis.R 中
if (params$analysis_type == "new_analysis") {
  formula <- bf(...)  # 新的模型公式
}
```

### 自定义可视化

```python
# 继承结果处理器类
class CustomResultProcessor(IntervalRegressionResultProcessor):
    def create_custom_plot(self, data):
        # 自定义绘图逻辑
        pass
```

## 版本信息

- **版本**: 1.0.0
- **创建日期**: 2024年
- **Python版本**: 3.7+
- **R版本**: 4.0+
- **主要依赖**: brms, pandas, matplotlib

## 许可证

本系统遵循项目整体许可证。

## 技术支持

如遇到问题或需要功能扩展，请：

1. 🔍 查看本README的故障排除章节
2. 🧪 运行 `demo_and_test.py` 诊断问题
3. 📝 查看生成的错误日志文件
4. 💬 联系项目维护者

---

📌 **重要提示**：本系统假设用户对区间回归和贝叶斯统计有基本了解。建议在使用前阅读相关统计学文献。