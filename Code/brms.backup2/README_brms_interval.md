# brms 区间 AAMR 分析模块

本模块实现使用 `brms` (Stan) 对 AAMR 区间估计进行多层线性回归（州随机截距），并保持与现有 LMM 输出结构一致的结果表。

## 目录
- `brms_interval_runner.py`: Python 调度脚本，负责分场景/癌症类型/ RUCC 分层拆分数据并调用 R 模型。
- `brms_interval_fit.R`: R 拟合脚本，对单批次数据构建区间删失模型并输出 JSON。

## 依赖环境
建议创建独立 Conda 环境 (示例名: `brms`)：

```bash
conda create -n brms python=3.11 -y
conda activate brms
# Python 侧依赖
pip install pandas pyyaml statsmodels
# R 侧 (请确保已安装 r-base)
# 在 R 中执行：install.packages(c('brms','jsonlite','dplyr','posterior'))
```

或使用 mamba 更快安装 r-base 与编译工具链。

## 数据来源
使用 `config.yaml` 中 `brms_analysis.data_file` 指定的区间数据：
`Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval.csv`
需包含列：`AAMR_lower, AAMR_upper, EQI, EQI_Air, EQI_Water, EQI_Land, EQI_Built, EQI_Social (可选), Smoking_Rate, State, EQI_Period, Time_Period, Cancer_Type, RUCC`

## 模型说明
- 区间删失响应：`cbind(AAMR_lower, AAMR_upper) | cens(cens_type)`
- 固定效应：EQI（或域分解） + 吸烟率
- 随机截距：`(1 | State)`
- Q1 (EQI=1) 为参照组
- 后验近似 p 值：通过双尾概率 `2 * min(P(beta>0), P(beta<0))`

## 调度运行示例

```bash
python Code/brms/brms_interval_runner.py --cancer-types C00_C97,C34 \
  --scenarios EQI0005_AAMR2006_2010,EQI0610_AAMR2016_2020 --apply-fdr
```

参数说明：
- `--cancer-types` 不指定则使用全部在 config.yaml 定义的癌症类型
- `--scenarios` 不指定则运行全部 4 个场景
- `--apply-fdr` 对所有系数 p 值执行 BH FDR 校正并在标记中增加 `†`
- `--dry-run` 仅打印将执行的 R 命令

## 输出
在 `Result/brms/` 生成：
- `brms_<ICD>_<timestamp>.csv` 每种癌症类型一份
- `brms_ALL_<timestamp>.csv` 汇总所有癌症类型

列结构保持与 LMM 一致：
`ICD_Code, EQI_Period, AAMR_Period, Lag, Model, Q1, Q2, Q3, Q4, Q5`

系数单元例：`-3.20(-5.10, -1.30)**`；FDR 校正会附加 `†`。

## 与原 LMM 差异
| 方面 | LMM (statsmodels) | brms 区间模型 |
|------|-------------------|---------------|
| 响应 | 点估计 AAMR | 区间 [lower, upper] |
| 不确定性 | 仅残差 | 区间 + 随机截距 |
| 检验 | Wald p 值 | 后验两尾概率近似 p |
| 拓展 | 需手动修改 | 可直接加非线性/交互/层次结构 |

## 后续扩展建议
1. 增加 RUCC × EQI 交互：在 R 公式中加入 `EQI:RUCC`（需确保 RUCC 作为因子）
2. 增加吸烟率随机斜率：`(1 + Smoking_Rate_std | State)`（样本量允许时）
3. 加入域联合的稀疏先验，降低共线性影响：使用 `horseshoe()` 或 `lasso()` prior
4. 输出后验概率矩阵：例如 `P(Q5 < Q1)` 等，便于方向性解释

## 诊断
首次运行可加入：
```r
pp_check(fit, nsamples = 50)
plot(fit)
```
若出现 `divergent transitions`：提高 `adapt_delta` 到 0.99 或增加迭代。

## 故障排查
- JSON 为空：可能样本量 < 50 或模型未收敛
- R 失败：检查 R 包是否完整安装（brms 依赖 rstan）
- FDR 未标记：可能没有显著 p 值或 p 列为空

## 轻量测试
```bash
python Code/brms/brms_interval_runner.py --cancer-types C00_C97 --scenarios EQI0005_AAMR2006_2010 --dry-run
```
确认命令后去掉 `--dry-run`。

---
维护者: 自动生成 (可根据需要扩展)
