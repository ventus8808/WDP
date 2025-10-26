# 区间删失混合效应模型介绍

## 概述
该脚本使用 **Stan**（通过 `cmdstanr`）实现了一个贝叶斯区间删失混合效应模型（Interval-Censored Mixed Effects Model），用于分析环境质量指数（EQI）与年龄调整死亡率（AAMR）的关系。该模型考虑了州级随机效应，并处理了区间删失数据（censored data）。

## 数据结构
- **观测数量**：$N$（样本数）
- **州数量**：$S$（州的数量）
- **响应变量**：
  - $y_{\text{lower}}$：AAMR 下界
  - $y_{\text{upper}}$：AAMR 上界
  - $\text{cens}$：删失指示符（0 = 精确观测，2 = 区间删失）
- **设计矩阵**：$X$（$N \times K$ 矩阵，包含协变量）
- **州索引**：$\text{state}$（每个观测对应的州索引）

## 模型公式

### 1. EQI 模型（Overall Model）
该模型使用整体 EQI 分位数作为主要协变量。

#### 均值模型
$$
\mu_i = \beta_0 + \beta_1 \cdot \text{Smoking\_Rate}_i + \sum_{q=2}^{5} \beta_{q} \cdot \mathbb{I}(\text{EQI}_{i} = q) + u_{\text{state}[i]}
$$

其中：
- $\beta_0$：截距
- $\beta_1$：吸烟率的系数
- $\beta_2, \beta_3, \beta_4, \beta_5$：EQI Q2-Q5 相对于 Q1 的系数
- $u_{\text{state}[i]}$：州级随机截距
- $\mathbb{I}(\cdot)$：指示函数

#### 似然函数
- 精确观测：$y_i \sim \mathcal{N}(\mu_i, \sigma^2)$
- 区间删失：$\log p(y_{\text{lower},i} < Y_i < y_{\text{upper},i}) = \log \left( \Phi\left(\frac{y_{\text{upper},i} - \mu_i}{\sigma}\right) - \Phi\left(\frac{y_{\text{lower},i} - \mu_i}{\sigma}\right) \right)$

### 2. EQI_domain 模型（Multi-domain Model）
该模型分别建模每个 EQI 域的分位数。

#### 均值模型
$$
\mu_i = \beta_0 + \beta_1 \cdot \text{Smoking\_Rate}_i + \sum_{d \in \{\text{Air, Water, Land, Built, Social}\}} \sum_{q=2}^{5} \beta_{d,q} \cdot \mathbb{I}(\text{EQI}_{d,i} = q) + u_{\text{state}[i]}
$$

其中：
- $\beta_0$：截距
- $\beta_1$：吸烟率的系数
- $\beta_{d,2}, \beta_{d,3}, \beta_{d,4}, \beta_{d,5}$：域 $d$ 的 EQI Q2-Q5 相对于 Q1 的系数
- $u_{\text{state}[i]}$：州级随机截距

#### 似然函数
同上，与 EQI 模型相同。

## 3. Delta 模型（变化模型）
该模型分析环境质量指数变化（EQI Change）与年龄调整死亡率变化（Delta AAMR）的关系，使用相同的区间删失混合效应框架。

### 均值模型
\[
\mu_i = \beta_0 + \beta_1 \cdot \Delta \text{Smoking\_Rate}_i + \sum_{c \in \{\text{Improved, Worsened}\}} \beta_c \cdot \mathbb{I}(\text{EQI_Change}_{i} = c) + u_{\text{state}[i]}
\]

其中：
- $\beta_0$：截距（基准变化）
- $\beta_1$：吸烟率变化的系数
- $\beta_{\text{Improved}}$：EQI 改善（Improved）相对于稳定（Stable）的系数
- $\beta_{\text{Worsened}}$：EQI 恶化（Worsened）相对于稳定（Stable）的系数
- $u_{\text{state}[i]}$：州级随机截距

### 多域版本
类似 EQI 模型的多域扩展，每个域独立建模变化类别。

### 似然函数
同上，处理区间删失的 Delta AAMR 数据。

## 共同参数
- 先验：$\beta \sim \mathcal{N}(0, 5^2)$, $\sigma \sim \text{Exponential}(1)$, $\sigma_u \sim \text{Exponential}(1)$
- 随机效应：$u_{\text{state}} \sim \mathcal{N}(0, \sigma_u^2)$

带似然函数的完整公式
$$

\log p(\beta, z, \sigma, \sigma_u \mid \text{data}) \propto

\sum_{i:\,\text{cens}_i=0} \log \varphi(y_i;\mu_i,\sigma)

+ \sum_{i:\,\text{cens}_i=2} \log\Big(\Phi\Big(\frac{y_i^{(u)}-\mu_i}{\sigma}\Big) - \Phi\Big(\frac{y_i^{(l)}-\mu_i}{\sigma}\Big)\Big)

+ \sum_j \log p(\beta_j) + \sum_s \log p(z_s) + \log p(\sigma) + \log p(\sigma_u)

$$