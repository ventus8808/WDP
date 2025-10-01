#!/usr/bin/env python3
"""
PyMC贝叶斯区间回归分析模块
=============================

基于PyMC和Bambi的贝叶斯多层区间回归分析

功能：
- 处理CDC WONDER区间删失数据
- 多层模型（州随机效应）
- 贝叶斯推断和不确定性量化
- EQI环境质量影响分析
"""

import numpy as np
import pandas as pd
import pymc as pm
# import bambi as bmb  # 暂时注释，使用纯PyMC
import arviz as az
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

# 设置绘图
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class PyMCIntervalRegression:
    """PyMC贝叶斯区间回归分析器"""
    
    def __init__(self):
        """初始化分析器"""
        self.project_root = Path(__file__).resolve().parents[2]
        self.code_dir = Path(__file__).parent
        self.results_dir = self.code_dir / "pymc_results"
        self.figures_dir = self.code_dir / "pymc_figures"
        
        # 创建输出目录
        self.results_dir.mkdir(exist_ok=True)
        self.figures_dir.mkdir(exist_ok=True)
        
        # 模型存储
        self.model = None
        self.trace = None
        self.data = None
        
        # 癌症类型映射
        self.cancer_names = {
            'C00_C97': '全部恶性肿瘤',
            'C34': '气管、支气管和肺癌',
            'C50': '乳腺癌',
            'C61': '前列腺癌',
            'C15_C26': '消化器官癌',
            'C18_C21': '结肠、直肠和肛门癌',
            'C25': '胰腺癌'
        }
    
    def load_and_prepare_data(self, 
                            cancer_types: List[str],
                            analysis_scenario: str = 'EQI0610_AAMR2016_2020') -> pd.DataFrame:
        """加载和准备分析数据"""
        print("🔄 加载和准备数据...")
        
        # 加载区间数据
        data_file = self.project_root / 'Data/df/EQI_LMM_Interval.csv'
        df = pd.read_csv(data_file)
        
        print(f"📊 原始数据: {df.shape}")
        
        # 筛选数据
        df = df[
            (df['Cancer_Type'].isin(cancer_types)) &
            (df['Analysis_Scenario'] == analysis_scenario)
        ].copy()
        
        print(f"📊 筛选后数据: {df.shape}")
        
        # 数据预处理
        df = self._preprocess_data(df)
        
        self.data = df
        return df
    
    def _preprocess_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """数据预处理"""
        print("  数据预处理...")
        
        # 确保必要的列存在
        required_cols = ['AAMR_lower', 'AAMR_upper', 'EQI', 'State', 'Cancer_Type', 'Smoking_Rate']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"缺少必要列: {missing_cols}")
        
        # 移除缺失值
        initial_rows = len(df)
        df = df.dropna(subset=required_cols)
        print(f"  移除缺失值: {initial_rows} → {len(df)} 行")
        
        # EQI五分位数（已预计算）
        df['EQI_quintile'] = df['EQI'].astype(int)
        
        # 标准化连续变量
        df['Smoking_Rate_std'] = (df['Smoking_Rate'] - df['Smoking_Rate'].mean()) / df['Smoking_Rate'].std()
        
        # 创建分类变量的索引
        df['State_idx'] = pd.Categorical(df['State']).codes
        df['Cancer_idx'] = pd.Categorical(df['Cancer_Type']).codes
        
        # 区间中点（用于初始值）
        df['AAMR_midpoint'] = (df['AAMR_lower'] + df['AAMR_upper']) / 2
        
        # 区间宽度
        df['interval_width'] = df['AAMR_upper'] - df['AAMR_lower']
        
        print(f"  ✅ 数据预处理完成")
        print(f"  - 州数量: {df['State'].nunique()}")
        print(f"  - 癌症类型: {df['Cancer_Type'].nunique()}")
        print(f"  - EQI五分位数分布: {dict(df['EQI_quintile'].value_counts().sort_index())}")
        
        return df
    
    def build_bayesian_model(self, 
                           include_random_effects: bool = True,
                           include_smoking: bool = True) -> pm.Model:
        """构建贝叶斯区间回归模型"""
        print("🏗️  构建贝叶斯模型...")
        
        if self.data is None:
            raise ValueError("请先加载数据")
        
        df = self.data
        n_states = df['State'].nunique()
        n_cancers = df['Cancer_Type'].nunique()
        
        with pm.Model() as model:
            print("  定义模型结构...")
            
            # === 先验分布 ===
            
            # 截距
            alpha = pm.Normal('alpha', mu=50, sigma=20)
            
            # EQI五分位数效应（相对于Q1）
            beta_eqi = pm.Normal('beta_eqi', mu=0, sigma=5, shape=4)  # Q2, Q3, Q4, Q5
            
            # 吸烟率效应
            if include_smoking:
                beta_smoking = pm.Normal('beta_smoking', mu=0, sigma=5)
            
            # 随机效应
            if include_random_effects:
                # 州随机截距
                sigma_state = pm.HalfNormal('sigma_state', sigma=5)
                alpha_state = pm.Normal('alpha_state', mu=0, sigma=sigma_state, shape=n_states)
                
                # 癌症类型随机截距
                sigma_cancer = pm.HalfNormal('sigma_cancer', sigma=10)
                alpha_cancer = pm.Normal('alpha_cancer', mu=0, sigma=sigma_cancer, shape=n_cancers)
            
            # 观测误差
            sigma_obs = pm.HalfNormal('sigma_obs', sigma=10)
            
            # === 线性预测子 ===
            
            # 创建EQI设计矩阵（相对于Q1）
            eqi_matrix = np.zeros((len(df), 4))
            for i, q in enumerate([2, 3, 4, 5]):
                eqi_matrix[:, i] = (df['EQI_quintile'] == q).astype(int)
            
            # 线性组合
            mu = alpha + pm.math.dot(eqi_matrix, beta_eqi)
            
            if include_smoking:
                mu += beta_smoking * df['Smoking_Rate_std'].values
            
            if include_random_effects:
                mu += alpha_state[df['State_idx'].values]
                mu += alpha_cancer[df['Cancer_idx'].values]
            
            # === 区间删失似然 ===
            
            # 对于每个观测值，定义区间删失的似然
            lower_bounds = df['AAMR_lower'].values
            upper_bounds = df['AAMR_upper'].values
            
            # 使用更简单的区间删失方法
            # 对于区间数据，使用Normal分布的CDF计算区间概率
            
            # 计算区间概率 P(lower < Y < upper) = Φ((upper-μ)/σ) - Φ((lower-μ)/σ)
            prob_lower = pm.math.erfc((lower_bounds - mu) / (sigma_obs * pm.math.sqrt(2))) / 2
            prob_upper = pm.math.erfc((upper_bounds - mu) / (sigma_obs * pm.math.sqrt(2))) / 2
            
            # 区间概率
            prob_interval = prob_lower - prob_upper
            
            # 确保概率为正且不为0
            prob_interval = pm.math.maximum(prob_interval, 1e-10)
            
            # 区间删失似然（对数概率）
            y_obs = pm.Potential('y_obs', pm.math.sum(pm.math.log(prob_interval)))
            
            print("  ✅ 模型构建完成")
            print(f"    - 参数总数: ~{alpha.size + beta_eqi.size + (beta_smoking.size if include_smoking else 0) + (alpha_state.size + alpha_cancer.size if include_random_effects else 0)}")
            print(f"    - 观测数量: {len(df)}")
        
        self.model = model
        return model
    
    def fit_model(self, 
                  draws: int = 2000,
                  tune: int = 1000,
                  chains: int = 4,
                  target_accept: float = 0.9) -> az.InferenceData:
        """拟合贝叶斯模型"""
        print("🎯 拟合贝叶斯模型...")
        
        if self.model is None:
            raise ValueError("请先构建模型")
        
        with self.model:
            print(f"  MCMC采样: {draws} draws × {chains} chains")
            print(f"  预热步数: {tune}")
            
            # MCMC采样
            trace = pm.sample(
                draws=draws,
                tune=tune,
                chains=chains,
                target_accept=target_accept,
                random_seed=42,
                return_inferencedata=True
            )
            
            print("  ✅ 模型拟合完成")
            
            # 模型诊断
            print("🔍 模型诊断:")
            rhat = az.rhat(trace)
            max_rhat = np.max([np.max(rhat[var].values) for var in rhat.data_vars])
            print(f"  - 最大 R̂: {max_rhat:.4f}")
            
            if max_rhat > 1.1:
                print("  ⚠️  收敛警告: 部分参数 R̂ > 1.1")
            else:
                print("  ✅ 收敛良好: 所有参数 R̂ < 1.1")
            
            # 有效样本量
            ess = az.ess(trace)
            min_ess = np.min([np.min(ess[var].values) for var in ess.data_vars])
            print(f"  - 最小有效样本量: {min_ess:.0f}")
            
        self.trace = trace
        return trace
    
    def extract_results(self) -> pd.DataFrame:
        """提取分析结果"""
        print("📊 提取分析结果...")
        
        if self.trace is None:
            raise ValueError("请先拟合模型")
        
        # 提取后验分布统计
        summary = az.summary(self.trace)
        
        # 重新格式化结果
        results = []
        
        # EQI效应
        for i, quintile in enumerate([2, 3, 4, 5]):
            param_name = f"beta_eqi[{i}]"
            if param_name in summary.index:
                row = summary.loc[param_name]
                results.append({
                    'Parameter': f'EQI_quintile_{quintile}',
                    'Mean': row['mean'],
                    'SD': row['sd'],
                    'HDI_3%': row['hdi_3%'],
                    'HDI_97%': row['hdi_97%'],
                    'ESS_bulk': row['ess_bulk'],
                    'ESS_tail': row['ess_tail'],
                    'R_hat': row['r_hat']
                })
        
        # 吸烟率效应
        if 'beta_smoking' in summary.index:
            row = summary.loc['beta_smoking']
            results.append({
                'Parameter': 'Smoking_Rate_std',
                'Mean': row['mean'],
                'SD': row['sd'],
                'HDI_3%': row['hdi_3%'],
                'HDI_97%': row['hdi_97%'],
                'ESS_bulk': row['ess_bulk'],
                'ESS_tail': row['ess_tail'],
                'R_hat': row['r_hat']
            })
        
        # 截距
        if 'alpha' in summary.index:
            row = summary.loc['alpha']
            results.append({
                'Parameter': 'Intercept',
                'Mean': row['mean'],
                'SD': row['sd'],
                'HDI_3%': row['hdi_3%'],
                'HDI_97%': row['hdi_97%'],
                'ESS_bulk': row['ess_bulk'],
                'ESS_tail': row['ess_tail'],
                'R_hat': row['r_hat']
            })
        
        results_df = pd.DataFrame(results)
        
        # 计算显著性（HDI不包含0）
        results_df['Significant'] = (
            (results_df['HDI_3%'] > 0) | (results_df['HDI_97%'] < 0)
        )
        
        print(f"✅ 结果提取完成: {len(results_df)} 个参数")
        
        return results_df
    
    def create_diagnostic_plots(self, save_path: Optional[str] = None) -> plt.Figure:
        """创建模型诊断图"""
        print("📈 创建诊断图...")
        
        if self.trace is None:
            raise ValueError("请先拟合模型")
        
        # 创建诊断图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('贝叶斯模型诊断图', fontsize=16, fontweight='bold')
        
        # 1. 轨迹图
        az.plot_trace(self.trace, var_names=['beta_eqi', 'beta_smoking'], axes=axes[:2])
        
        # 2. R-hat图
        ax = axes[1, 0]
        rhat_data = az.rhat(self.trace)
        rhat_values = []
        for var in rhat_data.data_vars:
            rhat_values.extend(rhat_data[var].values.flatten())
        
        ax.hist(rhat_values, bins=20, alpha=0.7, edgecolor='black')
        ax.axvline(x=1.1, color='red', linestyle='--', label='R̂ = 1.1')
        ax.set_xlabel('R̂')
        ax.set_ylabel('频数')
        ax.set_title('R̂ 分布')
        ax.legend()
        
        # 3. 有效样本量
        ax = axes[1, 1]
        ess_data = az.ess(self.trace)
        ess_values = []
        for var in ess_data.data_vars:
            ess_values.extend(ess_data[var].values.flatten())
        
        ax.hist(ess_values, bins=20, alpha=0.7, edgecolor='black')
        ax.set_xlabel('有效样本量')
        ax.set_ylabel('频数')
        ax.set_title('有效样本量分布')
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"📁 诊断图已保存: {save_path}")
        
        return fig
    
    def create_results_plot(self, results_df: pd.DataFrame, save_path: Optional[str] = None) -> plt.Figure:
        """创建结果图表"""
        print("📊 创建结果图表...")
        
        # 筛选EQI效应
        eqi_results = results_df[results_df['Parameter'].str.contains('EQI_quintile')].copy()
        
        if len(eqi_results) == 0:
            print("⚠️  没有EQI效应数据")
            return None
        
        # 提取五分位数
        eqi_results['Quintile'] = eqi_results['Parameter'].str.extract(r'(\d+)').astype(int)
        eqi_results = eqi_results.sort_values('Quintile')
        
        # 创建图表
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 绘制效应估计和置信区间
        colors = ['red' if sig else 'blue' for sig in eqi_results['Significant']]
        
        ax.errorbar(
            x=eqi_results['Quintile'],
            y=eqi_results['Mean'],
            yerr=[
                eqi_results['Mean'] - eqi_results['HDI_3%'],
                eqi_results['HDI_97%'] - eqi_results['Mean']
            ],
            fmt='o',
            capsize=8,
            capthick=3,
            linewidth=3,
            markersize=10,
            color='black',
            ecolor=colors[0] if len(colors) > 0 else 'blue'
        )
        
        # 为每个点添加颜色
        for i, (_, row) in enumerate(eqi_results.iterrows()):
            color = 'red' if row['Significant'] else 'blue'
            ax.scatter(row['Quintile'], row['Mean'], color=color, s=100, zorder=5)
        
        # 参考线
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        
        # 设置标签和标题
        ax.set_xlabel('EQI五分位数', fontsize=14)
        ax.set_ylabel('效应估计 (相对于Q1)', fontsize=14)
        ax.set_title('EQI环境质量对癌症死亡率的影响\n(贝叶斯区间回归结果)', fontsize=16, fontweight='bold')
        
        # 设置x轴刻度
        ax.set_xticks(sorted(eqi_results['Quintile']))
        ax.set_xticklabels([f'Q{q}' for q in sorted(eqi_results['Quintile'])])
        
        # 图例
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='显著 (HDI不含0)'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='不显著')
        ]
        ax.legend(handles=legend_elements, loc='best')
        
        # 网格
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"📁 结果图已保存: {save_path}")
        
        return fig
    
    def run_analysis(self, 
                   cancer_types: List[str],
                   analysis_scenario: str = 'EQI0610_AAMR2016_2020',
                   draws: int = 2000,
                   tune: int = 1000) -> Dict:
        """运行完整分析"""
        print("🚀 开始贝叶斯区间回归分析")
        print("=" * 60)
        
        results = {
            'success': False,
            'cancer_types': cancer_types,
            'analysis_scenario': analysis_scenario,
            'files_created': []
        }
        
        try:
            # 1. 数据准备
            data = self.load_and_prepare_data(cancer_types, analysis_scenario)
            
            # 2. 构建模型
            model = self.build_bayesian_model()
            
            # 3. 拟合模型
            trace = self.fit_model(draws=draws, tune=tune)
            
            # 4. 提取结果
            results_df = self.extract_results()
            
            # 5. 保存结果
            results_file = self.results_dir / f"bayesian_results_{'+'.join(cancer_types)}.csv"
            results_df.to_csv(results_file, index=False)
            results['files_created'].append(str(results_file))
            
            # 6. 创建图表
            # 诊断图
            diag_file = self.figures_dir / f"diagnostics_{'+'.join(cancer_types)}.png"
            self.create_diagnostic_plots(str(diag_file))
            results['files_created'].append(str(diag_file))
            
            # 结果图
            results_fig_file = self.figures_dir / f"results_{'+'.join(cancer_types)}.png"
            self.create_results_plot(results_df, str(results_fig_file))
            results['files_created'].append(str(results_fig_file))
            
            # 7. 生成报告
            report_file = self.results_dir / f"analysis_report_{'+'.join(cancer_types)}.txt"
            self._generate_report(results_df, data, report_file)
            results['files_created'].append(str(report_file))
            
            results['success'] = True
            results['results_summary'] = self._summarize_results(results_df)
            
            print("🎉 分析完成!")
            print(f"📁 生成文件: {len(results['files_created'])} 个")
            
        except Exception as e:
            print(f"❌ 分析失败: {e}")
            results['error'] = str(e)
            import traceback
            traceback.print_exc()
        
        return results
    
    def _generate_report(self, results_df: pd.DataFrame, data: pd.DataFrame, report_file: Path):
        """生成分析报告"""
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("贝叶斯区间回归分析报告\n")
            f.write("=" * 50 + "\n\n")
            
            # 数据摘要
            f.write("数据摘要:\n")
            f.write(f"- 观测数量: {len(data):,}\n")
            f.write(f"- 癌症类型: {', '.join(data['Cancer_Type'].unique())}\n")
            f.write(f"- 州数量: {data['State'].nunique()}\n")
            f.write(f"- EQI五分位数分布: {dict(data['EQI_quintile'].value_counts().sort_index())}\n\n")
            
            # 模型结果
            f.write("模型结果:\n")
            for _, row in results_df.iterrows():
                f.write(f"- {row['Parameter']}: {row['Mean']:.3f} [{row['HDI_3%']:.3f}, {row['HDI_97%']:.3f}]")
                if row['Significant']:
                    f.write(" *显著*")
                f.write("\n")
            
            f.write(f"\n显著效应数量: {results_df['Significant'].sum()}/{len(results_df)}\n")
        
        print(f"📄 分析报告已保存: {report_file}")
    
    def _summarize_results(self, results_df: pd.DataFrame) -> Dict:
        """结果摘要"""
        
        eqi_results = results_df[results_df['Parameter'].str.contains('EQI_quintile')]
        
        return {
            'total_parameters': len(results_df),
            'significant_parameters': int(results_df['Significant'].sum()),
            'eqi_effects': len(eqi_results),
            'significant_eqi_effects': int(eqi_results['Significant'].sum()),
            'convergence_good': bool((results_df['R_hat'] < 1.1).all())
        }


def main():
    """主函数"""
    print("🎯 PyMC贝叶斯区间回归分析")
    print("=" * 50)
    
    # 创建分析器
    analyzer = PyMCIntervalRegression()
    
    # 运行分析
    results = analyzer.run_analysis(
        cancer_types=['C34', 'C00_C97'],  # 肺癌和全部癌症
        draws=1000,  # 较少的采样数用于快速测试
        tune=500
    )
    
    if results['success']:
        summary = results['results_summary']
        print(f"\n📊 分析结果摘要:")
        print(f"  - 总参数数: {summary['total_parameters']}")
        print(f"  - 显著参数: {summary['significant_parameters']}")
        print(f"  - EQI效应: {summary['significant_eqi_effects']}/{summary['eqi_effects']} 显著")
        print(f"  - 收敛性: {'良好' if summary['convergence_good'] else '需要检查'}")
    else:
        print(f"❌ 分析失败: {results.get('error', '未知错误')}")


if __name__ == "__main__":
    main()