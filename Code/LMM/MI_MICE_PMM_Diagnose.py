"""
MICE+PMM 多重插补诊断分析模块

该模块提供MICE+PMM多重插补的全面诊断功能，包括：
- 收敛性诊断图
- 分布对比图  
- 散点图诊断
- 可视化分析

注意：此模块需要与 MI_MICE_PMM.py 配合使用
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message=".*downcasting behavior.*")
warnings.filterwarnings("ignore", message=".*Setting an item of incompatible dtype.*")
warnings.filterwarnings("ignore", message=".*labels.*boxplot.*")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging
from datetime import datetime
import json
import sys

# 设置英文字体和优化的绘图样式
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

# 设置优化的颜色主题
sns.set_style("whitegrid", {
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "grid.color": "#E5E5E5",
    "grid.alpha": 0.8
})

# 自定义配色方案
COLORS = {
    'observed': '#2E86AB',      # 深蓝色 - 观测值
    'imputed': '#A23B72',       # 深紫红色 - 插补值  
    'convergence': '#F18F01',   # 橙色 - 收敛线
    'trend': '#C73E1D',         # 红色 - 趋势线
    'grid': '#E5E5E5',          # 浅灰色 - 网格
    'text': '#2D3436'           # 深灰色 - 文字
}

logger = logging.getLogger(__name__)


class SimpleMICEDiagnostician:
    """
    简化的MICE+PMM多重插补诊断分析器
    """
    
    def __init__(self, final_dataset: pd.DataFrame, 
                 original_data: pd.DataFrame,
                 output_dir: Path = None):
        """
        初始化简化诊断分析器
        
        参数:
            final_dataset: 最终插补完成的数据集
            original_data: 原始数据
            output_dir: 诊断结果输出目录
        """
        self.final_dataset = final_dataset
        self.original_data = original_data
        
        # 设置输出目录
        if output_dir is None:
            output_dir = Path("/Users/ventus/Repository/WDP/Result/EQI_LMM_MI_Diagnose")
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def create_convergence_diagnostics(self) -> Dict[str, str]:
        """
        创建收敛性诊断图
        
        返回:
            收敛性诊断图文件路径字典
        """
        logger.info("=== Creating convergence diagnostics ===")
        
        if not self.convergence_history:
            logger.warning("No convergence history found, skipping convergence diagnostics")
            return {}
        
        diagnostic_files = {}
        
        # 为每个变量创建收敛性诊断图
        for var_name, history in self.convergence_history.items():
            if not history or len(history) < 2:
                continue
                
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            
            # 1. Convergence trace plot
            ax1 = axes[0, 0]
            iterations = range(1, len(history) + 1)
            ax1.plot(iterations, history, color=COLORS['convergence'], alpha=0.8, linewidth=2)
            ax1.set_title(f'{var_name} - Trace Plot', fontweight='bold', color=COLORS['text'])
            ax1.set_xlabel('Iteration')
            ax1.set_ylabel('Mean Value')
            ax1.grid(True, alpha=0.4, color=COLORS['grid'])
            
            # 2. Moving average trend
            ax2 = axes[0, 1]
            window_size = min(10, len(history) // 4)
            if window_size >= 2:
                moving_avg = pd.Series(history).rolling(window=window_size).mean()
                ax2.plot(iterations, history, color='lightgray', alpha=0.5, linewidth=1, label='Raw')
                ax2.plot(iterations, moving_avg, color=COLORS['trend'], linewidth=2.5, label=f'MA({window_size})')
                ax2.set_title(f'{var_name} - Convergence Trend', fontweight='bold', color=COLORS['text'])
                ax2.set_xlabel('Iteration')
                ax2.set_ylabel('Mean Value')
                ax2.legend(frameon=True, fancybox=True, shadow=True)
                ax2.grid(True, alpha=0.4, color=COLORS['grid'])
            
            # 3. Autocorrelation function
            ax3 = axes[1, 0]
            try:
                if len(history) > 10:
                    # Calculate autocorrelation
                    autocorr_lags = min(20, len(history) // 4)
                    
                    # Plot lag correlation
                    lags = range(1, min(autocorr_lags + 1, len(history)))
                    autocorrs = [pd.Series(history).autocorr(lag=lag) for lag in lags]
                    autocorrs = [ac for ac in autocorrs if not pd.isna(ac)]
                    
                    if autocorrs:
                        ax3.bar(lags[:len(autocorrs)], autocorrs, alpha=0.7, color=COLORS['imputed'])
                        ax3.axhline(y=0, color='black', linestyle='-', alpha=0.5)
                        ax3.set_title(f'{var_name} - Autocorrelation', fontweight='bold', color=COLORS['text'])
                        ax3.set_xlabel('Lag')
                        ax3.set_ylabel('Correlation')
                        ax3.grid(True, alpha=0.4, color=COLORS['grid'])
            except Exception as e:
                ax3.text(0.5, 0.5, 'Autocorrelation failed', transform=ax3.transAxes, 
                        ha='center', va='center', fontsize=12, color=COLORS['text'])
                logger.warning(f"自相关计算失败 {var_name}: {e}")
            
            # 4. Comparison curves (original vs smoothed)
            ax4 = axes[1, 1]
            
            # Plot original convergence trace and its smoothed version
            iterations = range(1, len(history) + 1)
            ax4.plot(iterations, history, color=COLORS['convergence'], alpha=0.6, linewidth=1.5, label='Original')
            
            # Apply smoothing (exponential moving average)
            alpha = 0.3  # smoothing parameter
            smoothed = [history[0]]
            for i in range(1, len(history)):
                smoothed.append(alpha * history[i] + (1 - alpha) * smoothed[i-1])
            
            ax4.plot(iterations, smoothed, color=COLORS['trend'], linewidth=2.5, label='Smoothed')
            
            # Add convergence band (±1 std)
            window_std = pd.Series(history).rolling(window=min(5, len(history)//3)).std().fillna(0)
            smoothed_series = pd.Series(smoothed)
            upper_band = smoothed_series + window_std
            lower_band = smoothed_series - window_std
            
            ax4.fill_between(iterations, upper_band, lower_band, 
                           color=COLORS['trend'], alpha=0.2, label='±1 SD band')
            
            ax4.set_title(f'{var_name} - Convergence Comparison', fontweight='bold', color=COLORS['text'])
            ax4.set_xlabel('Iteration')
            ax4.set_ylabel('Mean Value')
            ax4.legend(frameon=True, fancybox=True, shadow=True)
            ax4.grid(True, alpha=0.4, color=COLORS['grid'])
            
            plt.tight_layout()
            
            # 保存图像
            filename = f"convergence_{var_name.replace('/', '_')}.png"
            filepath = self.convergence_dir / filename
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            plt.close()
            
            diagnostic_files[var_name] = str(filepath)
            logger.info(f"收敛性诊断图已保存: {filename}")
        
        return diagnostic_files
    
    def create_distribution_diagnostics(self) -> Dict[str, str]:
        """
        创建插补前后分布对比诊断图
        
        返回:
            分布诊断图文件路径字典
        """
        logger.info("=== 创建分布对比诊断图 ===")
        
        if not self.imputed_datasets:
            return {}
        
        diagnostic_files = {}
        
        # 获取有缺失值的变量
        missing_vars = []
        for col in self.original_data.columns:
            if self.original_data[col].isna().any() and pd.api.types.is_numeric_dtype(self.original_data[col]):
                missing_vars.append(col)
        
        logger.info(f"发现 {len(missing_vars)} 个需要诊断的变量: {missing_vars[:5]}...")
        
        for var_name in missing_vars[:10]:  # 限制处理前10个变量避免过多图像
            try:
                fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                
                # 获取观测值和缺失值
                observed_mask = self.original_data[var_name].notna()
                observed_values = self.original_data.loc[observed_mask, var_name]
                
                if len(observed_values) == 0:
                    plt.close()
                    continue
                
                # 1. Original data distribution (observed only)
                ax1 = axes[0, 0]
                if len(observed_values) > 1:
                    ax1.hist(observed_values, bins=30, alpha=0.8, color=COLORS['observed'], 
                            density=True, label='Observed', edgecolor='white', linewidth=0.5)
                    ax1.set_title(f'{var_name} - Original Distribution', fontweight='bold', color=COLORS['text'])
                    ax1.set_xlabel('Value')
                    ax1.set_ylabel('Density')
                    ax1.legend(frameon=True, fancybox=True, shadow=True)
                    ax1.grid(True, alpha=0.4, color=COLORS['grid'])
                
                # 2-6. Distribution comparison for first 5 imputed datasets
                for i, df_imputed in enumerate(self.imputed_datasets[:5]):
                    if i < 5:  # Ensure not exceeding axes range
                        row = (i + 1) // 3
                        col = (i + 1) % 3
                        ax = axes[row, col]
                    
                    # Get imputed values
                    imputed_values = df_imputed.loc[~observed_mask, var_name]
                    
                    # Plot observed and imputed distributions
                    if len(observed_values) > 1:
                        ax.hist(observed_values, bins=20, alpha=0.7, color=COLORS['observed'], 
                               density=True, label='Observed', edgecolor='white', linewidth=0.5)
                    
                    if len(imputed_values) > 0:
                        ax.hist(imputed_values, bins=20, alpha=0.7, color=COLORS['imputed'], 
                               density=True, label='Imputed', edgecolor='white', linewidth=0.5)
                    
                    ax.set_title(f'Imputed Set {i+1}', fontweight='bold', color=COLORS['text'])
                    ax.set_xlabel('Value')
                    ax.set_ylabel('Density')
                    ax.legend(frameon=True, fancybox=True, shadow=True)
                    ax.grid(True, alpha=0.4, color=COLORS['grid'])
                
                plt.tight_layout()
                
                # 保存图像
                filename = f"distribution_{var_name.replace('/', '_')}.png"
                filepath = self.distribution_dir / filename
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                diagnostic_files[var_name] = str(filepath)
                logger.info(f"分布诊断图已保存: {filename}")
                
            except Exception as e:
                logger.error(f"创建 {var_name} 分布诊断图时出错: {e}")
                plt.close()
        
        return diagnostic_files
    
    def create_scatter_diagnostics(self) -> Dict[str, str]:
        """
        创建散点图诊断（检查变量间关系保持）
        
        返回:
            诊断图文件路径字典
        """
        logger.info("=== 创建散点图诊断 ===")
        
        if not self.imputed_datasets:
            return {}
        
        diagnostic_files = {}
        
        # 选择关键变量组合进行散点图分析
        scatter_pairs = [
            ('AAMR_C00_C97', 'EQI'),
            ('AAMR_C00_C97', 'SR'),
            ('AAMR_C34', 'EQI'),  # 肺癌 vs EQI
            ('SR', 'EQI')
        ]
        
        for var1, var2 in scatter_pairs:
            if var1 not in self.original_data.columns or var2 not in self.original_data.columns:
                continue
                
            try:
                fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                axes = axes.flatten()
                
                # Original observed data scatter plot
                ax = axes[0]
                obs_mask1 = self.original_data[var1].notna()
                obs_mask2 = self.original_data[var2].notna()
                obs_mask = obs_mask1 & obs_mask2
                
                if obs_mask.sum() > 0:
                    ax.scatter(self.original_data.loc[obs_mask, var2], 
                              self.original_data.loc[obs_mask, var1],
                              alpha=0.7, color=COLORS['observed'], s=25, edgecolors='white', linewidth=0.5)
                    ax.set_title('Original Observed Data', fontweight='bold', color=COLORS['text'])
                    ax.set_xlabel(var2)
                    ax.set_ylabel(var1)
                    ax.grid(True, alpha=0.4, color=COLORS['grid'])
                
                # Scatter plots for first 5 imputed datasets
                for i, df_imp in enumerate(self.imputed_datasets[:5]):
                    ax = axes[i + 1]
                    
                    # Distinguish observed and imputed points
                    obs_points = df_imp.loc[obs_mask, [var2, var1]]
                    
                    # Missing value imputed points
                    mis_mask1 = self.original_data[var1].isna()
                    mis_mask2 = self.original_data[var2].isna()
                    mis_mask = mis_mask1 | mis_mask2
                    
                    if obs_mask.sum() > 0:
                        ax.scatter(obs_points[var2], obs_points[var1], 
                                  alpha=0.7, color=COLORS['observed'], s=20, 
                                  label='Observed', edgecolors='white', linewidth=0.5)
                    
                    if mis_mask.sum() > 0:
                        imp_points = df_imp.loc[mis_mask, [var2, var1]]
                        ax.scatter(imp_points[var2], imp_points[var1],
                                  alpha=0.8, color=COLORS['imputed'], s=20, 
                                  label='Imputed', edgecolors='white', linewidth=0.5)
                    
                    ax.set_title(f'Imputed Dataset {i+1}', fontweight='bold', color=COLORS['text'])
                    ax.set_xlabel(var2)
                    ax.set_ylabel(var1)
                    ax.legend(frameon=True, fancybox=True, shadow=True)
                    ax.grid(True, alpha=0.4, color=COLORS['grid'])
                
                plt.tight_layout()
                
                # 保存图像
                filename = f"scatter_diagnostics_{var1}_vs_{var2}.png"
                filepath = self.scatter_dir / filename
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                diagnostic_files[f'scatter_{var1}_{var2}'] = str(filepath)
                logger.info(f"散点图诊断已保存: {filename}")
                
            except Exception as e:
                logger.error(f"创建 {var1} vs {var2} 散点图诊断时出错: {e}")
                plt.close()
        
        return diagnostic_files
    
    def create_boxplot_diagnostics(self) -> Dict[str, str]:
        """
        创建箱线图诊断（检查异常值和分布形状）
        
        返回:
            箱线图诊断文件路径字典
        """
        logger.info("=== 创建箱线图诊断 ===")
        
        if not self.imputed_datasets:
            return {}
        
        diagnostic_files = {}
        
        # 获取有缺失值的数值变量
        missing_vars = []
        for col in self.original_data.columns:
            if (self.original_data[col].isna().any() and 
                pd.api.types.is_numeric_dtype(self.original_data[col])):
                missing_vars.append(col)
        
        # 批量处理变量（每张图显示4个变量）
        vars_per_plot = 4
        for i in range(0, len(missing_vars), vars_per_plot):
            batch_vars = missing_vars[i:i+vars_per_plot]
            
            try:
                fig, axes = plt.subplots(2, 2, figsize=(16, 12))
                axes = axes.flatten()
                
                for j, var_name in enumerate(batch_vars):
                    if j >= len(axes):
                        break
                        
                    ax = axes[j]
                    
                    # 准备数据
                    box_data = []
                    box_labels = []
                    
                    # 原始观测值
                    observed_values = self.original_data[var_name].dropna()
                    if len(observed_values) > 0:
                        box_data.append(observed_values)
                        box_labels.append('观测值')
                    
                    # 各插补数据集的值
                    missing_mask = self.original_data[var_name].isna()
                    for k, df_imp in enumerate(self.imputed_datasets[:5]):
                        imputed_values = df_imp.loc[missing_mask, var_name].dropna()
                        if len(imputed_values) > 0:
                            # 确保数据是数值类型
                            try:
                                imputed_values = pd.to_numeric(imputed_values, errors='coerce').dropna()
                                if len(imputed_values) > 0:
                                    box_data.append(imputed_values)
                                    box_labels.append(f'插补{k+1}')
                            except (ValueError, TypeError):
                                continue
                    
                    # Draw boxplot
                    if box_data:
                        bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)
                        
                        # Set colors
                        colors = [COLORS['observed']] + [COLORS['imputed']] * (len(box_data) - 1)
                        for patch, color in zip(bp['boxes'], colors):
                            patch.set_facecolor(color)
                            patch.set_alpha(0.8)
                            patch.set_edgecolor('white')
                            patch.set_linewidth(1)
                        
                        # Style whiskers and median lines
                        for whisker in bp['whiskers']:
                            whisker.set_color(COLORS['text'])
                            whisker.set_linewidth(1.5)
                        for median in bp['medians']:
                            median.set_color('white')
                            median.set_linewidth(2)
                        
                        ax.set_title(f'{var_name} - Distribution Comparison', fontweight='bold', color=COLORS['text'])
                        ax.set_ylabel('Value')
                        ax.grid(True, alpha=0.4, color=COLORS['grid'])
                        ax.tick_params(axis='x', rotation=45)
                
                # 隐藏多余的子图
                for j in range(len(batch_vars), len(axes)):
                    axes[j].set_visible(False)
                
                plt.tight_layout()
                
                # 保存图像
                filename = f"boxplot_diagnostics_batch_{i//vars_per_plot + 1}.png"
                filepath = self.output_dir / filename
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                diagnostic_files[f'boxplot_batch_{i//vars_per_plot + 1}'] = str(filepath)
                logger.info(f"箱线图诊断已保存: {filename}")
                
            except Exception as e:
                logger.error(f"创建箱线图诊断时出错: {e}")
                plt.close()
        
        return diagnostic_files
    
    def create_qq_plot_diagnostics(self) -> Dict[str, str]:
        """
        创建Q-Q图诊断（检查数据分布的正态性）
        
        返回:
            Q-Q图诊断文件路径字典
        """
        logger.info("=== 创建Q-Q图诊断 ===")
        
        if not self.imputed_datasets:
            return {}
        
        diagnostic_files = {}
        
        # 选择关键AAMR变量进行Q-Q图分析
        aamr_vars = [col for col in self.original_data.columns 
                     if col.startswith('AAMR_') and self.original_data[col].isna().any()]
        
        # 限制处理数量避免过多图像
        for var_name in aamr_vars[:8]:
            try:
                fig, axes = plt.subplots(2, 3, figsize=(18, 12))
                axes = axes.flatten()
                
                # Original observed values Q-Q plot
                observed_values = self.original_data[var_name].dropna()
                if len(observed_values) > 10:
                    ax = axes[0]
                    stats.probplot(observed_values, dist="norm", plot=ax)
                    ax.set_title(f'{var_name} - Original Q-Q Plot', fontweight='bold', color=COLORS['text'])
                    ax.grid(True, alpha=0.4, color=COLORS['grid'])
                    # Style the Q-Q plot
                    ax.get_lines()[0].set_markerfacecolor(COLORS['observed'])
                    ax.get_lines()[0].set_markeredgecolor('white')
                    ax.get_lines()[0].set_markersize(6)
                    ax.get_lines()[1].set_color(COLORS['trend'])
                    ax.get_lines()[1].set_linewidth(2)
                
                # Q-Q plots for each imputed dataset
                missing_mask = self.original_data[var_name].isna()
                for i, df_imp in enumerate(self.imputed_datasets[:5]):
                    ax = axes[i + 1]
                    
                    # Complete data (observed + imputed values)
                    complete_values = df_imp[var_name].dropna()
                    
                    if len(complete_values) > 10:
                        stats.probplot(complete_values, dist="norm", plot=ax)
                        ax.set_title(f'Imputed Dataset {i+1} Q-Q Plot', fontweight='bold', color=COLORS['text'])
                        ax.grid(True, alpha=0.4, color=COLORS['grid'])
                        # Style the Q-Q plot
                        ax.get_lines()[0].set_markerfacecolor(COLORS['imputed'])
                        ax.get_lines()[0].set_markeredgecolor('white')
                        ax.get_lines()[0].set_markersize(6)
                        ax.get_lines()[1].set_color(COLORS['trend'])
                        ax.get_lines()[1].set_linewidth(2)
                
                plt.tight_layout()
                
                # 保存图像
                filename = f"qq_plot_{var_name.replace('/', '_')}.png"
                filepath = self.output_dir / filename
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                diagnostic_files[f'qq_{var_name}'] = str(filepath)
                logger.info(f"Q-Q图诊断已保存: {filename}")
                
            except Exception as e:
                logger.error(f"创建 {var_name} Q-Q图诊断时出错: {e}")
                plt.close()
        
        return diagnostic_files
    
    def run_simple_diagnostics(self) -> Dict[str, Any]:
        """
        运行简化的可视化诊断分析
        
        返回:
            简化的诊断结果字典
        """
        logger.info("开始运行简化的MICE+PMM诊断分析...")
        
        # 获取AAMR变量
        aamr_vars = [col for col in self.final_dataset.columns if col.startswith('AAMR_')]
        
        diagnostic_files = {}
        
        # 1. 插补效果对比图
        diagnostic_files['comparison'] = self.create_comparison_plot(aamr_vars)
        
        # 2. 缺失值统计图
        diagnostic_files['missing_stats'] = self.create_missing_stats_plot(aamr_vars)
        
        # 3. 数据分布对比图
        diagnostic_files['distribution'] = self.create_distribution_comparison(aamr_vars[:6])  # 只显示前6个
        
        diagnostics = {
            'diagnostic_files': diagnostic_files,
            'diagnostic_timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
            'final_dataset_shape': self.final_dataset.shape,
            'original_shape': self.original_data.shape,
            'n_aamr_vars': len(aamr_vars),
            'missing_before': self.original_data[aamr_vars].isnull().sum().sum(),
            'missing_after': self.final_dataset[aamr_vars].isnull().sum().sum()
        }
        
        # 保存诊断摘要
        summary_file = self.output_dir / "simple_diagnostic_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(diagnostics, f, ensure_ascii=False, indent=2, default=str)
        
        logger.info(f"简化诊断分析完成，结果保存至: {self.output_dir}")
        logger.info(f"- 诊断图表: {len(diagnostic_files)} 个")
        logger.info(f"- 插补前缺失值: {diagnostics['missing_before']}")
        logger.info(f"- 插补后缺失值: {diagnostics['missing_after']}")
        
        return diagnostics
    
    def create_comparison_plot(self, aamr_vars: List[str]) -> str:
        """创建插补前后对比图"""
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, var in enumerate(aamr_vars[:6]):
            ax = axes[i]
            
            # 原始观测值
            orig_obs = self.original_data[var].dropna()
            # 最终值（包含插补）
            final_vals = self.final_dataset[var]
            
            if len(orig_obs) > 0:
                ax.hist(orig_obs, bins=30, alpha=0.7, label='Original Observed', 
                       color=COLORS['observed'], density=True)
                ax.hist(final_vals, bins=30, alpha=0.7, label='After Imputation', 
                       color=COLORS['imputed'], density=True)
                ax.set_title(f'{var}', fontweight='bold', color=COLORS['text'])
                ax.legend()
                ax.grid(True, alpha=0.3, color=COLORS['grid'])
        
        # 隐藏多余的子图
        for i in range(len(aamr_vars), len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle('AAMR Variables: Before vs After Imputation', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        filepath = self.output_dir / 'aamr_imputation_comparison.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(filepath)
    
    def create_missing_stats_plot(self, aamr_vars: List[str]) -> str:
        """创建缺失值统计图"""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # 原始缺失值统计
        missing_orig = self.original_data[aamr_vars].isnull().sum()
        missing_orig.plot(kind='bar', ax=ax1, color=COLORS['observed'], alpha=0.8)
        ax1.set_title('Missing Values Before Imputation', fontweight='bold', color=COLORS['text'])
        ax1.set_ylabel('Number of Missing Values')
        ax1.tick_params(axis='x', rotation=45)
        ax1.grid(True, alpha=0.3, color=COLORS['grid'])
        
        # 插补后缺失值统计
        missing_final = self.final_dataset[aamr_vars].isnull().sum()
        missing_final.plot(kind='bar', ax=ax2, color=COLORS['imputed'], alpha=0.8)
        ax2.set_title('Missing Values After Imputation', fontweight='bold', color=COLORS['text'])
        ax2.set_ylabel('Number of Missing Values')
        ax2.tick_params(axis='x', rotation=45)
        ax2.grid(True, alpha=0.3, color=COLORS['grid'])
        
        plt.tight_layout()
        
        filepath = self.output_dir / 'missing_values_statistics.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(filepath)
    
    def create_distribution_comparison(self, aamr_vars: List[str]) -> str:
        """创建关键AAMR变量的分布对比"""
        n_vars = len(aamr_vars)
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, var in enumerate(aamr_vars):
            if i >= 6:  # 最多显示6个
                break
                
            ax = axes[i]
            
            # 原始完整观测
            orig_complete = self.original_data[var].dropna()
            # 最终完整数据
            final_complete = self.final_dataset[var].dropna()
            
            if len(orig_complete) > 0 and len(final_complete) > 0:
                # 绘制密度曲线
                orig_complete.hist(bins=30, alpha=0.6, label='Original Complete Data', 
                                 color=COLORS['observed'], density=True, ax=ax)
                final_complete.hist(bins=30, alpha=0.6, label='After Imputation Complete Data', 
                                  color=COLORS['imputed'], density=True, ax=ax)
                
                ax.set_title(f'{var}', fontweight='bold', color=COLORS['text'])
                ax.legend()
                ax.grid(True, alpha=0.3, color=COLORS['grid'])
                
                # 添加统计信息
                orig_mean = orig_complete.mean()
                final_mean = final_complete.mean()
                ax.axvline(orig_mean, color=COLORS['observed'], linestyle='--', alpha=0.8)
                ax.axvline(final_mean, color=COLORS['imputed'], linestyle='--', alpha=0.8)
                
                # 添加均值标注
                ax.text(0.02, 0.95, f'Original Mean: {orig_mean:.1f}', transform=ax.transAxes, 
                       fontsize=9, verticalalignment='top', color=COLORS['observed'])
                ax.text(0.02, 0.88, f'Imputed Mean: {final_mean:.1f}', transform=ax.transAxes, 
                       fontsize=9, verticalalignment='top', color=COLORS['imputed'])
        
        # 隐藏多余的子图
        for i in range(n_vars, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle('Key AAMR Variables Distribution Comparison', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        filepath = self.output_dir / 'aamr_distribution_comparison.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(filepath)


def main():
    """
    独立运行简化诊断分析的主函数
    """
    # 设置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    logger.info("简化MICE+PMM诊断分析模块运行")
    
    try:
        # 加载数据进行独立诊断
        original_path = "/Users/ventus/Repository/WDP/Data/df/EQI_LMM_MI_df.csv"
        final_path = "/Users/ventus/Repository/WDP/Data/df/EQI_LMM_MI_Imputed.csv"
        
        if Path(original_path).exists() and Path(final_path).exists():
            original_data = pd.read_csv(original_path)
            final_data = pd.read_csv(final_path)
            
            # 创建诊断分析器
            diagnostician = SimpleMICEDiagnostician(final_data, original_data)
            
            # 运行诊断
            results = diagnostician.run_simple_diagnostics()
            
            print("✅ 简化诊断分析完成!")
            print(f"📊 生成诊断图表: {len(results['diagnostic_files'])} 个")
            print(f"📁 输出目录: /Users/ventus/Repository/WDP/Result/EQI_LMM_MI_Diagnose")
            
        else:
            print("❌ 未找到必要的数据文件，请先运行插补流程")
            print("需要的文件:")
            print(f"  - 原始数据: {original_path}")
            print(f"  - 插补数据: {final_path}")
    
    except Exception as e:
        logger.error(f"诊断分析出错: {e}")
        print(f"❌ 诊断失败: {e}")
    
    print("\n该简化诊断模块提供以下功能：")
    print("- 插补前后AAMR分布对比")
    print("- 缺失值统计对比")
    print("- 关键变量分布分析")


if __name__ == "__main__":
    main()