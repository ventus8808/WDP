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

# 设置中文字体和绘图样式
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")
sns.set_palette("husl")

logger = logging.getLogger(__name__)


class MICEDiagnostician:
    """
    MICE+PMM多重插补诊断分析器
    """
    
    def __init__(self, imputed_datasets: List[pd.DataFrame], 
                 original_data: pd.DataFrame,
                 convergence_history: Dict[str, List[float]],
                 output_dir: Path):
        """
        初始化诊断分析器
        
        参数:
            imputed_datasets: 插补完成的数据集列表
            original_data: 原始数据
            convergence_history: 收敛历史记录
            output_dir: 诊断结果输出目录
        """
        self.imputed_datasets = imputed_datasets
        self.original_data = original_data
        self.convergence_history = convergence_history
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建诊断子目录
        self.convergence_dir = self.output_dir / "convergence"
        self.distribution_dir = self.output_dir / "distribution" 
        self.scatter_dir = self.output_dir / "scatter"
        
        for dir_path in [self.convergence_dir, self.distribution_dir, self.scatter_dir]:
            dir_path.mkdir(exist_ok=True)
    
    def create_convergence_diagnostics(self) -> Dict[str, str]:
        """
        创建收敛性诊断图
        
        返回:
            收敛性诊断图文件路径字典
        """
        logger.info("=== 创建收敛性诊断图 ===")
        
        if not self.convergence_history:
            logger.warning("无收敛历史记录，跳过收敛性诊断")
            return {}
        
        diagnostic_files = {}
        
        # 为每个变量创建收敛性诊断图
        for var_name, history in self.convergence_history.items():
            if not history or len(history) < 2:
                continue
                
            fig, axes = plt.subplots(2, 2, figsize=(15, 10))
            
            # 1. 链条迹线图
            ax1 = axes[0, 0]
            iterations = range(1, len(history) + 1)
            ax1.plot(iterations, history, 'b-', alpha=0.7, linewidth=1.5)
            ax1.set_title(f'{var_name} - 链条迹线图', fontweight='bold')
            ax1.set_xlabel('迭代次数')
            ax1.set_ylabel('均值')
            ax1.grid(True, alpha=0.3)
            
            # 2. 移动平均趋势
            ax2 = axes[0, 1]
            window_size = min(10, len(history) // 4)
            if window_size >= 2:
                moving_avg = pd.Series(history).rolling(window=window_size).mean()
                ax2.plot(iterations, history, 'lightgray', alpha=0.5, label='原始值')
                ax2.plot(iterations, moving_avg, 'red', linewidth=2, label=f'{window_size}点移动平均')
                ax2.set_title(f'{var_name} - 收敛趋势', fontweight='bold')
                ax2.set_xlabel('迭代次数')
                ax2.set_ylabel('均值')
                ax2.legend()
                ax2.grid(True, alpha=0.3)
            
            # 3. 自相关函数
            ax3 = axes[1, 0]
            try:
                if len(history) > 10:
                    # 计算自相关
                    autocorr_lags = min(20, len(history) // 4)
                    autocorr = pd.Series(history).autocorr(lag=1)
                    
                    # 绘制滞后图
                    lags = range(1, min(autocorr_lags + 1, len(history)))
                    autocorrs = [pd.Series(history).autocorr(lag=lag) for lag in lags]
                    autocorrs = [ac for ac in autocorrs if not pd.isna(ac)]
                    
                    if autocorrs:
                        ax3.bar(lags[:len(autocorrs)], autocorrs, alpha=0.7)
                        ax3.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                        ax3.set_title(f'{var_name} - 自相关函数', fontweight='bold')
                        ax3.set_xlabel('滞后阶数')
                        ax3.set_ylabel('自相关系数')
                        ax3.grid(True, alpha=0.3)
            except Exception as e:
                ax3.text(0.5, 0.5, '自相关计算失败', transform=ax3.transAxes, 
                        ha='center', va='center', fontsize=12)
                logger.warning(f"自相关计算失败 {var_name}: {e}")
            
            # 4. 收敛统计摘要
            ax4 = axes[1, 1]
            ax4.axis('off')
            
            # 计算收敛统计量
            last_10_pct = max(1, len(history) // 10)
            recent_values = history[-last_10_pct:]
            
            stats_text = [
                f'总迭代数: {len(history)}',
                f'最终值: {history[-1]:.6f}',
                f'最后{last_10_pct}次均值: {np.mean(recent_values):.6f}',
                f'最后{last_10_pct}次标准差: {np.std(recent_values):.6f}',
                f'变化范围: [{min(history):.6f}, {max(history):.6f}]',
                f'总体标准差: {np.std(history):.6f}'
            ]
            
            for i, text in enumerate(stats_text):
                ax4.text(0.05, 0.9 - i * 0.15, text, transform=ax4.transAxes,
                        fontsize=11, verticalalignment='top')
            
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
                
                # 1. 原始数据分布（仅观测值）
                ax1 = axes[0, 0]
                if len(observed_values) > 1:
                    ax1.hist(observed_values, bins=30, alpha=0.7, color='blue', 
                            density=True, label='观测值分布')
                    ax1.set_title(f'{var_name} - 原始观测值分布', fontweight='bold')
                    ax1.set_xlabel('数值')
                    ax1.set_ylabel('密度')
                    ax1.legend()
                    ax1.grid(True, alpha=0.3)
                
                # 2-6. 前5个插补数据集的分布对比
                for i, df_imputed in enumerate(self.imputed_datasets[:5]):
                    ax = axes[i // 3, (i % 3) + 1] if i < 3 else axes[1, i - 3]
                    
                    # 获取插补值
                    imputed_values = df_imputed.loc[~observed_mask, var_name]
                    
                    # 绘制观测值和插补值的分布
                    if len(observed_values) > 1:
                        ax.hist(observed_values, bins=20, alpha=0.6, color='blue', 
                               density=True, label='观测值')
                    
                    if len(imputed_values) > 0:
                        ax.hist(imputed_values, bins=20, alpha=0.6, color='red', 
                               density=True, label='插补值')
                    
                    ax.set_title(f'插补数据集 {i+1}', fontweight='bold')
                    ax.set_xlabel('数值')
                    ax.set_ylabel('密度')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                
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
                
                # 原始数据散点图
                ax = axes[0]
                obs_mask1 = self.original_data[var1].notna()
                obs_mask2 = self.original_data[var2].notna()
                obs_mask = obs_mask1 & obs_mask2
                
                if obs_mask.sum() > 0:
                    ax.scatter(self.original_data.loc[obs_mask, var2], 
                              self.original_data.loc[obs_mask, var1],
                              alpha=0.6, color='blue', s=20)
                    ax.set_title('原始观测数据', fontweight='bold')
                    ax.set_xlabel(var2)
                    ax.set_ylabel(var1)
                    ax.grid(True, alpha=0.3)
                
                # 前5个插补数据集的散点图
                for i, df_imp in enumerate(self.imputed_datasets[:5]):
                    ax = axes[i + 1]
                    
                    # 区分观测点和插补点
                    obs_points = df_imp.loc[obs_mask, [var2, var1]]
                    
                    # 缺失值插补点
                    mis_mask1 = self.original_data[var1].isna()
                    mis_mask2 = self.original_data[var2].isna()
                    mis_mask = mis_mask1 | mis_mask2
                    
                    if obs_mask.sum() > 0:
                        ax.scatter(obs_points[var2], obs_points[var1], 
                                  alpha=0.6, color='blue', s=15, label='观测值')
                    
                    if mis_mask.sum() > 0:
                        imp_points = df_imp.loc[mis_mask, [var2, var1]]
                        ax.scatter(imp_points[var2], imp_points[var1],
                                  alpha=0.8, color='red', s=15, label='插补值')
                    
                    ax.set_title(f'插补数据集 {i+1}', fontweight='bold')
                    ax.set_xlabel(var2)
                    ax.set_ylabel(var1)
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                
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
                        imputed_values = df_imp.loc[missing_mask, var_name]
                        if len(imputed_values) > 0:
                            box_data.append(imputed_values)
                            box_labels.append(f'插补{k+1}')
                    
                    # 绘制箱线图
                    if box_data:
                        bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)
                        
                        # 设置颜色
                        colors = ['lightblue'] + ['lightcoral'] * (len(box_data) - 1)
                        for patch, color in zip(bp['boxes'], colors):
                            patch.set_facecolor(color)
                            patch.set_alpha(0.7)
                        
                        ax.set_title(f'{var_name} - 分布对比', fontweight='bold')
                        ax.set_ylabel('数值')
                        ax.grid(True, alpha=0.3)
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
                
                # 原始观测值Q-Q图
                observed_values = self.original_data[var_name].dropna()
                if len(observed_values) > 10:
                    ax = axes[0]
                    stats.probplot(observed_values, dist="norm", plot=ax)
                    ax.set_title(f'{var_name} - 原始观测值 Q-Q图', fontweight='bold')
                    ax.grid(True, alpha=0.3)
                
                # 各插补数据集的Q-Q图
                missing_mask = self.original_data[var_name].isna()
                for i, df_imp in enumerate(self.imputed_datasets[:5]):
                    ax = axes[i + 1]
                    
                    # 完整数据（观测值 + 插补值）
                    complete_values = df_imp[var_name].dropna()
                    
                    if len(complete_values) > 10:
                        stats.probplot(complete_values, dist="norm", plot=ax)
                        ax.set_title(f'插补数据集 {i+1} Q-Q图', fontweight='bold')
                        ax.grid(True, alpha=0.3)
                
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
    
    def run_full_diagnostics(self) -> Dict[str, Any]:
        """
        运行完整的可视化诊断分析
        
        返回:
            完整的诊断结果字典
        """
        logger.info("开始运行完整的MICE+PMM诊断分析...")
        
        diagnostics = {
            'convergence_files': self.create_convergence_diagnostics(),
            'distribution_files': self.create_distribution_diagnostics(),
            'scatter_files': self.create_scatter_diagnostics(),
            'boxplot_files': self.create_boxplot_diagnostics(),
            'qq_plot_files': self.create_qq_plot_diagnostics(),
            'diagnostic_timestamp': datetime.now().strftime("%Y%m%d_%H%M%S"),
            'total_datasets': len(self.imputed_datasets),
            'original_shape': self.original_data.shape
        }
        
        # 保存完整诊断摘要
        summary_file = self.output_dir / "full_diagnostic_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(diagnostics, f, ensure_ascii=False, indent=2, default=str)
        
        logger.info(f"完整诊断分析完成，结果保存至: {self.output_dir}")
        logger.info(f"- 收敛性诊断图: {len(diagnostics['convergence_files'])} 个")
        logger.info(f"- 分布对比图: {len(diagnostics['distribution_files'])} 个")
        logger.info(f"- 散点图诊断: {len(diagnostics['scatter_files'])} 个")
        logger.info(f"- 箱线图诊断: {len(diagnostics['boxplot_files'])} 个")
        logger.info(f"- Q-Q图诊断: {len(diagnostics['qq_plot_files'])} 个")
        
        return diagnostics


def main():
    """
    独立运行诊断分析的主函数
    """
    # 设置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    logger.info("MICE+PMM诊断分析模块独立运行")
    
    # 这里可以添加从文件加载插补结果的代码
    # 例如：加载保存的插补数据集和原始数据进行诊断
    
    print("请配合 MI_MICE_PMM.py 使用本诊断模块")
    print("该模块提供以下诊断功能：")
    print("- 收敛性诊断图")
    print("- 插补前后分布对比图")
    print("- 变量间关系散点图")
    print("- 箱线图诊断")
    print("- Q-Q图正态性检验")


if __name__ == "__main__":
    main()