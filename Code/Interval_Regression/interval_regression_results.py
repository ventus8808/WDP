#!/usr/bin/env python3
"""
Interval Regression Results Processor
=====================================

处理和分析区间回归结果的模块

功能：
- 读取和处理R分析结果
- 生成结果汇总和比较
- 创建可视化图表
- 与传统LMM结果比较
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class IntervalRegressionResultProcessor:
    """区间回归结果处理器"""
    
    def __init__(self):
        """初始化结果处理器"""
        self.project_root = Path(__file__).resolve().parents[2]
        self.code_dir = Path(__file__).parent
        self.results_dir = self.code_dir / "results"
        self.figures_dir = self.code_dir / "figures"
        
        # 创建图表目录
        self.figures_dir.mkdir(exist_ok=True)
        
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
        
    def load_scenario_results(self, scenario_name: str) -> Optional[pd.DataFrame]:
        """加载场景分析结果"""
        
        # 查找合并结果文件
        combined_file = self.results_dir / f"{scenario_name}_combined_results.csv"
        
        if combined_file.exists():
            return pd.read_csv(combined_file)
        
        # 如果没有合并文件，查找单独的结果文件
        result_files = list(self.results_dir.glob(f"{scenario_name}_*_results.csv"))
        
        if not result_files:
            print(f"⚠️  未找到场景 {scenario_name} 的结果文件")
            return None
        
        # 合并单独的结果文件
        dfs = []
        for file_path in result_files:
            try:
                df = pd.read_csv(file_path)
                dfs.append(df)
            except Exception as e:
                print(f"❌ 读取文件 {file_path} 失败: {e}")
        
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            # 保存合并结果
            combined_df.to_csv(combined_file, index=False)
            return combined_df
        
        return None
    
    def extract_eqi_effects(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """提取EQI效应结果"""
        
        # 筛选EQI相关参数
        eqi_results = results_df[
            results_df['Parameter'].str.contains('EQI_quintile', na=False)
        ].copy()
        
        if len(eqi_results) == 0:
            return pd.DataFrame()
        
        # 解析参数名称
        eqi_results['Quintile'] = eqi_results['Parameter'].str.extract(r'EQI_quintile(\d+)')
        eqi_results['Quintile'] = eqi_results['Quintile'].astype(int)
        
        # 添加癌症中文名称
        eqi_results['Cancer_Name'] = eqi_results['Cancer_Type'].map(self.cancer_names)
        
        # 计算显著性
        eqi_results['Significant'] = (
            (eqi_results['Lower_CI'] > 0) | (eqi_results['Upper_CI'] < 0)
        )
        
        return eqi_results
    
    def create_eqi_effects_plot(self, 
                              eqi_results: pd.DataFrame,
                              scenario_name: str,
                              title_suffix: str = "") -> plt.Figure:
        """创建EQI效应图表"""
        
        if len(eqi_results) == 0:
            print("⚠️  没有EQI效应数据用于绘图")
            return None
        
        # 创建图表
        n_cancers = len(eqi_results['Cancer_Type'].unique())
        fig, axes = plt.subplots(
            nrows=max(1, (n_cancers + 1) // 2), 
            ncols=min(2, n_cancers),
            figsize=(15, 6 * max(1, (n_cancers + 1) // 2))
        )
        
        if n_cancers == 1:
            axes = [axes]
        elif n_cancers > 2:
            axes = axes.flatten()
        
        # 为每种癌症类型绘制子图
        cancer_types = sorted(eqi_results['Cancer_Type'].unique())
        
        for i, cancer_type in enumerate(cancer_types):
            if i >= len(axes):
                break
                
            ax = axes[i] if isinstance(axes, (list, np.ndarray)) else axes
            
            cancer_data = eqi_results[eqi_results['Cancer_Type'] == cancer_type]
            
            # 绘制效应估计和置信区间
            colors = ['red' if sig else 'blue' for sig in cancer_data['Significant']]
            
            ax.errorbar(
                x=cancer_data['Quintile'],
                y=cancer_data['Estimate'],
                yerr=[
                    cancer_data['Estimate'] - cancer_data['Lower_CI'],
                    cancer_data['Upper_CI'] - cancer_data['Estimate']
                ],
                fmt='o',
                capsize=5,
                capthick=2,
                linewidth=2,
                markersize=8,
                color='black',
                ecolor=colors[0] if len(colors) > 0 else 'blue'
            )
            
            # 添加参考线
            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            
            # 设置标题和标签
            cancer_name = self.cancer_names.get(cancer_type, cancer_type)
            ax.set_title(f'{cancer_name}', fontsize=14, fontweight='bold')
            ax.set_xlabel('EQI五分位数', fontsize=12)
            ax.set_ylabel('效应估计 (95% CI)', fontsize=12)
            
            # 设置x轴刻度
            ax.set_xticks(sorted(cancer_data['Quintile'].unique()))
            
            # 网格
            ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        if isinstance(axes, (list, np.ndarray)) and len(cancer_types) < len(axes):
            for i in range(len(cancer_types), len(axes)):
                axes[i].set_visible(False)
        
        # 整体标题
        title = f'区间回归分析 - EQI效应估计'
        if title_suffix:
            title += f' ({title_suffix})'
        
        fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
        
        plt.tight_layout()
        return fig
    
    def create_results_summary_table(self, 
                                   eqi_results: pd.DataFrame,
                                   scenario_name: str) -> pd.DataFrame:
        """创建结果摘要表"""
        
        if len(eqi_results) == 0:
            return pd.DataFrame()
        
        # 创建摘要表
        summary_rows = []
        
        for cancer_type in sorted(eqi_results['Cancer_Type'].unique()):
            cancer_data = eqi_results[eqi_results['Cancer_Type'] == cancer_type]
            
            for _, row in cancer_data.iterrows():
                summary_row = {
                    '场景': scenario_name,
                    '癌症类型': self.cancer_names.get(cancer_type, cancer_type),
                    'EQI五分位数': f"Q{row['Quintile']}",
                    '效应估计': f"{row['Estimate']:.4f}",
                    '95% CI下限': f"{row['Lower_CI']:.4f}",
                    '95% CI上限': f"{row['Upper_CI']:.4f}",
                    '显著性': '是' if row['Significant'] else '否',
                    'R_hat': f"{row['Rhat']:.3f}" if pd.notna(row['Rhat']) else 'N/A'
                }
                summary_rows.append(summary_row)
        
        return pd.DataFrame(summary_rows)
    
    def process_scenario_results(self, scenario_name: str) -> Dict:
        """处理单个场景的结果"""
        
        print(f"\n📊 处理场景结果: {scenario_name}")
        print("-" * 40)
        
        # 加载结果
        results_df = self.load_scenario_results(scenario_name)
        
        if results_df is None or len(results_df) == 0:
            return {
                'scenario_name': scenario_name,
                'success': False,
                'error': '无法加载结果数据'
            }
        
        # 提取EQI效应
        eqi_results = self.extract_eqi_effects(results_df)
        
        # 创建图表
        fig = None
        if len(eqi_results) > 0:
            fig = self.create_eqi_effects_plot(eqi_results, scenario_name)
            
            if fig is not None:
                # 保存图表
                fig_path = self.figures_dir / f"{scenario_name}_eqi_effects.png"
                fig.savefig(fig_path, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"📈 EQI效应图表已保存: {fig_path}")
        
        # 创建摘要表
        summary_table = self.create_results_summary_table(eqi_results, scenario_name)
        
        if len(summary_table) > 0:
            # 保存摘要表
            table_path = self.results_dir / f"{scenario_name}_summary_table.csv"
            summary_table.to_csv(table_path, index=False, encoding='utf-8-sig')
            print(f"📋 摘要表已保存: {table_path}")
        
        # 显示关键统计
        n_significant = len(eqi_results[eqi_results['Significant']]) if len(eqi_results) > 0 else 0
        n_total = len(eqi_results)
        
        print(f"✅ 处理完成:")
        print(f"  - 总参数数: {n_total}")
        print(f"  - 显著效应: {n_significant}")
        print(f"  - 显著比例: {n_significant/n_total*100:.1f}%" if n_total > 0 else "  - 显著比例: N/A")
        
        return {
            'scenario_name': scenario_name,
            'success': True,
            'n_parameters': n_total,
            'n_significant': n_significant,
            'significance_rate': n_significant/n_total if n_total > 0 else 0,
            'summary_table': summary_table.to_dict('records') if len(summary_table) > 0 else [],
            'files_created': [
                str(self.figures_dir / f"{scenario_name}_eqi_effects.png"),
                str(self.results_dir / f"{scenario_name}_summary_table.csv")
            ]
        }
    
    def process_all_scenarios(self) -> Dict:
        """处理所有场景的结果"""
        
        print("🔄 处理所有区间回归结果")
        print("=" * 50)
        
        # 查找所有结果文件
        result_files = list(self.results_dir.glob("*_combined_results.csv"))
        scenario_names = [f.stem.replace("_combined_results", "") for f in result_files]
        
        if not scenario_names:
            print("⚠️  未找到任何结果文件")
            return {'scenarios': [], 'summary': {}}
        
        print(f"📁 发现 {len(scenario_names)} 个场景: {', '.join(scenario_names)}")
        
        # 处理每个场景
        processed_scenarios = []
        
        for scenario_name in scenario_names:
            try:
                result = self.process_scenario_results(scenario_name)
                processed_scenarios.append(result)
            except Exception as e:
                print(f"❌ 处理场景 {scenario_name} 失败: {e}")
                processed_scenarios.append({
                    'scenario_name': scenario_name,
                    'success': False,
                    'error': str(e)
                })
        
        # 生成总体摘要
        successful_scenarios = [s for s in processed_scenarios if s.get('success', False)]
        
        total_summary = {
            'processed_time': datetime.now().isoformat(),
            'total_scenarios': len(processed_scenarios),
            'successful_scenarios': len(successful_scenarios),
            'total_parameters': sum(s.get('n_parameters', 0) for s in successful_scenarios),
            'total_significant': sum(s.get('n_significant', 0) for s in successful_scenarios),
            'overall_significance_rate': 0
        }
        
        if total_summary['total_parameters'] > 0:
            total_summary['overall_significance_rate'] = (
                total_summary['total_significant'] / total_summary['total_parameters']
            )
        
        # 保存总体摘要
        summary_data = {
            'processing_summary': total_summary,
            'scenarios': processed_scenarios
        }
        
        summary_file = self.results_dir / "processing_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary_data, f, indent=2, ensure_ascii=False)
        
        print(f"\n📊 处理完成摘要:")
        print(f"  - 总场景数: {total_summary['total_scenarios']}")
        print(f"  - 成功场景: {total_summary['successful_scenarios']}")
        print(f"  - 总参数数: {total_summary['total_parameters']}")
        print(f"  - 显著效应: {total_summary['total_significant']}")
        print(f"  - 总体显著率: {total_summary['overall_significance_rate']*100:.1f}%")
        
        print(f"\n💾 处理摘要已保存: {summary_file}")
        
        return summary_data
    
    def compare_with_lmm_results(self, 
                               interval_scenario: str,
                               lmm_results_path: Optional[str] = None) -> Optional[pd.DataFrame]:
        """与传统LMM结果比较"""
        
        print(f"\n🔄 比较区间回归与传统LMM结果")
        print("-" * 40)
        
        # 加载区间回归结果
        interval_df = self.load_scenario_results(interval_scenario)
        if interval_df is None:
            print("❌ 无法加载区间回归结果")
            return None
        
        # 如果没有提供LMM结果路径，尝试查找
        if lmm_results_path is None:
            lmm_dir = self.project_root / "Result" / "EQI_LMM"
            lmm_files = list(lmm_dir.glob("*.csv"))
            if lmm_files:
                lmm_results_path = lmm_files[0]  # 使用第一个找到的文件
        
        if lmm_results_path is None:
            print("⚠️  未找到LMM结果文件")
            return None
        
        try:
            lmm_df = pd.read_csv(lmm_results_path)
            print(f"📁 加载LMM结果: {lmm_results_path}")
        except Exception as e:
            print(f"❌ 加载LMM结果失败: {e}")
            return None
        
        # 进行比较分析（这里可以根据实际LMM结果格式调整）
        print("📊 比较分析功能待完善")
        print("  区间回归结果形状:", interval_df.shape)
        print("  LMM结果形状:", lmm_df.shape)
        
        return None


def main():
    """主函数"""
    print("📊 区间回归结果处理程序")
    print("=" * 50)
    
    try:
        # 创建处理器
        processor = IntervalRegressionResultProcessor()
        
        # 处理所有结果
        summary = processor.process_all_scenarios()
        
        print(f"\n🎉 结果处理完成!")
        
        return summary
        
    except Exception as e:
        print(f"❌ 结果处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()