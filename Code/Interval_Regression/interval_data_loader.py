#!/usr/bin/env python3
"""
Interval Regression Data Loader
===============================

加载和预处理区间回归分析数据的模块

功能：
- 加载EQI_LMM_Interval.csv数据
- 数据验证和清洗
- 为不同分析场景准备数据
- 导出R分析所需的格式
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
from typing import List, Optional, Dict, Tuple


class IntervalRegressionDataLoader:
    """区间回归数据加载器"""
    
    def __init__(self, data_file: Optional[Path] = None):
        """
        初始化数据加载器
        
        Parameters:
        -----------
        data_file : Path, optional
            数据文件路径，默认使用项目标准路径
        """
        self.project_root = Path(__file__).resolve().parents[2]
        
        if data_file is None:
            self.data_file = self.project_root / "Data" / "df" / "EQI_LMM_Interval.csv"
        else:
            self.data_file = data_file
            
        self.config = self._load_config()
        self.raw_data = None
        self.analysis_data = None
        
    def _load_config(self) -> Dict:
        """加载配置文件"""
        config_path = self.project_root / 'config.yaml'
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    
    def load_raw_data(self) -> pd.DataFrame:
        """加载原始数据"""
        if not self.data_file.exists():
            raise FileNotFoundError(f"数据文件不存在: {self.data_file}")
        
        print(f"📁 加载区间数据: {self.data_file}")
        self.raw_data = pd.read_csv(self.data_file)
        
        print(f"📊 原始数据: {self.raw_data.shape}")
        
        return self.raw_data
    
    def validate_data(self) -> Dict[str, bool]:
        """验证数据质量"""
        if self.raw_data is None:
            raise ValueError("请先加载数据")
        
        validation_results = {}
        
        # 检查必需列
        required_cols = [
            'COUNTY_FIPS', 'State', 'Cancer_Type', 'AAMR_lower', 'AAMR_upper',
            'EQI', 'Smoking_Rate', 'Analysis_Scenario'
        ]
        
        missing_cols = [col for col in required_cols if col not in self.raw_data.columns]
        validation_results['has_required_columns'] = len(missing_cols) == 0
        
        if missing_cols:
            print(f"⚠️  缺少必需列: {missing_cols}")
        
        # 检查区间有效性
        interval_valid = (self.raw_data['AAMR_lower'] <= self.raw_data['AAMR_upper']).all()
        validation_results['intervals_valid'] = interval_valid
        
        if not interval_valid:
            invalid_count = (self.raw_data['AAMR_lower'] > self.raw_data['AAMR_upper']).sum()
            print(f"⚠️  无效区间: {invalid_count} 个")
        
        # 检查缺失值
        essential_cols = ['COUNTY_FIPS', 'AAMR_lower', 'AAMR_upper', 'EQI']
        missing_counts = self.raw_data[essential_cols].isnull().sum()
        validation_results['no_missing_essential'] = missing_counts.sum() == 0
        
        if missing_counts.sum() > 0:
            print(f"⚠️  核心列缺失值: {missing_counts.to_dict()}")
        
        # 检查数据范围
        aamr_lower_range = (self.raw_data['AAMR_lower'] >= 0).all()
        aamr_upper_range = (self.raw_data['AAMR_upper'] >= 0).all()
        validation_results['valid_ranges'] = aamr_lower_range and aamr_upper_range
        
        print(f"✅ 数据验证: {sum(validation_results.values())}/{len(validation_results)} 项通过")
        
        return validation_results
    
    def prepare_analysis_data(self, 
                            cancer_types: Optional[List[str]] = None,
                            analysis_scenario: Optional[str] = None,
                            rucc_filter: Optional[List[int]] = None) -> pd.DataFrame:
        """
        准备分析数据
        
        Parameters:
        -----------
        cancer_types : List[str], optional
            要分析的癌症类型列表
        analysis_scenario : str, optional
            分析场景筛选
        rucc_filter : List[int], optional
            城乡分类筛选
            
        Returns:
        --------
        pd.DataFrame
            准备好的分析数据
        """
        if self.raw_data is None:
            self.load_raw_data()
        
        print("🔄 准备分析数据...")
        
        # 复制数据
        analysis_data = self.raw_data.copy()
        
        # 筛选癌症类型
        if cancer_types is not None:
            print(f"  筛选癌症类型: {cancer_types}")
            analysis_data = analysis_data[analysis_data['Cancer_Type'].isin(cancer_types)]
        
        # 筛选分析场景
        if analysis_scenario is not None:
            print(f"  筛选分析场景: {analysis_scenario}")
            analysis_data = analysis_data[analysis_data['Analysis_Scenario'] == analysis_scenario]
        
        # 城乡筛选
        if rucc_filter is not None:
            print(f"  筛选城乡类型: {rucc_filter}")
            analysis_data = analysis_data[analysis_data['RUCC'].isin(rucc_filter)]
        
        # 移除缺失值
        before_filter = len(analysis_data)
        analysis_data = analysis_data.dropna(subset=['AAMR_lower', 'AAMR_upper', 'EQI', 'Smoking_Rate'])
        after_filter = len(analysis_data)
        
        print(f"  过滤缺失值: {before_filter:,} → {after_filter:,} 行")
        
        # 创建EQI分类变量
        analysis_data = self._create_eqi_quintiles(analysis_data)
        
        # 标准化连续变量
        analysis_data = self._standardize_variables(analysis_data)
        
        self.analysis_data = analysis_data
        print(f"📊 分析数据: {analysis_data.shape}")
        
        return analysis_data
    
    def _create_eqi_quintiles(self, df: pd.DataFrame) -> pd.DataFrame:
        """使用现有的EQI五分位数"""
        print("  使用现有的EQI五分位数...")
        
        # EQI列已经是五分位数（1-5），直接使用
        df['EQI_quintile'] = df['EQI'].astype(int)
        
        # 各领域EQI也已经是五分位数
        eqi_domains = ['EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic']
        
        for domain in eqi_domains:
            if domain in df.columns:
                quintile_col = f"{domain}_quintile"
                df[quintile_col] = df[domain].astype(int)
        
        # 验证五分位数范围
        for col in ['EQI_quintile'] + [f"{domain}_quintile" for domain in eqi_domains if domain in df.columns]:
            if col in df.columns:
                unique_values = sorted(df[col].unique())
                expected_range = [1, 2, 3, 4, 5]
                if not set(unique_values).issubset(expected_range):
                    print(f"⚠️  {col} 值范围异常: {unique_values}")
        
        return df
    
    def _standardize_variables(self, df: pd.DataFrame) -> pd.DataFrame:
        """标准化连续变量"""
        print("  标准化连续变量...")
        
        # 需要标准化的变量
        continuous_vars = ['Smoking_Rate']
        
        for var in continuous_vars:
            if var in df.columns:
                standardized_col = f"{var}_std"
                df[standardized_col] = (df[var] - df[var].mean()) / df[var].std()
        
        return df
    
    def get_analysis_summary(self) -> Dict:
        """获取分析数据摘要"""
        if self.analysis_data is None:
            raise ValueError("请先准备分析数据")
        
        summary = {}
        
        # 基本统计
        summary['total_observations'] = len(self.analysis_data)
        summary['counties'] = self.analysis_data['COUNTY_FIPS'].nunique()
        summary['states'] = self.analysis_data['State'].nunique()
        summary['cancer_types'] = sorted(self.analysis_data['Cancer_Type'].unique())
        
        # 区间统计
        interval_width = self.analysis_data['AAMR_upper'] - self.analysis_data['AAMR_lower']
        summary['interval_stats'] = {
            'mean_width': interval_width.mean(),
            'median_width': interval_width.median(),
            'max_width': interval_width.max(),
            'zero_width_count': (interval_width == 0).sum()
        }
        
        # EQI分布
        if 'EQI_quintile' in self.analysis_data.columns:
            summary['eqi_distribution'] = self.analysis_data['EQI_quintile'].value_counts().to_dict()
        
        # RUCC分布
        if 'RUCC' in self.analysis_data.columns:
            summary['rucc_distribution'] = self.analysis_data['RUCC'].value_counts().to_dict()
        
        return summary
    
    def export_for_r_analysis(self, output_file: Optional[Path] = None) -> Path:
        """导出R分析格式的数据"""
        if self.analysis_data is None:
            raise ValueError("请先准备分析数据")
        
        if output_file is None:
            output_dir = Path(__file__).parent / "data"
            output_dir.mkdir(exist_ok=True)
            output_file = output_dir / "interval_analysis_data.csv"
        
        print(f"💾 导出R分析数据: {output_file}")
        
        # 选择R分析需要的列
        r_columns = [
            'COUNTY_FIPS', 'State', 'Cancer_Type', 
            'AAMR_lower', 'AAMR_upper',
            'EQI', 'EQI_quintile',
            'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic',
            'Smoking_Rate', 'Smoking_Rate_std',
            'RUCC', 'Analysis_Scenario'
        ]
        
        # 选择存在的列
        available_columns = [col for col in r_columns if col in self.analysis_data.columns]
        r_data = self.analysis_data[available_columns]
        
        # 保存数据
        r_data.to_csv(output_file, index=False)
        
        print(f"  📊 导出数据: {len(r_data):,} 行 × {len(r_data.columns)} 列")
        
        return output_file
    
    def print_data_summary(self):
        """打印数据摘要"""
        if self.analysis_data is None:
            print("❌ 尚未准备分析数据")
            return
        
        summary = self.get_analysis_summary()
        
        print("\n📊 分析数据摘要")
        print("=" * 40)
        print(f"总观测数: {summary['total_observations']:,}")
        print(f"县数量: {summary['counties']:,}")
        print(f"州数量: {summary['states']}")
        print(f"癌症类型: {len(summary['cancer_types'])}")
        
        print(f"\n📈 区间统计:")
        stats = summary['interval_stats']
        print(f"  平均宽度: {stats['mean_width']:.2f}")
        print(f"  中位数宽度: {stats['median_width']:.2f}")
        print(f"  最大宽度: {stats['max_width']:.2f}")
        print(f"  零宽度区间: {stats['zero_width_count']:,}")
        
        if 'eqi_distribution' in summary:
            print(f"\n🌍 EQI五分位数分布:")
            for quintile, count in sorted(summary['eqi_distribution'].items()):
                print(f"  {quintile}: {count:,}")
        
        if 'rucc_distribution' in summary:
            print(f"\n🏘️  城乡分类分布:")
            for rucc, count in sorted(summary['rucc_distribution'].items()):
                print(f"  RUCC {rucc}: {count:,}")


def main():
    """主函数 - 演示数据加载和准备"""
    print("🔄 区间回归数据加载器演示")
    print("=" * 50)
    
    try:
        # 1. 创建数据加载器
        loader = IntervalRegressionDataLoader()
        
        # 2. 加载数据
        raw_data = loader.load_raw_data()
        
        # 3. 验证数据
        validation = loader.validate_data()
        
        if not all(validation.values()):
            print("⚠️  数据验证未完全通过，请检查数据质量")
        
        # 4. 准备分析数据
        analysis_data = loader.prepare_analysis_data(
            cancer_types=['C00_C97', 'C34', 'C50'],  # 示例：总癌症、肺癌、乳腺癌
            analysis_scenario='EQI0610_AAMR2016_2020'
        )
        
        # 5. 打印摘要
        loader.print_data_summary()
        
        # 6. 导出R分析数据
        r_file = loader.export_for_r_analysis()
        
        print(f"\n✅ 数据准备完成!")
        print(f"📁 R分析文件: {r_file}")
        
    except Exception as e:
        print(f"❌ 处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()