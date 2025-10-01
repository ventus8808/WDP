#!/usr/bin/env python3
"""
Interval Data Processor - 区间数据处理器
=================================

将CDC WONDER数据转换为区间格式，充分利用原始数据中的置信区间信息：
- Reliable (≥20死亡): 使用95% CI区间 [Lower CI, Upper CI]
- Unreliable (10-19死亡): 使用95% CI区间 [Lower CI, Upper CI] 
- Suppressed (1-9死亡): 计算基于人口的AAMR区间 [AAMR_min, AAMR_max]
- Zero (0死亡): 精确区间 [0, 0]

输出格式：AAMR_lower, AAMR_upper 两列
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import re
from datetime import datetime


class IntervalDataProcessor:
    """区间数据处理器"""
    
    def __init__(self, use_2000_standard_pop=True):
        """
        初始化处理器
        
        Parameters:
        -----------
        use_2000_standard_pop : bool
            是否使用2000年美国标准人口进行年龄标准化
        """
        self.use_2000_standard_pop = use_2000_standard_pop
        self.config = self._load_config()
        
        # 2000年美国标准人口（每10万人口）
        self.std_pop_2000 = {
            '0-4': 6.765,
            '5-14': 14.465, 
            '15-24': 13.818,
            '25-34': 13.553,
            '35-44': 14.498,
            '45-54': 12.424,
            '55-64': 8.706,
            '65-74': 6.591,
            '75-84': 4.097,
            '85+': 1.409
        }
        
    def _load_config(self):
        """加载配置文件"""
        project_root = Path(__file__).resolve().parents[2]
        config_path = project_root / 'config.yaml'
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
            
    def _find_col_by_keywords(self, df: pd.DataFrame, include_keys, exclude_keys=None):
        """按关键词近似匹配列名（不区分大小写）"""
        include = [k.lower() for k in include_keys]
        exclude = [k.lower() for k in (exclude_keys or [])]
        for col in df.columns:
            s = str(col).lower()
            if all(k in s for k in include) and not any(k in s for k in exclude):
                return col
        return None
        
    def _extract_state_from_county(self, county_name: str) -> str:
        """从县名中提取州缩写，如 'Autauga County, AL'"""
        if pd.isna(county_name) or not isinstance(county_name, str):
            return ""
        if ", " in county_name:
            return county_name.split(", ")[-1].strip()
        return ""
        
    def _clean_dataframe(self, df):
        """清理数据，移除注释行和无效数据"""
        # 标准化列名
        cols_norm = {c: c.strip() for c in df.columns}
        df = df.rename(columns=cols_norm)
        
        # 删除Notes列（如存在）
        if 'Notes' in df.columns:
            df = df.drop(columns=['Notes'])
        
        # 移除County Code为空或nan的行（这些通常是注释）
        df_clean = df.dropna(subset=['County Code']).copy()
        
        # 确保County Code是数字格式
        df_clean = df_clean[pd.to_numeric(df_clean['County Code'], errors='coerce').notna()]
        
        # 移除包含Missing值的行
        df_clean = df_clean[
            (df_clean['Deaths'] != 'Missing') & 
            (df_clean['Population'] != 'Missing')
        ]
        
        return df_clean
        
    def _calculate_crude_rate_interval(self, deaths_min, deaths_max, population):
        """
        计算粗死亡率区间
        
        Parameters:
        -----------
        deaths_min : int
            最小死亡数
        deaths_max : int  
            最大死亡数
        population : int
            人口数
            
        Returns:
        --------
        tuple: (crude_rate_lower, crude_rate_upper)
        """
        if population <= 0:
            return (np.nan, np.nan)
            
        crude_lower = (deaths_min / population) * 100000
        crude_upper = (deaths_max / population) * 100000
        
        return (crude_lower, crude_upper)
        
    def _estimate_age_adjusted_interval_from_crude(self, crude_lower, crude_upper, 
                                                   age_adjustment_factor=1.0):
        """
        从粗死亡率区间估算年龄调整死亡率区间
        
        这是一个简化的方法。理想情况下应该有年龄别人口数据。
        
        Parameters:
        -----------
        crude_lower : float
            粗死亡率下限
        crude_upper : float
            粗死亡率上限
        age_adjustment_factor : float
            年龄调整因子（简化假设，实际应基于年龄结构）
            
        Returns:
        --------
        tuple: (aamr_lower, aamr_upper)
        """
        # 简化方法：假设年龄调整因子相对稳定
        # 在实际应用中，这应该基于该县/州的年龄结构进行更精确计算
        aamr_lower = crude_lower * age_adjustment_factor
        aamr_upper = crude_upper * age_adjustment_factor
        
        return (aamr_lower, aamr_upper)
        
    def process_single_file_to_intervals(self, file_path, icd_group):
        """
        处理单个文件，将AAMR转换为区间格式
        
        Parameters:
        -----------
        file_path : Path
            CDC原始数据文件路径
        icd_group : str
            ICD代码组，如 "C00-C97"
            
        Returns:
        --------
        pd.DataFrame: 包含COUNTY_FIPS, State, AAMR_lower_{icd}, AAMR_upper_{icd}的数据框
        """
        # 读取数据
        try:
            df = pd.read_csv(file_path, encoding='utf-8')
        except UnicodeDecodeError:
            df = pd.read_csv(file_path, encoding='latin-1')
            
        # 清理数据
        df_clean = self._clean_dataframe(df)
        
        # 标准化FIPS
        df_clean['COUNTY_FIPS'] = df_clean['County Code'].astype(int).astype(str).str.zfill(5)
        
        # 提取State列
        df_clean['State'] = df_clean['County'].apply(self._extract_state_from_county) if 'County' in df_clean.columns else ""
        
        # 提取Deaths和Population
        deaths_data = df_clean['Deaths']
        population_data = pd.to_numeric(df_clean['Population'], errors='coerce')
        
        # 查找AAMR相关列
        rate_col = self._find_col_by_keywords(
            df_clean,
            include_keys=['age', 'adjusted', 'rate'],
            exclude_keys=['lower', 'upper', 'confidence', 'interval', 'ci', 'standard', 'error']
        )
        
        lcl_col = self._find_col_by_keywords(
            df_clean, 
            ['age', 'adjusted', 'rate', 'lower', 'confidence']
        )
        
        ucl_col = self._find_col_by_keywords(
            df_clean, 
            ['age', 'adjusted', 'rate', 'upper', 'confidence']
        )
        
        # 提取AAMR数据
        aamr_point = pd.to_numeric(df_clean[rate_col], errors='coerce') if rate_col else None
        aamr_lower_ci = pd.to_numeric(df_clean[lcl_col], errors='coerce') if lcl_col else None
        aamr_upper_ci = pd.to_numeric(df_clean[ucl_col], errors='coerce') if ucl_col else None
        
        # 分类处理不同类型的数据
        deaths_numeric = pd.to_numeric(deaths_data, errors='coerce')
        deaths_text = deaths_data.astype(str).str.strip().str.lower()
        
        # 创建掩码
        zero_mask = deaths_numeric.eq(0)
        suppressed_mask = deaths_text.str.contains('suppress', na=False) | deaths_numeric.between(1, 9, inclusive='both')
        unreliable_mask = deaths_text.str.contains('unreliable', na=False) | deaths_numeric.between(10, 19, inclusive='both')
        reliable_mask = deaths_numeric.ge(20)
        
        # 初始化区间列
        n_rows = len(df_clean)
        aamr_lower = pd.Series([np.nan] * n_rows, dtype='Float64')
        aamr_upper = pd.Series([np.nan] * n_rows, dtype='Float64')
        
        # 1. Zero数据 (0死亡): [0, 0]
        aamr_lower = aamr_lower.mask(zero_mask, 0.0)
        aamr_upper = aamr_upper.mask(zero_mask, 0.0)
        
        # 2. Reliable数据 (≥20死亡): 使用95% CI
        if aamr_lower_ci is not None and aamr_upper_ci is not None:
            reliable_with_ci = reliable_mask & aamr_lower_ci.notna() & aamr_upper_ci.notna()
            aamr_lower = aamr_lower.mask(reliable_with_ci, aamr_lower_ci)
            aamr_upper = aamr_upper.mask(reliable_with_ci, aamr_upper_ci)
        
        # 如果没有CI但有点估计，则使用点估计作为区间
        if aamr_point is not None:
            reliable_no_ci = reliable_mask & (aamr_lower_ci.isna() | aamr_upper_ci.isna()) & aamr_point.notna()
            aamr_lower = aamr_lower.mask(reliable_no_ci, aamr_point)
            aamr_upper = aamr_upper.mask(reliable_no_ci, aamr_point)
        
        # 3. Unreliable数据 (10-19死亡): 使用95% CI（如果有的话）
        if aamr_lower_ci is not None and aamr_upper_ci is not None:
            unreliable_with_ci = unreliable_mask & aamr_lower_ci.notna() & aamr_upper_ci.notna()
            aamr_lower = aamr_lower.mask(unreliable_with_ci, aamr_lower_ci)
            aamr_upper = aamr_upper.mask(unreliable_with_ci, aamr_upper_ci)
        
        # 如果Unreliable没有CI，使用点估计
        if aamr_point is not None:
            unreliable_no_ci = unreliable_mask & (aamr_lower_ci.isna() | aamr_upper_ci.isna()) & aamr_point.notna()
            aamr_lower = aamr_lower.mask(unreliable_no_ci, aamr_point)
            aamr_upper = aamr_upper.mask(unreliable_no_ci, aamr_point)
        
        # 4. Suppressed数据 (1-9死亡): 基于人口计算粗死亡率区间
        suppressed_indices = suppressed_mask & population_data.notna() & (population_data > 0)
        if suppressed_indices.any():
            for idx in df_clean[suppressed_indices].index:
                pop = population_data.loc[idx]
                # 假设死亡数在1-9之间
                crude_lower, crude_upper = self._calculate_crude_rate_interval(1, 9, pop)
                
                # 简化的年龄调整：使用该州或全国的平均调整因子
                # 这里我们使用1.0作为简化，实际应该基于年龄结构
                age_adj_factor = 1.0  # 简化假设
                
                adj_lower, adj_upper = self._estimate_age_adjusted_interval_from_crude(
                    crude_lower, crude_upper, age_adj_factor
                )
                
                aamr_lower.loc[idx] = adj_lower
                aamr_upper.loc[idx] = adj_upper
        
        # 转换ICD代码格式
        icd_formatted = icd_group.replace('-', '_')
        
        # 构建结果数据框
        result = pd.DataFrame({
            'COUNTY_FIPS': df_clean['COUNTY_FIPS'].astype(str),
            'State': df_clean['State'].astype(str),
            f'AAMR_lower_{icd_formatted}': aamr_lower.round(2),
            f'AAMR_upper_{icd_formatted}': aamr_upper.round(2)
        })
        
        # 移除重复COUNTY_FIPS
        result = result.drop_duplicates(subset=['COUNTY_FIPS'], keep='first')
        
        return result
        
    def process_time_period_to_intervals(self, time_period, qualified_icd_groups):
        """
        处理指定时间周期的所有ICD分组，生成区间数据
        
        Parameters:
        -----------
        time_period : str
            时间周期，如 "2016-2020"
        qualified_icd_groups : list
            符合条件的ICD分组列表
            
        Returns:
        --------
        pd.DataFrame: 合并的区间数据
        """
        project_root = Path(__file__).resolve().parents[2]
        data_dir = project_root / "Data/Original/CDC WONDER EQI"
        
        # 查找该时间周期的所有文件
        pattern = f"*{time_period}*.csv"
        period_files = list(data_dir.glob(pattern))
        
        print(f"\n📅 处理时间周期: {time_period}")
        print(f"🔍 发现 {len(period_files)} 个文件")
        
        merged_df = None
        
        for file_path in period_files:
            # 提取ICD分组
            try:
                # 从文件名提取ICD分组
                filename = file_path.name.replace('.csv', '')
                # 移除时间周期部分
                icd_part = filename.replace(time_period, '').strip()
                # 清理可能的前缀
                icd_group = re.sub(r'^[^C]*', '', icd_part).strip()
                
                if not icd_group or icd_group not in qualified_icd_groups:
                    continue
                    
                print(f"  📊 处理 {icd_group}")
                
                # 处理文件
                interval_data = self.process_single_file_to_intervals(file_path, icd_group)
                
                # 合并数据
                if merged_df is None:
                    merged_df = interval_data.copy()
                    print(f"    🔵 基础数据: {len(merged_df):,} 记录")
                else:
                    # 获取新的AAMR列
                    new_cols = [col for col in interval_data.columns 
                               if col.startswith('AAMR_')]
                    merge_cols = ['COUNTY_FIPS', 'State'] + new_cols
                    
                    before_count = len(merged_df)
                    merged_df = merged_df.merge(
                        interval_data[merge_cols],
                        on=['COUNTY_FIPS', 'State'], 
                        how='outer'
                    )
                    print(f"    ➕ 合并后: {len(merged_df):,} 记录 (增加 {len(merged_df)-before_count:,})")
                    
            except Exception as e:
                print(f"    ⚠️  处理 {file_path.name} 失败: {e}")
                continue
        
        return merged_df
        
    def save_interval_data(self, interval_df, time_period):
        """保存区间数据"""
        project_root = Path(__file__).resolve().parents[2]
        output_dir = Path(self.config['data_directories']['processed']) / 'CDC'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成输出文件名
        time_formatted = time_period.replace('-', '_')
        output_filename = f"CDC_EQI_Interval_{time_formatted}.csv"
        output_path = output_dir / output_filename
        
        # 保存数据
        interval_df.to_csv(output_path, index=False)
        
        print(f"\n💾 保存区间数据: {output_filename}")
        print(f"  🧮 {len(interval_df):,} 记录 × {len(interval_df.columns)} 列")
        
        # 显示区间列统计
        interval_cols = [col for col in interval_df.columns if col.startswith('AAMR_')]
        lower_cols = [col for col in interval_cols if 'lower' in col]
        upper_cols = [col for col in interval_cols if 'upper' in col]
        
        print(f"  📊 区间列: {len(lower_cols)} 个lower列, {len(upper_cols)} 个upper列")
        
        # 检查区间完整性
        if lower_cols:
            sample_lower = lower_cols[0]
            sample_upper = sample_lower.replace('lower', 'upper')
            
            if sample_upper in interval_df.columns:
                valid_intervals = (
                    interval_df[sample_lower].notna() & 
                    interval_df[sample_upper].notna() &
                    (interval_df[sample_lower] <= interval_df[sample_upper])
                ).sum()
                
                print(f"  ✅ {sample_lower.replace('AAMR_lower_', '')} 有效区间: {valid_intervals:,}/{len(interval_df):,}")
        
        return output_path


def main():
    """主函数 - 处理所有时间周期的区间数据"""
    print("🔄 CDC区间数据处理器")
    print("=" * 50)
    
    processor = IntervalDataProcessor()
    
    # 从CDC_EQI_Merge.py获取符合条件的ICD分组
    # 这里我们使用已知的符合条件的分组
    qualified_icd_groups = [
        'C00-C97', 'C15-C26', 'C18-C21', 'C25', 'C30-C39', 'C34', 
        'C50', 'C51-C58', 'C60-C63', 'C61', 'C64-C68', 'C76-C80', 'C81-C96'
    ]
    
    # 处理所有时间周期
    time_periods = ['2006-2010', '2011-2015', '2016-2020']
    
    for time_period in time_periods:
        try:
            print(f"\n" + "="*50)
            print(f"📅 处理时间周期: {time_period}")
            print("="*50)
            
            # 处理该时间周期的区间数据
            interval_data = processor.process_time_period_to_intervals(
                time_period, qualified_icd_groups
            )
            
            if interval_data is not None and len(interval_data) > 0:
                # 保存数据
                output_path = processor.save_interval_data(interval_data, time_period)
                print(f"✅ {time_period} 处理完成")
            else:
                print(f"❌ {time_period} 无有效数据")
                
        except Exception as e:
            print(f"❌ {time_period} 处理失败: {e}")
            
    print(f"\n" + "="*70)
    print("🎉 所有时间周期的区间数据处理完成！")
    print("📊 可用于区间回归分析")
    print("="*70)


if __name__ == "__main__":
    main()