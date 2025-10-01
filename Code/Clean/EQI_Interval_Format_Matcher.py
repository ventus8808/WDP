#!/usr/bin/env python3
"""
EQI LMM Interval Data Formatter
==============================

将区间数据整理成与EQI_LMM_Delete_df.csv相同的格式，但AAMR改为区间形式

目标格式：
- COUNTY_FIPS, State, RUCC, EQI系列, Smoking_Rate
- Analysis_Scenario, Lag_Years, EQI_Period
- AAMR_lower, AAMR_upper (替代AAMR)
- Cancer_Type, Cancer_Description, State_FIPS
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml


def load_config():
    """加载配置文件"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / 'config.yaml'
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f), project_root


def get_cancer_description(cancer_type):
    """根据癌症类型获取描述"""
    cancer_descriptions = {
        'C00_C97': 'All Cancers',
        'C15_C26': 'Digestive System',
        'C18_C21': 'Colorectal',
        'C25': 'Pancreas',
        'C30_C39': 'Respiratory System',
        'C34': 'Lung and Bronchus',
        'C50': 'Breast',
        'C51_C58': 'Female Genital System',
        'C60_C63': 'Male Genital System',
        'C61': 'Prostate',
        'C64_C68': 'Urinary System',
        'C76_C80': 'Ill-defined Sites',
        'C81_C96': 'Lymphoid and Hematopoietic'
    }
    return cancer_descriptions.get(cancer_type, cancer_type)


def get_state_fips(state_abbr):
    """根据州缩写获取州FIPS代码"""
    state_fips_map = {
        'AL': 1, 'AK': 2, 'AZ': 4, 'AR': 5, 'CA': 6, 'CO': 8, 'CT': 9, 'DE': 10, 'DC': 11,
        'FL': 12, 'GA': 13, 'HI': 15, 'ID': 16, 'IL': 17, 'IN': 18, 'IA': 19, 'KS': 20,
        'KY': 21, 'LA': 22, 'ME': 23, 'MD': 24, 'MA': 25, 'MI': 26, 'MN': 27, 'MS': 28,
        'MO': 29, 'MT': 30, 'NE': 31, 'NV': 32, 'NH': 33, 'NJ': 34, 'NM': 35, 'NY': 36,
        'NC': 37, 'ND': 38, 'OH': 39, 'OK': 40, 'OR': 41, 'PA': 42, 'RI': 44, 'SC': 45,
        'SD': 46, 'TN': 47, 'TX': 48, 'UT': 49, 'VT': 50, 'VA': 51, 'WA': 53, 'WV': 54,
        'WI': 55, 'WY': 56
    }
    return state_fips_map.get(state_abbr, 0)


def create_formatted_interval_data():
    """创建格式化的区间数据"""
    
    print("🔄 格式化区间数据以匹配目标结构")
    print("=" * 50)
    
    # 加载当前区间数据
    project_root = Path(__file__).resolve().parents[2]
    current_file = project_root / "Data" / "df" / "EQI_LMM_Interval.csv"
    
    print(f"📁 加载区间数据: {current_file}")
    df = pd.read_csv(current_file)
    
    print(f"  原始数据: {df.shape}")
    
    # 重命名和调整列
    print("🔄 调整列名和格式...")
    formatted_df = df.copy()
    
    # 1. 重命名列
    formatted_df = formatted_df.rename(columns={
        'FIPS': 'COUNTY_FIPS',
        'SmokingRate': 'Smoking_Rate'
    })
    
    # 2. 添加缺失的列
    print("➕ 添加新列...")
    
    # Cancer_Description
    formatted_df['Cancer_Description'] = formatted_df['Cancer_Type'].apply(get_cancer_description)
    
    # State_FIPS
    formatted_df['State_FIPS'] = formatted_df['State'].apply(get_state_fips)
    
    # Analysis_Scenario - 基于2016-2020数据和EQI 0610
    formatted_df['Analysis_Scenario'] = 'EQI0610_AAMR2016_2020'
    
    # Lag_Years - 设为5年
    formatted_df['Lag_Years'] = 5
    
    # EQI_Period - 设为610 (表示2006-2010)
    formatted_df['EQI_Period'] = 610
    
    # 3. 重新排列列的顺序以匹配目标格式
    print("🔄 重新排列列...")
    
    target_columns = [
        'COUNTY_FIPS', 'State', 'RUCC', 'EQI', 'EQI_air', 'EQI_water', 'EQI_land', 
        'EQI_built', 'EQI_Sociodemographic', 'RUCC_EQI', 'RUCC_EQI_air', 'RUCC_EQI_water', 
        'RUCC_EQI_land', 'RUCC_EQI_built', 'RUCC_EQI_Sociodemographic', 'Smoking_Rate',
        'Analysis_Scenario', 'Lag_Years', 'EQI_Period', 
        'AAMR_lower', 'AAMR_upper',  # 区间列替代AAMR
        'Cancer_Type', 'Cancer_Description', 'State_FIPS'
    ]
    
    # 检查所有目标列是否存在
    missing_cols = [col for col in target_columns if col not in formatted_df.columns]
    if missing_cols:
        print(f"⚠️  缺少的列: {missing_cols}")
    
    # 选择存在的列
    available_cols = [col for col in target_columns if col in formatted_df.columns]
    formatted_df = formatted_df[available_cols]
    
    print(f"  格式化后数据: {formatted_df.shape}")
    
    return formatted_df


def save_formatted_data(formatted_df):
    """保存格式化数据"""
    
    project_root = Path(__file__).resolve().parents[2]
    output_dir = project_root / "Data" / "df"
    output_file = output_dir / "EQI_LMM_Interval_Formatted.csv"
    
    print(f"💾 保存格式化数据: {output_file}")
    formatted_df.to_csv(output_file, index=False)
    
    # 输出数据摘要
    print(f"\n📊 数据摘要:")
    print(f"  总行数: {len(formatted_df):,}")
    print(f"  总列数: {len(formatted_df.columns)}")
    print(f"  列名: {list(formatted_df.columns)}")
    
    # 验证与目标格式的匹配度
    target_file = project_root / "Data" / "df" / "EQI_LMM_Delete_df.csv"
    if target_file.exists():
        target_df = pd.read_csv(target_file)
        target_cols = set(target_df.columns)
        current_cols = set(formatted_df.columns)
        
        # 调整AAMR相关列的对比
        target_cols_adj = target_cols.copy()
        if 'AAMR' in target_cols_adj:
            target_cols_adj.remove('AAMR')
            target_cols_adj.add('AAMR_lower')
            target_cols_adj.add('AAMR_upper')
        
        print(f"\n🔍 与目标格式对比:")
        print(f"  匹配的列: {len(target_cols_adj & current_cols)}/{len(target_cols_adj)}")
        
        extra_in_target = target_cols_adj - current_cols
        extra_in_current = current_cols - target_cols_adj
        
        if extra_in_target:
            print(f"  目标格式中额外的列: {extra_in_target}")
        if extra_in_current:
            print(f"  当前数据中额外的列: {extra_in_current}")
        
        if not extra_in_target and not extra_in_current:
            print(f"  ✅ 完全匹配目标格式！")
    
    # 癌症类型统计
    print(f"\n📈 癌症类型统计:")
    cancer_counts = formatted_df['Cancer_Type'].value_counts().sort_index()
    for cancer, count in cancer_counts.items():
        desc = formatted_df[formatted_df['Cancer_Type']==cancer]['Cancer_Description'].iloc[0]
        print(f"  {cancer} ({desc}): {count:,}")
    
    # 区间质量检查
    if 'AAMR_lower' in formatted_df.columns and 'AAMR_upper' in formatted_df.columns:
        interval_width = formatted_df['AAMR_upper'] - formatted_df['AAMR_lower']
        print(f"\n📊 区间质量:")
        print(f"  有效区间: {(interval_width >= 0).sum():,}/{len(formatted_df):,}")
        print(f"  零宽度区间: {(interval_width == 0).sum():,}")
        print(f"  平均区间宽度: {interval_width.mean():.2f}")
        print(f"  最大区间宽度: {interval_width.max():.2f}")
    
    print(f"\n✅ 格式化完成!")
    print(f"📁 输出文件: {output_file}")
    
    return output_file


def main():
    """主函数"""
    print("🔄 EQI区间数据格式化工具")
    print("=" * 50)
    
    try:
        # 1. 创建格式化数据
        formatted_data = create_formatted_interval_data()
        
        # 2. 保存数据
        output_file = save_formatted_data(formatted_data)
        
        print(f"\n🎉 数据格式化完成!")
        print(f"📁 输出文件: {output_file}")
        
        # 3. 显示前几行作为验证
        print(f"\n🔍 数据样本验证:")
        print(formatted_data.head())
        
    except Exception as e:
        print(f"❌ 处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()