#!/usr/bin/env python3
"""
EQI LMM Interval Data Preparation
=================================

将区间数据整理成适合LMM分析的格式，保存到 /Users/ventus/Repository/WDP/Data/df/EQI_LMM_Interval.csv

输出格式：
- FIPS: 县FIPS代码  
- State: 州缩写
- Cancer_Type: 癌症类型 (C00_C97, C15_C26, etc.)
- AAMR_lower: AAMR下限
- AAMR_upper: AAMR上限
- EQI: 总EQI分数
- EQI_quintile: EQI五分位数
- SmokingRate: 吸烟率
- [其他协变量...]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import sys

def load_config():
    """加载配置文件"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / 'config.yaml'
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f), project_root

def load_interval_data(config, project_root, time_period="2016_2020"):
    """加载区间数据"""
    processed_dir = Path(config['data_directories']['processed'])
    interval_file = processed_dir / "CDC" / f"CDC_EQI_Interval_{time_period}.csv"
    
    if not interval_file.exists():
        raise FileNotFoundError(f"区间数据文件不存在: {interval_file}")
    
    print(f"📁 加载区间数据: {interval_file}")
    interval_df = pd.read_csv(interval_file)
    
    # 转换为长格式
    print("🔄 转换为长格式...")
    
    # 提取所有AAMR列对
    aamr_lower_cols = [col for col in interval_df.columns if col.startswith('AAMR_lower_')]
    aamr_upper_cols = [col for col in interval_df.columns if col.startswith('AAMR_upper_')]
    
    # 验证配对
    cancer_types = []
    for lower_col in aamr_lower_cols:
        cancer_type = lower_col.replace('AAMR_lower_', '')
        upper_col = f'AAMR_upper_{cancer_type}'
        if upper_col in aamr_upper_cols:
            cancer_types.append(cancer_type)
    
    print(f"  发现 {len(cancer_types)} 个癌症类型: {cancer_types}")
    
    # 转换为长格式
    long_data_list = []
    
    for cancer_type in cancer_types:
        lower_col = f'AAMR_lower_{cancer_type}'
        upper_col = f'AAMR_upper_{cancer_type}'
        
        cancer_data = interval_df[['COUNTY_FIPS', 'State', lower_col, upper_col]].copy()
        cancer_data = cancer_data.rename(columns={
            lower_col: 'AAMR_lower',
            upper_col: 'AAMR_upper'
        })
        
        # 添加癌症类型列
        cancer_data['Cancer_Type'] = cancer_type
        
        # 重命名FIPS列
        cancer_data = cancer_data.rename(columns={'COUNTY_FIPS': 'FIPS'})
        
        # 过滤有效数据
        cancer_data = cancer_data.dropna(subset=['AAMR_lower', 'AAMR_upper'])
        
        long_data_list.append(cancer_data)
    
    # 合并所有癌症类型
    long_interval_data = pd.concat(long_data_list, ignore_index=True)
    print(f"  长格式数据: {len(long_interval_data):,} 行")
    
    return long_interval_data

def load_covariates(config, project_root, eqi_period="0610"):
    """加载协变量数据"""
    processed_dir = Path(config['data_directories']['processed'])
    
    print("📊 加载协变量数据...")
    
    # 1. EQI数据
    eqi_file = processed_dir / "EQI" / f"EQI{eqi_period}.csv"
    if not eqi_file.exists():
        raise FileNotFoundError(f"EQI数据文件不存在: {eqi_file}")
    
    print(f"  EQI数据: {eqi_file}")
    eqi_data = pd.read_csv(eqi_file)
    
    # 2. 吸烟率数据
    smoking_file = processed_dir / "Smoking" / "County_Smoking.csv"
    if not smoking_file.exists():
        raise FileNotFoundError(f"吸烟数据文件不存在: {smoking_file}")
    
    print(f"  吸烟数据: {smoking_file}")
    smoking_data = pd.read_csv(smoking_file)
    
    # 3. 社会经济数据（如果存在）
    se_files = {
        'BEA': processed_dir / "Socioeconomic" / "BEA_Income.csv",
        'LAUS': processed_dir / "Socioeconomic" / "LAUS_Unemployment.csv",
        'SAIPE': processed_dir / "Socioeconomic" / "SAIPE_Poverty.csv",
        'Education': processed_dir / "Socioeconomic" / "USDA_ERS_Education.csv"
    }
    
    se_data_dict = {}
    for name, file_path in se_files.items():
        if file_path.exists():
            print(f"  {name}数据: {file_path}")
            se_data_dict[name] = pd.read_csv(file_path)
        else:
            print(f"  ⚠️  {name}数据不存在，跳过")
    
    return eqi_data, smoking_data, se_data_dict

def merge_analysis_data(interval_data, eqi_data, smoking_data, se_data_dict):
    """合并分析数据"""
    print("🔄 合并协变量数据...")
    
    # 确保FIPS列格式一致（字符串），并处理NaN值
    # 处理区间数据中的NaN FIPS
    interval_data = interval_data.dropna(subset=['FIPS'])
    interval_data['FIPS'] = interval_data['FIPS'].astype(int).astype(str).str.zfill(5)
    
    eqi_data['COUNTY_FIPS'] = eqi_data['COUNTY_FIPS'].astype(str).str.zfill(5)
    smoking_data['COUNTY_FIPS'] = smoking_data['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    print(f"  清理后区间数据: {len(interval_data):,} 行")
    
    # 合并EQI数据
    merged_data = interval_data.merge(
        eqi_data, 
        left_on='FIPS', 
        right_on='COUNTY_FIPS', 
        how='left'
    ).drop('COUNTY_FIPS', axis=1)
    
    print(f"  合并EQI后: {len(merged_data):,} 行")
    
    # 合并吸烟数据 (重命名列)
    smoking_data_renamed = smoking_data[['COUNTY_FIPS', 'SR_Total']].rename(columns={'SR_Total': 'SmokingRate'})
    merged_data = merged_data.merge(
        smoking_data_renamed, 
        left_on='FIPS', 
        right_on='COUNTY_FIPS', 
        how='left'
    ).drop('COUNTY_FIPS', axis=1)
    
    print(f"  合并吸烟率后: {len(merged_data):,} 行")
    
    # 合并其他社会经济数据
    for name, se_data in se_data_dict.items():
        if 'COUNTY_FIPS' in se_data.columns:
            se_data['COUNTY_FIPS'] = se_data['COUNTY_FIPS'].astype(str).str.zfill(5)
            
            # 选择关键列
            key_cols = ['COUNTY_FIPS']
            
            if name == 'BEA' and 'PerCapitaIncome' in se_data.columns:
                key_cols.append('PerCapitaIncome')
            elif name == 'LAUS' and 'UnemploymentRate' in se_data.columns:
                key_cols.append('UnemploymentRate')
            elif name == 'SAIPE' and 'PovertyRate' in se_data.columns:
                key_cols.append('PovertyRate')
            elif name == 'Education' and 'HighSchoolRate' in se_data.columns:
                key_cols.append('HighSchoolRate')
            
            if len(key_cols) > 1:
                merged_data = merged_data.merge(
                    se_data[key_cols], 
                    left_on='FIPS', 
                    right_on='COUNTY_FIPS', 
                    how='left'
                ).drop('COUNTY_FIPS', axis=1)
                
                print(f"  合并{name}后: {len(merged_data):,} 行")
    
    # 创建EQI五分位数 (先过滤有效EQI值)
    if 'EQI' in merged_data.columns:
        # 先过滤有效的EQI值
        valid_eqi_mask = merged_data['EQI'].notna()
        print(f"  EQI有效值: {valid_eqi_mask.sum():,}/{len(merged_data):,}")
        
        if valid_eqi_mask.sum() > 0:
            # 使用更简单的分位数方法
            quintiles = pd.qcut(
                merged_data.loc[valid_eqi_mask, 'EQI'], 
                q=5, 
                duplicates='drop'
            )
            merged_data.loc[valid_eqi_mask, 'EQI_quintile'] = quintiles.astype(str)
            
            print("  ✅ EQI五分位数已创建")
            print(f"    分位数分布: {quintiles.value_counts().sort_index().to_dict()}")
        else:
            print("  ⚠️  没有有效的EQI值")
    
    # 过滤完整数据
    essential_cols = ['FIPS', 'State', 'Cancer_Type', 'AAMR_lower', 'AAMR_upper', 'EQI', 'SmokingRate']
    before_filter = len(merged_data)
    merged_data = merged_data.dropna(subset=essential_cols)
    after_filter = len(merged_data)
    
    print(f"  过滤缺失值: {before_filter:,} → {after_filter:,} 行")
    
    return merged_data

def save_interval_analysis_data(merged_data, project_root):
    """保存区间分析数据"""
    output_dir = project_root / "Data" / "df"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = output_dir / "EQI_LMM_Interval.csv"
    
    print(f"💾 保存区间分析数据: {output_file}")
    merged_data.to_csv(output_file, index=False)
    
    # 输出数据摘要
    print(f"\n📊 数据摘要:")
    print(f"  总行数: {len(merged_data):,}")
    print(f"  总列数: {len(merged_data.columns)}")
    print(f"  癌症类型: {sorted(merged_data['Cancer_Type'].unique())}")
    print(f"  州数量: {len(merged_data['State'].unique())}")
    print(f"  县数量: {len(merged_data['FIPS'].unique())}")
    
    # EQI五分位数分布
    if 'EQI_quintile' in merged_data.columns:
        print(f"\n  EQI五分位数分布:")
        eqi_dist = merged_data['EQI_quintile'].value_counts().sort_index()
        for quintile, count in eqi_dist.items():
            print(f"    {quintile}: {count:,}")
    
    # 区间宽度统计
    merged_data['Interval_Width'] = merged_data['AAMR_upper'] - merged_data['AAMR_lower']
    print(f"\n  区间宽度统计:")
    print(f"    平均宽度: {merged_data['Interval_Width'].mean():.2f}")
    print(f"    中位数宽度: {merged_data['Interval_Width'].median():.2f}")
    print(f"    最大宽度: {merged_data['Interval_Width'].max():.2f}")
    print(f"    零宽度区间: {(merged_data['Interval_Width'] == 0).sum():,}")
    
    print(f"\n✅ 区间分析数据已保存到: {output_file}")
    return output_file

def main():
    """主函数"""
    print("🔄 EQI区间回归数据准备")
    print("=" * 50)
    
    try:
        # 1. 加载配置
        config, project_root = load_config()
        
        # 2. 加载区间数据
        interval_data = load_interval_data(config, project_root, "2016_2020")
        
        # 3. 加载协变量
        eqi_data, smoking_data, se_data_dict = load_covariates(config, project_root, "0610")
        
        # 4. 合并数据
        merged_data = merge_analysis_data(interval_data, eqi_data, smoking_data, se_data_dict)
        
        # 5. 保存数据
        output_file = save_interval_analysis_data(merged_data, project_root)
        
        print(f"\n🎉 数据准备完成!")
        print(f"📁 输出文件: {output_file}")
        
    except Exception as e:
        print(f"❌ 处理过程中出现错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()