#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 吸烟率数据清洗脚本
输出：1) County_Smoking_EQI.csv - EQI时期平均吸烟率
     2) County_Smoking.csv - 2010年各性别吸烟率
"""

import pandas as pd
from pathlib import Path

def create_state_mapping():
    """州名全称到缩写的映射"""
    return {
        'Alabama': 'AL', 'Alaska': 'AK', 'Arizona': 'AZ', 'Arkansas': 'AR', 'California': 'CA',
        'Colorado': 'CO', 'Connecticut': 'CT', 'Delaware': 'DE', 'Florida': 'FL', 'Georgia': 'GA',
        'Hawaii': 'HI', 'Idaho': 'ID', 'Illinois': 'IL', 'Indiana': 'IN', 'Iowa': 'IA',
        'Kansas': 'KS', 'Kentucky': 'KY', 'Louisiana': 'LA', 'Maine': 'ME', 'Maryland': 'MD',
        'Massachusetts': 'MA', 'Michigan': 'MI', 'Minnesota': 'MN', 'Mississippi': 'MS',
        'Missouri': 'MO', 'Montana': 'MT', 'Nebraska': 'NE', 'Nevada': 'NV', 'New Hampshire': 'NH',
        'New Jersey': 'NJ', 'New Mexico': 'NM', 'New York': 'NY', 'North Carolina': 'NC',
        'North Dakota': 'ND', 'Ohio': 'OH', 'Oklahoma': 'OK', 'Oregon': 'OR', 'Pennsylvania': 'PA',
        'Rhode Island': 'RI', 'South Carolina': 'SC', 'South Dakota': 'SD', 'Tennessee': 'TN',
        'Texas': 'TX', 'Utah': 'UT', 'Vermont': 'VT', 'Virginia': 'VA', 'Washington': 'WA',
        'West Virginia': 'WV', 'Wisconsin': 'WI', 'Wyoming': 'WY', 'District Of Columbia': 'DC'
    }

def create_fips_mapping(location_df):
    """创建FIPS映射表"""
    # 解析County列：从"Autauga County, AL"提取县名
    location_df['county_name'] = location_df['County'].str.split(',').str[0].str.strip()
    location_df['state_abbr'] = location_df['State'].str.strip()
    location_df['COUNTY_FIPS'] = location_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    # 创建映射字典：(州缩写, 县名) -> COUNTY_FIPS
    return {(row['state_abbr'], row['county_name']): row['COUNTY_FIPS'] 
            for _, row in location_df.iterrows()}

def process_data(county_data, fips_dict, state_mapping):
    """处理EQI时期和单年数据"""
    print("\n� 处理数据...")
    
    # 添加州缩写和FIPS映射
    county_data['state_abbr'] = county_data['state'].map(state_mapping)
    county_data['COUNTY_FIPS'] = county_data.apply(
        lambda row: fips_dict.get((row['state_abbr'], row['county']), None), axis=1
    )
    
    # 只保留匹配成功的县
    matched_data = county_data[county_data['COUNTY_FIPS'].notna()].copy()
    print(f"  ✅ 匹配成功: {len(matched_data)} 记录")
    
    # 1. EQI时期数据：计算2000-2005和2006-2010年平均值
    eqi_data = []
    periods = {'0005': (2000, 2005), '0610': (2006, 2010)}
    
    for period, (start, end) in periods.items():
        period_avg = matched_data[
            (matched_data['year'] >= start) & 
            (matched_data['year'] <= end) &
            (matched_data['sex'] == 'Both')
        ].groupby('COUNTY_FIPS')['total_mean'].mean().round(2).reset_index()
        period_avg.columns = ['COUNTY_FIPS', f'{period}_SR']
        eqi_data.append(period_avg)
    
    eqi_result = pd.merge(eqi_data[0], eqi_data[1], on='COUNTY_FIPS', how='outer')
    
    # 2. 2010年单年数据：按性别分组
    year_2010 = matched_data[matched_data['year'] == 2010]
    sex_mapping = {'Both': 'SR_Total', 'Males': 'SR_Male', 'Females': 'SR_Female'}
    
    single_data = []
    for sex, col_name in sex_mapping.items():
        sex_data = year_2010[year_2010['sex'] == sex].groupby('COUNTY_FIPS')['total_mean'].mean().reset_index()
        sex_data.columns = ['COUNTY_FIPS', col_name]
        single_data.append(sex_data)
    
    single_result = single_data[0]
    for data in single_data[1:]:
        single_result = pd.merge(single_result, data, on='COUNTY_FIPS', how='outer')
    
    print(f"  📅 EQI数据: {len(eqi_result)} 县")
    print(f"  � 2010数据: {len(single_result)} 县")
    
    return eqi_result, single_result

def main():
    """主函数"""
    print("🚬 IHME 吸烟率数据清洗")
    print("=" * 40)
    
    # 设置路径
    project_root = Path(__file__).resolve().parents[2]
    ihme_path = project_root / "Data/Original/IHME Smoking/IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012.csv"
    location_path = project_root / "Data/Processed/CDC/Location.csv"
    output_dir = project_root / "Data/Processed/Smoking"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"📂 处理: {ihme_path.name}")
    
    # 加载数据
    ihme_df = pd.read_csv(ihme_path)
    location_df = pd.read_csv(location_path)
    
    # 只保留县级数据
    county_data = ihme_df[ihme_df['county'].notna() & (ihme_df['county'] != '')].copy()
    print(f"🏘️  县级数据: {len(county_data):,} 行")
    
    # 创建映射
    state_mapping = create_state_mapping()
    fips_dict = create_fips_mapping(location_df)
    
    # 处理数据
    eqi_data, single_data = process_data(county_data, fips_dict, state_mapping)
    
    # 保存结果
    print("\n💾 保存结果...")
    eqi_data.to_csv(output_dir / "County_Smoking_EQI.csv", index=False)
    single_data.to_csv(output_dir / "County_Smoking.csv", index=False)
    
    print(f"  ✅ EQI数据: {len(eqi_data)} 县")
    print(f"  ✅ 2010数据: {len(single_data)} 县")
    print("\n🎯 完成！")

if __name__ == "__main__":
    main()