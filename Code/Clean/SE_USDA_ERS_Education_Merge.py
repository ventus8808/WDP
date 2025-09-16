#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USDA ERS 教育数据合并脚本
- 路径从 config.yaml 读取：data_sources.socioeconomic.usda_ers.original/processed
- 合并美国农业部经济研究局的教育数据
- 将长格式数据转换为宽表格式
"""

import sys
from pathlib import Path
import pandas as pd
import os
import yaml

def load_paths():
    """从config.yaml加载路径"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        sys.exit(f"ERROR: Configuration file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    try:
        usda_config = cfg["data_sources"]["socioeconomic"]["usda_ers"]
        input_rel = usda_config["original"]
        output_rel = usda_config["processed"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.socioeconomic.usda_ers")

    input_dir = (project_root / input_rel).resolve()
    output_dir = (project_root / output_rel).resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir

def merge_education_data(input_dir, output_dir):
    """合并教育数据"""
    print(f"开始合并教育数据文件夹: {input_dir}")
    
    # 读取数据
    print("读取原始数据...")
    input_file = input_dir / "Education2023.csv"
    
    if not input_file.exists():
        print(f"ERROR: 文件不存在: {input_file}")
        return None
    
    try:
        # 尝试不同编码
        for encoding in ['latin-1', 'utf-8', 'cp1252', 'iso-8859-1']:
            try:
                df = pd.read_csv(input_file, encoding=encoding)
                print(f"  ✅ 使用编码: {encoding}")
                break
            except UnicodeDecodeError:
                continue
        else:
            print("  ❌ 无法读取文件，尝试了所有编码格式")
            return None
    except Exception as e:
        print(f"  ❌ 读取文件失败: {e}")
        return None
    
    print(f"原始数据形状: {df.shape}")
    
    # 过滤县级数据（FIPS代码为5位数字）
    print("过滤县级数据...")
    df_county = df[df['FIPS Code'].astype(str).str.len() == 5].copy()
    print(f"县级数据形状: {df_county.shape}")
    
    # 标准化FIPS代码
    df_county['COUNTY_FIPS'] = df_county['FIPS Code'].astype(str).str.zfill(5)
    
    # 提取年份
    def extract_year(attribute):
        if ', ' in attribute:
            year_part = attribute.split(', ')[-1]
            if year_part.isdigit():
                return int(year_part)
            elif year_part == '2008-12':
                return 2010
            elif year_part == '2019-23':
                return 2021
        return None
    
    df_county['Year'] = df_county['Attribute'].apply(extract_year)
    
    # 过滤1999-2020年的数据
    df_filtered = df_county[df_county['Year'].between(1999, 2020)].copy()
    print(f"1999-2020年数据形状: {df_filtered.shape}")
    
    # 创建结果数据框
    result_data = []
    
    # 按县和年份分组处理
    for (fips, year), group in df_filtered.groupby(['COUNTY_FIPS', 'Year']):
        row_data = {'COUNTY_FIPS': fips, 'Year': year}
        
        # 提取教育水平数据
        for _, row in group.iterrows():
            attr = row['Attribute']
            value = row['Value']
            
            if 'percent' in attr.lower() or 'Percent' in attr:
                # 高中以下学历
                if any(keyword in attr.lower() for keyword in ['less than', 'not high school']):
                    row_data['Less_Than_High_School_Percent'] = value
                # 仅高中学历
                elif any(keyword in attr.lower() for keyword in ['high school diploma only', 'high school graduates (or equivalent)']):
                    row_data['High_School_Only_Percent'] = value
                # 部分大学学历
                elif any(keyword in attr.lower() for keyword in ['some college', 'associate degree']):
                    row_data['Some_College_Percent'] = value
                # 大学及以上学历
                elif any(keyword in attr.lower() for keyword in ['four years of college', "bachelor's degree", 'college or higher']):
                    row_data['College_Plus_Percent'] = value
        
        # 提取城乡分类代码
        rural_urban_rows = group[group['Attribute'].str.contains('Rural-urban Continuum Code', na=False)]
        if not rural_urban_rows.empty:
            row_data['Rural_Urban_Continuum_Code'] = rural_urban_rows.iloc[0]['Value']
        
        urban_influence_rows = group[group['Attribute'].str.contains('Urban Influence Code', na=False)]
        if not urban_influence_rows.empty:
            row_data['Urban_Influence_Code'] = urban_influence_rows.iloc[0]['Value']
        
        result_data.append(row_data)
    
    # 转换为DataFrame
    result_df = pd.DataFrame(result_data)
    
    print(f"处理后的数据形状: {result_df.shape}")
    
    # 使用线性插值填充缺失数据
    print("开始线性插值填充缺失数据...")
    
    # 获取所有县的FIPS代码
    all_counties = result_df['COUNTY_FIPS'].unique()
    all_years = list(range(1999, 2021))  # 1999-2020年
    
    # 创建完整的县-年份组合
    complete_index = pd.MultiIndex.from_product([all_counties, all_years], names=['COUNTY_FIPS', 'Year'])
    complete_df = pd.DataFrame(index=complete_index).reset_index()
    
    # 合并现有数据
    result_df = complete_df.merge(result_df, on=['COUNTY_FIPS', 'Year'], how='left')
    
    # 定义需要插值的数值列
    numeric_columns = [
        'Less_Than_High_School_Percent',
        'High_School_Only_Percent', 
        'Some_College_Percent',
        'College_Plus_Percent'
    ]
    
    # 对每个县进行线性插值
    interpolated_data = []
    for county_fips in all_counties:
        county_data = result_df[result_df['COUNTY_FIPS'] == county_fips].copy()
        county_data = county_data.sort_values('Year')
        
        # 对数值列进行线性插值
        for col in numeric_columns:
            county_data[col] = county_data[col].interpolate(method='linear', limit_direction='both')
        
        # 对分类变量使用前向填充，如果列存在
        if 'Rural_Urban_Continuum_Code' in county_data.columns:
            county_data['Rural_Urban_Continuum_Code'] = county_data['Rural_Urban_Continuum_Code'].fillna(method='ffill').fillna(method='bfill')
        
        if 'Urban_Influence_Code' in county_data.columns:
            county_data['Urban_Influence_Code'] = county_data['Urban_Influence_Code'].fillna(method='ffill').fillna(method='bfill')
        
        interpolated_data.append(county_data)
    
    # 合并所有县的数据
    result_df = pd.concat(interpolated_data, ignore_index=True)
    
    # 排序
    result_df = result_df.sort_values(['COUNTY_FIPS', 'Year'])
    
    print(f"插值完成，总共 {len(result_df)} 行")
    print(f"覆盖年份: 1999-2020")
    print(f"县数量: {result_df['COUNTY_FIPS'].nunique()}")
    
    # 四舍五入保留两位小数
    for col in numeric_columns:
        if col in result_df.columns:
            result_df[col] = result_df[col].round(2)

    # 保存合并后的数据
    output_file = output_dir / "Education.csv"
    result_df.to_csv(output_file, index=False)
    print(f"数据已保存到: {output_file}")
    
    return result_df

def main():
    """主函数"""
    # 加载路径
    input_dir, output_dir = load_paths()
    
    # 合并数据
    merged_data = merge_education_data(input_dir, output_dir)
    
    if merged_data is not None:
        print("\n教育数据合并成功!")
        print(f"数据形状: {merged_data.shape}")
        print(f"年份范围: {merged_data['Year'].min()} - {merged_data['Year'].max()}")
        print(f"县数量: {merged_data['COUNTY_FIPS'].nunique()}")
    else:
        print("\n教育数据合并失败!")

if __name__ == "__main__":
    main()
