#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LAUS 失业率数据合并脚本
- 路径从 config.yaml 读取：data_sources.socioeconomic.laus.original/processed
- 合并 BLS 的劳动力统计和失业率数据
- 处理1999-2019年的失业率数据
"""

import sys
from pathlib import Path
import pandas as pd
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
        laus_config = cfg["data_sources"]["socioeconomic"]["laus"]
        input_rel = laus_config["original"]
        output_rel = laus_config["processed"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.socioeconomic.laus")

    input_dir = (project_root / input_rel).resolve()
    output_dir = (project_root / output_rel).resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir

def merge_laus_data(input_dir, output_dir):
    """合并 LAUS 数据"""
    print(f"开始合并 LAUS 数据文件夹: {input_dir}")
    
    # 获取1999-2019年的文件列表
    years = range(1999, 2020)
    all_data = []
    
    for year in years:
        # 修复文件命名规则：1999年是99，2000-2009年是00-09，2010+是10-24
        if year == 1999:
            year_str = "99"
        elif year < 2010:
            year_str = f"0{year % 100}"
        else:
            year_str = str(year % 100)
            
        file_path = input_dir / f"laucnty{year_str}.xlsx"
        
        if not file_path.exists():
            print(f"警告: 文件 {file_path} 不存在，跳过")
            continue
            
        print(f"处理 {year} 年数据...")
        
        try:
            # 读取Excel文件，使用第2行作为列名
            df = pd.read_excel(file_path, header=1)
            
            # 添加年份列
            df['Year'] = year
            
            # 标准化FIPS代码 - 处理浮点数格式和缺失值
            state_fips = df['State FIPS Code'].fillna(0).astype(float).astype(int).astype(str).str.zfill(2)
            county_fips = df['County FIPS Code'].fillna(0).astype(float).astype(int).astype(str).str.zfill(3)
            df['COUNTY_FIPS'] = state_fips + county_fips
            
            # 选择需要的列并重命名
            df_clean = df[['COUNTY_FIPS', 'Year', 'Labor Force', 'Employed', 'Unemployed', 'Unemployment Rate (%)']].copy()
            df_clean.columns = ['COUNTY_FIPS', 'Year', 'Labor_Force', 'Employed', 'Unemployed', 'Unemployment_Rate']
            
            # 数据类型转换
            numeric_cols = ['Labor_Force', 'Employed', 'Unemployed', 'Unemployment_Rate']
            for col in numeric_cols:
                if col in ['Labor_Force', 'Employed', 'Unemployed']:
                    # 劳动力、就业、失业人数转换为整数
                    df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce').astype('Int64')
                else:
                    # 失业率保持浮点数
                    df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
            
            # 过滤掉无效的FIPS代码（通常以'CN'开头的是汇总行，或者FIPS为00000）
            df_clean = df_clean[~df_clean['COUNTY_FIPS'].str.startswith('CN', na=False)]
            df_clean = df_clean[df_clean['COUNTY_FIPS'] != '00000']
            
            all_data.append(df_clean)
            print(f"  ✅ 成功处理，{len(df_clean)} 行")
            
        except Exception as e:
            print(f"  ❌ 处理 {year} 年数据时出错: {e}")
            continue
    
    if not all_data:
        print("错误: 没有成功读取任何数据文件")
        return None
    
    # 合并所有年份的数据
    print("合并所有年份数据...")
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # 排序
    combined_df = combined_df.sort_values(['COUNTY_FIPS', 'Year'])
    
    print(f"合并完成，总共 {len(combined_df)} 行")
    
    # 保存合并后的数据
    output_file = output_dir / "Unemployment.csv"
    combined_df.to_csv(output_file, index=False)
    print(f"数据已保存到: {output_file}")
    
    return combined_df

def main():
    """主函数"""
    # 加载路径
    input_dir, output_dir = load_paths()
    
    # 合并数据
    merged_data = merge_laus_data(input_dir, output_dir)
    
    if merged_data is not None:
        print("\nLAUS 数据合并成功!")
        print(f"数据形状: {merged_data.shape}")
        print(f"年份范围: {merged_data['Year'].min()} - {merged_data['Year'].max()}")
        print(f"县数量: {merged_data['COUNTY_FIPS'].nunique()}")
    else:
        print("\nLAUS 数据合并失败!")

if __name__ == "__main__":
    main()
