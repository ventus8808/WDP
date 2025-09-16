#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
处理CACES LUR数据
删除经纬度和州信息，将污染物数据整理成宽表格式
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
        lur_config = cfg["data_sources"]["air_pollution"]
        input_rel = lur_config["original"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.environmental.lur")

    input_dir = (project_root / input_rel).resolve()
    # The output path is explicitly defined
    output_path = (project_root / "Data/Processed/Environmental/Air_Pollution.csv").resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    # Ensure the output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return input_dir, output_path

def process_lur_data(input_dir, output_file):
    """处理LUR数据，转换为宽表格式"""
    print("开始处理CACES LUR数据...")

    # Define input files
    o3_file = input_dir / "1999-2019 O3.csv"
    other_pollutants_file = input_dir / "1999-2020 CO SO2 NO2 PM10 PM25.csv"

    if not o3_file.exists() or not other_pollutants_file.exists():
        sys.exit(f"ERROR: Input data not found in {input_dir}")

    # 1. 处理O3数据
    print("处理O3数据...")
    o3_df = pd.read_csv(o3_file)
    o3_df = o3_df.drop(columns=['state_abbr', 'lat', 'lon'])
    o3_df = o3_df.rename(columns={'fips': 'COUNTY_FIPS', 'year': 'Year', 'pred_wght': 'O3'})
    o3_df['COUNTY_FIPS'] = o3_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    o3_df = o3_df[o3_df['pollutant'] == 'o3'].copy()
    o3_df = o3_df.drop(columns=['pollutant'])
    print(f"O3数据: {len(o3_df)} 行")

    # 2. 处理其他污染物数据
    print("处理其他污染物数据...")
    other_df = pd.read_csv(other_pollutants_file)
    other_df = other_df.drop(columns=['state_abbr', 'lat', 'lon'])
    other_df = other_df.rename(columns={'fips': 'COUNTY_FIPS', 'year': 'Year', 'pred_wght': 'Value'})
    other_df['COUNTY_FIPS'] = other_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    other_df = other_df[other_df['pollutant'] != 'pollutant'].copy()
    other_df = other_df[other_df['Year'] <= 2019].copy()

    other_df_wide = other_df.pivot_table(index=['COUNTY_FIPS', 'Year'], columns='pollutant', values='Value', aggfunc='first').reset_index()
    other_df_wide = other_df_wide.rename(columns={'co': 'CO', 'so2': 'SO2', 'no2': 'NO2', 'pm10': 'PM10', 'pm25': 'PM25'})
    print(f"其他污染物数据: {len(other_df_wide)} 行")

    # 3. 合并所有污染物数据
    print("合并所有污染物数据...")
    combined_df = pd.merge(o3_df, other_df_wide, on=['COUNTY_FIPS', 'Year'], how='outer')

    # 4. 重新排列列的顺序并排序
    column_order = ['COUNTY_FIPS', 'Year', 'O3', 'CO', 'SO2', 'NO2', 'PM10', 'PM25']
    combined_df = combined_df[column_order]
    combined_df = combined_df.sort_values(['COUNTY_FIPS', 'Year']).reset_index(drop=True)

    # 5. 保存数据
    combined_df.to_csv(output_file, index=False, encoding='utf-8')

    print(f"\n处理完成！")
    print(f"最终数据形状: {combined_df.shape}")
    print(f"数据已保存到: {output_file}")

def main():
    """主函数"""
    input_dir, output_file = load_paths()
    process_lur_data(input_dir, output_file)

if __name__ == "__main__":
    main()
