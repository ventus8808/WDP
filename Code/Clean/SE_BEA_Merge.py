#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BEA 经济数据合并脚本
- 路径从 config.yaml 读取：data_sources.socioeconomic.bea.original/processed
- 合并经济分析局的经济指标数据
- 处理CAINC1（个人收入）和CAGDP1（GDP）数据
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import glob
import os
import warnings
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
        bea_config = cfg["data_sources"]["socioeconomic"]["bea"]
        input_rel = bea_config["original"]
        output_rel = bea_config["processed"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.socioeconomic.bea")

    input_dir = (project_root / input_rel).resolve()
    output_dir = (project_root / output_rel).resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir

def load_and_process_cainc1_data(input_dir):
    """Load and process CAINC1 (Personal Income) data"""
    print("Loading CAINC1 (Personal Income) data...")
    
    cainc1_files = glob.glob(str(input_dir / "CAINC1/CAINC1_*_1969_2023.csv"))
    all_data = []
    
    for file_path in cainc1_files:
        try:
            state = file_path.split('/')[-1].split('_')[1]
            if len(state) != 2 or state in ['MS', 'MSA', 'PORT', 'CSA', 'MDIV', 'MIC']:
                continue
                
            print(f"  Processing {state}...")
            
            for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                try:
                    df = pd.read_csv(file_path, encoding=encoding)
                    break
                except UnicodeDecodeError:
                    continue
            else:
                print(f"    Could not read {file_path} with any encoding")
                continue
            
            df['GeoFIPS'] = df['GeoFIPS'].str.strip().str.replace('"', '').str.replace(' ', '')
            df = df[df['GeoFIPS'].str.len() == 5]
            
            if df['LineCode'].dtype == 'object':
                df['LineCode'] = pd.to_numeric(df['LineCode'], errors='coerce')
            
            df = df[df['LineCode'].isin([1, 2, 3])]
            
            year_cols = [str(year) for year in range(1999, 2020)]
            available_years = [col for col in year_cols if col in df.columns]
            
            if len(available_years) == 0:
                print(f"    No year columns found for {state}")
                continue
                
            df = df[['GeoFIPS', 'GeoName', 'LineCode', 'Description'] + available_years]
            all_data.append(df)
            
        except Exception as e:
            print(f"  Error processing {file_path}: {e}")
            continue
    
    if not all_data:
        raise ValueError("No CAINC1 data loaded")
    
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"  Loaded {len(combined_df)} records from CAINC1")
    
    return combined_df

def load_and_process_cagdp1_data(input_dir):
    """Load and process CAGDP1 (GDP) data"""
    print("Loading CAGDP1 (GDP) data...")
    
    cagdp1_files = glob.glob(str(input_dir / "CAGDP1/CAGDP1_*_2001_2023.csv"))
    all_data = []
    
    for file_path in cagdp1_files:
        try:
            state = file_path.split('/')[-1].split('_')[1]
            if len(state) != 2 or state in ['MS', 'MSA', 'PORT', 'CSA', 'MDIV', 'MIC']:
                continue
                
            print(f"  Processing {state}...")
            
            for encoding in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                try:
                    df = pd.read_csv(file_path, encoding=encoding)
                    break
                except UnicodeDecodeError:
                    continue
            else:
                print(f"    Could not read {file_path} with any encoding")
                continue
            
            df['GeoFIPS'] = df['GeoFIPS'].str.strip().str.replace('"', '').str.replace(' ', '')
            df = df[df['GeoFIPS'].str.len() == 5]
            
            if df['LineCode'].dtype == 'object':
                df['LineCode'] = pd.to_numeric(df['LineCode'], errors='coerce')
            
            df = df[df['LineCode'].isin([1, 2])]
            
            year_cols = [str(year) for year in range(2001, 2020)]
            available_years = [col for col in year_cols if col in df.columns]
            
            if len(available_years) == 0:
                print(f"    No year columns found for {state}")
                continue
                
            df = df[['GeoFIPS', 'GeoName', 'LineCode', 'Description'] + available_years]
            all_data.append(df)
            
        except Exception as e:
            print(f"  Error processing {file_path}: {e}")
            continue
    
    if not all_data:
        raise ValueError("No CAGDP1 data loaded")
    
    combined_df = pd.concat(all_data, ignore_index=True)
    print(f"  Loaded {len(combined_df)} records from CAGDP1")
    
    return combined_df

def reshape_cainc1_to_long(cainc1_df):
    """Reshape CAINC1 data from wide to long format"""
    print("Reshaping CAINC1 data to long format...")
    
    if len(cainc1_df) == 0:
        raise ValueError("No data to reshape")
    
    id_vars = ['GeoFIPS', 'GeoName', 'LineCode', 'Description']
    value_vars = [col for col in cainc1_df.columns if col.isdigit() and 1999 <= int(col) <= 2019]
    
    if not value_vars:
        raise ValueError("No year columns found in CAINC1 data")
    
    long_df = cainc1_df.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name='Year',
        value_name='Value'
    )
    
    long_df['Year'] = long_df['Year'].astype(int)
    
    pivot_df = long_df.pivot_table(
        index=['GeoFIPS', 'GeoName', 'Year'],
        columns='LineCode',
        values='Value',
        aggfunc='first'
    ).reset_index()
    
    expected_columns = ['GeoFIPS', 'GeoName', 'Year', 1, 2, 3]
    for col in expected_columns:
        if col not in pivot_df.columns:
            pivot_df[col] = np.nan
    
    pivot_df.columns = ['GeoFIPS', 'GeoName', 'Year', 'Personal_Income', 'Population', 'Per_Capita_Income']
    
    pivot_df['Personal_Income'] = pd.to_numeric(pivot_df['Personal_Income'], errors='coerce')
    pivot_df['Population'] = pd.to_numeric(pivot_df['Population'], errors='coerce')
    pivot_df['Per_Capita_Income'] = pd.to_numeric(pivot_df['Per_Capita_Income'], errors='coerce')
    
    print(f"  Reshaped to {len(pivot_df)} county-year records")
    return pivot_df

def reshape_cagdp1_to_long(cagdp1_df):
    """Reshape CAGDP1 data from wide to long format"""
    print("Reshaping CAGDP1 data to long format...")
    
    if len(cagdp1_df) == 0:
        raise ValueError("No data to reshape")
    
    id_vars = ['GeoFIPS', 'GeoName', 'LineCode', 'Description']
    value_vars = [col for col in cagdp1_df.columns if col.isdigit() and 2001 <= int(col) <= 2019]
    
    if not value_vars:
        raise ValueError("No year columns found in CAGDP1 data")
    
    long_df = cagdp1_df.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name='Year',
        value_name='Value'
    )
    
    long_df['Year'] = long_df['Year'].astype(int)
    
    pivot_df = long_df.pivot_table(
        index=['GeoFIPS', 'GeoName', 'Year'],
        columns='LineCode',
        values='Value',
        aggfunc='first'
    ).reset_index()
    
    expected_columns = ['GeoFIPS', 'GeoName', 'Year', 1, 2]
    for col in expected_columns:
        if col not in pivot_df.columns:
            pivot_df[col] = np.nan
    
    pivot_df.columns = ['GeoFIPS', 'GeoName', 'Year', 'Real_GDP', 'Current_Dollar_GDP']
    
    pivot_df['Real_GDP'] = pd.to_numeric(pivot_df['Real_GDP'], errors='coerce')
    pivot_df['Current_Dollar_GDP'] = pd.to_numeric(pivot_df['Current_Dollar_GDP'], errors='coerce')
    
    print(f"  Reshaped to {len(pivot_df)} county-year records")
    return pivot_df

def merge_economic_data(cainc1_long, cagdp1_long):
    """Merge personal income and GDP data"""
    print("Merging economic data...")
    
    merged_df = pd.merge(
        cainc1_long,
        cagdp1_long,
        on=['GeoFIPS', 'GeoName', 'Year'],
        how='left'
    )
    
    merged_df['State_FIPS'] = merged_df['GeoFIPS'].str[:2]
    merged_df['County_FIPS'] = merged_df['GeoFIPS'].str[2:]
    merged_df['County'] = merged_df['GeoName'].str.replace(r', [A-Z]{2}$', '', regex=True)
    
    final_columns = [
        'GeoFIPS', 'County', 'State_FIPS', 'County_FIPS', 'Year',
        'Personal_Income', 'Population', 'Per_Capita_Income',
        'Real_GDP', 'Current_Dollar_GDP'
    ]
    
    merged_df = merged_df[final_columns]
    merged_df = merged_df.sort_values(['GeoFIPS', 'Year'])
    
    print(f"  Merged data contains {len(merged_df)} records")
    return merged_df

def clean_and_validate_data(merged_df):
    """Clean and validate the merged economic data"""
    print("Cleaning and validating data...")
    
    initial_count = len(merged_df)
    merged_df = merged_df.dropna(subset=['GeoFIPS'])
    
    # Filter for valid 5-digit FIPS codes and remove state/national summaries
    merged_df = merged_df[merged_df['GeoFIPS'].str.match(r'^\d{5}$')]
    merged_df = merged_df[~merged_df['GeoFIPS'].str.endswith('000')]

    merged_df['GeoFIPS'] = merged_df['GeoFIPS'].astype(str).str.zfill(5)
    merged_df['State_FIPS'] = merged_df['State_FIPS'].astype(str).str.zfill(2)
    merged_df['County_FIPS'] = merged_df['County_FIPS'].astype(str).str.zfill(3)
    
    merged_df['Personal_Income'] = merged_df['Personal_Income'].fillna(0)
    merged_df['Population'] = merged_df['Population'].fillna(0)
    merged_df['Per_Capita_Income'] = merged_df['Per_Capita_Income'].fillna(0)
    
    merged_df = merged_df[merged_df['Population'] > 0]
    
    final_count = len(merged_df)
    print(f"  Data cleaning: {initial_count} -> {final_count} records")
    
    return merged_df

def optimize_economic_data(final_df):
    """Optimize economic data for statistical analysis"""
    print("\n=== Optimizing Economic Data for Statistical Analysis ===")
    
    df = final_df.copy()
    
    print("1. Creating 5-digit FIPS code...")
    df['COUNTY_FIPS'] = df['State_FIPS'] + df['County_FIPS']
    
    print("2. Creating Total GDP in 10K USD...")
    df['Total_GDP_10K_USD'] = df['Real_GDP'] / 10
    df['Total_GDP_10K_USD'] = df['Total_GDP_10K_USD'].round(2)
    
    print("3. Converting data types...")
    df['Population'] = df['Population'].astype(int)
    df['Per_Capita_Income'] = df['Per_Capita_Income'].astype(int)
    
    print("4. Finalizing optimized dataset...")
    
    final_cols = ['COUNTY_FIPS', 'Year', 'Population', 'Total_GDP_10K_USD', 'Per_Capita_Income']
    
    existing_cols = [col for col in final_cols if col in df.columns]
    df = df[existing_cols]
    
    print(f"   Final dataset shape: {df.shape}")
    print(f"   Variables: {list(df.columns)}")
    
    return df

def interpolate_gdp_data(merged_df):
    """Interpolate missing GDP data for 1999 and 2000"""
    print("Interpolating missing GDP data for 1999-2000...")
    
    gdp_cols = ['Real_GDP', 'Current_Dollar_GDP']
    
    merged_df = merged_df.sort_values(['GeoFIPS', 'Year'])
    
    merged_df[gdp_cols] = merged_df.groupby('GeoFIPS')[gdp_cols].transform(
        lambda x: x.interpolate(method='linear', limit_direction='backward', limit=2)
    )
    
    print(f"  Interpolation complete.")
    return merged_df

def merge_bea_data(input_dir, output_dir):
    """合并 BEA 数据"""
    print(f"开始合并 BEA 数据文件夹: {input_dir}")
    
    try:
        cainc1_df = load_and_process_cainc1_data(input_dir)
        cagdp1_df = load_and_process_cagdp1_data(input_dir)
        
        cainc1_long = reshape_cainc1_to_long(cainc1_df)
        cagdp1_long = reshape_cagdp1_to_long(cagdp1_df)
        
        merged_df = merge_economic_data(cainc1_long, cagdp1_long)
        interpolated_df = interpolate_gdp_data(merged_df)
        cleaned_df = clean_and_validate_data(interpolated_df)
        
        optimized_df = optimize_economic_data(cleaned_df)
        
        output_file = output_dir / "GDP.csv"
        optimized_df.to_csv(output_file, index=False)
        
        print(f"\nOutput saved to: {output_file}")
        
        return optimized_df
        
    except Exception as e:
        print(f"\nError during processing: {e}")
        return None

def main():
    """主函数"""
    input_dir, output_dir = load_paths()
    
    merged_data = merge_bea_data(input_dir, output_dir)
    
    if merged_data is not None:
        print("\nBEA 数据合并成功!")
        print(f"数据形状: {merged_data.shape}")
        print(f"年份范围: {merged_data['Year'].min()} - {merged_data['Year'].max()}")
        print(f"县数量: {merged_data['COUNTY_FIPS'].nunique()}")
    else:
        print("\nBEA 数据合并失败!")

if __name__ == "__main__":
    main()
