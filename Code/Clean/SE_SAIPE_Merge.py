#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAIPE 贫困和收入数据合并脚本
- 路径从 config.yaml 读取：data_sources.socioeconomic.saipe.original/processed
- 处理1999-2019年的贫困和收入估计数据
- 支持.dat文件（1999-2002）和Excel文件（2003-2019）
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
        saipe_config = cfg["data_sources"]["socioeconomic"]["saipe"]
        input_rel = saipe_config["original"]
        output_rel = saipe_config["processed"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.socioeconomic.saipe")

    input_dir = (project_root / input_rel).resolve()
    output_dir = (project_root / output_rel).resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir

def parse_dat_file(file_path, year):
    """解析固定宽度格式的.dat文件（1999-2002年）"""
    data = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) < 50:
                continue
                
            try:
                # 解析FIPS代码
                state_fips = line[0:2].strip()
                county_fips = line[2:5].strip()
                
                # 跳过州级汇总数据
                if county_fips == '0':
                    continue
                
                if not state_fips.isdigit() or not county_fips.isdigit():
                    continue
                
                # 构建完整的FIPS代码
                full_fips = f"{int(state_fips):02d}{int(county_fips):03d}"
                
                # 解析数据字段
                parts = line[5:].split()
                
                if len(parts) < 10:
                    continue
                
                # 查找贫困率和收入数据
                poverty_rate = None
                median_income = None
                
                # 贫困率通常在10-30%范围内
                for i, part in enumerate(parts[:len(parts)//2]):
                    try:
                        rate = float(part)
                        if 5 <= rate <= 50:
                            poverty_rate = rate
                            break
                    except ValueError:
                        continue
                
                # 收入数据通常是大数字
                for part in parts[len(parts)//2:]:
                    try:
                        clean_part = part.replace(',', '')
                        income = int(clean_part)
                        if 10000 <= income <= 200000:
                            median_income = income
                            break
                    except (ValueError, TypeError):
                        continue
                
                if poverty_rate is not None and median_income is not None:
                    data.append({
                        'COUNTY_FIPS': full_fips,
                        'Year': year,
                        'Poverty_Percent_All_Ages': poverty_rate,
                        'Median_Household_Income': median_income
                    })
                    
            except Exception as e:
                continue
    
    return pd.DataFrame(data)

def process_excel_file(file_path, year):
    """智能处理Excel文件，适应不同年份的不同格式"""
    try:
        # 根据年份确定header行和列名映射
        if year in [2003, 2004]:
            header_row = 1
            state_col = 'State FIPS'
            county_col = 'County FIPS'
            poverty_col = 'Poverty Percent All Ages'
            income_col = 'Median Household Income'
        elif year in [2005, 2006, 2007, 2008, 2009, 2010, 2011]:
            header_row = 2
            state_col = 'State FIPS'
            county_col = 'County FIPS'
            poverty_col = 'Poverty Percent All Ages'
            income_col = 'Median Household Income'
        elif year in [2012]:
            header_row = 2
            state_col = 'State FIPS Code'
            county_col = 'County FIPS Code'
            poverty_col = 'Poverty Percent, All Ages'
            income_col = 'Median Household Income'
        elif year in [2013, 2014, 2015, 2016, 2017, 2018, 2019]:
            header_row = 3
            state_col = 'State FIPS Code'
            county_col = 'County FIPS Code'
            poverty_col = 'Poverty Percent, All Ages'
            income_col = 'Median Household Income'
        else:
            print(f"    - 未知年份格式: {year}")
            return None
        
        # 读取Excel文件
        df = pd.read_excel(file_path, header=header_row)
        
        # 检查必要的列是否存在
        required_cols = [state_col, county_col, poverty_col, income_col]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"    - 缺少必要的列: {missing_cols}")
            return None
        
        print(f"    - 使用header={header_row}成功识别列名")
        
        # 添加年份列
        df['Year'] = year
        
        # 标准化FIPS代码
        state_fips = df[state_col].astype(str).str.replace('nan', '0').str.replace('None', '0')
        county_fips = df[county_col].astype(str).str.replace('nan', '0').str.replace('None', '0')
        
        # 清理FIPS代码
        state_fips = state_fips.str.extract(r'(\d+)')[0].fillna('0').astype(int).astype(str).str.zfill(2)
        county_fips = county_fips.str.extract(r'(\d+)')[0].fillna('0').astype(int).astype(str).str.zfill(3)
        
        df['COUNTY_FIPS'] = state_fips + county_fips
        
        # 选择需要的列并重命名
        df_clean = df[['COUNTY_FIPS', 'Year', poverty_col, income_col]].copy()
        df_clean.columns = ['COUNTY_FIPS', 'Year', 'Poverty_Percent_All_Ages', 'Median_Household_Income']
        
        # 过滤掉州级汇总数据
        df_clean = df_clean[df_clean['COUNTY_FIPS'].str[-3:] != '000']
        
        # 数据类型转换
        df_clean['Poverty_Percent_All_Ages'] = pd.to_numeric(df_clean['Poverty_Percent_All_Ages'], errors='coerce')
        df_clean['Median_Household_Income'] = pd.to_numeric(df_clean['Median_Household_Income'], errors='coerce')
        
        # 过滤掉无效数据
        df_clean = df_clean[
            (df_clean['Poverty_Percent_All_Ages'].notna()) & 
            (df_clean['Median_Household_Income'].notna()) &
            (df_clean['Poverty_Percent_All_Ages'] >= 0) &
            (df_clean['Poverty_Percent_All_Ages'] <= 100) &
            (df_clean['Median_Household_Income'] > 0)
        ]
        
        # 去除重复记录
        df_clean = df_clean.drop_duplicates(subset=['COUNTY_FIPS', 'Year'])
        
        return df_clean
        
    except Exception as e:
        print(f"    - 处理Excel文件失败: {e}")
        return None

def merge_saipe_data(input_dir, output_dir):
    """合并 SAIPE 数据"""
    print(f"开始合并 SAIPE 数据文件夹: {input_dir}")
    
    # 获取1999-2019年的文件列表
    years = range(1999, 2020)
    all_data = []
    
    for year in years:
        print(f"处理 {year} 年数据...")
        
        # 处理.dat文件（1999-2002年）
        if year <= 2002:
            file_path = input_dir / f"est{str(year)[-2:]}all.dat"
            if file_path.exists():
                try:
                    df = parse_dat_file(file_path, year)
                    if not df.empty:
                        all_data.append(df)
                        print(f"  - 成功解析 {len(df)} 条记录")
                    else:
                        print(f"  - 警告: 没有解析到有效数据")
                except Exception as e:
                    print(f"  - 错误: 解析.dat文件失败 - {e}")
            else:
                print(f"  - 警告: 文件 {file_path} 不存在")
            continue
        else:
            # .xls文件（Excel格式）
            file_path = input_dir / f"est{str(year)[-2:]}all.xls"
            if file_path.exists():
                df_clean = process_excel_file(file_path, year)
                if df_clean is not None and not df_clean.empty:
                    all_data.append(df_clean)
                    print(f"  - 成功解析 {len(df_clean)} 条记录")
                else:
                    print(f"  - 警告: 没有解析到有效数据")
            else:
                print(f"  - 警告: 文件 {file_path} 不存在")
    
    if not all_data:
        print("错误: 没有成功读取任何数据文件")
        return None
    
    # 合并所有年份的数据
    print("合并所有年份数据...")
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # 去除重复记录
    print("去除重复记录...")
    combined_df = combined_df.drop_duplicates(subset=['COUNTY_FIPS', 'Year'])
    
    # 排序
    combined_df = combined_df.sort_values(['COUNTY_FIPS', 'Year'])
    
    print(f"合并完成，总共 {len(combined_df)} 行")
    
    # 保存合并后的数据
    output_file = output_dir / "Poverty_Income.csv"
    combined_df.to_csv(output_file, index=False)
    print(f"数据已保存到: {output_file}")
    
    return combined_df

def main():
    """主函数"""
    # 加载路径
    input_dir, output_dir = load_paths()
    
    # 合并数据
    merged_data = merge_saipe_data(input_dir, output_dir)
    
    if merged_data is not None:
        print("\nSAIPE 数据合并成功!")
        print(f"数据形状: {merged_data.shape}")
        print(f"年份范围: {merged_data['Year'].min()} - {merged_data['Year'].max()}")
        print(f"县数量: {merged_data['COUNTY_FIPS'].nunique()}")
    else:
        print("\nSAIPE 数据合并失败!")

if __name__ == "__main__":
    main()
