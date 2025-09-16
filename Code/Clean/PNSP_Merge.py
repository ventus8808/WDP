#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PNSP 农药数据合并与重塑脚本
- 合并1999-2019年所有PNSP数据文件
- 使用mapping.csv进行化合物映射和分类
- 输出格式：COUNTY_FIPS | Year | cat1_min/avg/max | ... | chem1_min/avg/max | ...
- 路径从config.yaml读取：data_sources.usgs_pnsp.original/processed
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

try:
    import yaml
except Exception as e:
    print("ERROR: 需要PyYAML。请安装: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

def load_paths():
    """从config.yaml加载路径"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        print(f"ERROR: 未找到配置文件: {config_path}", file=sys.stderr)
        sys.exit(1)
    
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    
    ds = (cfg.get("data_sources") or {}).get("usgs_pnsp") or {}
    original_rel = ds.get("original")
    processed_rel = ds.get("processed")
    
    if not original_rel or not processed_rel:
        print("ERROR: config.yaml缺少data_sources.usgs_pnsp.original/processed", file=sys.stderr)
        sys.exit(1)
    
    input_dir = (project_root / original_rel).resolve()
    output_dir = (project_root / processed_rel).resolve()
    
    if not input_dir.exists():
        print(f"ERROR: 输入目录不存在: {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir

def load_mapping():
    """加载化合物映射表"""
    project_root = Path(__file__).resolve().parents[2]
    mapping_path = project_root / "Data/Processed/Pesticide/mapping.csv"
    
    if not mapping_path.exists():
        print(f"ERROR: 映射文件不存在: {mapping_path}", file=sys.stderr)
        sys.exit(1)
    
    mapping = pd.read_csv(mapping_path)
    print(f"加载映射表: {len(mapping)}个化合物，{mapping['category1_id'].nunique()}个类别")
    return mapping

def merge_pnsp_data(input_dir):
    """合并所有PNSP数据文件"""
    print("开始合并PNSP数据...")
    
    all_data = []
    
    # 处理单独年份文件 (1999-2012, 2018-2019)
    single_years = list(range(1999, 2013)) + list(range(2018, 2020))
    
    for year in single_years:
        # 尝试不同的文件名格式
        possible_files = [
            f"EPest.county.estimates.{year}.txt",
            f"EPest_county_estimates_{year}.txt"
        ]
        
        file_found = False
        for filename in possible_files:
            file_path = input_dir / filename
            if file_path.exists():
                print(f"处理文件: {filename}")
                try:
                    # 读取数据，自动检测分隔符
                    df = pd.read_csv(file_path, sep=None, engine='python', encoding='utf-8')
                    
                    # 标准化列名
                    df.columns = [col.strip() for col in df.columns]
                    
                    # 检查必需列
                    required_cols = ['COMPOUND', 'YEAR', 'STATE_FIPS_CODE', 'COUNTY_FIPS_CODE', 'EPEST_LOW_KG', 'EPEST_HIGH_KG']
                    if not all(col in df.columns for col in required_cols):
                        print(f"  ❌ 缺少必需列: {[col for col in required_cols if col not in df.columns]}")
                        continue
                    
                    # 只保留需要的列
                    df = df[required_cols].copy()
                    
                    # 标准化FIPS代码
                    df['COUNTY_FIPS'] = (df['STATE_FIPS_CODE'].astype(str).str.zfill(2) + 
                                        df['COUNTY_FIPS_CODE'].astype(str).str.zfill(3))
                    
                    # 确保数值列
                    df['EPEST_LOW_KG'] = pd.to_numeric(df['EPEST_LOW_KG'], errors='coerce').fillna(0)
                    df['EPEST_HIGH_KG'] = pd.to_numeric(df['EPEST_HIGH_KG'], errors='coerce').fillna(0)
                    df['YEAR'] = pd.to_numeric(df['YEAR'], errors='coerce').fillna(year).astype(int)
                    
                    # 过滤有效数据
                    df = df[(df['EPEST_LOW_KG'] >= 0) & (df['EPEST_HIGH_KG'] >= 0)]
                    
                    all_data.append(df)
                    print(f"  ✅ 成功处理，{len(df)}行")
                    file_found = True
                    break
                    
                except Exception as e:
                    print(f"  ❌ 错误: {e}")
                    continue
        
        if not file_found:
            print(f"⚠️  未找到{year}年的数据文件")
    
    # 处理2013-2017合并文件
    combined_file = input_dir / "EPest_county_estimates_2013_2017_v2.txt"
    if combined_file.exists():
        print(f"处理合并文件: {combined_file.name}")
        try:
            # 读取数据，自动检测分隔符
            df = pd.read_csv(combined_file, sep=None, engine='python', encoding='utf-8')
            
            # 标准化列名
            df.columns = [col.strip() for col in df.columns]
            
            # 检查必需列
            required_cols = ['COMPOUND', 'YEAR', 'STATE_FIPS_CODE', 'COUNTY_FIPS_CODE', 'EPEST_LOW_KG', 'EPEST_HIGH_KG']
            if not all(col in df.columns for col in required_cols):
                print(f"  ❌ 缺少必需列: {[col for col in required_cols if col not in df.columns]}")
            else:
                # 只保留需要的列
                df = df[required_cols].copy()
                
                # 标准化FIPS代码
                df['COUNTY_FIPS'] = (df['STATE_FIPS_CODE'].astype(str).str.zfill(2) + 
                                    df['COUNTY_FIPS_CODE'].astype(str).str.zfill(3))
                
                # 确保数值列
                df['EPEST_LOW_KG'] = pd.to_numeric(df['EPEST_LOW_KG'], errors='coerce').fillna(0)
                df['EPEST_HIGH_KG'] = pd.to_numeric(df['EPEST_HIGH_KG'], errors='coerce').fillna(0)
                df['YEAR'] = pd.to_numeric(df['YEAR'], errors='coerce').fillna(2013).astype(int)
                
                # 过滤有效数据
                df = df[(df['EPEST_LOW_KG'] >= 0) & (df['EPEST_HIGH_KG'] >= 0)]
                
                all_data.append(df)
                print(f"  ✅ 成功处理，{len(df)}行")
        except Exception as e:
            print(f"  ❌ 错误: {e}")
    
    if not all_data:
        print("ERROR: 没有成功处理任何文件")
        sys.exit(1)
    
    # 合并所有数据
    print("\n合并所有数据...")
    merged_df = pd.concat(all_data, ignore_index=True)
    print(f"合并完成，总共{len(merged_df)}行")
    
    return merged_df

def reshape_data(merged_df, mapping):
    """重塑数据为宽格式"""
    print("开始重塑数据...")
    
    # 合并映射表
    merged_df = merged_df.merge(mapping[['compound_name', 'compound_id', 'category1_id']], 
                                left_on='COMPOUND', right_on='compound_name', how='left')
    
    # 检查未匹配的化合物
    unmatched = merged_df[merged_df['compound_id'].isna()]['COMPOUND'].unique()
    if len(unmatched) > 0:
        print(f"⚠️  发现{len(unmatched)}个未匹配的化合物: {unmatched[:5]}...")
        print(f"   总未匹配化合物数: {len(unmatched)}")
    
    # 只保留有映射的数据
    original_rows = len(merged_df)
    merged_df = merged_df.dropna(subset=['compound_id']).copy()
    matched_rows = len(merged_df)
    print(f"原始数据: {original_rows}行，匹配后: {matched_rows}行")
    
    merged_df['compound_id'] = merged_df['compound_id'].astype(int)
    merged_df['category1_id'] = merged_df['category1_id'].astype(int)
    
    # 检查匹配后的化合物和类别数量
    unique_compounds = merged_df['compound_id'].nunique()
    unique_categories = merged_df['category1_id'].nunique()
    print(f"匹配的化合物数: {unique_compounds}, 类别数: {unique_categories}")
    
    print(f"有效数据: {len(merged_df)}行")
    
    # 按类别聚合（将同一类别下所有化合物的用量相加）
    print("计算类别级别汇总（同一类别下所有化合物用量相加）...")
    category_agg = merged_df.groupby(['COUNTY_FIPS', 'YEAR', 'category1_id']).agg({
        'EPEST_LOW_KG': 'sum',
        'EPEST_HIGH_KG': 'sum'
    }).reset_index()
    
    # 计算平均值
    category_agg['EPEST_AVG_KG'] = (category_agg['EPEST_LOW_KG'] + category_agg['EPEST_HIGH_KG']) / 2
    
    # 重塑类别数据为宽格式
    category_wide = category_agg.pivot_table(
        index=['COUNTY_FIPS', 'YEAR'],
        columns='category1_id',
        values=['EPEST_LOW_KG', 'EPEST_AVG_KG', 'EPEST_HIGH_KG'],
        fill_value=0
    ).reset_index()
    
    # 重命名列 - 修正列名格式
    new_cols = []
    for col in category_wide.columns:
        if col[1] != '':  # 多层索引的数据列
            if 'LOW' in col[0]:
                suffix = 'min'
            elif 'HIGH' in col[0]:
                suffix = 'max'
            elif 'AVG' in col[0]:
                suffix = 'avg'
            else:
                suffix = col[0].split('_')[1].lower()
            new_cols.append(f"cat{col[1]}_{suffix}")
        else:  # 单层索引的标识列
            new_cols.append(col[0])
    category_wide.columns = new_cols
    
    # 按化合物聚合
    print("计算化合物级别数据...")
    compound_agg = merged_df.groupby(['COUNTY_FIPS', 'YEAR', 'compound_id']).agg({
        'EPEST_LOW_KG': 'sum',
        'EPEST_HIGH_KG': 'sum'
    }).reset_index()
    
    # 计算平均值
    compound_agg['EPEST_AVG_KG'] = (compound_agg['EPEST_LOW_KG'] + compound_agg['EPEST_HIGH_KG']) / 2
    
    # 重塑化合物数据为宽格式
    compound_wide = compound_agg.pivot_table(
        index=['COUNTY_FIPS', 'YEAR'],
        columns='compound_id',
        values=['EPEST_LOW_KG', 'EPEST_AVG_KG', 'EPEST_HIGH_KG'],
        fill_value=0
    ).reset_index()
    
    # 重命名列 - 修正列名格式
    new_cols = []
    for col in compound_wide.columns:
        if col[1] != '':  # 多层索引的数据列
            if 'LOW' in col[0]:
                suffix = 'min'
            elif 'HIGH' in col[0]:
                suffix = 'max'
            elif 'AVG' in col[0]:
                suffix = 'avg'
            else:
                suffix = col[0].split('_')[1].lower()
            new_cols.append(f"chem{col[1]}_{suffix}")
        else:  # 单层索引的标识列
            new_cols.append(col[0])
    compound_wide.columns = new_cols
    
    # 合并类别和化合物数据
    print("合并类别和化合物数据...")
    final_df = category_wide.merge(compound_wide, on=['COUNTY_FIPS', 'YEAR'], how='outer')
    
    # 填充缺失值并保留一位小数
    final_df = final_df.fillna(0)
    
    # 对所有数值列保留一位小数
    numeric_cols = [col for col in final_df.columns if col not in ['COUNTY_FIPS', 'YEAR']]
    for col in numeric_cols:
        final_df[col] = final_df[col].round(1)
    
    print(f"重塑完成，最终数据: {len(final_df)}行，{len(final_df.columns)}列")
    return final_df

def main():
    """主函数"""
    print("=== PNSP数据合并与重塑脚本 ===\n")
    
    # 加载路径和映射
    input_dir, output_dir = load_paths()
    mapping = load_mapping()
    
    # 合并数据
    merged_df = merge_pnsp_data(input_dir)
    
    # 重塑数据
    final_df = reshape_data(merged_df, mapping)
    
    # 保存结果
    output_file = output_dir / "PNSP.csv"
    final_df.to_csv(output_file, index=False)
    print(f"\n数据已保存到: {output_file}")
    
    # 输出统计信息
    print(f"\n=== 处理完成 ===")
    print(f"输出数据形状: {final_df.shape}")
    print(f"县数量: {final_df['COUNTY_FIPS'].nunique()}")
    print(f"年份范围: {final_df['YEAR'].min()} - {final_df['YEAR'].max()}")
    print(f"类别列数: {len([col for col in final_df.columns if col.startswith('cat')])}")
    print(f"化合物列数: {len([col for col in final_df.columns if col.startswith('chem')])}")

if __name__ == "__main__":
    main()
