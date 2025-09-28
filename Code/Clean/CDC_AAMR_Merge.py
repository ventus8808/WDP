#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDC WONDER 死亡与人口数据清洗脚本 (精简版)

目标:
- 为PyMC等贝叶斯模型准备数据。
- 仅提取 FIPS, Year, Population, 和 Deaths。
- 核心功能: 将'Suppressed'的死亡人数正确处理为区间审查数据 [1, 9]。

输出列:
- COUNTY_FIPS: 5位县级编码
- Year: 年份
- Population: 总人口
- Deaths_Observed: 观测到的死亡人数 (0 或 >=10), Suppressed时为NA
- Deaths_Censored_Lower: 审查数据下界 (Suppressed时为1)
- Deaths_Censored_Upper: 审查数据上界 (Suppressed时为9)
- Deaths_Type: 数据类型标识 (observed/censored/missing)
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml

# ============ 手动指定疾病子目录（相对输入根目录，例："C81-C96"） ============
MANUAL_ICD_GROUP = "C81-C96"
# =====================================================================

# --- 路径配置 (与之前脚本相同) ---
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

if not CONFIG_PATH.exists():
    print(f"ERROR: 未找到配置文件: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

# 获取数据源配置
ds = (cfg.get("data_sources") or {}).get("cdc_wonder") or {}
aamr_base_rel = ds.get("aamr_original")

# 获取输出目录配置
dd = (cfg.get("data_directories") or {})
processed_rel = dd.get("processed")

if not aamr_base_rel or not processed_rel:
    print("ERROR: config.yaml 缺少必要配置项", file=sys.stderr)
    sys.exit(1)

_base_dir = (PROJECT_ROOT / aamr_base_rel).resolve()
_output_dir = (PROJECT_ROOT / processed_rel).resolve()
_input_dir = (_base_dir / MANUAL_ICD_GROUP.strip()).resolve()

def merge_and_clean_data(folder_path):
    """合并并清洗CDC数据，只保留核心列并处理审查数据"""
    print(f"开始处理文件夹: {folder_path}")

    # 仅需这几列即可
    expected_columns = ['Year', 'County Code', 'Deaths', 'Population']
    all_data = []

    for file_path in sorted(folder_path.glob("*.csv")):
        print(f"  处理文件: {file_path.name}")
        
        try:
            # 只读取需要的列
            df = pd.read_csv(file_path, usecols=lambda c: c in expected_columns, encoding='latin1')
            df = df.dropna(how='all').copy()
        except Exception as e:
            print(f"    警告: 读取文件失败: {e}")
            continue

        if df.empty:
            print(f"    跳过空文件: {file_path.name}")
            continue

        # 标准化年份并过滤
        df['Year'] = pd.to_numeric(df['Year'], errors='coerce')
        df = df.dropna(subset=['Year'])
        df['Year'] = df['Year'].astype(int)
        df = df[(df['Year'] >= 1999) & (df['Year'] <= 2020)].copy()

        if df.empty:
            print(f"    文件无有效年份数据: {file_path.name}")
            continue
            
        # 标准化FIPS
        fips_num = pd.to_numeric(df['County Code'], errors='coerce')
        df = df[fips_num.notna()].copy()
        
        if df.empty:
            print(f"    文件无有效FIPS数据: {file_path.name}")
            continue
        
        # 重新获取有效的fips_num（因为df已经被过滤）
        fips_num = pd.to_numeric(df['County Code'], errors='coerce')
        df['COUNTY_FIPS'] = fips_num.astype(int).astype(str).str.zfill(5)

        # --- 核心逻辑：处理死亡人数的审查数据 ---
        
        # 强制转换为字符串以便处理 "Suppressed"
        deaths_str = df['Deaths'].astype(str).str.strip().str.lower()
        
        # 将原始数值转换为数字，非数字（如 "Suppressed"）变为 NaN
        deaths_numeric = pd.to_numeric(df['Deaths'], errors='coerce')
        
        # 1. 创建观测值列 (Deaths_Observed)
        # 只有当原始值是数字时，我们才保留它
        df['Deaths_Observed'] = deaths_numeric.astype('Int64')

        # 2. 创建审查区间列 (Deaths_Censored_Lower, Deaths_Censored_Upper)
        # 仅当原始值为 "Suppressed" 时，定义区间 [1, 9]
        is_suppressed = deaths_str == 'suppressed'
        df['Deaths_Censored_Lower'] = pd.NA
        df['Deaths_Censored_Upper'] = pd.NA
        df.loc[is_suppressed, 'Deaths_Censored_Lower'] = 1
        df.loc[is_suppressed, 'Deaths_Censored_Upper'] = 9
        
        # 转换为整数类型
        df['Deaths_Censored_Lower'] = df['Deaths_Censored_Lower'].astype('Int64')
        df['Deaths_Censored_Upper'] = df['Deaths_Censored_Upper'].astype('Int64')
        
        # 3. 创建数据类型标识列
        df['Deaths_Type'] = 'observed'  # 默认为观测值
        df.loc[is_suppressed, 'Deaths_Type'] = 'censored'  # 审查数据
        df.loc[deaths_numeric.isna() & ~is_suppressed, 'Deaths_Type'] = 'missing'  # 真正缺失

        # 处理人口数据
        df['Population'] = pd.to_numeric(df['Population'], errors='coerce').astype('Int64')
        
        # 保存有效记录数
        valid_count = len(df)
        print(f"    ✅ 成功处理，{valid_count} 行数据")
        
        all_data.append(df)

    if not all_data:
        print("❌ 没有处理任何有效数据。")
        return None

    print("\n合并所有数据...")
    merged_df = pd.concat(all_data, ignore_index=True)

    # 去重：如果同一年同一个县有多个记录，优先保留有确切死亡人数的记录
    print("处理重复记录...")
    merged_df['priority_score'] = merged_df['Deaths_Observed'].notna().astype(int)
    merged_df = merged_df.sort_values(
        ['COUNTY_FIPS', 'Year', 'priority_score'], 
        ascending=[True, True, False]
    )
    
    before_dedup = len(merged_df)
    merged_df = merged_df.drop_duplicates(subset=['COUNTY_FIPS', 'Year'], keep='first')
    after_dedup = len(merged_df)
    
    if before_dedup != after_dedup:
        print(f"  去重完成: {before_dedup} -> {after_dedup} 条记录")

    # 整理最终输出的列
    final_cols = [
        'COUNTY_FIPS', 'Year', 'Population', 
        'Deaths_Observed', 'Deaths_Censored_Lower', 'Deaths_Censored_Upper',
        'Deaths_Type'
    ]
    final_df = merged_df[final_cols].copy()
    
    # 确保所有数值列都是整数类型
    int_columns = ['Year', 'Population', 'Deaths_Observed', 'Deaths_Censored_Lower', 'Deaths_Censored_Upper']
    for col in int_columns:
        if col in final_df.columns:
            final_df[col] = final_df[col].astype('Int64')

    print(f"✅ 数据处理完成，总共 {len(final_df)} 条记录。")
    return final_df

def print_data_summary(df):
    """打印数据摘要报告"""
    print("\n" + "="*60)
    print("数据摘要报告")
    print("="*60)
    
    print(f"数据形状: {df.shape}")
    print(f"年份范围: {df['Year'].min()} - {df['Year'].max()}")
    print(f"县数量: {df['COUNTY_FIPS'].nunique()}")
    print(f"年份×县组合数: {len(df)}")
    
    # 审查数据统计
    print(f"\n数据类型统计:")
    type_counts = df['Deaths_Type'].value_counts()
    for dtype, count in type_counts.items():
        percentage = count / len(df) * 100
        print(f"  {dtype}: {count} ({percentage:.1f}%)")
    
    # 缺失统计
    print(f"\n缺失统计:")
    for col in ['Deaths_Observed', 'Population']:
        missing = df[col].isna().sum()
        percentage = missing / len(df) * 100
        print(f"  {col} 缺失: {missing} ({percentage:.1f}%)")

def main():
    """主函数"""
    if not _input_dir.exists():
        print(f"ERROR: 输入目录不存在: {_input_dir}", file=sys.stderr)
        sys.exit(1)
        
    print(f"开始处理CDC数据: {MANUAL_ICD_GROUP}")
    print(f"输入目录: {_input_dir}")
    
    final_data = merge_and_clean_data(_input_dir)

    if final_data is not None:
        # 确保输出目录存在
        cdc_output_dir = _output_dir / "CDC"
        cdc_output_dir.mkdir(parents=True, exist_ok=True)
        
        output_file = cdc_output_dir / f"AAMR_{MANUAL_ICD_GROUP}.csv"
        final_data.to_csv(output_file, index=False)
        print(f"\n数据已保存到: {output_file}")
        
        # 打印数据摘要
        print_data_summary(final_data)
        
        # 打印数据预览
        print(f"\n数据预览 (前5行):")
        print(final_data.head())
        
        # 审查数据示例
        censored_data = final_data[final_data['Deaths_Type'] == 'censored']
        if len(censored_data) > 0:
            print(f"\n审查数据示例 ({len(censored_data)}条):")
            print(censored_data[['COUNTY_FIPS', 'Year', 'Deaths_Type', 'Deaths_Censored_Lower', 'Deaths_Censored_Upper']].head())
        else:
            print("\n无审查数据")
        
        print(f"\n✅ CDC精简数据清理完成!")

if __name__ == "__main__":
    main()