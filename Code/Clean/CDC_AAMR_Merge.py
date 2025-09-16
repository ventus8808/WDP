#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDC AAMR（年龄调整死亡率）清洗脚本（精简版）
- 读取指定ICD子目录下的所有CSV
- 仅保留列：Year, County, County Code→COUNTY_FIPS, Deaths, Population,
  Crude Rate→CMR, Crude Rate Standard Error→CMR_SE,
  Age Adjusted Rate→AAMR, Age Adjusted Rate Standard Error→AAMR_SE
- 年份限定 1999–2020，同一县同一年多条记录优先保留SE非缺失
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Import YAML configuration
try:
    import yaml  # type: ignore
except Exception:
    print("ERROR: 需要 PyYAML。请安装: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

# ============ 手动指定疾病子目录（相对输入根目录，例："C81-C96"） ============
MANUAL_ICD_GROUP = "C81-C96"
# =====================================================================

# Parse configuration and compute paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

if not CONFIG_PATH.exists():
    print(f"ERROR: 未找到配置文件: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

# 读取路径配置
ds = (cfg.get("data_sources") or {}).get("cdc_wonder") or {}
aamr_base_rel = ds.get("aamr_original")
processed_rel = ds.get("processed")

if not aamr_base_rel or not processed_rel:
    print("ERROR: config.yaml 缺少 data_sources.cdc_wonder.aamr_original 或 processed", file=sys.stderr)
    sys.exit(1)

_base_dir = (PROJECT_ROOT / aamr_base_rel).resolve()
_output_dir = (PROJECT_ROOT / processed_rel).resolve()

if not MANUAL_ICD_GROUP.strip():
    print("ERROR: 请在脚本顶部设置 MANUAL_ICD_GROUP，例如 'C81-C96'", file=sys.stderr)
    sys.exit(1)

_input_dir = (_base_dir / MANUAL_ICD_GROUP.strip()).resolve()

if not _base_dir.exists():
    print(f"ERROR: 输入根目录不存在: {_base_dir}", file=sys.stderr)
    sys.exit(1)

if not _input_dir.exists():
    print(f"ERROR: 疾病目录不存在: {_input_dir}", file=sys.stderr)
    sys.exit(1)

FOLDER_PATH = _input_dir

def to_numeric(series):
    return pd.to_numeric(series, errors='coerce')

def is_valid_row(row):
    year_val = str(row.get('Year', ''))
    if not year_val.isdigit() or len(year_val) != 4:
        return False
    county_code = str(row.get('County Code', ''))
    if len(county_code) < 1:
        return False
    return True

def merge_aamr_data(folder_path):
    """合并AAMR数据（精简列）"""
    print(f"开始合并AAMR文件夹: {folder_path}")

    expected_columns = ['Year','County','County Code','Deaths','Population','Crude Rate','Crude Rate Standard Error','Age Adjusted Rate','Age Adjusted Rate Standard Error']

    all_data = []

    # Process all CSV files
    for file_path in sorted(folder_path.glob("*.csv")):
        filename = file_path.name
        print(f"处理AAMR文件: {filename}")
        df = pd.read_csv(file_path, encoding='latin1')

        # Remove completely empty rows
        df = df.dropna(how='all')

        # Check for required columns
        missing_cols = [col for col in expected_columns if col not in df.columns]
        if missing_cols:
            raise ValueError(f"缺少必需列: {missing_cols}")

        # Keep only expected columns and rename
        df = df[expected_columns].copy()

        # 年份标准化与过滤
        initial_len = len(df)
        df['Year'] = pd.to_numeric(df['Year'], errors='coerce')
        df = df.dropna(subset=['Year']).copy()
        df['Year'] = df['Year'].astype(int)
        df = df[(df['Year'] >= 1999) & (df['Year'] <= 2020)].copy()

        if len(df) == 0:
            print(f"  ⚠️ 过滤后无有效数据")
            continue

        # 去除非县级记录：County必须包含州缩写分隔符", "
        df = df[df['County'].astype(str).str.contains(', ')].copy()

        # FIPS 标准化（仅保留县级，去除全国/总计等无FIPS记录）
        fips_num = pd.to_numeric(df['County Code'], errors='coerce')
        df = df[~fips_num.isna()].copy()
        df['County Code'] = fips_num.astype(int).astype(str).str.zfill(5)

        # 数值化（保持缺失为NaN）
        df['Deaths'] = to_numeric(df['Deaths']).astype('Int64')
        df['Population'] = to_numeric(df['Population']).astype('Int64')
        df['Crude Rate'] = to_numeric(df['Crude Rate'])
        df['Crude Rate Standard Error'] = to_numeric(df['Crude Rate Standard Error'])
        df['Age Adjusted Rate'] = to_numeric(df['Age Adjusted Rate'])
        df['Age Adjusted Rate Standard Error'] = to_numeric(df['Age Adjusted Rate Standard Error'])

        all_data.append(df)
        print(f"  ✅ 成功处理，{len(df)} 行有效数据（移除 {initial_len - len(df)} 行无效年份/空值）")

    if not all_data:
        print("没有成功处理任何文件")
        return None

    # Combine all data
    print("\n合并所有AAMR数据...")
    merged_df = pd.concat(all_data, ignore_index=True)

    # 去重：优先保留SE非缺失，其次Deaths非缺失
    print("处理重叠年份数据...")
    initial_count = len(merged_df)
    merged_df['se_score'] = (~merged_df['Age Adjusted Rate Standard Error'].isna()).astype(int) * 2 + (~merged_df['Crude Rate Standard Error'].isna()).astype(int) * 1
    merged_df['death_score'] = (~merged_df['Deaths'].isna()).astype(int)
    merged_df = merged_df.sort_values(['County Code','Year','se_score','death_score'], ascending=[True,True,False,False])
    merged_df = merged_df.drop_duplicates(subset=['County Code','Year'], keep='first')
    merged_df = merged_df.drop(columns=['se_score','death_score'])

    final_count = len(merged_df)
    print(f"  去重完成: {initial_count} -> {final_count} 条记录")

    # 重命名并整理列
    final_df = merged_df.rename(columns={
        'County Code': 'COUNTY_FIPS',
        'Crude Rate': 'CMR',
        'Crude Rate Standard Error': 'CMR_SE',
        'Age Adjusted Rate': 'AAMR',
        'Age Adjusted Rate Standard Error': 'AAMR_SE'
    }).copy()

    final_df = final_df[['COUNTY_FIPS','Year','County','Deaths','Population','CMR','CMR_SE','AAMR','AAMR_SE']]

    print(f"AAMR数据合并完成，总共 {len(final_df)} 行")

    # Save merged data
    output_dir = _output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"AAMR_{MANUAL_ICD_GROUP}.csv"
    final_df.to_csv(output_file, index=False)
    print(f"AAMR数据已保存到: {output_file}")

    return final_df

def print_data_summary(df):
    """Print comprehensive data summary"""
    if df is None or len(df) == 0:
        print("无数据可显示摘要")
        return

    print(f"\n{'='*60}")
    print(f"AAMR 数据摘要报告")
    print(f"{'='*60}")
    print(f"数据形状: {df.shape}")
    print(f"年份范围: {df['Year'].min()} - {df['Year'].max()}")
    print(f"县数量: {df['COUNTY_FIPS'].nunique()}")
    print(f"年份×县组合数: {len(df)}")
    print(f"\n缺失统计:")
    for col in ['Deaths','Population','CMR','CMR_SE','AAMR','AAMR_SE']:
        missing = df[col].isna().sum()
        print(f"  {col} 缺失: {missing} ({missing/len(df)*100:.1f}%)")

def main():
    """主函数"""
    print(f"开始处理AAMR数据: {MANUAL_ICD_GROUP}")
    print(f"输入目录: {FOLDER_PATH}")

    # Merge AAMR data
    merged_data = merge_aamr_data(FOLDER_PATH)

    if merged_data is not None:
        print_data_summary(merged_data)
        print(f"\n✅ AAMR数据合并成功!")
    else:
        print(f"\n❌ AAMR数据合并失败!")

if __name__ == "__main__":
    main()
