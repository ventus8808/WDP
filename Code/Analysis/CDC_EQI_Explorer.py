#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDC EQI 数据探索性分析脚本

目标:
- 批量读取所有CDC WONDER EQI数据文件
- 统计不同数据类型的分布情况
- 分析数据质量和覆盖程度
- 为后续GLMM建模提供数据洞察

输出:
- 详细的数据质量报告
- 各类数据类型的统计摘要
- ICD分组和时间趋势分析
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import re
from collections import defaultdict, Counter

# ============ 路径配置 ============
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

if not CONFIG_PATH.exists():
    print(f"ERROR: 未找到配置文件: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

# 获取数据源配置
ds = (cfg.get("data_sources") or {}).get("cdc_wonder") or {}
eqi_base_rel = ds.get("eqi_original") or "Data/Original/CDC WONDER EQI"

_base_dir = (PROJECT_ROOT / eqi_base_rel).resolve()

def classify_deaths_data(deaths_value):
    """
    分类deaths数据类型
    """
    if pd.isna(deaths_value):
        return 'missing'
    
    deaths_str = str(deaths_value).strip().lower()
    
    # 检查是否为Suppressed
    if 'suppress' in deaths_str:
        return 'suppressed'
    
    # 检查是否为Unreliable
    if 'unreliable' in deaths_str:
        return 'unreliable'
    
    # 尝试转换为数字
    try:
        deaths_num = float(deaths_value)
        if deaths_num == 0:
            return 'zero'
        elif 1 <= deaths_num <= 9:
            return 'small_count'  # 这些应该被suppress，但如果出现说明数据有问题
        elif 10 <= deaths_num < 20:
            return 'medium_count'
        elif deaths_num >= 20:
            return 'large_count'
        else:
            return 'invalid_number'
    except (ValueError, TypeError):
        return 'non_numeric'

def extract_metadata_from_filename(filename):
    """
    从文件名提取ICD分组和时间周期信息
    """
    # 示例: "2006-2010 C81-C96.csv"
    
    # 提取时间周期
    time_pattern = r'(\d{4}-\d{4})'
    time_match = re.search(time_pattern, filename)
    time_period = time_match.group(1) if time_match else 'Unknown'
    
    # 提取ICD分组 - 更精确的模式
    # 去掉时间部分和.csv后缀，然后提取ICD
    name_without_time = re.sub(r'\d{4}-\d{4}\s+', '', filename)
    name_without_ext = name_without_time.replace('.csv', '')
    
    # 直接使用剩余部分作为ICD分组
    icd_group = name_without_ext.strip()
    
    return time_period, icd_group

def analyze_single_file(file_path):
    """
    分析单个CDC EQI文件
    """
    print(f"  分析文件: {file_path.name}")
    
    # 从文件名提取元数据
    time_period, icd_group = extract_metadata_from_filename(file_path.name)
    
    try:
        # 先读取第一行来确定列名
        with open(file_path, 'r', encoding='latin1') as f:
            first_line = f.readline().strip()
        
        # 检查是否有Notes列，如果有则跳过第一行
        if first_line.startswith('"Notes"'):
            df = pd.read_csv(file_path, encoding='latin1')
            # 删除第一行（Notes行）
            df = df.iloc[1:].reset_index(drop=True)
        else:
            df = pd.read_csv(file_path, encoding='latin1')
        
        if df.empty:
            return {
                'file': file_path.name,
                'time_period': time_period,
                'icd_group': icd_group,
                'total_rows': 0,
                'status': 'empty'
            }
        
        # 标准化列名
        if 'County Code' in df.columns:
            df['FIPS'] = df['County Code']
        elif 'County_Code' in df.columns:
            df['FIPS'] = df['County_Code']

        # 分析Deaths列
        deaths_classification = df['Deaths'].apply(classify_deaths_data)
        deaths_stats = Counter(deaths_classification)

        # 数值化以支持基于规则的计数
        deaths_numeric = pd.to_numeric(df['Deaths'], errors='coerce')
        total_rows = len(df)

        # 基础布尔掩码
        missing_mask = deaths_classification == 'missing'
        zero_mask = deaths_numeric == 0
        suppressed_mask = deaths_classification == 'suppressed'
        small_1_9_mask = (deaths_numeric >= 1) & (deaths_numeric <= 9)
        surpassed_mask = suppressed_mask | small_1_9_mask
        unreliable_text_mask = deaths_classification == 'unreliable'
        unreliable_10_19_mask = (deaths_numeric >= 10) & (deaths_numeric < 20)  # 采用 <20 边界
        unreliable_mask = unreliable_text_mask | unreliable_10_19_mask
        completed_mask = deaths_numeric >= 20

        # 互斥优先级：Missing > Zero > Surpassed > Unreliable > Completed
        avail = ~missing_mask
        zero_only = zero_mask & avail
        avail = avail & (~zero_only)
        surpassed_only = surpassed_mask & avail
        avail = avail & (~surpassed_only)
        unreliable_only = unreliable_mask & avail
        avail = avail & (~unreliable_only)
        completed_only = completed_mask & avail

        missing_count = int(missing_mask.sum())
        zero_count = int(zero_only.sum())
        surpassed_count = int(surpassed_only.sum())
        unreliable_count = int(unreliable_only.sum())
        completed_count = int(completed_only.sum())
        
        # 分析Population列
        pop_missing = df['Population'].isna().sum()
        pop_numeric = pd.to_numeric(df['Population'], errors='coerce').notna().sum()
        
        # 分析FIPS数据
        fips_missing = df['FIPS'].isna().sum() if 'FIPS' in df.columns else len(df)
        fips_valid = 0
        if 'FIPS' in df.columns:
            fips_numeric = pd.to_numeric(df['FIPS'], errors='coerce')
            fips_valid = ((fips_numeric >= 1001) & (fips_numeric <= 56045)).sum()
        
        return {
            'file': file_path.name,
            'time_period': time_period,
            'icd_group': icd_group,
            'total_rows': total_rows,
            'deaths_stats': dict(deaths_stats),
            'missing_count': missing_count,
            'zero_count': zero_count,
            'surpassed_count': surpassed_count,
            'unreliable_count': unreliable_count,
            'completed_count': completed_count,
            'population_missing': pop_missing,
            'population_valid': pop_numeric,
            'fips_missing': fips_missing,
            'fips_valid': fips_valid,
            'status': 'success'
        }
        
    except Exception as e:
        print(f"    ❌ 读取失败: {e}")
        return {
            'file': file_path.name,
            'time_period': time_period,
            'icd_group': icd_group,
            'total_rows': 0,
            'status': f'error: {e}',
            'deaths_stats': {}
        }

def analyze_all_files():
    """
    分析所有CDC EQI文件
    """
    if not _base_dir.exists():
        print(f"ERROR: 目录不存在: {_base_dir}")
        return None
    
    print(f"开始分析CDC EQI数据...")
    print(f"数据目录: {_base_dir}")
    
    # 获取所有CSV文件
    csv_files = list(_base_dir.glob("*.csv"))
    if not csv_files:
        print("❌ 未找到任何CSV文件")
        return None
    
    print(f"发现 {len(csv_files)} 个CSV文件")
    
    # 分析每个文件
    all_results = []
    for file_path in sorted(csv_files):
        result = analyze_single_file(file_path)
        all_results.append(result)
    
    return all_results

def create_clean_summary_table(results):
    """
    创建清晰的汇总表格
    """
    # 准备表格数据
    table_data = []
    
    for result in results:
        if result['status'] == 'success':
            deaths_stats = result.get('deaths_stats', {})
            
            row = {
                'file': result['file'],
                'time_period': result.get('time_period', 'Unknown'),
                'icd_group': result.get('icd_group', 'Unknown'),
                'total_rows': result.get('total_rows', 0),
                # 互斥五类：Missing/Zero/Surpassed/Unreliable/Completed
                'missing': int(result.get('missing_count', deaths_stats.get('missing', 0) or 0)),
                'surpassed': int(result.get('surpassed_count', deaths_stats.get('suppressed', 0) or 0)),
                'zero': int(result.get('zero_count', deaths_stats.get('zero', 0) or 0)),
                'unreliable': int(result.get('unreliable_count', deaths_stats.get('unreliable', 0) or 0)),
                'completed': int(result.get('completed_count', 0)),
                # 可用性（用于报告）：Deaths≥10 = Unreliable + Completed
                'usable': int(result.get('unreliable_count', 0)) + int(result.get('completed_count', 0)),
                'pop_missing': result.get('population_missing', 0),
                'fips_valid': result.get('fips_valid', 0)
            }
            
            # 计算百分比
            total = result.get('total_rows', 1)  # 避免除零
            row['usable_pct'] = round(row['usable'] / total * 100, 1) if total > 0 else 0
            row['suppressed_pct'] = round(row.get('surpassed', 0) / total * 100, 1) if total > 0 else 0
            
            table_data.append(row)
    
    return pd.DataFrame(table_data)

def generate_summary_report(results):
    """
    生成综合分析报告
    """
    print("\n" + "="*80)
    print("CDC EQI 数据质量分析报告")
    print("="*80)
    
    # 创建清晰的汇总表
    df_summary = create_clean_summary_table(results)
    
    if df_summary.empty:
        print("❌ 没有成功处理的数据")
        return
    
    # 基本统计
    print(f"\n📊 基本统计:")
    print(f"  总文件数: {len(results)}")
    print(f"  成功处理: {len(df_summary)}")
    print(f"  总记录数: {df_summary['total_rows'].sum():,}")
    
    # 按时间周期统计
    print(f"\n⏰ 时间周期分布:")
    time_stats = df_summary.groupby('time_period').agg({
        'total_rows': 'sum',
        'usable': 'sum',
        'surpassed': 'sum',
        'usable_pct': 'mean'
    }).round(1)
    
    for period, row in time_stats.iterrows():
        print(f"  {period}: {row['total_rows']:,} 记录, 可用: {row['usable_pct']:.1f}%")
    
    # 按ICD分组统计（排序显示）
    print(f"\n🏷️  ICD分组数据质量排名:")
    icd_stats = df_summary.groupby('icd_group').agg({
        'total_rows': 'sum',
        'usable': 'sum',
        'surpassed': 'sum',
        'missing': 'sum',
        'usable_pct': 'mean'
    }).round(1)
    
    # 按可用数据百分比排序
    icd_stats_sorted = icd_stats.sort_values('usable_pct', ascending=False)
    
    print(f"{'ICD分组':<12} {'总记录':<8} {'可用':<8} {'可用率':<8} {'抑制':<8} {'缺失':<8}")
    print("-" * 60)
    
    for icd, row in icd_stats_sorted.iterrows():
        status_icon = "🟢" if row['usable_pct'] >= 60 else "🟡" if row['usable_pct'] >= 40 else "🔴"
    print(f"{status_icon} {icd:<10} {row['total_rows']:>8.0f} {row['usable']:>8.0f} {row['usable_pct']:>6.1f}% {row['surpassed']:>8.0f} {row['missing']:>8.0f}")
    
    # 总体统计
    total_usable = df_summary['usable'].sum()
    total_records = df_summary['total_rows'].sum()
    overall_usable_pct = total_usable / total_records * 100 if total_records > 0 else 0
    
    print(f"\n🎯 GLMM建模评估:")
    print(f"  总可用数据 (Deaths≥10): {total_usable:,} / {total_records:,} ({overall_usable_pct:.1f}%)")
    
    if overall_usable_pct >= 50:
        print("  ✅ 数据质量良好，适合GLMM建模")
    elif overall_usable_pct >= 30:
        print("  ⚠️  数据质量中等，建议谨慎建模")
    else:
        print("  ❌ 数据质量较差，可能需要合并ICD分组")
    
    # 建议
    print(f"\n💡 建模建议:")
    high_quality = icd_stats_sorted[icd_stats_sorted['usable_pct'] >= 60].index.tolist()
    medium_quality = icd_stats_sorted[(icd_stats_sorted['usable_pct'] >= 40) & (icd_stats_sorted['usable_pct'] < 60)].index.tolist()
    low_quality = icd_stats_sorted[icd_stats_sorted['usable_pct'] < 40].index.tolist()
    
    if high_quality:
        print(f"  🟢 高质量 (≥60%): {', '.join(high_quality)} - 单独建模")
    if medium_quality:
        print(f"  🟡 中等质量 (40-60%): {', '.join(medium_quality)} - 可建模但需注意")
    if low_quality:
        print(f"  🔴 低质量 (<40%): {', '.join(low_quality)} - 建议合并或谨慎处理")
    
    return df_summary

def create_ordered_summary_table(df_summary):
    """
    按指定ICD顺序重新排序汇总表
    """
    # 定义ICD分组的优先级顺序（按您要求的顺序）
    icd_order = [
        'C00-C97',   # 总体
        'C00-C14',   # Lip, oral cavity and pharynx
        'C15-C26',   # Digestive organs
        'C30-C39',   # Respiratory and intrathoracic organs
        'C40-C41',   # Bone and articular cartilage
        'C43-C44',   # Skin
        'C45-C49',   # Mesothelial and soft tissue
        'C50',       # Breast
        'C51-C58',   # Female genital organs
        'C60-C63',   # Male genital organs
        'C64-C68',   # Urinary tract
        'C69-C72',   # Eye, brain and other parts of CNS
        'C73-C75',   # Thyroid and other endocrine glands
        'C76-C80',   # Ill-defined, secondary and unspecified sites
        'C81-C96',   # Lymphoid, hematopoietic and related tissue
        'C97',       # Independent (primary) multiple sites
        # 子分组
        'C18-C21',   # Colorectal Cancer
        'C25',       # Malignant neoplasms of pancreas
        'C34',       # Malignant neoplasms of bronchus and lung
        'C61',       # Malignant neoplasms of prostate
        'C91-C95'    # Leukemia
    ]
    
    # 时间周期顺序
    time_order = ['2006-2010', '2011-2015', '2016-2020']
    
    # 创建排序键
    def get_sort_key(row):
        icd = row['icd_group']
        time = row['time_period']
        
        # 查找ICD在顺序中的位置
        icd_idx = 999  # 默认值，放在最后
        for i, order_icd in enumerate(icd_order):
            if icd == order_icd:
                icd_idx = i
                break
        
        # 查找时间在顺序中的位置
        time_idx = 999
        for i, order_time in enumerate(time_order):
            if time == order_time:
                time_idx = i
                break
        
        return (icd_idx, time_idx)
    
    # 应用排序
    df_sorted = df_summary.copy()
    df_sorted['sort_key'] = df_sorted.apply(get_sort_key, axis=1)
    df_sorted = df_sorted.sort_values('sort_key').drop('sort_key', axis=1)
    
    return df_sorted

def main():
    """
    主函数
    """
    print("🔍 CDC EQI 数据探索性分析")
    print("-" * 50)
    
    # 分析所有文件
    results = analyze_all_files()
    
    if results:
        # 生成报告
        df_summary = generate_summary_report(results)
        
        # 保存清晰的汇总表到CSV
        output_dir = PROJECT_ROOT / "Result" / "Tables"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存清晰格式的汇总表
        if df_summary is not None and not df_summary.empty:
            # 按指定顺序排序
            df_ordered = create_ordered_summary_table(df_summary)
            
            summary_file = output_dir / "CDC_EQI_Clean_Summary.csv"
            
            # 重新构建导出数据，满足用户指定格式
            def fmt_count_pct(count, total):
                pct = 0 if total == 0 else round(count / total * 100)
                return f"{int(count)}({pct}%)"

            export_rows = []
            for _, r in df_ordered.iterrows():
                total = int(r.get('total_rows', 0) or 0)
                export_rows.append({
                    'File': r.get('file', ''),
                    'Time_Period': r.get('time_period', ''),
                    'ICD_Group': r.get('icd_group', ''),
                    'Total_Records': total,
                    'Valid_FIPS': int(r.get('fips_valid', 0) or 0),
                    'Surpassed(%)': fmt_count_pct(int(r.get('surpassed', 0) or 0), total),
                    'Unreliable(%)': fmt_count_pct(int(r.get('unreliable', 0) or 0), total),
                    'Zero(%)': fmt_count_pct(int(r.get('zero', 0) or 0), total),
                    'Missing(%)': fmt_count_pct(int(r.get('missing', 0) or 0), total),
                    'Completed(%)': fmt_count_pct(int(r.get('completed', 0) or 0), total)
                })

            df_clean = pd.DataFrame(export_rows, columns=[
                'File', 'Time_Period', 'ICD_Group', 'Total_Records',
                'Valid_FIPS', 'Surpassed(%)', 'Unreliable(%)', 'Zero(%)', 'Missing(%)', 'Completed(%)'
            ])

            df_clean.to_csv(summary_file, index=False)
            print(f"\n📄 清晰汇总表已保存到: {summary_file}")
        
        print(f"\n✅ 分析完成！")
    else:
        print(f"\n❌ 分析失败，请检查数据目录和文件格式")

if __name__ == "__main__":
    main()