#!/usr/bin/env python3
"""
CDC EQI Data Merge - 最终版本（扩展 AAMR 输出）
================================================

改进点：
- 删除 Notes 注释列
- 新增 AAMR 点估计解析/反推与规则化输出（Zero/Surpassed/Unreliable/Completed）
- 每个时间段输出两个文件：
    - CDC_EQI_Death_YYYY_YYYY.csv（原 Death 合并表，改名）
    - CDC_EQI_AAMR_YYYY_YYYY.csv（仅 AAMR 点估计列 AAMR_{ICD}）
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import re
from datetime import datetime

# 配置
SUPPRESSION_THRESHOLD = 40.0

def load_config():
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / 'config.yaml'
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def extract_metadata_from_filename(filename: str):
    base_name = filename.replace('.csv', '')
    time_match = re.search(r'(\d{4}-\d{4})', base_name)
    if not time_match:
        raise ValueError(f"无法从文件名提取时间周期: {filename}")
    time_period = time_match.group(1)
    icd_group = base_name.replace(time_period, '').strip()
    return time_period, icd_group

def _extract_state_from_county(county_name: str) -> str:
    """Extract state abbreviation from county name like 'Autauga County, AL'"""
    if pd.isna(county_name) or not isinstance(county_name, str):
        return ""
    if ", " in county_name:
        return county_name.split(", ")[-1].strip()
    return ""

def clean_dataframe(df):
    """清理数据，移除注释行和无效数据"""
    # 删除 Notes 注释列（如存在）
    cols_norm = {c: c.strip() for c in df.columns}
    df = df.rename(columns=cols_norm)
    if 'Notes' in df.columns:
        df = df.drop(columns=['Notes'])
    
    # 移除County Code为空或nan的行（这些通常是注释）
    df_clean = df.dropna(subset=['County Code']).copy()
    
    # 确保County Code是数字格式
    df_clean = df_clean[pd.to_numeric(df_clean['County Code'], errors='coerce').notna()]
    
    # 移除包含Missing值的行
    df_clean = df_clean[
        (df_clean['Deaths'] != 'Missing') & 
        (df_clean['Population'] != 'Missing')
    ]
    
    return df_clean

def _find_col_by_keywords(df: pd.DataFrame, include_keys, exclude_keys=None):
    """按关键词近似匹配列名（不区分大小写）"""
    include = [k.lower() for k in include_keys]
    exclude = [k.lower() for k in (exclude_keys or [])]
    for col in df.columns:
        s = str(col).lower()
        if all(k in s for k in include) and not any(k in s for k in exclude):
            return col
    return None

def extract_aamr_components(df: pd.DataFrame):
    """提取/反推 AAMR 点估计所需组件（点估计、LCL、UCL、SE）"""
    # 直接点估计列：尽量排除区间/SE
    rate_col = _find_col_by_keywords(
        df,
        include_keys=['age', 'adjusted', 'rate'],
        exclude_keys=['lower', 'upper', 'confidence', 'interval', 'ci', 'standard', 'error', 'se']
    )
    # 兼容常见写法
    if rate_col is None:
        for alias in [
            'Age Adjusted Rate', 'Age-Adjusted Rate', 'Age adjusted rate',
            'Age-adjusted Rate', 'Age Adjusted Rate (per 100,000)'
        ]:
            if alias in df.columns:
                rate_col = alias
                break

    lcl_col = _find_col_by_keywords(df, ['age', 'adjusted', 'rate', 'lower', 'confidence'])
    if lcl_col is None:
        for alias in [
            'Age Adjusted Rate Lower 95% Confidence Interval',
            'Lower 95% CI for Age Adjusted Rate'
        ]:
            if alias in df.columns:
                lcl_col = alias
                break

    ucl_col = _find_col_by_keywords(df, ['age', 'adjusted', 'rate', 'upper', 'confidence'])
    if ucl_col is None:
        for alias in [
            'Age Adjusted Rate Upper 95% Confidence Interval',
            'Upper 95% CI for Age Adjusted Rate'
        ]:
            if alias in df.columns:
                ucl_col = alias
                break

    se_col = _find_col_by_keywords(df, ['age', 'adjusted', 'rate', 'standard', 'error'])
    if se_col is None:
        for alias in [
            'Age Adjusted Rate Standard Error',
            'Standard Error for Age Adjusted Rate'
        ]:
            if alias in df.columns:
                se_col = alias
                break

    rate = pd.to_numeric(df[rate_col], errors='coerce') if rate_col else pd.Series([pd.NA] * len(df), dtype='Float64')
    lcl = pd.to_numeric(df[lcl_col], errors='coerce') if lcl_col else pd.Series([pd.NA] * len(df), dtype='Float64')
    ucl = pd.to_numeric(df[ucl_col], errors='coerce') if ucl_col else pd.Series([pd.NA] * len(df), dtype='Float64')
    se = pd.to_numeric(df[se_col], errors='coerce') if se_col else pd.Series([pd.NA] * len(df), dtype='Float64')

    # 反推点估计
    estimate = rate.copy()
    # 优先 (L+U)/2
    mask_need = estimate.isna() & lcl.notna() & ucl.notna()
    estimate = estimate.mask(mask_need, (lcl + ucl) / 2)
    # 次选 L + 1.96*SE
    mask_need = estimate.isna() & lcl.notna() & se.notna()
    estimate = estimate.mask(mask_need, lcl + 1.96 * se)
    # 再选 U - 1.96*SE
    mask_need = estimate.isna() & ucl.notna() & se.notna()
    estimate = estimate.mask(mask_need, ucl - 1.96 * se)

    return estimate.astype('Float64')

def analyze_suppression_rate(df: pd.DataFrame):
    df_clean = clean_dataframe(df)
    total_records = len(df_clean)
    suppressed_count = (df_clean['Deaths'] == 'Suppressed').sum()
    return (suppressed_count / total_records) * 100 if total_records > 0 else 0

def process_deaths_column(deaths_series: pd.Series):
    def convert_death_value(value):
        if pd.isna(value):
            return pd.NA  # 使用pandas的NA
        
        value_str = str(value).strip()
        if value_str in ['Suppressed', 'Unreliable', '']:
            return pd.NA  # Suppressed/Unreliable转为NA，Missing已在清理阶段删除
        elif value_str == '0':
            return 0
        else:
            try:
                return int(float(value_str))  # 先转float再转int，处理可能的浮点格式
            except (ValueError, TypeError):
                return pd.NA
    
    # 应用转换并确保是nullable integer类型
    result = deaths_series.apply(convert_death_value)
    return result.astype('Int64')  # pandas的nullable integer类型

def process_single_file(file_path, icd_group):
    """处理单个文件，返回 Death 与 AAMR 两套数据（按互斥口径处理 AAMR）"""
    # 读取数据
    try:
        df = pd.read_csv(file_path, encoding='utf-8')
    except UnicodeDecodeError:
        df = pd.read_csv(file_path, encoding='latin-1')
    
    # 清理数据
    df_clean = clean_dataframe(df)
    
    # 标准化FIPS - 确保是字符串格式
    df_clean['COUNTY_FIPS'] = df_clean['County Code'].astype(int).astype(str).str.zfill(5)
    
    # 提取State列
    df_clean['State'] = df_clean['County'].apply(_extract_state_from_county) if 'County' in df_clean.columns else ""
    
    # 处理Deaths列（整数类型），同时构建互斥口径掩码
    deaths_processed = process_deaths_column(df_clean['Deaths'])
    deaths_numeric = pd.to_numeric(df_clean['Deaths'], errors='coerce')
    deaths_text = df_clean['Deaths'].astype(str).str.strip().str.lower()

    zero_mask = deaths_numeric.eq(0)
    suppressed_mask = deaths_text.str.contains('suppress', na=False) | deaths_numeric.between(1, 9, inclusive='both')
    unreliable_mask = deaths_text.str.contains('unreliable', na=False) | deaths_numeric.between(10, 19, inclusive='both')
    completed_mask = deaths_numeric.ge(20)
    # 缺失：无法转为数字且非上述文本（清洗阶段已去除'Missing'字符串）
    missing_mask = deaths_numeric.isna() & ~suppressed_mask & ~unreliable_mask & ~zero_mask
    
    # 处理Population - 确保是整数类型，Missing已在清理阶段删除
    population_processed = df_clean['Population'].astype('int64')

    # AAMR 原始/反推点估计（未套口径）
    aamr_estimate = extract_aamr_components(df_clean)
    # 粗死亡率（每10万），避免除零
    with np.errstate(divide='ignore', invalid='ignore'):
        crude_rate = (deaths_processed.astype('Float64') / population_processed.astype('Float64')) * 100000
    # 按互斥规则设置 AAMR 点估计
    aamr_final = aamr_estimate.copy()
    aamr_final = aamr_final.mask(suppressed_mask, pd.NA)            # Surpassed → NA
    aamr_final = aamr_final.mask(unreliable_mask, crude_rate)        # Unreliable → 粗率
    aamr_final = aamr_final.mask(zero_mask, 0.0)                     # Zero → 0
    # Completed → 保留 aamr_estimate；Missing → NA（保持原值）
    # 统一保留两位小数显示，但内部保留浮点
    aamr_final = aamr_final.astype('Float64')
    
    # 转换ICD代码格式：将"-"替换为"_"
    icd_formatted = icd_group.replace('-', '_')
    
    # 返回标准格式（Deaths 合并行）
    result_death = pd.DataFrame({
        'COUNTY_FIPS': df_clean['COUNTY_FIPS'].astype(str),  # 确保COUNTY_FIPS是字符串
        'State': df_clean['State'].astype(str),
        f'Deaths_{icd_formatted}': deaths_processed,
        'Population': population_processed
    })
    
    # 移除重复COUNTY_FIPS（保留第一个）
    result_death = result_death.drop_duplicates(subset=['COUNTY_FIPS'], keep='first')

    # AAMR 合并行（仅点估计，不含区间/SE）
    result_aamr = pd.DataFrame({
        'COUNTY_FIPS': df_clean['COUNTY_FIPS'].astype(str),
        'State': df_clean['State'].astype(str),
        f'AAMR_{icd_formatted}': aamr_final.round(2)
    })
    result_aamr = result_aamr.drop_duplicates(subset=['COUNTY_FIPS'], keep='first')
    
    return result_death, result_aamr

def reorder_columns(df, qualified_icd_groups):
    """按指定顺序重新排列列"""
    
    # 定义ICD分组的优先顺序（按您之前发给我的顺序）
    icd_order = [
        'C00_C97',   # 所有恶性肿瘤
        'C00_C14',   # 口腔咽部
        'C15_C26',   # 消化系统
        'C18_C21',   # 结直肠癌
        'C25',       # 胰腺癌
        'C30_C39',   # 呼吸系统
        'C34',       # 肺癌
        'C40_C41',   # 骨肉瘤
        'C43_C44',   # 皮肤癌
        'C45_C49',   # 结缔组织
        'C50',       # 乳腺癌
        'C51_C58',   # 女性生殖器
        'C60_C63',   # 男性生殖器
        'C61',       # 前列腺癌
        'C64_C68',   # 泌尿系统
        'C69_C72',   # 眼脑神经
        'C73_C75',   # 内分泌腺
        'C76_C80',   # 未定部位
        'C81_C96',   # 淋巴造血系统
        'C91_C95',   # 白血病
        'C97'        # 多发性肿瘤
    ]
    
    # 构建最终列顺序
    columns_ordered = ['COUNTY_FIPS', 'State', 'Population']
    
    # 按指定顺序添加Deaths列
    qualified_formatted = [icd.replace('-', '_') for icd in qualified_icd_groups]
    
    for icd in icd_order:
        deaths_col = f'Deaths_{icd}'
        if deaths_col in df.columns:
            columns_ordered.append(deaths_col)
    
    # 添加任何遗漏的Deaths列
    for col in df.columns:
        if col.startswith('Deaths_') and col not in columns_ordered:
            columns_ordered.append(col)
    
    return df[columns_ordered]

def reorder_aamr_columns(df, qualified_icd_groups):
    """AAMR 列重排：COUNTY_FIPS + State + AAMR_* 按既定 ICD 顺序"""
    icd_order = [
        'C00_C97','C00_C14','C15_C26','C18_C21','C25','C30_C39','C34','C40_C41','C43_C44',
        'C45_C49','C50','C51_C58','C60_C63','C61','C64_C68','C69_C72','C73_C75','C76_C80',
        'C81_C96','C91_C95','C97'
    ]
    columns_ordered = ['COUNTY_FIPS', 'State']
    for icd in icd_order:
        col = f'AAMR_{icd}'
        if col in df.columns:
            columns_ordered.append(col)
    # 附加遗漏的 AAMR 列
    for col in df.columns:
        if col.startswith('AAMR_') and col not in columns_ordered:
            columns_ordered.append(col)
    return df[columns_ordered]

def main():
    print("🧬 CDC EQI 数据合并处理")
    print("=" * 40)
    
    # 设置路径
    config = load_config()
    project_root = Path(__file__).resolve().parents[2]
    data_dir = project_root / "Data/Original/CDC WONDER EQI"
    output_dir = Path(config['data_directories']['processed']) / 'CDC'
    
    print(f"📂 数据源: {data_dir}")
    print(f"📂 输出: {output_dir}")
    
    # 获取所有CSV文件并按时间周期分组
    csv_files = list(data_dir.glob("*.csv"))
    
    # 按时间周期分组
    time_periods = {}
    for file_path in csv_files:
        try:
            time_period, icd_group = extract_metadata_from_filename(file_path.name)
            if time_period not in time_periods:
                time_periods[time_period] = []
            time_periods[time_period].append((icd_group, file_path))
        except Exception as e:
            print(f"⚠️  跳过文件 {file_path.name}: {e}")
    
    print(f"\n🗓️  发现时间周期: {list(time_periods.keys())}")
    
    # 分析所有时间周期的抑制率，找出一致符合条件的ICD分组
    print(f"\n🔍 筛选符合条件的ICD分组...")
    
    icd_suppression_rates = {}  # {icd_group: {time_period: rate}}
    
    for time_period, period_files in time_periods.items():
        for icd_group, file_path in period_files:
            try:
                # 读取并分析
                try:
                    df = pd.read_csv(file_path, encoding='utf-8')
                except UnicodeDecodeError:
                    df = pd.read_csv(file_path, encoding='latin-1')
                
                suppression_rate = analyze_suppression_rate(df)
                
                if icd_group not in icd_suppression_rates:
                    icd_suppression_rates[icd_group] = {}
                icd_suppression_rates[icd_group][time_period] = suppression_rate
                
            except Exception:
                continue
    
    # 筛选在所有时间周期都符合条件的ICD分组
    qualified_icd_groups = []
    
    for icd_group, rates_by_period in icd_suppression_rates.items():
        # 检查是否在所有时间周期都有数据且都符合条件
        if len(rates_by_period) == len(time_periods):  # 在所有时间周期都存在
            all_qualify = all(rate <= SUPPRESSION_THRESHOLD for rate in rates_by_period.values())
            
            if all_qualify:
                qualified_icd_groups.append(icd_group)
    
    if not qualified_icd_groups:
        print("❌ 没有ICD分组在所有时间周期都符合条件！")
        return
    
    # 处理每个时间周期
    for time_period in sorted(time_periods.keys()):
        print(f"\n" + "="*50)
        print(f"📅 处理时间周期: {time_period}")
        print("="*50)
        
        period_files = time_periods[time_period]
        
        # 只处理符合条件的ICD分组
        print(f"🔍 筛选符合条件的文件:")
        qualified_files = []
        
        for icd_group, file_path in period_files:
            if icd_group in qualified_icd_groups:
                qualified_files.append((icd_group, file_path))
                print(f"  ✅ {icd_group}")
            else:
                print(f"  ❌ {icd_group} (不符合全周期条件)")
        
        print(f"\n🎯 {len(qualified_files)}/{len(period_files)} 个ICD分组符合条件")
        
        if not qualified_files:
            print("❌ 没有符合条件的数据，跳过此时间周期")
            continue
        
        # 合并数据
        print(f"\n📊 合并数据:")
        
        merged_death_df = None
        merged_aamr_df = None
        for i, (icd_group, file_path) in enumerate(qualified_files):
            print(f"  处理 {i+1}/{len(qualified_files)}: {icd_group}")
            
            try:
                current_death_df, current_aamr_df = process_single_file(file_path, icd_group)
                
                # 合并 Death 表
                if merged_death_df is None:
                    merged_death_df = current_death_df.copy()
                    print(f"    🔵 基础死亡数据: {len(merged_death_df):,} 记录")
                else:
                    deaths_col = [col for col in current_death_df.columns if col.startswith('Deaths_')][0]
                    before_count = len(merged_death_df)
                    merged_death_df = merged_death_df.merge(
                        current_death_df[['COUNTY_FIPS', 'State', deaths_col]],
                        on=['COUNTY_FIPS', 'State'], how='outer'
                    )
                    print(f"    ➕ 死亡合并后: {len(merged_death_df):,} 记录 (增加 {len(merged_death_df)-before_count:,})")

                # 合并 AAMR 表
                if merged_aamr_df is None:
                    merged_aamr_df = current_aamr_df.copy()
                    print(f"    🔵 基础AAMR数据: {len(merged_aamr_df):,} 记录")
                else:
                    aamr_col = [col for col in current_aamr_df.columns if col.startswith('AAMR_')][0]
                    before_count = len(merged_aamr_df)
                    merged_aamr_df = merged_aamr_df.merge(
                        current_aamr_df[['COUNTY_FIPS', 'State', aamr_col]],
                        on=['COUNTY_FIPS', 'State'], how='outer'
                    )
                    print(f"    ➕ AAMR合并后: {len(merged_aamr_df):,} 记录 (增加 {len(merged_aamr_df)-before_count:,})")
                    
            except Exception as e:
                print(f"    ⚠️  处理失败: {e}")
        
        if merged_death_df is None or merged_aamr_df is None:
            print("❌ 合并失败")
            continue
        
        # 重新排序列
        print(f"\n🔄 重新排序列...")
        merged_death_df = reorder_columns(merged_death_df, qualified_icd_groups)
        merged_aamr_df = reorder_aamr_columns(merged_aamr_df, qualified_icd_groups)
        
        # 保存数据
        print(f"\n💾 保存 {time_period} 数据:")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_filename_death = f"CDC_EQI_Death_{time_period.replace('-', '_')}.csv"
        output_filename_aamr = f"CDC_EQI_AAMR_{time_period.replace('-', '_')}.csv"
        (output_dir / output_filename_death).write_text('', encoding='utf-8') if False else None
        merged_death_df.to_csv(output_dir / output_filename_death, index=False)
        merged_aamr_df.to_csv(output_dir / output_filename_aamr, index=False)
        
        print(f"  📁 {output_filename_death}")
        print(f"  � {output_filename_aamr}")
        print(f"  🧮 Death: {len(merged_death_df):,} 记录 × {len(merged_death_df.columns)} 列")
        print(f"  🧮 AAMR : {len(merged_aamr_df):,} 记录 × {len(merged_aamr_df.columns)} 列")
        

    
    print(f"\n" + "="*70)
    print("✅ 所有时间周期处理完成！")
    print("🎯 数据已准备就绪，可用于GLMM分析")
    print("="*70)

if __name__ == "__main__":
    main()