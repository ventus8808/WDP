#!/usr/bin/env python3
"""
EQI × CDC WONDER — AAMR point estimates (long table)
Outputs a single long CSV with rows per county × time_period × ICD × Lag_Years,
including RUCC/EQI covariates mapped by lagged EQI period.

Output path: config.yaml eqi_aamr_outputs.base_dir + eqi_aamr_outputs.eqi_aamr_point
"""
import sys
from pathlib import Path
import re
import pandas as pd
import numpy as np
import yaml

# ---------------------- Config & Paths ----------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / 'config.yaml'

with CONFIG_PATH.open('r', encoding='utf-8') as f:
    CFG = yaml.safe_load(f)

CDC_EQI_SRC_DIR = PROJECT_ROOT / CFG['data_sources']['cdc_wonder']['eqi_original']
# 使用config中的路径配置
DF_OUTPUT = PROJECT_ROOT / CFG.get('eqi_aamr_outputs', {}).get('base_dir', 'Data/Processed/df_EQI_AAMR') / CFG.get('eqi_aamr_outputs', {}).get('eqi_aamr_point', 'EQI_AAMR_Point.csv')
EQI_DIR = PROJECT_ROOT / CFG['data_sources']['epa_eqi']['processed']
SUPPRESSION_THRESHOLD = float(CFG['data_sources']['cdc_wonder'].get('eqi_suppression_threshold', 40.0))
SMOKING_PATH = PROJECT_ROOT / CFG['data_directories']['processed'] / 'Smoking' / 'County_Smoking.csv'

# ---------------------- Helpers ----------------------

def _extract_state_from_county(county_name: str) -> str:
    if pd.isna(county_name) or not isinstance(county_name, str):
        return ""
    if ", " in county_name:
        return county_name.split(", ")[-1].strip()
    return ""

def _clean_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    cols_norm = {c: c.strip() for c in df.columns}
    df = df.rename(columns=cols_norm)
    if 'Notes' in df.columns:
        df = df.drop(columns=['Notes'])
    df_clean = df.dropna(subset=['County Code']).copy()
    df_clean = df_clean[pd.to_numeric(df_clean['County Code'], errors='coerce').notna()]
    df_clean = df_clean[(df_clean['Deaths'] != 'Missing') & (df_clean['Population'] != 'Missing')]
    return df_clean

def _find_col_by_keywords(df: pd.DataFrame, include_keys, exclude_keys=None):
    include = [k.lower() for k in include_keys]
    exclude = [k.lower() for k in (exclude_keys or [])]
    for col in df.columns:
        s = str(col).lower()
        if all(k in s for k in include) and not any(k in s for k in exclude):
            return col
    return None

def _extract_aamr_point(df: pd.DataFrame) -> pd.Series:
    rate_col = _find_col_by_keywords(df, ['age','adjusted','rate'], ['lower','upper','confidence','interval','ci','standard','error','se'])
    if rate_col is None:
        for alias in ['Age Adjusted Rate','Age-Adjusted Rate','Age adjusted rate','Age-adjusted Rate','Age Adjusted Rate (per 100,000)']:
            if alias in df.columns:
                rate_col = alias
                break
    lcl_col = _find_col_by_keywords(df, ['age','adjusted','rate','lower','confidence'])
    if lcl_col is None:
        for alias in ['Age Adjusted Rate Lower 95% Confidence Interval','Lower 95% CI for Age Adjusted Rate']:
            if alias in df.columns:
                lcl_col = alias
                break
    ucl_col = _find_col_by_keywords(df, ['age','adjusted','rate','upper','confidence'])
    if ucl_col is None:
        for alias in ['Age Adjusted Rate Upper 95% Confidence Interval','Upper 95% CI for Age Adjusted Rate']:
            if alias in df.columns:
                ucl_col = alias
                break
    se_col = _find_col_by_keywords(df, ['age','adjusted','rate','standard','error'])
    if se_col is None:
        for alias in ['Age Adjusted Rate Standard Error','Standard Error for Age Adjusted Rate']:
            if alias in df.columns:
                se_col = alias
                break
    rate = pd.to_numeric(df[rate_col], errors='coerce') if rate_col else pd.Series([pd.NA]*len(df), dtype='Float64')
    lcl = pd.to_numeric(df[lcl_col], errors='coerce') if lcl_col else pd.Series([pd.NA]*len(df), dtype='Float64')
    ucl = pd.to_numeric(df[ucl_col], errors='coerce') if ucl_col else pd.Series([pd.NA]*len(df), dtype='Float64')
    se  = pd.to_numeric(df[se_col], errors='coerce') if se_col else pd.Series([pd.NA]*len(df), dtype='Float64')
    est = rate.copy()
    mask = est.isna() & lcl.notna() & ucl.notna()
    est = est.mask(mask, (lcl+ucl)/2)
    mask = est.isna() & lcl.notna() & se.notna()
    est = est.mask(mask, lcl + 1.96*se)
    mask = est.isna() & ucl.notna() & se.notna()
    est = est.mask(mask, ucl - 1.96*se)
    return est.astype('Float64')

def _process_deaths_series(s: pd.Series) -> pd.Series:
    def conv(v):
        if pd.isna(v):
            return pd.NA
        t = str(v).strip()
        if t in ['Suppressed','Unreliable','']:
            return pd.NA
        if t == '0':
            return 0
        try:
            return int(float(t))
        except Exception:
            return pd.NA
    return s.apply(conv).astype('Int64')

EQI_COLS = ['RUCC','EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
            'RUCC_EQI','RUCC_EQI_Air','RUCC_EQI_Water','RUCC_EQI_Land','RUCC_EQI_Built','RUCC_EQI_Social']

def _extract_meta_from_name(name: str):
    base = name.replace('.csv','')
    m = re.search(r'(\d{4}-\d{4})', base)
    if not m:
        raise ValueError(f'无法解析时间段: {name}')
    period = m.group(1)
    icd = base.replace(period,'').strip()
    return period, icd

def _suppression_rate(df_raw: pd.DataFrame) -> float:
    dfc = _clean_dataframe(df_raw)
    tot = len(dfc)
    sup = (dfc['Deaths'] == 'Suppressed').sum()
    return (sup/tot*100.0) if tot>0 else 0.0

def _load_eqi_by_period() -> dict:
    d = {}
    # 列名映射字典：旧列名 -> 新列名
    column_mapping = {
        'EQI_air': 'EQI_Air',
        'EQI_water': 'EQI_Water', 
        'EQI_land': 'EQI_Land',
        'EQI_built': 'EQI_Built',
        'EQI_Sociodemographic': 'EQI_Social',
        'RUCC_EQI_air': 'RUCC_EQI_Air',
        'RUCC_EQI_water': 'RUCC_EQI_Water',
        'RUCC_EQI_land': 'RUCC_EQI_Land', 
        'RUCC_EQI_built': 'RUCC_EQI_Built',
        'RUCC_EQI_Sociodemographic': 'RUCC_EQI_Social'
    }
    
    fp_0005 = EQI_DIR / 'EQI0005.csv'
    fp_0610 = EQI_DIR / 'EQI0610.csv'
    if fp_0005.exists():
        df = pd.read_csv(fp_0005)
        df['COUNTY_FIPS'] = df['COUNTY_FIPS'].astype(str).str.zfill(5)
        
        # 检测并映射列名（如果需要的话）
        columns_to_rename = {}
        for old_col, new_col in column_mapping.items():
            if old_col in df.columns and new_col not in df.columns:
                columns_to_rename[old_col] = new_col
        
        if columns_to_rename:
            df = df.rename(columns=columns_to_rename)
            print(f"  📝 EQI0005: 映射列名 {list(columns_to_rename.keys())} -> {list(columns_to_rename.values())}")
        
        d['0005'] = df
    if fp_0610.exists():
        df = pd.read_csv(fp_0610)
        df['COUNTY_FIPS'] = df['COUNTY_FIPS'].astype(str).str.zfill(5)
        
        # 检测并映射列名（如果需要的话）
        columns_to_rename = {}
        for old_col, new_col in column_mapping.items():
            if old_col in df.columns and new_col not in df.columns:
                columns_to_rename[old_col] = new_col
        
        if columns_to_rename:
            df = df.rename(columns=columns_to_rename)
            print(f"  📝 EQI0610: 映射列名 {list(columns_to_rename.keys())} -> {list(columns_to_rename.values())}")
        
        d['0610'] = df
    return d

def _load_smoking() -> pd.DataFrame | None:
    if SMOKING_PATH.exists():
        df = pd.read_csv(SMOKING_PATH)
        if 'COUNTY_FIPS' in df.columns and 'SR_Total' in df.columns:
            df = df[['COUNTY_FIPS','SR_Total']].copy()
            df['COUNTY_FIPS'] = df['COUNTY_FIPS'].astype(str).str.zfill(5)
            df = df.rename(columns={'SR_Total':'Smoking_Rate'})
            return df
    return None

def _map_eqi_period(time_period: str, lag_years: int):
    # Based on available EQI datasets, return corresponding EQI year ranges
    # EQI 2000-2005 + 5年滞后 → AAMR 2006-2010
    # EQI 2000-2005 + 10年滞后 → AAMR 2011-2015  
    # EQI 2006-2010 + 5年滞后 → AAMR 2011-2015
    # EQI 2006-2010 + 10年滞后 → AAMR 2016-2020
    if time_period == '2006-2010':
        return '2000-2005' if lag_years == 5 else None
    if time_period == '2011-2015':
        return '2006-2010' if lag_years == 5 else '2000-2005' if lag_years == 10 else None
    if time_period == '2016-2020':
        return '2006-2010' if lag_years == 10 else None
    return None

def _map_eqi_code(time_period: str, lag_years: int):
    # Helper function to get EQI file codes for loading data
    # EQI 0005(2000-2005) + 5年滞后 → AAMR 2006-2010
    # EQI 0005(2000-2005) + 10年滞后 → AAMR 2011-2015  
    # EQI 0610(2006-2010) + 5年滞后 → AAMR 2011-2015
    # EQI 0610(2006-2010) + 10年滞后 → AAMR 2016-2020
    if time_period == '2006-2010':
        return '0005' if lag_years == 5 else None
    if time_period == '2011-2015':
        return '0610' if lag_years == 5 else '0005' if lag_years == 10 else None
    if time_period == '2016-2020':
        return '0610' if lag_years == 10 else None
    return None

# ---------------------- Main build ----------------------

def main():
    print('🧮 Building AAMR point long table with RUCC/EQI...')
    files = sorted(CDC_EQI_SRC_DIR.glob('*.csv'))
    if not files:
        print(f'⚠️ No files under {CDC_EQI_SRC_DIR}')
        sys.exit(0)

    # Discover time periods and suppression-qualified ICDs
    grouped = {}
    for fp in files:
        try:
            period, icd = _extract_meta_from_name(fp.name)
        except Exception:
            continue
        grouped.setdefault(period, []).append((icd, fp))

    icd_rates = {}
    for period, pairs in grouped.items():
        for icd, fp in pairs:
            try:
                try:
                    df = pd.read_csv(fp, encoding='utf-8')
                except UnicodeDecodeError:
                    df = pd.read_csv(fp, encoding='latin-1')
                r = _suppression_rate(df)
                icd_rates.setdefault(icd, {})[period] = r
            except Exception:
                pass

    qualified_icds = [icd for icd, d in icd_rates.items() if len(d)==len(grouped) and all(v<=SUPPRESSION_THRESHOLD for v in d.values())]
    if not qualified_icds:
        print('⚠️ No ICDs pass suppression threshold across all periods. Proceeding with intersection empty → no output.')
        qualified_icds = []

    eqi_dict = _load_eqi_by_period()
    smoking_df = _load_smoking()

    rows = []
    for period in sorted(grouped.keys()):
        pairs = [(icd, fp) for icd, fp in grouped[period] if (not qualified_icds) or (icd in qualified_icds)]
        if not pairs:
            continue
        print(f'📅 {period}: {len(pairs)} ICDs')
        for icd, fp in pairs:
            try:
                try:
                    df = pd.read_csv(fp, encoding='utf-8')
                except UnicodeDecodeError:
                    df = pd.read_csv(fp, encoding='latin-1')
                dfc = _clean_dataframe(df)
                dfc['COUNTY_FIPS'] = dfc['County Code'].astype(int).astype(str).str.zfill(5)
                dfc['State'] = dfc.get('County','').apply(_extract_state_from_county) if 'County' in dfc.columns else ''

                deaths_num = pd.to_numeric(dfc['Deaths'], errors='coerce')
                deaths_txt = dfc['Deaths'].astype(str).str.strip().str.lower()
                zero_mask = deaths_num.eq(0)
                suppressed_mask = deaths_txt.str.contains('suppress', na=False) | deaths_num.between(1,9, inclusive='both')
                unreliable_mask = deaths_txt.str.contains('unreliable', na=False) | deaths_num.between(10,19, inclusive='both')
                # completed_mask = deaths_num.ge(20)

                deaths_processed = _process_deaths_series(dfc['Deaths'])
                pop = dfc['Population'].astype('int64')
                with np.errstate(divide='ignore', invalid='ignore'):
                    crude = (deaths_processed.astype('Float64') / pop.astype('Float64')) * 100000

                est = _extract_aamr_point(dfc)
                est = est.mask(suppressed_mask, pd.NA)
                est = est.mask(unreliable_mask, crude)
                est = est.mask(zero_mask, 0.0)
                est = est.astype('Float64').round(2)

                icd_fmt = icd.replace('-', '_')
                base = pd.DataFrame({
                    'COUNTY_FIPS': dfc['COUNTY_FIPS'].astype(str),
                    'State': dfc['State'].astype(str),
                    'Time_Period': period,
                    'Cancer_Type': icd_fmt,
                    'AAMR': est,
                })
                # Deduplicate
                base = base.drop_duplicates(subset=['COUNTY_FIPS'], keep='first')

                for lag in (5, 10):
                    eqi_period = _map_eqi_period(period, lag)
                    eqi_code = _map_eqi_code(period, lag)
                    out = base.copy()
                    out['Lag_Years'] = lag
                    out['EQI_Period'] = eqi_period if eqi_period is not None else pd.NA
                    if eqi_code and eqi_code in eqi_dict:
                        eqidf = eqi_dict[eqi_code][['COUNTY_FIPS'] + EQI_COLS].copy()
                        out = out.merge(eqidf, on='COUNTY_FIPS', how='left')
                    else:
                        # ensure columns exist with NaN
                        for c in EQI_COLS:
                            out[c] = pd.NA
                    # merge Smoking_Rate
                    if smoking_df is not None:
                        out = out.merge(smoking_df, on='COUNTY_FIPS', how='left')
                    else:
                        out['Smoking_Rate'] = pd.NA
                    rows.append(out)
            except Exception as e:
                print(f'  ⚠️ {fp.name} failed: {e}')
                continue

    if not rows:
        print('⚠️ No rows to write.')
        DF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=['COUNTY_FIPS','State','Time_Period','Lag_Years','EQI_Period','Cancer_Type','AAMR']+EQI_COLS).to_csv(DF_OUTPUT, index=False)
        print(f'💾 Wrote empty skeleton to {DF_OUTPUT}')
        return

    final = pd.concat(rows, ignore_index=True)
    # Order columns
    first_cols = ['COUNTY_FIPS','State','Time_Period','Lag_Years','EQI_Period','Cancer_Type','AAMR','Smoking_Rate']
    ordered = first_cols + [c for c in EQI_COLS if c in final.columns]
    final = final[ordered]

    # Cast RUCC/EQI quintiles to nullable int to avoid 1.0/2.0 formatting
    for c in EQI_COLS:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors='coerce').astype('Int64')

    DF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(DF_OUTPUT, index=False)
    print(f'💾 Saved {len(final):,} rows to {DF_OUTPUT}')

if __name__ == '__main__':
    main()