#!/usr/bin/env python3
"""
EQI × CDC WONDER — Deaths (long table)
Outputs a long CSV with rows per county × time_period × ICD × Lag_Years,
including Population, Deaths (nullable int), and RUCC/EQI covariates.

Output path: config.yaml df_outputs.eqi_deaths (default Data/df/EQI_Deaths.csv)
"""
import sys
from pathlib import Path
import re
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / 'config.yaml'

with CONFIG_PATH.open('r', encoding='utf-8') as f:
    CFG = yaml.safe_load(f)

CDC_EQI_SRC_DIR = PROJECT_ROOT / CFG['data_sources']['cdc_wonder']['eqi_original']
DF_OUTPUT = PROJECT_ROOT / CFG.get('df_outputs', {}).get('eqi_deaths', 'Data/df/EQI_Deaths.csv')
EQI_DIR = PROJECT_ROOT / CFG['data_sources']['epa_eqi']['processed']
SUPPRESSION_THRESHOLD = float(CFG['data_sources']['cdc_wonder'].get('eqi_suppression_threshold', 40.0))
SMOKING_PATH = PROJECT_ROOT / CFG['data_directories']['processed'] / 'Smoking' / 'County_Smoking.csv'

EQI_COLS = ['RUCC','EQI','EQI_air','EQI_water','EQI_land','EQI_built','EQI_Sociodemographic',
            'RUCC_EQI','RUCC_EQI_air','RUCC_EQI_water','RUCC_EQI_land','RUCC_EQI_built','RUCC_EQI_Sociodemographic']

# --------------- Helpers ---------------

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
    for code in ('0005','0610'):
        fp = EQI_DIR / f'EQI{code}.csv'
        if fp.exists():
            t = pd.read_csv(fp)
            t['COUNTY_FIPS'] = t['COUNTY_FIPS'].astype(str).str.zfill(5)
            d[code] = t
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
    if time_period == '2006-2010':
        return '0005' if lag_years == 5 else None
    if time_period == '2011-2015':
        return '0610' if lag_years == 5 else '0005'
    if time_period == '2016-2020':
        return None if lag_years == 5 else '0610'
    return None

# --------------- Main ---------------

def main():
    print('🧮 Building Deaths long table with RUCC/EQI...')
    files = sorted(CDC_EQI_SRC_DIR.glob('*.csv'))
    if not files:
        print(f'⚠️ No files under {CDC_EQI_SRC_DIR}')
        sys.exit(0)

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
                icd_rates.setdefault(icd, {})[period] = _suppression_rate(df)
            except Exception:
                pass

    qualified_icds = [icd for icd, d in icd_rates.items() if len(d)==len(grouped) and all(v<=SUPPRESSION_THRESHOLD for v in d.values())]
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

                deaths = _process_deaths_series(dfc['Deaths'])
                pop = pd.to_numeric(dfc['Population'], errors='coerce').astype('Int64')

                icd_fmt = icd.replace('-', '_')
                base = pd.DataFrame({
                    'COUNTY_FIPS': dfc['COUNTY_FIPS'].astype(str),
                    'State': dfc['State'].astype(str),
                    'Time_Period': period,
                    'Cancer_Type': icd_fmt,
                    'Population': pop,
                    'Deaths': deaths,
                })
                base = base.drop_duplicates(subset=['COUNTY_FIPS'], keep='first')

                for lag in (5, 10):
                    eqi_code = _map_eqi_period(period, lag)
                    out = base.copy()
                    out['Lag_Years'] = lag
                    out['EQI_Period'] = eqi_code if eqi_code is not None else pd.NA
                    if eqi_code and eqi_code in eqi_dict:
                        eqidf = eqi_dict[eqi_code][['COUNTY_FIPS'] + EQI_COLS].copy()
                        out = out.merge(eqidf, on='COUNTY_FIPS', how='left')
                    else:
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
        pd.DataFrame(columns=['COUNTY_FIPS','State','Time_Period','Lag_Years','EQI_Period','Cancer_Type','Population','Deaths']+EQI_COLS).to_csv(DF_OUTPUT, index=False)
        print(f'💾 Wrote empty skeleton to {DF_OUTPUT}')
        return

    final = pd.concat(rows, ignore_index=True)
    first = ['COUNTY_FIPS','State','Time_Period','Lag_Years','EQI_Period','Cancer_Type','Population','Deaths','Smoking_Rate']
    final = final[first + [c for c in EQI_COLS if c in final.columns]]

    # Cast RUCC/EQI quintiles to nullable int to avoid 1.0/2.0 formatting
    for c in EQI_COLS:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors='coerce').astype('Int64')

    DF_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(DF_OUTPUT, index=False)
    print(f'💾 Saved {len(final):,} rows to {DF_OUTPUT}')

if __name__ == '__main__':
    main()
