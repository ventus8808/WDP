#!/usr/bin/env python3
"""
EQI × CDC WONDER — AAMR deltas (change table)
Computes changes for TWO lag scenarios:
  - Lag 5: AAMR 2006-2010 → 2011-2015 (EQI 2000-2005 → 2006-2010)
  - Lag 10: AAMR 2011-2015 → 2016-2020 (EQI 2000-2005 → 2006-2010, same EQI change)

Output: Single CSV with Lag column to distinguish 5-year vs 10-year lag effects.
Output path: config.yaml eqi_aamr_outputs.base_dir + eqi_aamr_outputs.eqi_aamr_delta
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / 'config.yaml'

with CONFIG_PATH.open('r', encoding='utf-8') as f:
    CFG = yaml.safe_load(f)

# Define period pairs for TWO lag scenarios
# Lag 5: Compare AAMR periods with 5-year lag
LAG5_AAMR_P1 = '2006-2010'
LAG5_AAMR_P2 = '2011-2015'

# Lag 10: Compare AAMR periods with 10-year lag
LAG10_AAMR_P1 = '2011-2015'
LAG10_AAMR_P2 = '2016-2020'

# Input and output paths
INTERVAL_PATH = PROJECT_ROOT / CFG.get('eqi_aamr_outputs', {}).get('base_dir', 'Data/Processed/df_EQI_AAMR') / CFG.get('eqi_aamr_outputs', {}).get('eqi_aamr_interval', 'EQI_AAMR_Interval.csv')
SMOKING_EQI_PATH = PROJECT_ROOT / 'Data/Processed/Smoking/County_Smoking_EQI.csv'
DELTA_OUTPUT = PROJECT_ROOT / CFG.get('eqi_aamr_outputs', {}).get('base_dir', 'Data/Processed/df_EQI_AAMR') / CFG.get('eqi_aamr_outputs', {}).get('eqi_aamr_delta', 'EQI_AAMR_Delta.csv')

# EQI columns that need delta calculation (excluding RUCC which is constant)
DELTA_EQI_COLS = ['EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
                  'RUCC_EQI','RUCC_EQI_Air','RUCC_EQI_Water','RUCC_EQI_Land','RUCC_EQI_Built','RUCC_EQI_Social']

# ---------------- Helpers ----------------

def _load_interval_data() -> pd.DataFrame:
    """Load the interval data CSV."""
    if not INTERVAL_PATH.exists():
        print(f'⚠️ Input file not found: {INTERVAL_PATH}')
        sys.exit(1)
    df = pd.read_csv(INTERVAL_PATH)
    return df

def _filter_periods(df: pd.DataFrame, period: str) -> pd.DataFrame:
    """Filter dataframe for a specific time period."""
    return df[df['Time_Period'] == period].copy()

def _compute_deltas(df_p1: pd.DataFrame, df_p2: pd.DataFrame, lag_value: int) -> pd.DataFrame:
    """
    Compute delta values between period 1 and period 2 for a specific lag.
    
    Args:
        df_p1: DataFrame for earlier AAMR period
        df_p2: DataFrame for later AAMR period  
        lag_value: Lag years (5 or 10) to filter and label results
    """
    # Filter to specific lag and merge
    df_p1_lag = df_p1[df_p1['Lag_Years'] == lag_value].copy()
    df_p2_lag = df_p2[df_p2['Lag_Years'] == lag_value].copy()
    
    # Merge on COUNTY_FIPS and Cancer_Type only (Lag_Years is already filtered)
    merged = pd.merge(
        df_p1_lag, df_p2_lag, 
        on=['COUNTY_FIPS', 'Cancer_Type'], 
        how='inner', 
        suffixes=('_P1', '_P2')
    )
    
    if merged.empty:
        print(f'⚠️ No matching counties/cancers for Lag={lag_value}.')
        return pd.DataFrame()
    
    # Load smoking data and compute delta
    smoking_df = pd.read_csv(SMOKING_EQI_PATH, dtype={'COUNTY_FIPS': str})
    smoking_df['COUNTY_FIPS'] = smoking_df['COUNTY_FIPS'].str.zfill(5)
    smoking_df['delta_Smoking_Rate'] = smoking_df['0610_SR'] - smoking_df['0005_SR']
    
    # Merge smoking delta
    merged['COUNTY_FIPS'] = merged['COUNTY_FIPS'].astype(str)
    merged = pd.merge(merged, smoking_df[['COUNTY_FIPS', 'delta_Smoking_Rate']], on='COUNTY_FIPS', how='left')
    
    # Base columns
    result = pd.DataFrame({
        'COUNTY_FIPS': merged['COUNTY_FIPS'],
        'State': merged['State_P2'],  # Use P2 state
        'RUCC': merged['RUCC_P2'].astype('Int64'),  # Ensure INT
        'Cancer_Type': merged['Cancer_Type'],
        'Lag': lag_value,  # Constant lag value for this analysis
    })
    
    # Delta calculations for AAMR
    result['delta_AAMR_lower'] = (merged['AAMR_lower_P2'] - merged['AAMR_upper_P1']).round(2)
    result['delta_AAMR_upper'] = (merged['AAMR_upper_P2'] - merged['AAMR_lower_P1']).round(2)
    result['delta_Smoking_Rate'] = merged['delta_Smoking_Rate']
    
    # Delta for EQI quintiles (both P1 and P2 should have same EQI period for same lag)
    for col in DELTA_EQI_COLS:
        delta_col = f'delta_{col}_quintile'
        result[delta_col] = merged[f'{col}_P2'] - merged[f'{col}_P1']
    
    return result

def _add_categories(df: pd.DataFrame) -> pd.DataFrame:
    """Add category columns based on delta values."""
    def _categorize(delta):
        if pd.isna(delta):
            return pd.NA
        elif delta < 0:
            return 'Improved'
        elif delta == 0:
            return 'Stable'
        else:
            return 'Worsened'
    
    # Category for EQI
    df['EQI_Change_Category'] = df['delta_EQI_quintile'].apply(_categorize)
    
    # Categories for other dimensions
    for col in DELTA_EQI_COLS[1:]:  # Skip EQI itself
        delta_col = f'delta_{col}_quintile'
        if col.startswith('RUCC_EQI_'):
            cat_col = f'{col}_Change_Category'
        else:
            cat_col = f'{col.replace("EQI_", "")}_Change_Category'
        df[cat_col] = df[delta_col].apply(_categorize)
    
    return df

def _validate_output(df: pd.DataFrame) -> None:
    """Validate output data and print summary."""
    if df.empty:
        print('⚠️ Output dataframe is empty.')
        return
    
    total_rows = len(df)
    missing_smoking = df['delta_Smoking_Rate'].isna().sum() / total_rows * 100
    
    # Category distribution
    eqi_cat_dist = df['EQI_Change_Category'].value_counts(normalize=True, dropna=False) * 100
    
    print(f'✅ EQI_AAMR_Delta.csv successfully created.')
    print(f'   Total rows: {total_rows}')
    print(f'   Missing delta_Smoking_Rate rate: {missing_smoking:.1f}%')
    print(f'   Distribution of EQI_Change_Category:')
    for cat, pct in eqi_cat_dist.items():
        print(f'     - {cat}: {pct:.1f}%')

# ---------------- Main ----------------

def main():
    print('🧮 Computing AAMR deltas for TWO lag scenarios...')
    print('='*70)
    
    # Load interval data
    df = _load_interval_data()
    print(f'📊 Loaded {len(df):,} rows from {INTERVAL_PATH}')
    
    all_deltas = []
    
    # =========================================================================
    # LAG 5 ANALYSIS
    # =========================================================================
    print('\n' + '='*70)
    print('🔬 LAG 5 ANALYSIS')
    print('='*70)
    print(f'EQI Change: 2000-2005 → 2006-2010')
    print(f'AAMR Change: {LAG5_AAMR_P1} → {LAG5_AAMR_P2}')
    
    df_lag5_p1 = _filter_periods(df, LAG5_AAMR_P1)
    df_lag5_p2 = _filter_periods(df, LAG5_AAMR_P2)
    print(f'📅 Period 1 ({LAG5_AAMR_P1}): {len(df_lag5_p1)} rows')
    print(f'📅 Period 2 ({LAG5_AAMR_P2}): {len(df_lag5_p2)} rows')
    
    delta_lag5 = _compute_deltas(df_lag5_p1, df_lag5_p2, lag_value=5)
    
    if not delta_lag5.empty:
        delta_lag5 = _add_categories(delta_lag5)
        delta_lag5['delta_Smoking_Rate'] = delta_lag5['delta_Smoking_Rate'].fillna(0)
        
        initial_rows = len(delta_lag5)
        delta_lag5 = delta_lag5.dropna()
        dropped_rows = initial_rows - len(delta_lag5)
        if dropped_rows > 0:
            print(f'🧹 Dropped {dropped_rows} rows with missing values.')
        
        all_deltas.append(delta_lag5)
        print(f'✅ Lag 5 delta: {len(delta_lag5):,} complete observations')
    else:
        print('⚠️ No Lag 5 deltas computed.')
    
    # =========================================================================
    # LAG 10 ANALYSIS
    # =========================================================================
    print('\n' + '='*70)
    print('🔬 LAG 10 ANALYSIS')
    print('='*70)
    print(f'EQI Change: 2000-2005 → 2006-2010 (same as Lag 5)')
    print(f'AAMR Change: {LAG10_AAMR_P1} → {LAG10_AAMR_P2}')
    
    df_lag10_p1 = _filter_periods(df, LAG10_AAMR_P1)
    df_lag10_p2 = _filter_periods(df, LAG10_AAMR_P2)
    print(f'📅 Period 1 ({LAG10_AAMR_P1}): {len(df_lag10_p1)} rows')
    print(f'📅 Period 2 ({LAG10_AAMR_P2}): {len(df_lag10_p2)} rows')
    
    delta_lag10 = _compute_deltas(df_lag10_p1, df_lag10_p2, lag_value=10)
    
    if not delta_lag10.empty:
        delta_lag10 = _add_categories(delta_lag10)
        delta_lag10['delta_Smoking_Rate'] = delta_lag10['delta_Smoking_Rate'].fillna(0)
        
        initial_rows = len(delta_lag10)
        delta_lag10 = delta_lag10.dropna()
        dropped_rows = initial_rows - len(delta_lag10)
        if dropped_rows > 0:
            print(f'🧹 Dropped {dropped_rows} rows with missing values.')
        
        all_deltas.append(delta_lag10)
        print(f'✅ Lag 10 delta: {len(delta_lag10):,} complete observations')
    else:
        print('⚠️ No Lag 10 deltas computed.')
    
    # =========================================================================
    # COMBINE AND SAVE
    # =========================================================================
    if not all_deltas:
        print('\n⚠️ No deltas to save.')
        return
    
    print('\n' + '='*70)
    print('💾 COMBINING AND SAVING RESULTS')
    print('='*70)
    
    # Concatenate both lag analyses
    final_df = pd.concat(all_deltas, ignore_index=True)
    
    # Validate and print summary
    print(f'\n📊 Final combined dataset:')
    print(f'   Total rows: {len(final_df):,}')
    print(f'   Lag 5 rows: {len(final_df[final_df["Lag"] == 5]):,}')
    print(f'   Lag 10 rows: {len(final_df[final_df["Lag"] == 10]):,}')
    
    # Print category distribution by lag
    for lag in [5, 10]:
        lag_subset = final_df[final_df['Lag'] == lag]
        if len(lag_subset) > 0:
            eqi_cat_dist = lag_subset['EQI_Change_Category'].value_counts(normalize=True) * 100
            print(f'\n   Lag {lag} - EQI_Change_Category distribution:')
            for cat, pct in eqi_cat_dist.items():
                print(f'     - {cat}: {pct:.1f}%')
    
    # Save combined output
    DELTA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(DELTA_OUTPUT, index=False)
    print(f'\n✅ Successfully saved combined delta analysis to:')
    print(f'   {DELTA_OUTPUT}')
    print('='*70)

if __name__ == '__main__':
    main()