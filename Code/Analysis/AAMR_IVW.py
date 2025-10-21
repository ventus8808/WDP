import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import yaml

def calculate_pooled_aamr_meta_analysis(df, group_col=None):
    """
    使用固定效应meta-analysis方法
    这是最标准的合并率的方法
    """
    
    def pool_group(data):
        data = data.copy()
        
        # 权重 = 1/SE^2
        data['weight'] = 1 / (data['SE'] ** 2)
        
        # 固定效应估计
        pooled_aamr = np.sum(data['AAMR'] * data['weight']) / np.sum(data['weight'])
        
        # 标准误
        pooled_se = np.sqrt(1 / np.sum(data['weight']))
        
        # 95% CI
        ci_lower = pooled_aamr - 1.96 * pooled_se
        ci_upper = pooled_aamr + 1.96 * pooled_se
        
        # Cochran's Q检验异质性
        Q = np.sum(data['weight'] * (data['AAMR'] - pooled_aamr) ** 2)
        df_q = len(data) - 1
        p_heterogeneity = 1 - stats.chi2.cdf(Q, df_q) if df_q > 0 else None
        
        # I² 统计量
        I2 = max(0, (Q - df_q) / Q * 100) if Q > 0 else 0
        
        return pd.Series({
            'AAMR': pooled_aamr,
            'SE': pooled_se,
            'CI_Lower': ci_lower,
            'CI_Upper': ci_upper,
            'Population': data['Population'].sum(),
            'N_Counties': len(data),
            'Q_statistic': Q,
            'P_heterogeneity': p_heterogeneity,
            'I2': I2
        })
    
    if group_col is None:
        result = pool_group(df)
        return result.to_dict()
    else:
        result = df.groupby(group_col).apply(pool_group).reset_index()
        return result

def main():
    """
    Generate National_AAMR.csv, RUCC_AAMR.csv, merge them, and output to /Result/Table.
    Additionally, output AAMR_Top5.csv for Period 2006-2010 with selected cancer types.
    """
    # Define project root relative to this script's location
    # WDP/Code/Analysis/AAMR_IVW.py -> WDP/
    project_root = Path(__file__).resolve().parents[2]
    
    # Load the project configuration file (for consistency, though not used here)
    config_path = project_root / 'config.yaml'
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Warning: Configuration file not found at {config_path}")
        config = {}
    
    # Define paths
    eqi_dir = project_root / 'Data' / 'Original' / 'CDC WONDER EQI'
    rucc_file = project_root / 'Data' / 'Processed' / 'EQI' / 'EQI0005.csv'
    national_dir = project_root / 'Data' / 'Original' / 'CDC WONDER EQI AAMR TOTAL'
    output_dir = project_root / 'Result' / 'Tables'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Generate National_AAMR.csv
    print("Generating National_AAMR.csv...")
    national_csv_files = list(national_dir.glob('*.csv'))
    national_dfs = []
    for file_path in national_csv_files:
        filename = file_path.stem  # e.g., 'C00-C97 AAMR 2006-2010'
        parts = filename.split(' AAMR ')
        cancer_type = parts[0]
        period = parts[1]
        
        df = pd.read_csv(file_path, encoding='latin1')
        df = df[df['Notes'] == 'Total']
        df['Cancer_Type'] = cancer_type
        df['Period'] = period
        relevant_columns = ['Deaths', 'Population', 'Crude Rate', 'Crude Rate Standard Error', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error', 'Cancer_Type', 'Period']
        df = df[relevant_columns]
        national_dfs.append(df)
    
    national_df = pd.concat(national_dfs, ignore_index=True)
    national_df = national_df.dropna(subset=['Deaths'])
    national_df['Deaths'] = national_df['Deaths'].astype(int)
    national_df['Population'] = national_df['Population'].astype(int)
    national_df['Crude Rate'] = national_df['Crude Rate'].round(2)
    national_df['Crude Rate Standard Error'] = national_df['Crude Rate Standard Error'].round(2)
    national_df['Age Adjusted Rate'] = national_df['Age Adjusted Rate'].round(2)
    national_df['Age Adjusted Rate Standard Error'] = national_df['Age Adjusted Rate Standard Error'].round(2)
    column_order = ['Cancer_Type', 'Period', 'Deaths', 'Population', 'Crude Rate', 'Crude Rate Standard Error', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error']
    national_df = national_df[column_order]
    national_df = national_df.sort_values(by=['Cancer_Type', 'Period'])
    national_path = output_dir / 'AAMR_National.csv'
    national_df.to_csv(national_path, index=False, columns=column_order)
    print(f"AAMR_National.csv saved to {national_path}")
    
    # Step 2: Generate RUCC_AAMR.csv
    print("Generating RUCC_AAMR.csv...")
    rucc_df = pd.read_csv(rucc_file)
    rucc_df = rucc_df[['COUNTY_FIPS', 'RUCC']]
    
    csv_files = list(eqi_dir.glob('*.csv'))
    dfs = []
    for file_path in csv_files:
        filename = file_path.stem
        parts = filename.split(' ', 1)
        period = parts[0]
        cancer_type = parts[1]
        
        df = pd.read_csv(file_path, encoding='latin1')
        df['Deaths'] = pd.to_numeric(df['Deaths'], errors='coerce')
        df['Population'] = pd.to_numeric(df['Population'], errors='coerce')
        df = df[(df['Deaths'] == 0) | (df['Deaths'] >= 10)]
        df = df.dropna(subset=['Deaths'])
        
        df['Age Adjusted Rate'] = pd.to_numeric(df['Age Adjusted Rate'], errors='coerce')
        df['Age Adjusted Rate Lower 95% Confidence Interval'] = pd.to_numeric(df['Age Adjusted Rate Lower 95% Confidence Interval'], errors='coerce')
        df['Age Adjusted Rate Upper 95% Confidence Interval'] = pd.to_numeric(df['Age Adjusted Rate Upper 95% Confidence Interval'], errors='coerce')
        mask_unreliable = (df['Deaths'].between(10, 20)) & df['Age Adjusted Rate'].isna()
        if mask_unreliable.any():
            df.loc[mask_unreliable, 'Age Adjusted Rate'] = (
                df.loc[mask_unreliable, 'Age Adjusted Rate Lower 95% Confidence Interval'] +
                df.loc[mask_unreliable, 'Age Adjusted Rate Upper 95% Confidence Interval']
            ) / 2
        
        df['Cancer_Type'] = cancer_type
        df['Period'] = period
        relevant_columns = [
            'County Code', 'Deaths', 'Population', 
            'Crude Rate',
            'Age Adjusted Rate', 'Age Adjusted Rate Standard Error',
            'Cancer_Type', 'Period'
        ]
        df = df[relevant_columns]
        dfs.append(df)
    
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.merge(rucc_df, left_on='County Code', right_on='COUNTY_FIPS', how='left')
    combined_df = combined_df.dropna(subset=['RUCC'])
    combined_df['RUCC'] = combined_df['RUCC'].astype(int)
    
    results = []
    for (cancer_type, period), group in combined_df.groupby(['Cancer_Type', 'Period']):
        total_pop = group['Population'].sum()
        total_deaths = int(group['Deaths'].sum())
        
        cr_results = {}
        for r in [1, 2, 3, 4]:
            rucc_group = group[group['RUCC'] == r]
            if not rucc_group.empty:
                cr_values = pd.to_numeric(rucc_group['Crude Rate'], errors='coerce')
                if not cr_values.isna().all():
                    mean_cr = cr_values.mean()
                    if len(cr_values.dropna()) > 1:
                        se_cr = cr_values.std() / (len(cr_values.dropna()) ** 0.5)
                    else:
                        se_cr = 0.0
                else:
                    total_deaths_rucc = rucc_group['Deaths'].sum()
                    total_pop_rucc = rucc_group['Population'].sum()
                    if total_pop_rucc > 0:
                        mean_cr = (total_deaths_rucc / total_pop_rucc) * 100000
                        se_cr = (total_deaths_rucc ** 0.5) / total_pop_rucc * 100000
                    else:
                        mean_cr = float('nan')
                        se_cr = float('nan')
                if not pd.isna(mean_cr):
                    cr_results[r] = f"{mean_cr:.2f} ± {se_cr:.2f}"
                else:
                    cr_results[r] = 'N/A'
            else:
                cr_results[r] = 'N/A'
        
        aamr_data = group[['RUCC', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error', 'Population']].copy()
        aamr_data = aamr_data.rename(columns={'Age Adjusted Rate': 'AAMR', 'Age Adjusted Rate Standard Error': 'SE'})
        aamr_data['SE'] = pd.to_numeric(aamr_data['SE'], errors='coerce')
        aamr_data = aamr_data.dropna(subset=['AAMR', 'SE'])
        aamr_pooled = calculate_pooled_aamr_meta_analysis(aamr_data, group_col='RUCC')
        
        row = {
            'Cancer_Type': cancer_type,
            'Period': period,
            'Population': total_pop,
            'Deaths': total_deaths
        }
        for r in [1, 2, 3, 4]:
            row[f'RUCC{r}_CR'] = cr_results[r]
            
            aamr_row = aamr_pooled[aamr_pooled['RUCC'] == r]
            if not aamr_row.empty:
                mean_aamr = aamr_row['AAMR'].values[0]
                se_aamr = aamr_row['SE'].values[0]
                row[f'RUCC{r}_AAMR'] = f"{mean_aamr:.2f} ± {se_aamr:.2f}"
            else:
                row[f'RUCC{r}_AAMR'] = 'N/A'
        
        results.append(row)
    
    rucc_output_df = pd.DataFrame(results)
    rucc_output_df = rucc_output_df.sort_values(by=['Cancer_Type', 'Period'])
    rucc_path = output_dir / 'AAMR_RUCC.csv'
    rucc_output_df.to_csv(rucc_path, index=False)
    print(f"AAMR_RUCC.csv saved to {rucc_path}")
    
    # Step 3: Merge National and RUCC
    print("Merging National and RUCC data...")
    rucc_df_renamed = rucc_output_df.rename(columns={
        'Deaths': 'RUCC_Deaths',
        'Population': 'RUCC_Population'
    })
    
    merged_df = pd.merge(
        national_df, 
        rucc_df_renamed, 
        on=['Cancer_Type', 'Period'], 
        how='left'
    )
    print(f"Merged dataframe has {len(merged_df)} rows")
    
    merged_df['RUCC_Population'] = merged_df['RUCC_Population'].fillna(0).astype(int)
    merged_df['RUCC_Deaths'] = merged_df['RUCC_Deaths'].fillna(0).astype(int)
    print("Converted RUCC_Population and RUCC_Deaths to integers")
    
    merged_df['D_Population'] = ((merged_df['RUCC_Population'] - merged_df['Population']) / merged_df['Population'] * 100).round(2)
    merged_df['D_Death'] = ((merged_df['RUCC_Deaths'] - merged_df['Deaths']) / merged_df['Deaths'] * 100).round(2)
    print("Calculated D_Population and D_Death percentage differences")
    
    merged_df['AAMR'] = merged_df.apply(
        lambda row: f"{row['Age Adjusted Rate']} ± {row['Age Adjusted Rate Standard Error']}", 
        axis=1
    )
    
    merged_df['CR'] = merged_df.apply(
        lambda row: f"{row['Crude Rate']} ± {row['Crude Rate Standard Error']}", 
        axis=1
    )
    
    merged_df = merged_df.drop(columns=[
        'Age Adjusted Rate', 'Age Adjusted Rate Standard Error',
        'Crude Rate', 'Crude Rate Standard Error'
    ])
    print("Merged rate columns into formatted strings")
    
    column_order = [
        'Cancer_Type', 'Period', 'Population', 'RUCC_Population', 'D_Population', 'Deaths', 'RUCC_Deaths', 'D_Death',
        'AAMR', 'RUCC1_AAMR', 'RUCC2_AAMR', 'RUCC3_AAMR', 'RUCC4_AAMR',
        'CR', 'RUCC1_CR', 'RUCC2_CR', 'RUCC3_CR', 'RUCC4_CR'
    ]
    merged_df = merged_df[column_order]
    print("Reordered columns as specified")
    
    merged_output_path = output_dir / 'AAMR_IVW.csv'
    merged_df.to_csv(merged_output_path, index=False)
    print(f"Merged data saved to {merged_output_path}")
    
    # Step 4: Generate AAMR_Top5.csv
    print("Generating AAMR_Top5.csv...")
    top5_cancers = ['C00-C97', 'C34', 'C18-C21', 'C50', 'C25', 'C61']
    outcomes = ['All-site cancer', 'Lung cancer', 'Colorectal cancer', 'Breast cancer', 'Pancreatic cancer', 'Prostate cancer']
    outcome_map = dict(zip(top5_cancers, outcomes))
    
    top5_df = merged_df[(merged_df['Period'] == '2006-2010') & (merged_df['Cancer_Type'].isin(top5_cancers))].copy()
    top5_df['Outcome'] = top5_df['Cancer_Type'].map(outcome_map)
    top5_df = top5_df[['Cancer_Type', 'Outcome', 'Deaths', 'AAMR', 'RUCC1_AAMR', 'RUCC2_AAMR', 'RUCC3_AAMR', 'RUCC4_AAMR']]
    top5_df = top5_df.rename(columns={
        'AAMR': 'National',
        'RUCC1_AAMR': 'Metropolitan Urbanized (RUCC1)',
        'RUCC2_AAMR': 'Nonmetropolitan Urbanized (RUCC2)',
        'RUCC3_AAMR': 'Less Urbanized (RUCC3)',
        'RUCC4_AAMR': 'Thinly Populated (RUCC4)'
    })
    top5_df = top5_df.set_index('Cancer_Type').reindex(top5_cancers).reset_index()
    top5_path = output_dir / 'AAMR_Top5.csv'
    top5_df.to_csv(top5_path, index=False)
    print(f"AAMR_Top5.csv saved to {top5_path}")
    
    # Print summary
    print("\nMerge Summary:")
    print(f"- National rows: {len(national_df)}")
    print(f"- RUCC rows: {len(rucc_output_df)}")
    print(f"- Merged rows: {len(merged_df)}")
    print(f"- Top5 rows: {len(top5_df)}")

if __name__ == "__main__":
    main()