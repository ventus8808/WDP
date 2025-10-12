import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

# Define paths
eqi_dir = Path('/Users/ventus/Repository/WDP/Data/Original/CDC WONDER EQI')
rucc_file = Path('/Users/ventus/Repository/WDP/Data/Processed/EQI/EQI0005.csv')
output_dir = Path('/Users/ventus/Repository/WDP/Result/Total_AAMR')
output_file = output_dir / 'RUCC_AAMR.csv'

# Ensure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

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

# Read RUCC data
rucc_df = pd.read_csv(rucc_file)
rucc_df = rucc_df[['COUNTY_FIPS', 'RUCC']]

# List all EQI CSV files
csv_files = list(eqi_dir.glob('*.csv'))

# Initialize list for DataFrames
dfs = []

for file_path in csv_files:
    # Extract period and cancer_type from filename, e.g., '2006-2010 C15-C26.csv'
    filename = file_path.stem
    parts = filename.split(' ', 1)
    period = parts[0]
    cancer_type = parts[1]
    
    # Read CSV
    df = pd.read_csv(file_path, encoding='latin1')
    
    # Convert Deaths to numeric, coerce errors to NaN
    df['Deaths'] = pd.to_numeric(df['Deaths'], errors='coerce')
    df['Population'] = pd.to_numeric(df['Population'], errors='coerce')
    
    # Filter out suppressed data: Deaths [1,9] and NaN
    df = df[(df['Deaths'] == 0) | (df['Deaths'] >= 10)]
    df = df.dropna(subset=['Deaths'])
    
    # Fix AAMR for Deaths [10,20]: if Age Adjusted Rate is not numeric, calculate from CI
    df['Age Adjusted Rate'] = pd.to_numeric(df['Age Adjusted Rate'], errors='coerce')
    df['Age Adjusted Rate Lower 95% Confidence Interval'] = pd.to_numeric(df['Age Adjusted Rate Lower 95% Confidence Interval'], errors='coerce')
    df['Age Adjusted Rate Upper 95% Confidence Interval'] = pd.to_numeric(df['Age Adjusted Rate Upper 95% Confidence Interval'], errors='coerce')
    mask_unreliable = (df['Deaths'].between(10, 20)) & df['Age Adjusted Rate'].isna()
    if mask_unreliable.any():
        df.loc[mask_unreliable, 'Age Adjusted Rate'] = (
            df.loc[mask_unreliable, 'Age Adjusted Rate Lower 95% Confidence Interval'] +
            df.loc[mask_unreliable, 'Age Adjusted Rate Upper 95% Confidence Interval']
        ) / 2
    
    # Add columns
    df['Cancer_Type'] = cancer_type
    df['Period'] = period
    
    # Select relevant columns
    relevant_columns = [
        'County Code', 'Deaths', 'Population', 
        'Crude Rate',
        'Age Adjusted Rate', 'Age Adjusted Rate Standard Error',
        'Cancer_Type', 'Period'
    ]
    df = df[relevant_columns]
    
    dfs.append(df)

# Concatenate all
combined_df = pd.concat(dfs, ignore_index=True)

# Merge with RUCC
combined_df = combined_df.merge(rucc_df, left_on='County Code', right_on='COUNTY_FIPS', how='left')

# Drop rows without RUCC
combined_df = combined_df.dropna(subset=['RUCC'])

# Convert RUCC to int
combined_df['RUCC'] = combined_df['RUCC'].astype(int)

# Group by Cancer_Type, Period, RUCC and calculate pooled rates
results = []
for (cancer_type, period), group in combined_df.groupby(['Cancer_Type', 'Period']):
    total_pop = group['Population'].sum()
    total_deaths = int(group['Deaths'].sum())
    
    # Calculate pooled CR (mean, since no SE available)
    cr_results = {}
    for r in [1, 2, 3, 4]:
        rucc_group = group[group['RUCC'] == r]
        if not rucc_group.empty:
            # Convert Crude Rate to numeric, ignoring 'Unreliable'
            cr_values = pd.to_numeric(rucc_group['Crude Rate'], errors='coerce')
            if not cr_values.isna().all():
                mean_cr = cr_values.mean()
                if len(cr_values.dropna()) > 1:
                    se_cr = cr_values.std() / (len(cr_values.dropna()) ** 0.5)
                else:
                    se_cr = 0.0
            else:
                # Fallback to aggregated deaths / population * 100000
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
    
    # Prepare data for AAMR
    aamr_data = group[['RUCC', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error', 'Population']].copy()
    aamr_data = aamr_data.rename(columns={'Age Adjusted Rate': 'AAMR', 'Age Adjusted Rate Standard Error': 'SE'})
    aamr_data['SE'] = pd.to_numeric(aamr_data['SE'], errors='coerce')
    aamr_data = aamr_data.dropna(subset=['AAMR', 'SE'])
    aamr_pooled = calculate_pooled_aamr_meta_analysis(aamr_data, group_col='RUCC')
    
    # Build row
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

# Create output DataFrame
output_df = pd.DataFrame(results)

# Sort
output_df = output_df.sort_values(by=['Cancer_Type', 'Period'])

# Save
output_df.to_csv(output_file, index=False)

print(f"Data processed and saved to {output_file}")