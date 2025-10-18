import pandas as pd
import yaml
import numpy as np
from pathlib import Path

def main():
    """
    Process Climate_Zone.csv with specific modifications and filtering.
    """
    # Find project root and load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Define input and output paths
    input_path = project_root / config['eqi_aamr_outputs']['base_dir'] / "Climate_Zone.csv"
    output_path = project_root / config['eqi_aamr_outputs']['base_dir'] / "Climate_Zone_Processed.csv"
    
    # Read the dataframe
    print("Reading Climate_Zone.csv...")
    df = pd.read_csv(input_path, dtype={'COUNTY_FIPS': str})
    
    # Group by COUNTY_FIPS and select the best koppen row
    def select_best_koppen(group):
        if (group['koppen_coverage'] == 1.0).any():
            return group[group['koppen_coverage'] == 1.0].iloc[0]
        else:
            return group.nlargest(1, 'koppen_coverage').iloc[0]
    
    df = df.groupby('COUNTY_FIPS').apply(select_best_koppen).reset_index(drop=True)
    
    print(f"After selecting best koppen: {df.shape}")
    
    # Modify doe_zone_number
    print("Modifying doe_zone_number...")
    df['doe_zone_number'] = df['doe_zone_number'].replace({1: 2, 8: 7})
    
    # Modify koppen_major
    print("Modifying koppen_major...")
    df['koppen_major'] = df['koppen_major'].replace({'A': 'B', 'E': 'D'})
    
    # For doe_zone_code and koppen_code, keep only those with >30 samples
    print("Filtering doe_zone_code and koppen_code...")
    
    # Get counts for doe_zone_code
    doe_counts = df.groupby('doe_zone_code')['COUNTY_FIPS'].nunique()
    valid_doe = doe_counts[doe_counts > 30].index.tolist()
    
    # Get counts for koppen_code
    koppen_counts = df.groupby('koppen_code')['COUNTY_FIPS'].nunique()
    valid_koppen = koppen_counts[koppen_counts > 30].index.tolist()
    
    # Apply filtering
    df['doe_zone_code'] = df['doe_zone_code'].where(df['doe_zone_code'].isin(valid_doe), np.nan)
    df['koppen_code'] = df['koppen_code'].where(df['koppen_code'].isin(valid_koppen), np.nan)
    
    # Ensure RUCC and doe_zone_number are integers
    df['RUCC'] = pd.to_numeric(df['RUCC'], errors='coerce').astype('Int64')
    df['doe_zone_number'] = pd.to_numeric(df['doe_zone_number'], errors='coerce').astype('Int64')
    
    # Extract numeric codes from Census_Region and Census_Division
    print("Extracting numeric codes...")
    df['census_region'] = df['Census_Region'].str.extract(r'Census Region (\d+):').astype('Int64')
    df['census_division'] = df['Census_Division'].str.extract(r'Division (\d+):').astype('Int64')
    
    # Rename columns
    column_renames = {
        'RUCC': 'rucc',
        'doe_zone_number': 'doe_major',
        'doe_zone_code': 'doe_code'
    }
    df = df.rename(columns=column_renames)
    
    # Select final columns
    final_columns = [
        'COUNTY_FIPS',
        'census_region',
        'census_division', 
        'rucc',
        'koppen_code',
        'koppen_major',
        'koppen_coverage',
        'doe_major',
        'doe_code'
    ]
    df = df[final_columns]
    
    print(f"Final columns: {list(df.columns)}")
    print(f"DataFrame shape: {df.shape}")
    print("Sample data:")
    print(df.head(2))
    
    # Save to CSV
    print(f"Saving to {output_path}...")
    df.to_csv(output_path, index=False)
    print("Done.")