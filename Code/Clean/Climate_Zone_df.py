import pandas as pd
import yaml
from pathlib import Path

def main():
    """
    Merge EQI_AAMR_Interval.csv with Climate_Zone_Processed.csv after removing RUCC* columns.

    Reads EQI_AAMR_Interval.csv, removes RUCC_EQI* columns, merges with Climate_Zone_Processed.csv
    on COUNTY_FIPS, outputs to EQI_AAMR_Interval_Climate.csv
    """
    # Find project root and load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Define input paths
    eqi_interval_path = project_root / config['eqi_aamr_outputs']['base_dir'] / config['eqi_aamr_outputs']['eqi_aamr_interval']
    climate_path = project_root / config['eqi_aamr_outputs']['base_dir'] / "Climate_Zone_Processed.csv"
    
    # Define output path
    output_path = project_root / config['eqi_aamr_outputs']['base_dir'] / "EQI_AAMR_Interval_Climate.csv"
    
    # Read EQI_AAMR_Interval.csv
    print("Reading EQI_AAMR_Interval.csv...")
    eqi_df = pd.read_csv(eqi_interval_path, dtype={'COUNTY_FIPS': str})
    
    # Remove RUCC* columns (those starting with RUCC_)
    rucc_cols = [col for col in eqi_df.columns if col.startswith('RUCC_')]
    print(f"Removing columns: {rucc_cols}")
    eqi_df = eqi_df.drop(columns=rucc_cols)
    
    # Read Climate_Zone_Processed.csv
    print("Reading Climate_Zone_Processed.csv...")
    climate_df = pd.read_csv(climate_path, dtype={'COUNTY_FIPS': str})
    print(f"Climate columns: {list(climate_df.columns)}")
    
    # Merge on COUNTY_FIPS
    print("Merging dataframes...")
    merged_df = eqi_df.merge(climate_df, on='COUNTY_FIPS', how='left')
    
    # Convert specified columns to integers
    columns_to_int = ['RUCC', 'EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social', 'doe_major', 'census_region', 'census_division', 'rucc']
    for col in columns_to_int:
        merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').astype('Int64')
    
    # Sort by COUNTY_FIPS, State, etc. for consistency
    sort_cols = ['COUNTY_FIPS', 'State', 'EQI_Period', 'Time_Period', 'Lag_Years', 'Cancer_Type']
    merged_df = merged_df.sort_values(sort_cols).reset_index(drop=True)
    
    # Save to CSV
    print(f"Saving to {output_path}...")
    merged_df.to_csv(output_path, index=False)
    print("Done.")

if __name__ == '__main__':
    main()