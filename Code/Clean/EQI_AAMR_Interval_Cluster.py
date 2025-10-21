# Code/Clean/EQI_AAMR_Interval_Cluster.py

import pandas as pd
import os
import yaml

def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config.yaml')
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def merge_clusters_with_eqi_aamr():
    """Merge cluster IDs with EQI_AAMR_Interval data"""
    print("Merging cluster IDs with EQI_AAMR_Interval data...")

    # Configuration
    config = get_config()
    processed_path = config['data_directories']['processed']
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    cluster_file = os.path.join(project_root, 'Result', 'Cluster_Visualization', 'EQI_Clusters_All_K.csv')
    k_values = [3, 4, 5, 6]  # Include clusters for k=3 to 6

    # Input files
    eqi_aamr_file = os.path.join(processed_path, 'df_EQI_AAMR', 'EQI_AAMR_Interval.csv')
    output_file = os.path.join(processed_path, 'df_EQI_AAMR', 'EQI_AAMR_Interval_Clustered.csv')

    # Load data
    print(f"Loading EQI_AAMR data from: {eqi_aamr_file}")
    eqi_aamr_df = pd.read_csv(eqi_aamr_file)

    print(f"Loading cluster data from: {cluster_file}")
    cluster_df = pd.read_csv(cluster_file)
    
    # Select the cluster columns for k=3 to 6
    cluster_cols = [f'cluster_{k}' for k in k_values]
    missing_cols = [col for col in cluster_cols if col not in cluster_df.columns]
    if missing_cols:
        raise ValueError(f"Cluster columns {missing_cols} not found in cluster data")
    
    cluster_df = cluster_df[['COUNTY_FIPS'] + cluster_cols]

    # Check required columns
    if 'COUNTY_FIPS' not in eqi_aamr_df.columns:
        raise ValueError("COUNTY_FIPS column not found in EQI_AAMR data")
    if 'COUNTY_FIPS' not in cluster_df.columns or not any(col.startswith('cluster_') for col in cluster_df.columns):
        raise ValueError("COUNTY_FIPS and cluster columns required in cluster data")

    print(f"EQI_AAMR data shape: {eqi_aamr_df.shape}")
    print(f"Cluster data shape: {cluster_df.shape}")

    # Merge on COUNTY_FIPS
    # Keep all rows from EQI_AAMR, add cluster columns where available
    merged_df = eqi_aamr_df.merge(cluster_df, on='COUNTY_FIPS', how='left')

    print(f"Merged data shape: {merged_df.shape}")
    for col in cluster_cols:
        unique_vals = sorted(merged_df[col].dropna().unique())
        print(f"Unique {col}: {unique_vals}")

    # Check for missing clusters
    for col in cluster_cols:
        missing = merged_df[col].isna().sum()
        if missing > 0:
            print(f"Warning: {missing} rows have missing {col} assignments")

    # Order columns to match original format
    first = ['COUNTY_FIPS','State','EQI_Period','Time_Period','Lag_Years','Cancer_Type','AAMR_lower','AAMR_upper','Smoking_Rate']
    eqi_cols = ['RUCC','EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
                'RUCC_EQI','RUCC_EQI_Air','RUCC_EQI_Water','RUCC_EQI_Land','RUCC_EQI_Built','RUCC_EQI_Social']
    ordered = first + [c for c in eqi_cols if c in merged_df.columns] + cluster_cols
    merged_df = merged_df[ordered]

    # Cast RUCC/EQI quintiles to nullable int to avoid 1.0/2.0 formatting
    for c in eqi_cols:
        if c in merged_df.columns:
            merged_df[c] = pd.to_numeric(merged_df[c], errors='coerce').astype('Int64')

    # Cast cluster columns to nullable int
    for col in cluster_cols:
        merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce').astype('Int64')

    # Ensure COUNTY_FIPS is properly formatted as 5-digit string
    merged_df['COUNTY_FIPS'] = merged_df['COUNTY_FIPS'].astype(str).str.zfill(5)

    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save merged data
    merged_df.to_csv(output_file, index=False)
    print(f"Successfully saved merged data to: {output_file}")

    # Summary statistics
    cancer_types = merged_df['Cancer_Type'].unique()
    print(f"Number of cancer types: {len(cancer_types)}")
    print(f"Cancer types: {', '.join(sorted(cancer_types))}")

    for col in cluster_cols:
        cluster_counts = merged_df.groupby(col).size()
        print(f"{col} distribution:")
        for cluster, count in cluster_counts.items():
            print(f"  Cluster {cluster}: {count} records")

if __name__ == '__main__':
    merge_clusters_with_eqi_aamr()