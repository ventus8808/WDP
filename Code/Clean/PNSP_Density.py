# -*- coding: utf-8 -*-
"""
WDP Pesticide Exposure Density Calculation Script

This script calculates pesticide exposure density (kg/km²) by normalizing
pesticide application weight by the agricultural area of each county.

Input:
- Processed pesticide application data (kg/county/year)
- Processed land cover data (agricultural area in km²)

Output:
- A new CSV file containing pesticide exposure density values.
"""

import pandas as pd
import numpy as np
import os
import yaml
from pathlib import Path

def load_config(config_path):
    """Loads the YAML configuration file."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at: {config_path}")
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def calculate_density(config):
    """
    Loads pesticide and land use data, calculates pesticide density,
    and saves the result to a new CSV file.
    """
    # Corrected path access based on the actual config.yaml structure
    base_path = Path(config['base_paths']['wdp_root'])
    
    # Define input paths from config using the correct keys
    pesticide_processed_dir = config['data_sources']['usgs_pnsp']['processed']
    environmental_processed_dir = config['data_sources']['gee']['processed']

    pesticide_path = base_path / pesticide_processed_dir / 'PNSP.csv'
    land_use_path = base_path / environmental_processed_dir / 'NLCD_JRC.csv'
    
    # Define output path
    output_dir = base_path / pesticide_processed_dir
    output_path = output_dir / 'PNSP_Density.csv'
    
    print("Starting pesticide density calculation...")
    print(f"Loading pesticide weights from: {pesticide_path}")
    print(f"Loading land use data from: {land_use_path}")

    # Load datasets
    try:
        pesticide_df = pd.read_csv(pesticide_path)
        land_use_df = pd.read_csv(land_use_path)
    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
        return

    # --- Data Preparation ---
    # Standardize column names to handle 'Year' vs 'YEAR' inconsistency
    if 'YEAR' in pesticide_df.columns and 'Year' not in pesticide_df.columns:
        pesticide_df.rename(columns={'YEAR': 'Year'}, inplace=True)
    if 'YEAR' in land_use_df.columns and 'Year' not in land_use_df.columns:
        land_use_df.rename(columns={'YEAR': 'Year'}, inplace=True)

    # Ensure key columns are of the same type
    pesticide_df['COUNTY_FIPS'] = pesticide_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    land_use_df['COUNTY_FIPS'] = land_use_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    # Select only necessary columns from land use data
    ag_area_df = land_use_df[['COUNTY_FIPS', 'Year', 'nlcd_agriculture_km2']]

    # Merge pesticide data with agricultural area data
    merged_df = pd.merge(pesticide_df, ag_area_df, on=['COUNTY_FIPS', 'Year'], how='left')
    print(f"Successfully merged {len(merged_df)} records.")

    # --- Density Calculation ---
    # Identify all pesticide columns (starting with 'cat' or 'chem')
    pesticide_cols = [col for col in merged_df.columns if col.startswith('cat') or col.startswith('chem')]
    
    # Handle cases where agricultural area is zero or missing to avoid division by zero
    # Replace area <= 0 with NaN, so that division results in NaN
    merged_df['nlcd_agriculture_km2_safe'] = merged_df['nlcd_agriculture_km2'].replace(0, np.nan)

    print(f"Calculating density for {len(pesticide_cols)} pesticide variables...")
    # Calculate density
    for col in pesticide_cols:
        merged_df[col] = merged_df[col] / merged_df['nlcd_agriculture_km2_safe']

    # --- Finalize and Save ---
    # Create the final density dataframe
    density_df = merged_df[['COUNTY_FIPS', 'Year'] + pesticide_cols]
    
    # Fill any resulting NaN/inf values with 0, assuming no density if area is 0 or no application
    density_df.fillna(0, inplace=True)
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save to CSV
    density_df.to_csv(output_path, index=False)
    print(f"✓ Pesticide density data successfully saved to: {output_path}")
    print(f"  - Shape of the output file: {density_df.shape}")


if __name__ == '__main__':
    # Assuming the script is run from the root of the project directory
    # or has access to the main config file.
    # We construct a relative path to the config file.
    try:
        config_path = 'config.yaml'
        config = load_config(config_path)
        calculate_density(config)
    except FileNotFoundError:
        # Fallback for different execution context
        try:
            # Assuming execution from Code/Clean
            config_path = '../../config.yaml'
            config = load_config(config_path)
            calculate_density(config)
        except Exception as e:
            print(f"An error occurred: {e}")
            print("Please ensure you are running this script from the project root or 'Code/Clean' directory,")
            print("and the 'config.yaml' file is correctly located.")

