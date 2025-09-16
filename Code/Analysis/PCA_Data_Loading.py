#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PCA Data Loading Module

This module loads and merges all processed data files required for PCA analysis.
It handles the integration of CDC location/urbanization data, socioeconomic data, 
and environmental data to create a master dataset for PCA.

Author: WDP Analysis Team (with Qoder AI)
Date: 2024-09-04
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# --- Configuration ---
# 使用相对于项目根目录的路径
PROJECT_ROOT = Path(__file__).resolve().parents[2]

def aggregate_population_to_total(pop_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate the detailed population structure to get total population by county-year."""
    print("Aggregating population data...")
    total_pop = pop_df.groupby(['COUNTY_FIPS', 'Year'])['Population'].sum().reset_index()
    total_pop.columns = ['COUNTY_FIPS', 'Year', 'Total_Population']
    return total_pop

def load_all_processed_data() -> pd.DataFrame:
    """
    Load and merge all processed data files to create a master dataset for PCA.
    
    Returns:
        pd.DataFrame: Master dataset with all variables needed for PCA
    """
    print("\n" + "="*60)
    print("LOADING DATA FOR PCA ANALYSIS")
    print("="*60)
    
    # Define data paths
    processed_cdc = PROJECT_ROOT / "Data/Processed/CDC"
    processed_socio = PROJECT_ROOT / "Data/Processed/Socioeconomic"
    processed_env = PROJECT_ROOT / "Data/Processed/Environmental"
    
    # 1. Load Location Data (baseline with COUNTY_FIPS)
    print("Loading CDC location data...")
    location_df = pd.read_csv(processed_cdc / "Location.csv")
    location_df = location_df[['COUNTY_FIPS', 'County']].drop_duplicates()
    # Clean and standardize COUNTY_FIPS in location data
    location_df = location_df[~location_df['COUNTY_FIPS'].str.contains('nan|000', na=False)]  # type: ignore
    location_df['COUNTY_FIPS'] = location_df['COUNTY_FIPS'].astype(str).str.zfill(5)  # type: ignore
    
    # 2. Load Urbanization Data (contains Year information)
    print("Loading urbanization data...")
    urban_df = pd.read_csv(processed_cdc / "Urbanization.csv")
    urban_df = urban_df[['COUNTY_FIPS', 'Year', 'Urbanization_Code']].copy()
    # Standardize COUNTY_FIPS format in urbanization data
    urban_df['COUNTY_FIPS'] = urban_df['COUNTY_FIPS'].astype(str).str.zfill(5)  # type: ignore
    
    # Start with urbanization as the base (has COUNTY_FIPS and Year)
    master_df = urban_df.copy()
    
    # 3. Load Socioeconomic Data
    print("Loading socioeconomic data...")
    
    # Poverty and Income
    poverty_df = pd.read_csv(processed_socio / "Poverty_Income.csv")
    poverty_df['COUNTY_FIPS'] = poverty_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    master_df = master_df.merge(poverty_df, on=['COUNTY_FIPS', 'Year'], how='left')
    
    # Unemployment
    unemployment_df = pd.read_csv(processed_socio / "Unemployment.csv")
    unemployment_df = unemployment_df[['COUNTY_FIPS', 'Year', 'Unemployment_Rate']].copy()
    unemployment_df['COUNTY_FIPS'] = unemployment_df['COUNTY_FIPS'].astype(str).str.zfill(5)  # type: ignore
    master_df = master_df.merge(unemployment_df, on=['COUNTY_FIPS', 'Year'], how='left')
    
    # Education
    education_df = pd.read_csv(processed_socio / "Education.csv")
    education_df['COUNTY_FIPS'] = education_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    master_df = master_df.merge(education_df, on=['COUNTY_FIPS', 'Year'], how='left')
    
    # BEA Economic Data (Per Capita Income)
    try:
        gdp_df = pd.read_csv(processed_socio / "GDP.csv")
        gdp_df['COUNTY_FIPS'] = gdp_df['COUNTY_FIPS'].astype(str).str.zfill(5)
        # Assuming GDP has Per_Capita_Income column
        if 'Per_Capita_Income' in gdp_df.columns:
            gdp_subset = gdp_df[['COUNTY_FIPS', 'Year', 'Per_Capita_Income']].copy()
            master_df = master_df.merge(gdp_subset, on=['COUNTY_FIPS', 'Year'], how='left')
        else:
            print("  Warning: Per_Capita_Income not found in GDP data")
    except FileNotFoundError:
        print("  Warning: GDP.csv not found, continuing without Per_Capita_Income")
    
    # Population Data (aggregate to total)
    try:
        population_df = pd.read_csv(processed_socio / "Population_Structure.csv")
        population_df['COUNTY_FIPS'] = population_df['COUNTY_FIPS'].astype(str).str.zfill(5)
        total_pop_df = aggregate_population_to_total(population_df)
        master_df = master_df.merge(total_pop_df, on=['COUNTY_FIPS', 'Year'], how='left')
    except FileNotFoundError:
        print("  Warning: Population_Structure.csv not found")
    
    # 4. Load Environmental Data
    print("Loading environmental data...")
    
    # NLDAS Climate Data
    try:
        nldas_df = pd.read_csv(processed_env / "NLDAS.csv")
        nldas_df['COUNTY_FIPS'] = nldas_df['COUNTY_FIPS'].astype(str).str.zfill(5)
        # Select only the climate variables we need for PCA based on user preference
        climate_vars = [
            'COUNTY_FIPS', 'year', 'tas_mean_annual', 'prcp_sum_annual', 
            'wind_mean_annual', 'rh_mean_annual', 'swrad_mean_annual', 'potevap_sum_annual',
            'tas_mean_DJF', 'tas_mean_MAM', 'tas_mean_JJA', 'tas_mean_SON',
            'prcp_sum_DJF', 'prcp_sum_MAM', 'prcp_sum_JJA', 'prcp_sum_SON'
        ]
        
        # Check which columns actually exist
        available_climate_vars = ['COUNTY_FIPS']
        for var in climate_vars[1:]:  # Skip COUNTY_FIPS
            if var in nldas_df.columns:
                available_climate_vars.append(var)
            elif var == 'year' and 'Year' in nldas_df.columns:
                available_climate_vars.append('Year')
            else:
                print(f"  Warning: {var} not found in NLDAS data")
        
        nldas_subset = nldas_df[available_climate_vars].copy()
        
        # Rename 'year' to 'Year' if needed for consistency
        if 'year' in nldas_subset.columns:
            nldas_subset.rename(columns={'year': 'Year'}, inplace=True)  # type: ignore
        
        master_df = master_df.merge(nldas_subset, on=['COUNTY_FIPS', 'Year'], how='left')
        
    except FileNotFoundError:
        print("  Warning: NLDAS.csv not found")
    
    # Add location information
    print("Adding location information...")
    master_df = master_df.merge(location_df, on=['COUNTY_FIPS'], how='left')  # type: ignore
    
    # 5. Data Quality and Preprocessing
    print("Performing data quality checks...")
    
    # Convert COUNTY_FIPS to string and ensure 5-digit format
    master_df['COUNTY_FIPS'] = master_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    # Remove rows with invalid FIPS codes
    master_df = master_df[~master_df['COUNTY_FIPS'].str.contains('nan|000')]
    
    # Filter to reasonable year range (1999-2020)
    master_df = master_df[(master_df['Year'] >= 1999) & (master_df['Year'] <= 2020)]
    
    # Add Per_Capita_Income from GDP data if not already present
    if 'Per_Capita_Income' not in master_df.columns and 'Median_Household_Income' in master_df.columns:  # type: ignore
        # Use Median_Household_Income as a proxy for Per_Capita_Income if GDP data isn't available
        master_df['Per_Capita_Income'] = master_df['Median_Household_Income']
    
    # Print summary statistics
    print(f"\nData loading complete!")
    print(f"Final dataset shape: {master_df.shape}")
    print(f"Counties: {master_df['COUNTY_FIPS'].nunique()}")  # type: ignore
    print(f"Years: {sorted(master_df['Year'].unique())}")  # type: ignore
    print(f"Available variables: {list(master_df.columns)}")  # type: ignore
    
    # Check for missing data in key variables
    key_vars = ['Poverty_Percent_All_Ages', 'Unemployment_Rate', 'Less_Than_High_School_Percent',
                'College_Plus_Percent', 'Median_Household_Income', 'Urbanization_Code']
    
    print(f"\nMissing data summary for key variables:")
    for var in key_vars:
        if var in master_df.columns:  # type: ignore
            missing_pct = (master_df[var].isna().sum() / len(master_df)) * 100  # type: ignore
            print(f"  {var}: {missing_pct:.1f}% missing")
    
    return master_df  # type: ignore