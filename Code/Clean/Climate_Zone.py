import pandas as pd
import yaml
from pathlib import Path

def main():
    """
    Combine climate zone data with census regions and RUCC codes by county FIPS.

    Reads:
    - Location.csv for census regions/divisions
    - koeppengeigerUScounty.txt for Köppen-Geiger climate classes
    - climate_zones.csv for IECC/BA climate zones (fixes FIPS by combining state+county)
    - EQI0005.csv for RUCC codes

    Outputs merged dataframe to Data/Processed/df_EQI_AAMR/Climate_Zone.csv
    with COUNTY_FIPS as string.
    """
    # Find project root and load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Define input paths
    location_path = project_root / config['data_sources']['cdc_wonder']['location_output_file']
    urbanization_path = project_root / config['data_sources']['cdc_wonder']['urbanization_output_file']
    koeppen_path = project_root / "Data/Original/Climate Zone/koeppengeigerUScounty.txt"
    climate_zones_path = project_root / "Data/Original/Climate Zone/climate_zones.csv"
    eqi_path = project_root / config['data_sources']['epa_eqi']['processed'] / "EQI0005.csv"
    
    # Define output path
    output_dir = project_root / config['eqi_aamr_outputs']['base_dir']
    output_path = output_dir / "Climate_Zone.csv"
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read Location.csv for census data
    print("Reading Location.csv...")
    location_df = pd.read_csv(location_path, dtype={'COUNTY_FIPS': str})
    location_df = location_df[['COUNTY_FIPS', 'Census_Region', 'Census_Division']]
    
    # Read Köppen-Geiger data
    print("Reading koeppengeigerUScounty.txt...")
    koeppen_df = pd.read_csv(koeppen_path, sep='\t', dtype={'FIPS': str})
    koeppen_df = koeppen_df[['FIPS', 'CLS', 'PROP']]
    koeppen_df.rename(columns={'FIPS': 'COUNTY_FIPS'}, inplace=True)
    
    # Read climate_zones.csv and fix FIPS
    print("Reading climate_zones.csv...")
    climate_df = pd.read_csv(climate_zones_path, dtype={'State FIPS': str, 'County FIPS': str})
    # Combine State FIPS + County FIPS to make 5-digit FIPS
    climate_df['COUNTY_FIPS'] = climate_df['State FIPS'] + climate_df['County FIPS'].str.zfill(3)
    climate_df = climate_df[['COUNTY_FIPS', 'IECC Climate Zone', 'IECC Moisture Regime', 'BA Climate Zone']]
    
    # Read EQI0005.csv for RUCC
    print("Reading EQI0005.csv...")
    eqi_df = pd.read_csv(eqi_path, dtype={'COUNTY_FIPS': str})
    eqi_df = eqi_df[['COUNTY_FIPS', 'RUCC']]
    eqi_df['RUCC'] = pd.to_numeric(eqi_df['RUCC'], errors='coerce').astype('Int64')
    
    # Merge all dataframes on COUNTY_FIPS
    print("Merging dataframes...")
    merged_df = location_df.merge(koeppen_df, on='COUNTY_FIPS', how='outer')
    merged_df = merged_df.merge(climate_df, on='COUNTY_FIPS', how='outer')
    merged_df = merged_df.merge(eqi_df, on='COUNTY_FIPS', how='outer')
    
    # Rename columns
    column_renames = {
        'CLS': 'koppen_code',
        'PROP': 'koppen_coverage',
        'IECC Climate Zone': 'doe_zone_number',
        'IECC Moisture Regime': 'doe_moisture',
        'BA Climate Zone': 'doe_zone_name'
    }
    merged_df = merged_df.rename(columns=column_renames)
    
    # Convert doe_zone_number to int (remove .0)
    merged_df['doe_zone_number'] = pd.to_numeric(merged_df['doe_zone_number'], errors='coerce').astype('Int64')
    
    # Add new columns
    # koppen_major: first letter of koppen_code
    merged_df['koppen_major'] = merged_df['koppen_code'].str[0]
    
    # koppen_major_name: map to names
    koppen_name_map = {
        'A': 'Tropical',
        'B': 'Dry', 
        'C': 'Temperate',
        'D': 'Continental',
        'E': 'Polar'
    }
    merged_df['koppen_major_name'] = merged_df['koppen_major'].map(koppen_name_map)
    
    # doe_zone_code: combine number and moisture
    merged_df['doe_zone_code'] = merged_df.apply(
        lambda row: f"{row['doe_zone_number']}{row['doe_moisture']}" 
        if pd.notna(row['doe_moisture']) and row['doe_moisture'] != ''
        else str(row['doe_zone_number']) if pd.notna(row['doe_zone_number']) else None,
        axis=1
    )
    
    # Reorder columns as per specification
    column_order = [
        'COUNTY_FIPS',
        'Census_Region', 
        'Census_Division',
        'RUCC',
        'koppen_code',
        'koppen_major',
        'koppen_major_name', 
        'koppen_coverage',
        'doe_zone_number',
        'doe_moisture',
        'doe_zone_code',
        'doe_zone_name'
    ]
    merged_df = merged_df[column_order]
    
    # Ensure COUNTY_FIPS is string
    merged_df['COUNTY_FIPS'] = merged_df['COUNTY_FIPS'].astype(str)
    
    # Sort by COUNTY_FIPS
    merged_df = merged_df.sort_values('COUNTY_FIPS').reset_index(drop=True)
    
    # Save to CSV
    print(f"Saving to {output_path}...")
    merged_df.to_csv(output_path, index=False)
    print("Done.")

if __name__ == '__main__':
    main()