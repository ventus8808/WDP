import pandas as pd
import yaml
from pathlib import Path

def main():
    """
    Merge National_AAMR.csv and RUCC_AAMR.csv with left join on Cancer_Type and Period.
    Rename Deaths and Population in RUCC file to RUCC_Deaths and RUCC_Population.
    """
    # Define project root relative to this script's location
    # WDP/Code/Analysis/Total_AAMR.py -> WDP/
    project_root = Path(__file__).resolve().parents[2]
    
    # Load the project configuration file (for consistency, though not used here)
    config_path = project_root / 'config.yaml'
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Warning: Configuration file not found at {config_path}")
        config = {}
    
    # Define file paths
    national_path = project_root / 'Result' / 'Total_AAMR' / 'National_AAMR.csv'
    rucc_path = project_root / 'Result' / 'Total_AAMR' / 'RUCC_AAMR.csv'
    output_path = project_root / 'Result' / 'Total_AAMR' / 'Total_AAMR_Merged.csv'
    
    # Load dataframes
    print("Loading National_AAMR.csv...")
    national_df = pd.read_csv(national_path)
    print(f"Loaded {len(national_df)} rows from National_AAMR.csv")
    
    print("Loading RUCC_AAMR.csv...")
    rucc_df = pd.read_csv(rucc_path)
    print(f"Loaded {len(rucc_df)} rows from RUCC_AAMR.csv")
    
    # Rename columns in RUCC dataframe
    rucc_df = rucc_df.rename(columns={
        'Deaths': 'RUCC_Deaths',
        'Population': 'RUCC_Population'
    })
    print("Renamed Deaths -> RUCC_Deaths and Population -> RUCC_Population in RUCC dataframe")
    
    # Perform left join on Cancer_Type and Period
    merged_df = pd.merge(
        national_df, 
        rucc_df, 
        on=['Cancer_Type', 'Period'], 
        how='left'
    )
    print(f"Merged dataframe has {len(merged_df)} rows")
    
    # Convert RUCC_Population and RUCC_Deaths to integers
    merged_df['RUCC_Population'] = merged_df['RUCC_Population'].fillna(0).astype(int)
    merged_df['RUCC_Deaths'] = merged_df['RUCC_Deaths'].fillna(0).astype(int)
    print("Converted RUCC_Population and RUCC_Deaths to integers")
    
    # Calculate difference percentages
    merged_df['D_Population'] = ((merged_df['RUCC_Population'] - merged_df['Population']) / merged_df['Population'] * 100).round(2)
    merged_df['D_Death'] = ((merged_df['RUCC_Deaths'] - merged_df['Deaths']) / merged_df['Deaths'] * 100).round(2)
    print("Calculated D_Population and D_Death percentage differences")
    
    # Merge Age Adjusted Rate columns
    merged_df['AAMR'] = merged_df.apply(
        lambda row: f"{row['Age Adjusted Rate']} ± {row['Age Adjusted Rate Standard Error']}", 
        axis=1
    )
    
    # Merge Crude Rate columns  
    merged_df['CR'] = merged_df.apply(
        lambda row: f"{row['Crude Rate']} ± {row['Crude Rate Standard Error']}", 
        axis=1
    )
    
    # Drop the original separate columns
    merged_df = merged_df.drop(columns=[
        'Age Adjusted Rate', 'Age Adjusted Rate Standard Error',
        'Crude Rate', 'Crude Rate Standard Error'
    ])
    print("Merged rate columns into formatted strings")
    
    # Reorder columns
    column_order = [
        'Cancer_Type', 'Period', 'Population', 'RUCC_Population', 'D_Population', 'Deaths', 'RUCC_Deaths', 'D_Death',
        'AAMR', 'RUCC1_AAMR', 'RUCC2_AAMR', 'RUCC3_AAMR', 'RUCC4_AAMR',
        'CR', 'RUCC1_CR', 'RUCC2_CR', 'RUCC3_CR', 'RUCC4_CR'
    ]
    merged_df = merged_df[column_order]
    print("Reordered columns as specified")
    
    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save merged dataframe
    merged_df.to_csv(output_path, index=False)
    print(f"Merged data saved to {output_path}")
    
    # Print summary
    print("\nMerge Summary:")
    print(f"- National rows: {len(national_df)}")
    print(f"- RUCC rows: {len(rucc_df)}")
    print(f"- Merged rows: {len(merged_df)}")
    print(f"- Columns in merged: {list(merged_df.columns)}")

if __name__ == "__main__":
    main()