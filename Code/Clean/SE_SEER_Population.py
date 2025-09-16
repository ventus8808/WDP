import pandas as pd
import yaml
from pathlib import Path
import numpy as np

def load_config():
    """Load configuration from config.yaml"""
    # The script is in WDP/Code/Clean, so config.yaml is 2 levels up
    config_path = Path(__file__).resolve().parents[2] / 'config.yaml'
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def parse_seer_population_record(line):
    """
    Parse a single SEER population record according to the file dictionary.
    
    Fixed length ASCII text records (26 bytes)
    """
    # Check for minimum length to avoid errors on empty or malformed lines
    if len(line.strip()) < 26:
        return None
    
    # Variable Name and Values	Start Column	Length	Data Type
    # Year (1969, 1970, 1971...)	1	4	numeric
    # State postal abbreviation	5	2	character
    # State FIPS code	7	2	numeric
    # County FIPS code	9	3	numeric
    # Race	14	1	numeric
    # Origin	15	1	numeric
    # Sex	16	1	numeric
    # Age	17	2	numeric
    # Population	19	8	numeric
    record = {
        'Year': int(line[0:4]),
        'State_FIPS': line[6:8].strip(),
        'County_FIPS': line[8:11].strip(),
        'Race': int(line[13:14]),
        'Origin': int(line[14:15]),
        'Sex': int(line[15:16]),
        'Age': int(line[16:18]),
        'Population': int(line[18:26])
    }
    
    return record

def create_county_fips(state_fips, county_fips):
    """Create 5-digit county FIPS code"""
    # Handle special cases for hurricane evacuees (9-filled)
    if state_fips == '99' or county_fips == '999':
        return None
    
    try:
        # Ensure proper formatting
        state_str = state_fips.zfill(2)
        county_str = county_fips.zfill(3)
        return state_str + county_str
    except:
        return None

def process_seer_population_data(input_file, output_file):
    """
    Process SEER population data from fixed-width format to CSV.
    
    Data Mappings:
    
    Race (1990+ data):
        1 = White
        2 = Black
        3 = American Indian/Alaska Native
        4 = Asian or Pacific Islander
        
    Origin (Applicable to 1990+ data):
        0 = Non-Hispanic
        1 = Hispanic
        9 = Not applicable
        
    Sex:
        1 = Male
        2 = Female
        
    Age (Single age data):
        00 = 0 years
        01 = 1 year
        ...
        89 = 89 years
        90 = 90+ years
    """
    print(f"Processing SEER population data from {input_file}")
    print(f"Output will be saved to {output_file}")
    
    records = []
    line_count = 0
    error_count = 0
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                line_count += 1
                
                if line_count % 1000000 == 0:
                    print(f"Processed {line_count:,} lines, {len(records):,} valid records")
                
                try:
                    record = parse_seer_population_record(line)
                    if record is None:
                        error_count += 1
                        continue
                    
                    # Create county FIPS
                    county_fips = create_county_fips(record['State_FIPS'], record['County_FIPS'])
                    if county_fips is None:
                        continue  # Skip hurricane evacuee records and invalid FIPS
                    
                    # Keep original codes, no decoding
                    processed_record = {
                        'COUNTY_FIPS': county_fips,
                        'Year': record['Year'],
                        'Race': record['Race'],
                        'Origin': record['Origin'],
                        'Sex': record['Sex'],
                        'Age': record['Age'],
                        'Population': record['Population']
                    }
                    
                    records.append(processed_record)
                    
                except Exception as e:
                    error_count += 1
                    if error_count <= 10:  # Only print first 10 errors
                        print(f"Error processing line {line_count}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading file: {e}")
        return False
    
    print(f"\nProcessing complete:")
    print(f"Total lines processed: {line_count:,}")
    print(f"Valid records: {len(records):,}")
    print(f"Errors/skipped: {error_count:,}")
    
    if not records:
        print("No valid records found!")
        return False
    
    # Convert to DataFrame
    print("Converting to DataFrame...")
    df = pd.DataFrame(records)
    
    # Filter for years 1999-2020
    print(f"Filtering data to years 1999-2020...")
    df = df[(df['Year'] >= 1999) & (df['Year'] <= 2020)]
    print(f"Filtered data contains {len(df):,} records.")
    
    # Data quality checks
    print(f"\nData quality summary:")
    print(f"Years covered: {df['Year'].min()} - {df['Year'].max()}")
    print(f"Unique counties: {df['COUNTY_FIPS'].nunique():,}")
    print(f"Race codes: {np.sort(df['Race'].unique())}")
    print(f"Origin codes: {np.sort(df['Origin'].unique())}")
    print(f"Sex codes: {np.sort(df['Sex'].unique())}")
    print(f"Age codes: {np.sort(df['Age'].unique())}")
    print(f"Total population (1999-2020): {df['Population'].sum():,}")
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Save to CSV
    print(f"\nSaving to {output_file}...")
    df.to_csv(output_file, index=False)
    print("File saved successfully!")
    
    return True

def main():
    """Main function to process SEER population data"""
    try:
        # Load configuration
        config = load_config()
        project_root = Path(__file__).resolve().parents[2]
        
        # Define input file path
        input_file = project_root / "Data/Original/SEER Population/us.1990_2023.singleages.through89.90plus.adjusted.txt"
        
        # Define output file path using config
        processed_dir = config['data_directories']['processed']
        output_file = project_root / processed_dir / "Socioeconomic/Population_Structure.csv"
        
        # Check if input file exists
        if not input_file.exists():
            print(f"Error: Input file not found at {input_file}")
            return False
        
        # Process the data
        success = process_seer_population_data(input_file, output_file)
        
        if success:
            print(f"\n✅ SEER population data processing completed successfully!")
            print(f"Output saved to: {output_file}")
        else:
            print(f"\n❌ SEER population data processing failed!")
        
        return success
        
    except Exception as e:
        print(f"Error in main(): {e}")
        return False

if __name__ == "__main__":
    main()
