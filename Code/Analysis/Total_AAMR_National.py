import pandas as pd
import os
from pathlib import Path

# Define paths
data_dir = Path('/Users/ventus/Repository/WDP/Data/Original/CDC WONDER EQI AAMR TOTAL')
output_dir = Path('/Users/ventus/Repository/WDP/Result/Total_AAMR')
output_file = output_dir / 'National_AAMR.csv'

# Ensure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

# List all CSV files
csv_files = list(data_dir.glob('*.csv'))

# Initialize an empty list to hold DataFrames
dfs = []

for file_path in csv_files:
    # Extract cancer type and period from filename
    filename = file_path.stem  # e.g., 'C00-C97 AAMR 2006-2010'
    parts = filename.split(' AAMR ')
    cancer_type = parts[0]
    period = parts[1]
    
    # Read the CSV file with proper encoding
    df = pd.read_csv(file_path, encoding='latin1')
    
    # Filter to keep only the "Total" row
    df = df[df['Notes'] == 'Total']
    
    # Add cancer type and period columns
    df['Cancer_Type'] = cancer_type
    df['Period'] = period
    
    # Select relevant columns
    relevant_columns = ['Deaths', 'Population', 'Crude Rate', 'Crude Rate Standard Error', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error', 'Cancer_Type', 'Period']
    df = df[relevant_columns]
    
    # Append to list
    dfs.append(df)

# Concatenate all DataFrames
combined_df = pd.concat(dfs, ignore_index=True)

# Drop rows with NaN in Deaths (to remove comment rows)
combined_df = combined_df.dropna(subset=['Deaths'])

# Convert Deaths and Population to integers
combined_df['Deaths'] = combined_df['Deaths'].astype(int)
combined_df['Population'] = combined_df['Population'].astype(int)

# Round the rates to 2 decimal places
combined_df['Crude Rate'] = combined_df['Crude Rate'].round(2)
combined_df['Crude Rate Standard Error'] = combined_df['Crude Rate Standard Error'].round(2)
combined_df['Age Adjusted Rate'] = combined_df['Age Adjusted Rate'].round(2)
combined_df['Age Adjusted Rate Standard Error'] = combined_df['Age Adjusted Rate Standard Error'].round(2)

# Reorder columns
column_order = ['Cancer_Type', 'Period', 'Deaths', 'Population', 'Crude Rate', 'Crude Rate Standard Error', 'Age Adjusted Rate', 'Age Adjusted Rate Standard Error']
combined_df = combined_df[column_order]

# Sort by Cancer_Type and Period
combined_df = combined_df.sort_values(by=['Cancer_Type', 'Period'])

# Save to CSV
combined_df.to_csv(output_file, index=False, columns=column_order)

print(f"Data combined and saved to {output_file}")