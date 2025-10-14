import pandas as pd
import os

# Define the base directory
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Define the result directory
result_dir = os.path.join(base_dir, 'Result', 'brms')

# Define the output directory
output_dir = os.path.join(base_dir, 'Result', 'brms_extract')
os.makedirs(output_dir, exist_ok=True)

# List of files to read for top5
top5_files = ['C00_C97_brms.csv', 'C34_brms.csv', 'C18_C21_brms.csv', 'C50_brms.csv', 'C61_brms.csv','C15_C26_brms.csv']

# Get all brms files
all_files = [f for f in os.listdir(result_dir) if f.endswith('_brms.csv')]

# Mapping
icd_mapping = {
    'C00_C97': 'All-site Cancer',
    'C15_C26': 'Digestive System Cancer',
    'C18_C21': 'Colorectal Cancer',
    'C22': 'Liver Cancer',
    'C25': 'Pancreatic Cancer',
    'C30_C39': 'Respiratory System Cancer',
    'C34': 'Lung Cancer',
    'C40_C41': 'Bone and Articular Cartilage Cancer',
    'C43_C44': 'Skin Cancer',
    'C45_C49': 'Mesothelial and Soft Tissue Cancer',
    'C50': 'Breast Cancer',
    'C51_C58': 'Female Genital Organs Cancer',
    'C60_C63': 'Male Genital Organs Cancer',
    'C61': 'Prostate Cancer',
    'C64_C68': 'Urinary Tract Cancer',
    'C64': 'Kidney Cancer',
    'C67': 'Bladder Cancer',
    'C69_C72': 'Brain and Central Nervous System Cancer',
    'C73_C75': 'Thyroid and Endocrine Glands Cancer',
    'C76_C80': 'Ill-defined and Secondary Sites Cancer',
    'C81_C96': 'Lymphoid and Hematopoietic Cancer',
    'C82_C86': 'Non-Hodgkin Lymphoma',
    'C91_C95': 'Leukemia',
    'C97': 'Multiple Primary Sites Cancer',
}


# Function to process files
def process_files(files, output_name, period='2000_2005'):
    dfs = []
    for file in files:
        file_path = os.path.join(result_dir, file)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            # Filter rows where Model is 'EQI' and EQI_Period is the specified period
            filtered_df = df[(df['Model'] == 'EQI') & (df['EQI_Period'] == period)]
            dfs.append(filtered_df)
        else:
            print(f"File {file_path} does not exist.")
    
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        # Add Chinese name column
        combined_df['Chinese_Name'] = combined_df['ICD_Code'].map(icd_mapping).fillna('未知')
        # Drop AAMR_Period and Model columns (keep EQI_Period)
        combined_df = combined_df.drop(columns=['AAMR_Period', 'Model'])
        # Convert Lag to int
        combined_df['Lag'] = combined_df['Lag'].astype(int)
        # Rename Chinese_Name to type
        combined_df = combined_df.rename(columns={'Chinese_Name': 'type'})
        # Rename type to Outcome
        combined_df = combined_df.rename(columns={'type': 'Outcome'})
        # Reorder columns: ICD_Code, Outcome, EQI_Period, Lag, Q1, Q2, Q3, Q4, Q5
        combined_df = combined_df[['ICD_Code', 'Outcome', 'EQI_Period', 'Lag', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']]
        # Sort by ICD_Code, EQI_Period, Lag for non-Top5, or by custom order for Top5
        if 'Top5' in output_name:
            order_mapping = {'C00_C97': 0, 'C34': 1, 'C18_C21': 2, 'C50': 3, 'C61': 4, 'C15_C26': 5}
            combined_df['order'] = combined_df['ICD_Code'].map(order_mapping)
            combined_df = combined_df.sort_values(by=['order', 'EQI_Period', 'Lag'])
            combined_df = combined_df.drop(columns=['order'])
        else:
            combined_df = combined_df.sort_values(by=['ICD_Code', 'EQI_Period', 'Lag'])
        # Save to CSV
        output_path = os.path.join(output_dir, output_name)
        combined_df.to_csv(output_path, index=False)
        print(f"Extracted data saved to {output_path}")
    else:
        print(f"No data to save for {output_name}.")

# Process 0005 top5
process_files(top5_files, 'Overview_0005_Top5.csv', '2000_2005')

# Process 0005 all
process_files(all_files, 'Overview_0005_All.csv', '2000_2005')

# Process 0610 top5
process_files(top5_files, 'Overview_0610_Top5.csv', '2006_2010')

# Process 0610 all
process_files(all_files, 'Overview_0610_All.csv', '2006_2010')

# Process combined top5 (both periods)
dfs_top5 = []
for period in ['2000_2005', '2006_2010']:
    dfs = []
    for file in top5_files:
        file_path = os.path.join(result_dir, file)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            filtered_df = df[(df['Model'] == 'EQI') & (df['EQI_Period'] == period)]
            dfs.append(filtered_df)
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        dfs_top5.append(combined_df)

if dfs_top5:
    combined_top5 = pd.concat(dfs_top5, ignore_index=True)
    combined_top5['Chinese_Name'] = combined_top5['ICD_Code'].map(icd_mapping).fillna('未知')
    combined_top5 = combined_top5.drop(columns=['AAMR_Period', 'Model'])
    combined_top5['Lag'] = combined_top5['Lag'].astype(int)
    combined_top5 = combined_top5.rename(columns={'Chinese_Name': 'Outcome'})
    combined_top5 = combined_top5[['ICD_Code', 'Outcome', 'EQI_Period', 'Lag', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']]
    # Sort by custom order for Top5
    order_mapping = {'C00_C97': 0, 'C34': 1, 'C18_C21': 2, 'C50': 3, 'C61': 4, 'C15_C26': 5}
    combined_top5['order'] = combined_top5['ICD_Code'].map(order_mapping)
    combined_top5 = combined_top5.sort_values(by=['order', 'EQI_Period', 'Lag'])
    combined_top5 = combined_top5.drop(columns=['order'])
    output_path = os.path.join(output_dir, 'Overview_Top5.csv')
    combined_top5.to_csv(output_path, index=False)
    print(f"Extracted data saved to {output_path}")

# Process combined all (both periods)
dfs_all = []
for period in ['2000_2005', '2006_2010']:
    dfs = []
    for file in all_files:
        file_path = os.path.join(result_dir, file)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            filtered_df = df[(df['Model'] == 'EQI') & (df['EQI_Period'] == period)]
            dfs.append(filtered_df)
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        dfs_all.append(combined_df)

if dfs_all:
    combined_all = pd.concat(dfs_all, ignore_index=True)
    combined_all['Chinese_Name'] = combined_all['ICD_Code'].map(icd_mapping).fillna('未知')
    combined_all = combined_all.drop(columns=['AAMR_Period', 'Model'])
    combined_all['Lag'] = combined_all['Lag'].astype(int)
    combined_all = combined_all.rename(columns={'Chinese_Name': 'Outcome'})
    combined_all = combined_all[['ICD_Code', 'Outcome', 'EQI_Period', 'Lag', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']]
    combined_all = combined_all.sort_values(by=['ICD_Code', 'EQI_Period', 'Lag'])
    output_path = os.path.join(output_dir, 'Overview_All.csv')
    combined_all.to_csv(output_path, index=False)
    print(f"Extracted data saved to {output_path}")