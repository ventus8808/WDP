import pandas as pd
import os

# Define the base directory
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Define the result directory
result_dir = os.path.join(base_dir, 'Result', 'brms')

# Define the output directory
output_dir = os.path.join(base_dir, 'Result', 'brms_heatmap')
os.makedirs(output_dir, exist_ok=True)

# List of files to read for top5
top5_files = ['C00_C97_brms.csv', 'C34_brms.csv', 'C18_C21_brms.csv', 'C50_brms.csv', 'C25_brms.csv', 'C61_brms.csv']

# Models to extract
models = ['EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social', 'RUCC1_EQI', 'RUCC2_EQI', 'RUCC3_EQI', 'RUCC4_EQI', 'RUCC5_EQI', 'RUCC6_EQI', 'RUCC7_EQI', 'RUCC8_EQI', 'RUCC9_EQI']

# Periods
periods = ['2000_2005', '2006_2010']

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

# Function to process files for a specific model
def process_model(model, period='2000_2005'):
    dfs = []
    for file in top5_files:
        file_path = os.path.join(result_dir, file)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            # Filter rows where Model is the specified model and EQI_Period is the specified period
            filtered_df = df[(df['Model'] == model) & (df['EQI_Period'] == period)]
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
        # Rename Chinese_Name to Outcome
        combined_df = combined_df.rename(columns={'Chinese_Name': 'Outcome'})
        # Add Model column
        combined_df['Model'] = model
        # Reorder columns: Model, ICD_Code, Outcome, EQI_Period, Lag, Q1, Q2, Q3, Q4, Q5
        combined_df = combined_df[['Model', 'ICD_Code', 'Outcome', 'EQI_Period', 'Lag', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']]
        # Sort by custom order for Top5
        order_mapping = {'C00_C97': 0, 'C34': 1, 'C18_C21': 2, 'C50': 3, 'C25': 4, 'C61': 5}
        combined_df['order'] = combined_df['ICD_Code'].map(order_mapping)
        combined_df = combined_df.sort_values(by=['order', 'EQI_Period', 'Lag'])
        combined_df = combined_df.drop(columns=['order'])
        return combined_df
    else:
        return None

# Collect all data
all_dfs = []
for model in models:
    for period in periods:
        df = process_model(model, period)
        if df is not None:
            all_dfs.append(df)

if all_dfs:
    final_df = pd.concat(all_dfs, ignore_index=True)
    output_path = os.path.join(output_dir, 'brms_heatmap.csv')
    final_df.to_csv(output_path, index=False)
    print(f"All extracted data saved to {output_path}")
    
    # Save separate files for each period
    for period in periods:
        df_period = final_df[final_df['EQI_Period'] == period]
        output_path_period = os.path.join(output_dir, f'brms_heatmap_{period}.csv')
        df_period.to_csv(output_path_period, index=False)
        print(f"Data for {period} saved to {output_path_period}")
else:
    print("No data to save.")