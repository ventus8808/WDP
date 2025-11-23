import os

import pandas as pd
import yaml

# Define the base directory
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load config.yaml
config_path = os.path.join(base_dir, "config.yaml")
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

icd_mapping = config["brms_analysis"]["icd_mapping"]

# Define the result directory
result_dir = os.path.join(base_dir, "Result", "brms")

# Define the output directory
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

# List of files to read for top5 cancer
top5_cancer_files = [
    "C00_C97_main.csv",
    "C34_main.csv",
    "C18_C21_main.csv",
    "C50_main.csv",
    "C25_main.csv",
    "C61_main.csv",
]

# List of files to read for top5 NDD
top5_ndd_files = [
    "G20_G30_G12.2_F01_F03_main.csv",
    "G20_main.csv",
    "G30_main.csv",
    "G12.2_main.csv",
    "F01_main.csv",
    "F03_main.csv",
]

# Get all main files
all_files = [f for f in os.listdir(result_dir) if f.endswith("_main.csv")]

# icd_mapping loaded from config.yaml


# Function to process files for combined periods
def process_combined_files(files, output_name, order_mapping=None):
    dfs = []
    for period in ["2000_2005", "2006_2010"]:
        for file in files:
            file_path = os.path.join(result_dir, file)
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                # Filter rows where Model is 'EQI' and EQI_Period is the specified period
                filtered_df = df[(df["Model"] == "EQI") & (df["EQI_Period"] == period)]
                dfs.append(filtered_df)
            else:
                print(f"File {file_path} does not exist.")

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        # Add Outcome column
        combined_df["Outcome"] = (
            combined_df["ICD_Code"].map(icd_mapping).fillna("Unknown")
        )
        # Drop AAMR_Period and Model columns (keep EQI_Period)
        combined_df = combined_df.drop(columns=["AAMR_Period", "Model"])
        # Convert Lag to int
        combined_df["Lag"] = combined_df["Lag"].astype(int)
        # Reorder columns: ICD_Code, Outcome, EQI_Period, Lag, Q1, Q2, Q3, Q4, Q5
        combined_df = combined_df[
            ["ICD_Code", "Outcome", "EQI_Period", "Lag", "Q1", "Q2", "Q3", "Q4", "Q5"]
        ]
        # Sort
        if order_mapping:
            combined_df["order"] = combined_df["ICD_Code"].map(order_mapping)
            combined_df = combined_df.sort_values(by=["order", "EQI_Period", "Lag"])
            combined_df = combined_df.drop(columns=["order"])
        else:
            combined_df = combined_df.sort_values(by=["ICD_Code", "EQI_Period", "Lag"])
        # Save to CSV
        output_path = os.path.join(output_dir, output_name)
        combined_df.to_csv(output_path, index=False)
        print(f"Extracted data saved to {output_path}")
    else:
        print(f"No data to save for {output_name}.")


# Process brms_Main.csv (all results)
process_combined_files(all_files, "brms_Main.csv")

# Process brms_Top5_Cancer.csv
cancer_order = {"C00_C97": 0, "C34": 1, "C18_C21": 2, "C50": 3, "C25": 4, "C61": 5}
process_combined_files(top5_cancer_files, "brms_Top5_Cancer.csv", cancer_order)

# Process brms_Top5_NDD.csv
ndd_order = {
    "G20_G30_G12.2_F01_F03": 0,
    "G20": 1,
    "G30": 2,
    "G12.2": 3,
    "F01": 4,
    "F03": 5,
}
process_combined_files(top5_ndd_files, "brms_Top5_NDD.csv", ndd_order)
