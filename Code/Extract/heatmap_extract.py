import os
from pathlib import Path

import pandas as pd
import yaml

# Define the base directory
base_dir = Path(__file__).resolve().parents[2]

# Load config.yaml
config_path = base_dir / "config.yaml"
with open(config_path, "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

# Get ICD mapping from config
icd_mapping = config.get("brms_analysis", {}).get("icd_mapping", {})

# Define the result directory
result_dir = base_dir / "Result" / "brms"

# Define the output directory
output_dir = base_dir / "Result" / "brms_heatmap"
output_dir.mkdir(parents=True, exist_ok=True)

# List of ICD codes for top5 (updated file naming convention: {ICD}_main.csv)
top5_icds = ["C00_C97", "C34", "C18_C21", "C50", "C25", "C61"]

# Models to extract (only main EQI models, no RUCC stratification)
models = [
    "EQI",
    "EQI_Air",
    "EQI_Water",
    "EQI_Land",
    "EQI_Built",
    "EQI_Social",
]

# Periods (EQI periods from config: 0005 = 2000_2005, 0610 = 2006_2010)
periods = ["2000_2005", "2006_2010"]


# Function to process files for a specific model
def process_model(model, period="2000_2005"):
    dfs = []
    for icd_code in top5_icds:
        # New file naming convention: {ICD}_main.csv
        file_path = result_dir / f"{icd_code}_main.csv"

        if file_path.exists():
            df = pd.read_csv(file_path)
            # Filter rows where Model is the specified model and EQI_Period is the specified period
            filtered_df = df[(df["Model"] == model) & (df["EQI_Period"] == period)]
            dfs.append(filtered_df)
        else:
            print(f"File {file_path} does not exist.")

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        # Add outcome name column using config mapping
        combined_df["Outcome"] = (
            combined_df["ICD_Code"].map(icd_mapping).fillna("Unknown")
        )
        # Drop AAMR_Period and Model columns (keep EQI_Period)
        combined_df = combined_df.drop(columns=["AAMR_Period", "Model"])
        # Convert Lag to int
        combined_df["Lag"] = combined_df["Lag"].astype(int)
        # Add Model column
        combined_df["Model"] = model
        # Reorder columns: Model, ICD_Code, Outcome, EQI_Period, Lag, Q1, Q2, Q3, Q4, Q5
        combined_df = combined_df[
            [
                "Model",
                "ICD_Code",
                "Outcome",
                "EQI_Period",
                "Lag",
                "Q1",
                "Q2",
                "Q3",
                "Q4",
                "Q5",
            ]
        ]
        # Sort by custom order for Top5
        order_mapping = {
            "C00_C97": 0,
            "C34": 1,
            "C18_C21": 2,
            "C50": 3,
            "C25": 4,
            "C61": 5,
        }
        combined_df["order"] = combined_df["ICD_Code"].map(order_mapping)
        combined_df = combined_df.sort_values(by=["order", "EQI_Period", "Lag"])
        combined_df = combined_df.drop(columns=["order"])
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
    output_path = output_dir / "Top5Cancer.csv"
    final_df.to_csv(output_path, index=False)
    print(f"All extracted data saved to {output_path}")
else:
    print("No data to save.")
