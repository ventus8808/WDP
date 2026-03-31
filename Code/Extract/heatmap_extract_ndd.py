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

# Define NDD mapping
ndd_mapping = {
    "G10": "HD",
    "G20": "PD",
    "G30": "AD",
    "G12.2": "ALS",
    "F01": "VD",
    "G30_F01_F03": "Dementia",
    "G20_G30_G12.2_F01_F03": "NDD",
}

# Define the result directory
result_dir = base_dir / "Result" / "brms"

# Define the output directory
output_dir = base_dir / "Result" / "brms_heatmap"
output_dir.mkdir(parents=True, exist_ok=True)

# List of NDD codes
ndd_codes = ["G20_G30_G12.2_F01_F03", "G30_F01_F03", "G10", "G20", "G30", "G12.2", "F01"]

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
    for ndd_code in ndd_codes:
        # New file naming convention: {NDD_CODE}_main.csv
        file_path = result_dir / f"{ndd_code}_main.csv"

        if file_path.exists():
            df = pd.read_csv(file_path)
            # Filter rows where Model is the specified model and EQI_Period is the specified period
            filtered_df = df[(df["Model"] == model) & (df["EQI_Period"] == period)]
            dfs.append(filtered_df)
        else:
            print(f"File {file_path} does not exist.")

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        # Add outcome name column using NDD mapping
        combined_df["Outcome"] = (
            combined_df["ICD_Code"].map(ndd_mapping).fillna("Unknown")
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
        # Sort by custom order for NDD
        order_mapping = {
            "G20_G30_G12.2_F01_F03": 0,
            "G30_F01_F03": 1,
            "G10": 2,
            "G20": 3,
            "G30": 4,
            "G12.2": 5,
            "F01": 6,
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

model_short = {
    "EQI": "EQI",
    "EQI_Air": "Air",
    "EQI_Water": "Water",
    "EQI_Land": "Land",
    "EQI_Built": "Built",
    "EQI_Social": "Social",
}

if all_dfs:
    final_df = pd.concat(all_dfs, ignore_index=True)
    output_path = output_dir / "NDD.csv"
    final_df.to_csv(output_path, index=False)
    print(f"All extracted data saved to {output_path}")

    # Build formatted version: disease name as section header rows
    # Only 2000_2005 period; Model and Lag as separate columns; Q1 fixed at 0
    outcome_order = ["NDD", "Dementia", "AD", "PD", "VD", "ALS", "HD"]
    base_df = final_df[final_df["EQI_Period"] == "2000_2005"].copy()
    formatted_rows = []
    for outcome in outcome_order:
        subset = base_df[base_df["Outcome"] == outcome].copy()
        if subset.empty:
            continue
        formatted_rows.append({"Model": outcome})
        for _, row in subset.iterrows():
            formatted_rows.append({
                "Model": model_short.get(row["Model"], row["Model"]),
                "Lag": str(int(row["Lag"])),
                "Q1": "0",
                "Q2": row["Q2"],
                "Q3": row["Q3"],
                "Q4": row["Q4"],
                "Q5": row["Q5"],
            })

    formatted_df = pd.DataFrame(formatted_rows)
    fmt_path = output_dir / "NDD_formatted.csv"
    formatted_df.to_csv(fmt_path, index=False)
    print(f"Formatted data saved to {fmt_path}")
else:
    print("No data to save.")
