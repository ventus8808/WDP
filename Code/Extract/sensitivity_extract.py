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

# Define directories
sensitivity_dir = os.path.join(base_dir, "Result", "brms_Sensativity")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

# Prefixes and their ICD codes
PREFIXES = {
    "NDD": "G20_G30_G12.2_F01_F03",
    "Cancer": "C00_C97",
}

MODEL_GROUPS = ["Baseline_Full", "Insurance", "Doctor", "Forest", "Monitoring"]

# Sensitivity model display order
SENSITIVITY_MODEL_ORDER = [
    "Baseline",
    "Full",
    "+Insurance",
    "+Doctor",
    "+Forest",
    "+Monitoring",
]


def load_and_combine(file_pattern, prefixes, model_groups):
    """Load split CSV files and combine them."""
    dfs = []
    for prefix in prefixes:
        for group in model_groups:
            fname = f"{prefix}_{file_pattern}_{group}.csv"
            fpath = os.path.join(sensitivity_dir, fname)
            if os.path.exists(fpath):
                df = pd.read_csv(fpath)
                if not df.empty:
                    dfs.append(df)
            else:
                print(f"Warning: {fname} not found")
    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def process_eqi_table(df, output_name):
    """Process EQI quintile results into output table."""
    if df.empty:
        print(f"No data for {output_name}")
        return

    df["Outcome"] = df["ICD_Code"].map(icd_mapping).fillna("Unknown")
    df["Lag"] = df["Lag"].astype(int)

    # Order by sensitivity model
    df["model_order"] = df["Sensitivity_Model"].map(
        {m: i for i, m in enumerate(SENSITIVITY_MODEL_ORDER)}
    )

    # Order by ICD code (NDD first, then Cancer)
    icd_order = {v: i for i, v in enumerate(PREFIXES.values())}
    df["icd_order"] = df["ICD_Code"].map(icd_order)

    df = df.sort_values(by=["icd_order", "AAMR_Period", "Lag", "model_order"])
    df = df.drop(columns=["model_order", "icd_order"])

    # Select output columns
    cols = ["ICD_Code", "Outcome", "AAMR_Period", "Lag", "Sensitivity_Model",
            "Q1", "Q2", "Q3", "Q4", "Q5"]
    df = df[cols]

    output_path = os.path.join(output_dir, output_name)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path} ({len(df)} rows)")


def process_covar_table(df, output_name):
    """Process covariate coefficient results into output table."""
    if df.empty:
        print(f"No data for {output_name}")
        return

    df["Outcome"] = df["ICD_Code"].map(icd_mapping).fillna("Unknown")
    df["Lag"] = df["Lag"].astype(int)

    df["model_order"] = df["Sensitivity_Model"].map(
        {m: i for i, m in enumerate(SENSITIVITY_MODEL_ORDER)}
    )
    icd_order = {v: i for i, v in enumerate(PREFIXES.values())}
    df["icd_order"] = df["ICD_Code"].map(icd_order)

    df = df.sort_values(by=["icd_order", "AAMR_Period", "Lag", "model_order"])
    df = df.drop(columns=["model_order", "icd_order"])

    cols = ["ICD_Code", "Outcome", "AAMR_Period", "Lag", "Sensitivity_Model",
            "uninsured_rate_2005", "physician_density_per100k",
            "forest_coverage", "site_number_mean"]
    # Only keep columns that exist
    cols = [c for c in cols if c in df.columns]
    df = df[cols]

    output_path = os.path.join(output_dir, output_name)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path} ({len(df)} rows)")


# --- Main ---
prefixes = list(PREFIXES.keys())

# Combine EQI results
eqi_df = load_and_combine("sensitivity_EQI", prefixes, MODEL_GROUPS)
process_eqi_table(eqi_df, "brms_Sensitivity_EQI.csv")

# Combine Covar results
covar_df = load_and_combine("sensitivity_Covar", prefixes, MODEL_GROUPS)
process_covar_table(covar_df, "brms_Sensitivity_Covar.csv")
