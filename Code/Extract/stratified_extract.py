import os

import pandas as pd
import yaml

# Define the base directory
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Load config.yaml for icd_mapping (though not used here, but in case)
config_path = os.path.join(base_dir, "config.yaml")
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

icd_mapping = config["brms_analysis"]["icd_mapping"]

# Define directories
stratified_dir = os.path.join(base_dir, "Result", "brms_stratified")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

# Strata mapping
strata_mapping = {
    "Male": "Sex(Male)",
    "Female": "Sex(Female)",
    "Asian": "Race(Asian)",
    "Black": "Race(Black)",
    "White": "Race(White)",
    "Indian": "Race(Indian)",
}


# Function to parse ICD_Code
def parse_icd_strata(icd_str):
    if icd_str.startswith("C00_C97_"):
        icd_code = "C00_C97"
        strata_suffix = icd_str[len("C00_C97_") :]
    elif icd_str.startswith("NDD_"):
        icd_code = "NDD"
        strata_suffix = icd_str[len("NDD_") :]
    else:
        raise ValueError(f"Unexpected ICD_Code format: {icd_str}")
    strata = strata_mapping.get(strata_suffix, strata_suffix)
    return icd_code, strata


# Strata order for sorting
strata_order = [
    "Sex(Male)",
    "Sex(Female)",
    "Race(White)",
    "Race(Black)",
    "Race(Asian)",
    "Race(Indian)",
]


# Model order for sorting
model_order = [
    "Stratified_EQI",
    "Stratified_Air",
    "Stratified_Water",
    "Stratified_Land",
    "Stratified_Built",
    "Stratified_Social",
]


def strata_sort_key(strata):
    return strata_order.index(strata) if strata in strata_order else len(strata_order)


def model_sort_key(model):
    return model_order.index(model) if model in model_order else len(model_order)


# Collect data
cancer_data = []
ndd_data = []

# Get all csv files
files = [f for f in os.listdir(stratified_dir) if f.endswith("_cmdstan.csv")]

for file in files:
    file_path = os.path.join(stratified_dir, file)
    df = pd.read_csv(file_path)
    for _, row in df.iterrows():
        icd_code, strata = parse_icd_strata(row["ICD_Code"])
        new_row = {
            "Strata": strata,
            "Model": row["Model"],
            "EQI_Period": row["EQI_Period"],
            "AAMR_Period": row["AAMR_Period"],
            "Lag": row["Lag"],
            "Q1": row["Q1"],
            "Q2": row["Q2"],
            "Q3": row["Q3"],
            "Q4": row["Q4"],
            "Q5": row["Q5"],
        }
        if icd_code == "C00_C97":
            cancer_data.append(new_row)
        elif icd_code == "NDD":
            ndd_data.append(new_row)

# Create DataFrames
cancer_df = pd.DataFrame(cancer_data)
ndd_df = pd.DataFrame(ndd_data)


# Filter data: keep only 2000_2005 EQI_Period and second occurrence of each model combination
def filter_second_occurrence(df):
    """Keep only rows where EQI_Period is 2000_2005 and take the second occurrence of each model combination"""
    # Filter to keep only 2000_2005 EQI period
    df = df[df["EQI_Period"] == "2000_2005"].copy()

    # Group by (Strata, Model, EQI_Period, AAMR_Period) and keep only the second occurrence
    filtered_rows = []
    for (strata, model, eqi_period, aamr_period), group in df.groupby(
        ["Strata", "Model", "EQI_Period", "AAMR_Period"]
    ):
        if len(group) >= 2:
            # Keep the second occurrence
            filtered_rows.append(group.iloc[1])
        elif len(group) == 1:
            # If there's only one occurrence, skip it (we want the second one)
            continue

    return pd.DataFrame(filtered_rows)


# Function to format results with significance stars
def format_result(q_val):
    """Format Q values as estimate(CI_lower,CI_upper) with significance stars"""
    if pd.isna(q_val) or q_val == 0.0:
        return "0.0"

    # Parse the value - assuming format like "12.02 (9.82, 14.29)***"
    val_str = str(q_val)

    # If already formatted with parentheses, return as is
    if "(" in val_str:
        return val_str

    # Otherwise just return the numeric value
    return val_str


# Function to create formatted output
def create_formatted_output(df, output_path):
    """Create grouped output by strata with custom formatting"""

    # Filter to keep only second occurrences
    df = filter_second_occurrence(df)

    if len(df) == 0:
        print(f"Warning: No data to write to {output_path}")
        return

    # Add sort keys
    df["strata_order"] = df["Strata"].apply(strata_sort_key)
    df["model_order"] = df["Model"].apply(model_sort_key)

    # Sort by strata first, then model order, then AAMR_Period
    df = df.sort_values(by=["strata_order", "model_order", "AAMR_Period", "Lag"])

    # Open file for writing
    with open(output_path, "w") as f:
        # Write header
        f.write("Model\tLag\tQ1\tQ2\tQ3\tQ4\tQ5\n")

        # Group by strata
        for strata in strata_order:
            strata_df = df[df["Strata"] == strata]
            if len(strata_df) == 0:
                continue

            # Write strata header
            f.write(f"{strata}\n")

            # Write each row
            for _, row in strata_df.iterrows():
                # Remove "Stratified_" prefix from model name
                model_name = row["Model"].replace("Stratified_", "")

                line = (
                    f"{model_name}\t"
                    f"{row['Lag']}\t"
                    f"{format_result(row['Q1'])}\t"
                    f"{format_result(row['Q2'])}\t"
                    f"{format_result(row['Q3'])}\t"
                    f"{format_result(row['Q4'])}\t"
                    f"{format_result(row['Q5'])}\n"
                )
                f.write(line)


# Output formatted files
cancer_output_path = os.path.join(output_dir, "brms_stratified_Cancer(2000-2005).csv")
create_formatted_output(cancer_df, cancer_output_path)

ndd_output_path = os.path.join(output_dir, "brms_stratified_NDD(2000-2005).csv")
create_formatted_output(ndd_df, ndd_output_path)

print(f"Cancer stratified data saved to {cancer_output_path}")
print(f"NDD stratified data saved to {ndd_output_path}")
