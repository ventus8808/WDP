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


# Function to format results (strip significance stars)
def format_result(q_val):
    """Format Q values as estimate(CI_lower,CI_upper), no significance stars."""
    if pd.isna(q_val) or q_val == 0.0:
        return "0.0"
    val_str = str(q_val).rstrip("*")
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
cancer_output_path = os.path.join(output_dir, "brms_stratified_Cancer_MRD(2000-2005).csv")
create_formatted_output(cancer_df, cancer_output_path)

ndd_output_path = os.path.join(output_dir, "brms_stratified_NDD_MRD(2000-2005).csv")
create_formatted_output(ndd_df, ndd_output_path)

print(f"Cancer stratified MRD saved to {cancer_output_path}")
print(f"NDD stratified MRD saved to {ndd_output_path}")

# ============================================================================
# RUCC Stratification (from brms directory)
# ============================================================================

rucc_dir = os.path.join(base_dir, "Result", "brms")

# RUCC strata mapping
rucc_strata_mapping = {
    "RUCC1": "RUCC: Metropolitan urbanized",
    "RUCC2": "RUCC: Non-metropolitan urbanized",
    "RUCC3": "RUCC: Less urbanized",
    "RUCC4": "RUCC: Thinly populated",
}

# RUCC model order
rucc_model_order = [
    "EQI",
    "Air",
    "Water",
    "Land",
    "Built",
    "Social",
]

rucc_strata_order = [
    "RUCC: Metropolitan urbanized",
    "RUCC: Non-metropolitan urbanized",
    "RUCC: Less urbanized",
    "RUCC: Thinly populated",
]


def rucc_strata_sort_key(strata):
    return (
        rucc_strata_order.index(strata)
        if strata in rucc_strata_order
        else len(rucc_strata_order)
    )


def rucc_model_sort_key(model):
    return (
        rucc_model_order.index(model)
        if model in rucc_model_order
        else len(rucc_model_order)
    )


def parse_rucc_model(model_str):
    """Parse RUCC model string to extract strata and model name"""
    # Models look like: RUCC1_EQI, RUCC1_EQI_Air, RUCC2_EQI_Water, etc.
    for rucc_code, strata_name in rucc_strata_mapping.items():
        if model_str.startswith(rucc_code + "_"):
            suffix = model_str[len(rucc_code) + 1 :]  # Remove "RUCC1_" prefix
            # suffix is like "EQI", "EQI_Air", "EQI_Water", etc.
            if suffix == "EQI":
                model_name = "EQI"
            elif suffix.startswith("EQI_"):
                model_name = suffix[
                    4:
                ]  # Remove "EQI_" prefix to get "Air", "Water", etc.
            else:
                continue
            return strata_name, model_name
    return None, None


def filter_rucc_data(df, icd_code):
    """Filter RUCC data for specific ICD code"""
    # Filter to keep only 2000_2005 EQI period and specific ICD code
    df = df[(df["EQI_Period"] == "2000_2005") & (df["ICD_Code"] == icd_code)].copy()

    # Filter for RUCC models only
    rucc_rows = []
    for _, row in df.iterrows():
        model = row["Model"]
        strata, model_name = parse_rucc_model(model)
        if strata is not None:
            new_row = {
                "Strata": strata,
                "Model": model_name,
                "EQI_Period": row["EQI_Period"],
                "AAMR_Period": row["AAMR_Period"],
                "Lag": row["Lag"],
                "Q1": row["Q1"],
                "Q2": row["Q2"],
                "Q3": row["Q3"],
                "Q4": row["Q4"],
                "Q5": row["Q5"],
            }
            rucc_rows.append(new_row)

    return pd.DataFrame(rucc_rows)


def create_rucc_formatted_output(df, output_path):
    """Create grouped output by RUCC strata with custom formatting"""
    if len(df) == 0:
        print(f"Warning: No RUCC data to write to {output_path}")
        return

    # Add sort keys
    df["strata_order"] = df["Strata"].apply(rucc_strata_sort_key)
    df["model_order"] = df["Model"].apply(rucc_model_sort_key)

    # Sort by strata first, then model order, then AAMR_Period
    df = df.sort_values(by=["strata_order", "model_order", "AAMR_Period", "Lag"])

    # Append to file
    with open(output_path, "a") as f:
        # Group by strata
        for strata in rucc_strata_order:
            strata_df = df[df["Strata"] == strata]
            if len(strata_df) == 0:
                continue

            # Write strata header
            f.write(f"{strata}\n")

            # Write each row
            for _, row in strata_df.iterrows():
                line = (
                    f"{row['Model']}\t"
                    f"{row['Lag']}\t"
                    f"{format_result(row['Q1'])}\t"
                    f"{format_result(row['Q2'])}\t"
                    f"{format_result(row['Q3'])}\t"
                    f"{format_result(row['Q4'])}\t"
                    f"{format_result(row['Q5'])}\n"
                )
                f.write(line)


# Process RUCC data for Cancer (C00_C97)
cancer_rucc_file = os.path.join(rucc_dir, "C00_C97_main.csv")
if os.path.exists(cancer_rucc_file):
    rucc_df = pd.read_csv(cancer_rucc_file)
    cancer_rucc_df = filter_rucc_data(rucc_df, "C00_C97")
    create_rucc_formatted_output(cancer_rucc_df, cancer_output_path)
    print(f"Cancer RUCC data appended to {cancer_output_path}")

# Process RUCC data for NDD (G20_G30_G12.2_F01_F03)
ndd_rucc_file = os.path.join(rucc_dir, "G20_G30_G12.2_F01_F03_main.csv")
if os.path.exists(ndd_rucc_file):
    rucc_df = pd.read_csv(ndd_rucc_file)
    ndd_rucc_df = filter_rucc_data(rucc_df, "G20_G30_G12.2_F01_F03")
    create_rucc_formatted_output(ndd_rucc_df, ndd_output_path)
    print(f"NDD RUCC data appended to {ndd_output_path}")

# ============================================================================
# Census Region Stratification (from brms_Climate directory)
# ============================================================================

climate_dir = os.path.join(base_dir, "Result", "brms_Climate")

# Census Region strata mapping
census_region_strata_mapping = {
    "census_region_1": "Census Region: Northeast",
    "census_region_2": "Census Region: Midwest",
    "census_region_3": "Census Region: South",
    "census_region_4": "Census Region: West",
}

census_region_strata_order = [
    "Census Region: Northeast",
    "Census Region: Midwest",
    "Census Region: South",
    "Census Region: West",
]


def census_region_strata_sort_key(strata):
    return (
        census_region_strata_order.index(strata)
        if strata in census_region_strata_order
        else len(census_region_strata_order)
    )


def parse_census_region_model(model_str):
    """Parse Census Region model string to extract strata and model name"""
    # Models look like: census_region_1_EQI, census_region_1_EQI_Air, etc.
    for region_code, strata_name in census_region_strata_mapping.items():
        if model_str.startswith(region_code + "_"):
            suffix = model_str[
                len(region_code) + 1 :
            ]  # Remove "census_region_1_" prefix
            # suffix is like "EQI", "EQI_Air", "EQI_Water", etc.
            if suffix == "EQI":
                model_name = "EQI"
            elif suffix.startswith("EQI_"):
                model_name = suffix[
                    4:
                ]  # Remove "EQI_" prefix to get "Air", "Water", etc.
            else:
                continue
            return strata_name, model_name
    return None, None


def filter_census_region_data(df, icd_code):
    """Filter Census Region data for specific ICD code"""
    # Filter to keep only 2000_2005 EQI period and specific ICD code
    df = df[(df["EQI_Period"] == "2000_2005") & (df["ICD_Code"] == icd_code)].copy()

    # Filter for Census Region models only
    census_rows = []
    for _, row in df.iterrows():
        model = row["Model"]
        strata, model_name = parse_census_region_model(model)
        if strata is not None:
            new_row = {
                "Strata": strata,
                "Model": model_name,
                "EQI_Period": row["EQI_Period"],
                "AAMR_Period": row["AAMR_Period"],
                "Lag": row["Lag"],
                "Q1": row["Q1"],
                "Q2": row["Q2"],
                "Q3": row["Q3"],
                "Q4": row["Q4"],
                "Q5": row["Q5"],
            }
            census_rows.append(new_row)

    return pd.DataFrame(census_rows)


def create_census_region_formatted_output(df, output_path):
    """Create grouped output by Census Region strata with custom formatting"""
    if len(df) == 0:
        print(f"Warning: No Census Region data to write to {output_path}")
        return

    # Add sort keys
    df["strata_order"] = df["Strata"].apply(census_region_strata_sort_key)
    df["model_order"] = df["Model"].apply(
        rucc_model_sort_key
    )  # Reuse rucc_model_sort_key

    # Sort by strata first, then model order, then AAMR_Period
    df = df.sort_values(by=["strata_order", "model_order", "AAMR_Period", "Lag"])

    # Append to file
    with open(output_path, "a") as f:
        # Group by strata
        for strata in census_region_strata_order:
            strata_df = df[df["Strata"] == strata]
            if len(strata_df) == 0:
                continue

            # Write strata header
            f.write(f"{strata}\n")

            # Write each row
            for _, row in strata_df.iterrows():
                line = (
                    f"{row['Model']}\t"
                    f"{row['Lag']}\t"
                    f"{format_result(row['Q1'])}\t"
                    f"{format_result(row['Q2'])}\t"
                    f"{format_result(row['Q3'])}\t"
                    f"{format_result(row['Q4'])}\t"
                    f"{format_result(row['Q5'])}\n"
                )
                f.write(line)


# Process Census Region data for Cancer (C00_C97)
cancer_census_file = os.path.join(climate_dir, "C00_C97_census_region.csv")
if os.path.exists(cancer_census_file):
    census_df = pd.read_csv(cancer_census_file)
    cancer_census_df = filter_census_region_data(census_df, "C00_C97")
    create_census_region_formatted_output(cancer_census_df, cancer_output_path)
    print(f"Cancer Census Region data appended to {cancer_output_path}")

# Process Census Region data for NDD (G20_G30_G12.2_F01_F03)
ndd_census_file = os.path.join(climate_dir, "G20_G30_G12.2_F01_F03_census_region.csv")
if os.path.exists(ndd_census_file):
    census_df = pd.read_csv(ndd_census_file)
    ndd_census_df = filter_census_region_data(census_df, "G20_G30_G12.2_F01_F03")
    create_census_region_formatted_output(ndd_census_df, ndd_output_path)
    print(f"NDD Census Region data appended to {ndd_output_path}")

# ============================================================================
# Köppen-Geiger Climate Zone Stratification (from brms_Climate directory)
# ============================================================================

# Köppen-Geiger strata mapping
koppen_strata_mapping = {
    "koppen_major_B": "Köppen-Geiger Climate Zone: Dry",
    "koppen_major_C": "Köppen-Geiger Climate Zone: Temperate",
    "koppen_major_D": "Köppen-Geiger Climate Zone: Continental",
}

koppen_strata_order = [
    "Köppen-Geiger Climate Zone: Dry",
    "Köppen-Geiger Climate Zone: Temperate",
    "Köppen-Geiger Climate Zone: Continental",
]


def koppen_strata_sort_key(strata):
    return (
        koppen_strata_order.index(strata)
        if strata in koppen_strata_order
        else len(koppen_strata_order)
    )


def parse_koppen_model(model_str):
    """Parse Köppen-Geiger model string to extract strata and model name"""
    # Models look like: koppen_major_B_EQI, koppen_major_B_EQI_Air, etc.
    for koppen_code, strata_name in koppen_strata_mapping.items():
        if model_str.startswith(koppen_code + "_"):
            suffix = model_str[
                len(koppen_code) + 1 :
            ]  # Remove "koppen_major_B_" prefix
            # suffix is like "EQI", "EQI_Air", "EQI_Water", etc.
            if suffix == "EQI":
                model_name = "EQI"
            elif suffix.startswith("EQI_"):
                model_name = suffix[
                    4:
                ]  # Remove "EQI_" prefix to get "Air", "Water", etc.
            else:
                continue
            return strata_name, model_name
    return None, None


def filter_koppen_data(df, icd_code):
    """Filter Köppen-Geiger data for specific ICD code"""
    # Filter to keep only 2000_2005 EQI period and specific ICD code
    df = df[(df["EQI_Period"] == "2000_2005") & (df["ICD_Code"] == icd_code)].copy()

    # Filter for Köppen-Geiger models only
    koppen_rows = []
    for _, row in df.iterrows():
        model = row["Model"]
        strata, model_name = parse_koppen_model(model)
        if strata is not None:
            new_row = {
                "Strata": strata,
                "Model": model_name,
                "EQI_Period": row["EQI_Period"],
                "AAMR_Period": row["AAMR_Period"],
                "Lag": row["Lag"],
                "Q1": row["Q1"],
                "Q2": row["Q2"],
                "Q3": row["Q3"],
                "Q4": row["Q4"],
                "Q5": row["Q5"],
            }
            koppen_rows.append(new_row)

    return pd.DataFrame(koppen_rows)


def create_koppen_formatted_output(df, output_path):
    """Create grouped output by Köppen-Geiger strata with custom formatting"""
    if len(df) == 0:
        print(f"Warning: No Köppen-Geiger data to write to {output_path}")
        return

    # Add sort keys
    df["strata_order"] = df["Strata"].apply(koppen_strata_sort_key)
    df["model_order"] = df["Model"].apply(
        rucc_model_sort_key
    )  # Reuse rucc_model_sort_key

    # Sort by strata first, then model order, then AAMR_Period
    df = df.sort_values(by=["strata_order", "model_order", "AAMR_Period", "Lag"])

    # Append to file
    with open(output_path, "a") as f:
        # Group by strata
        for strata in koppen_strata_order:
            strata_df = df[df["Strata"] == strata]
            if len(strata_df) == 0:
                continue

            # Write strata header
            f.write(f"{strata}\n")

            # Write each row
            for _, row in strata_df.iterrows():
                line = (
                    f"{row['Model']}\t"
                    f"{row['Lag']}\t"
                    f"{format_result(row['Q1'])}\t"
                    f"{format_result(row['Q2'])}\t"
                    f"{format_result(row['Q3'])}\t"
                    f"{format_result(row['Q4'])}\t"
                    f"{format_result(row['Q5'])}\n"
                )
                f.write(line)


# Process Köppen-Geiger data for Cancer (C00_C97)
cancer_koppen_file = os.path.join(climate_dir, "C00_C97_koppen_major.csv")
if os.path.exists(cancer_koppen_file):
    koppen_df = pd.read_csv(cancer_koppen_file)
    cancer_koppen_df = filter_koppen_data(koppen_df, "C00_C97")
    create_koppen_formatted_output(cancer_koppen_df, cancer_output_path)
    print(f"Cancer Köppen-Geiger data appended to {cancer_output_path}")

# Process Köppen-Geiger data for NDD (G20_G30_G12.2_F01_F03)
ndd_koppen_file = os.path.join(climate_dir, "G20_G30_G12.2_F01_F03_koppen_major.csv")
if os.path.exists(ndd_koppen_file):
    koppen_df = pd.read_csv(ndd_koppen_file)
    ndd_koppen_df = filter_koppen_data(koppen_df, "G20_G30_G12.2_F01_F03")
    create_koppen_formatted_output(ndd_koppen_df, ndd_output_path)
    print(f"NDD Köppen-Geiger data appended to {ndd_output_path}")

# ============================================================================
# County Economic Typology Stratification (from brms_Typology_LandUse directory)
# ============================================================================

typology_dir = os.path.join(base_dir, "Result", "brms_Typology_LandUse")

# Typology strata mapping
typology_strata_mapping = {
    "Typology_Farming": "County Economic Typology: Farming",
    "Typology_Mining": "County Economic Typology: Mining",
    "Typology_Manufacturing": "County Economic Typology: Manufacturing",
    "Typology_Government": "County Economic Typology: Government",
    "Typology_Services": "County Economic Typology: Services",
    "Typology_Nonspecialized": "County Economic Typology: Nonspecialized",
}

typology_strata_order = [
    "County Economic Typology: Farming",
    "County Economic Typology: Mining",
    "County Economic Typology: Manufacturing",
    "County Economic Typology: Government",
    "County Economic Typology: Services",
    "County Economic Typology: Nonspecialized",
]


def typology_strata_sort_key(strata):
    return (
        typology_strata_order.index(strata)
        if strata in typology_strata_order
        else len(typology_strata_order)
    )


def parse_typology_model(model_str):
    """Parse Typology model string to extract strata and model name"""
    # Models look like: Typology_Farming_EQI, Typology_Farming_EQI_Air, etc.
    for typology_code, strata_name in typology_strata_mapping.items():
        if model_str.startswith(typology_code + "_"):
            suffix = model_str[
                len(typology_code) + 1 :
            ]  # Remove "Typology_Farming_" prefix
            # suffix is like "EQI", "EQI_Air", "EQI_Water", etc.
            if suffix == "EQI":
                model_name = "EQI"
            elif suffix.startswith("EQI_"):
                model_name = suffix[
                    4:
                ]  # Remove "EQI_" prefix to get "Air", "Water", etc.
            else:
                continue
            return strata_name, model_name
    return None, None


def filter_typology_data(df, icd_code):
    """Filter Typology data for specific ICD code"""
    # Filter to keep only 2000_2005 EQI period and specific ICD code
    df = df[(df["EQI_Period"] == "2000_2005") & (df["ICD_Code"] == icd_code)].copy()

    # Filter for Typology models only
    typology_rows = []
    for _, row in df.iterrows():
        model = row["Model"]
        strata, model_name = parse_typology_model(model)
        if strata is not None:
            new_row = {
                "Strata": strata,
                "Model": model_name,
                "EQI_Period": row["EQI_Period"],
                "AAMR_Period": row["AAMR_Period"],
                "Lag": row["Lag"],
                "Q1": row["Q1"],
                "Q2": row["Q2"],
                "Q3": row["Q3"],
                "Q4": row["Q4"],
                "Q5": row["Q5"],
            }
            typology_rows.append(new_row)

    return pd.DataFrame(typology_rows)


def create_typology_formatted_output(df, output_path):
    """Create grouped output by Typology strata with custom formatting"""
    if len(df) == 0:
        print(f"Warning: No Typology data to write to {output_path}")
        return

    # Add sort keys
    df["strata_order"] = df["Strata"].apply(typology_strata_sort_key)
    df["model_order"] = df["Model"].apply(
        rucc_model_sort_key
    )  # Reuse rucc_model_sort_key

    # Sort by strata first, then model order, then AAMR_Period
    df = df.sort_values(by=["strata_order", "model_order", "AAMR_Period", "Lag"])

    # Append to file
    with open(output_path, "a") as f:
        # Group by strata
        for strata in typology_strata_order:
            strata_df = df[df["Strata"] == strata]
            if len(strata_df) == 0:
                continue

            # Write strata header
            f.write(f"{strata}\n")

            # Write each row
            for _, row in strata_df.iterrows():
                line = (
                    f"{row['Model']}\t"
                    f"{row['Lag']}\t"
                    f"{format_result(row['Q1'])}\t"
                    f"{format_result(row['Q2'])}\t"
                    f"{format_result(row['Q3'])}\t"
                    f"{format_result(row['Q4'])}\t"
                    f"{format_result(row['Q5'])}\n"
                )
                f.write(line)


# Process Typology data for NDD (G20_G30_G12.2_F01_F03) only
ndd_typology_file = os.path.join(typology_dir, "G20_G30_G12.2_F01_F03_Typology.csv")
if os.path.exists(ndd_typology_file):
    typology_df = pd.read_csv(ndd_typology_file)
    ndd_typology_df = filter_typology_data(typology_df, "G20_G30_G12.2_F01_F03")
    create_typology_formatted_output(ndd_typology_df, ndd_output_path)
    print(f"NDD Typology data appended to {ndd_output_path}")

# ============================================================================
# Stratified MRR — 2 combined output files (Cancer / NDD)
# Sources: Result/brms_Stratified_Typo_LandUse_MRR/
# Strata order: Sex → Race → RUCC → Census Region → Köppen → Typology
# Format: mean(lower,upper) 2 dp, no asterisks, Q1 = "1.00", EQI 2000-2005 only
# ============================================================================

mrr_dir = os.path.join(base_dir, "Result", "brms_Stratified_Typo_LandUse_MRR")

# Targets: (base ICD code, output file name)
mrr_targets = [
    ("C00_C97",                 "brms_stratified_Cancer_MRR(2000-2005).csv"),
    ("G20_G30_G12.2_F01_F03",  "brms_stratified_NDD_MRR(2000-2005).csv"),
]

# Ordered strata sections: (display label, Model-column prefix in MRR CSV)
# Each entry: (section_header, {model_value: display_label}, [ordered display labels])
# Sex/race: stratum is in ICD_Code suffix, model is domain (matches MRD format)
# e.g. ICD_Code="NDD_Asian", Model="Stratified_EQI"
mrr_sex_race_icd_map = {
    "Male":   "Sex(Male)",
    "Female": "Sex(Female)",
    "White":  "Race(White)",
    "Black":  "Race(Black)",
    "Asian":  "Race(Asian)",
    "Indian": "Race(Indian)",
}
mrr_sex_race_order = [
    "Sex(Male)", "Sex(Female)",
    "Race(White)", "Race(Black)", "Race(Asian)", "Race(Indian)",
]

# Typology: stratum is in Model prefix (e.g. "Typology_Farming_EQI"), ICD_Code = base cancer
mrr_typology_strata_prefixes = [
    "Typology_Farming",
    "Typology_Mining",
    "Typology_Manufacturing",
    "Typology_Government",
    "Typology_Services",
    "Typology_Nonspecialized",
]
mrr_typology_strata_labels = [
    "County Economic Typology: Farming",
    "County Economic Typology: Mining",
    "County Economic Typology: Manufacturing",
    "County Economic Typology: Government",
    "County Economic Typology: Services",
    "County Economic Typology: Nonspecialized",
]
mrr_typology_order = mrr_typology_strata_labels

# Domain model display order (suffix in Model column after stratum prefix)
mrr_domain_order = ["EQI", "Air", "Water", "Land", "Built", "Social"]


def format_mrr_cell(row):
    """Format MRR as mean(lower,upper) with 2 decimal places, no asterisks."""
    try:
        mean  = float(row["MRR_mean"])
        lower = float(row["MRR_lower"])
        upper = float(row["MRR_upper"])
        return f"{mean:.2f}({lower:.2f},{upper:.2f})"
    except (ValueError, TypeError):
        return ""


def load_mrr_file(filepath):
    """Load long-format MRR CSV, filter EQI 2000-2005, format cells, pivot wide."""
    if not os.path.exists(filepath):
        return pd.DataFrame()
    df = pd.read_csv(filepath)
    df = df[df["EQI_Period"] == "2000_2005"].copy()
    if df.empty:
        return pd.DataFrame()
    df["cell"] = df.apply(format_mrr_cell, axis=1)
    pivot = df.pivot_table(
        index=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "Model"],
        columns="Quintile",
        values="cell",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
        if q not in pivot.columns:
            pivot[q] = ""
    pivot["Q1"] = "1.00"
    return pivot


def write_sex_race_mrr(pivot_df, base_icd, fh):
    """
    Write sex/race MRR section.
    ICD_Code encodes stratum (e.g. NDD_Asian); Model encodes domain (Stratified_EQI, Stratified_Air…).
    Groups by stratum (strata_order), then within each stratum rows by domain order × Lag.
    """
    if pivot_df.empty:
        return

    rows = []
    for _, row in pivot_df.iterrows():
        icd = row["ICD_Code"]
        # Extract stratum suffix: everything after base_icd_ prefix
        suffix = icd[len(base_icd) + 1:] if icd.startswith(base_icd + "_") else None
        if suffix is None:
            continue
        strata_label = mrr_sex_race_icd_map.get(suffix)
        if strata_label is None:
            continue
        # Extract domain from Model: "Stratified_EQI" → "EQI", "Stratified_Air" → "Air"
        model_val = row["Model"]
        domain = model_val[len("Stratified_"):] if model_val.startswith("Stratified_") else None
        if domain not in mrr_domain_order:
            continue
        rows.append({
            "Strata": strata_label,
            "Domain": domain,
            "Lag":    row["Lag"],
            "Q1":     row.get("Q1", ""),
            "Q2":     row.get("Q2", ""),
            "Q3":     row.get("Q3", ""),
            "Q4":     row.get("Q4", ""),
            "Q5":     row.get("Q5", ""),
        })

    if not rows:
        return

    result_df = pd.DataFrame(rows)

    def _strata_key(s):
        return mrr_sex_race_order.index(s) if s in mrr_sex_race_order else len(mrr_sex_race_order)
    def _domain_key(d):
        return mrr_domain_order.index(d) if d in mrr_domain_order else len(mrr_domain_order)

    result_df["_sk"] = result_df["Strata"].apply(_strata_key)
    result_df["_dk"] = result_df["Domain"].apply(_domain_key)
    result_df = result_df.sort_values(["_sk", "_dk", "Lag"])

    for strata in mrr_sex_race_order:
        sub = result_df[result_df["Strata"] == strata]
        if sub.empty:
            continue
        fh.write(f"{strata}\n")
        for _, r in sub.iterrows():
            fh.write(
                f"{r['Domain']}\t{r['Lag']}\t{r['Q1']}\t{r['Q2']}\t"
                f"{r['Q3']}\t{r['Q4']}\t{r['Q5']}\n"
            )


def write_typology_mrr(pivot_df, fh):
    """
    Write Typology MRR section.
    ICD_Code = base cancer; Model = "Typology_Farming_EQI", "Typology_Farming_Air" …
    Groups by stratum, then domain order × Lag.
    """
    if pivot_df.empty:
        return

    rows = []
    for _, row in pivot_df.iterrows():
        model_val = row["Model"]
        strata_label = None
        domain = None
        for prefix, label in zip(mrr_typology_strata_prefixes, mrr_typology_strata_labels):
            if model_val.startswith(prefix + "_"):
                domain = model_val[len(prefix) + 1:]  # e.g. "EQI", "Air"
                strata_label = label
                break
        if strata_label is None or domain not in mrr_domain_order:
            continue
        rows.append({
            "Strata": strata_label,
            "Domain": domain,
            "Lag":    row["Lag"],
            "Q1":     row.get("Q1", ""),
            "Q2":     row.get("Q2", ""),
            "Q3":     row.get("Q3", ""),
            "Q4":     row.get("Q4", ""),
            "Q5":     row.get("Q5", ""),
        })

    if not rows:
        return

    result_df = pd.DataFrame(rows)

    def _strata_key(s):
        return mrr_typology_order.index(s) if s in mrr_typology_order else len(mrr_typology_order)
    def _domain_key(d):
        return mrr_domain_order.index(d) if d in mrr_domain_order else len(mrr_domain_order)

    result_df["_sk"] = result_df["Strata"].apply(_strata_key)
    result_df["_dk"] = result_df["Domain"].apply(_domain_key)
    result_df = result_df.sort_values(["_sk", "_dk", "Lag"])

    for strata in mrr_typology_order:
        sub = result_df[result_df["Strata"] == strata]
        if sub.empty:
            continue
        fh.write(f"{strata}\n")
        for _, r in sub.iterrows():
            fh.write(
                f"{r['Domain']}\t{r['Lag']}\t{r['Q1']}\t{r['Q2']}\t"
                f"{r['Q3']}\t{r['Q4']}\t{r['Q5']}\n"
            )


for cancer_icd, out_filename in mrr_targets:
    out_path = os.path.join(output_dir, out_filename)
    with open(out_path, "w") as fh:
        fh.write("Model\tLag\tQ1\tQ2\tQ3\tQ4\tQ5\n")

        # 1. Sex/race strata
        strat_pivot = load_mrr_file(os.path.join(mrr_dir, f"{cancer_icd}_Stratified_MRR.csv"))
        write_sex_race_mrr(strat_pivot, cancer_icd, fh)

        # 2. Typology strata
        typo_pivot = load_mrr_file(os.path.join(mrr_dir, f"{cancer_icd}_Typology_MRR.csv"))
        write_typology_mrr(typo_pivot, fh)

    print(f"MRR table written: {out_path}")
