# Code/Clean/USDA_ERS_CountyTypology.py

import os
from pathlib import Path

import pandas as pd
import yaml


def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_county_typology_data(input_file):
    """Load USDA ERS County Typology data from Excel file"""
    print(f"Loading data from: {input_file}")

    try:
        # Try reading with xlrd engine for old .xls files
        df = pd.read_excel(input_file)
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        print("Please ensure xlrd is installed: pip install xlrd")
        raise

    print(f"Loaded {df.shape[0]} rows, {df.shape[1]} columns")
    return df


def clean_county_typology_data(df):
    """Clean and standardize county typology data"""
    print("\nCleaning data...")

    # Check required columns exist
    required_cols = ["FIPSTXT", "econdep"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Select only FIPSTXT and econdep columns
    df_clean = df[["FIPSTXT", "econdep"]].copy()

    # Rename FIPSTXT to COUNTY_FIPS
    df_clean = df_clean.rename(columns={"FIPSTXT": "COUNTY_FIPS"})

    # Convert FIPS to 5-digit zero-padded string
    df_clean["COUNTY_FIPS"] = df_clean["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Validate econdep values (should be 1-6)
    valid_econdep = [1, 2, 3, 4, 5, 6]
    invalid_econdep = df_clean[~df_clean["econdep"].isin(valid_econdep)]
    if len(invalid_econdep) > 0:
        print(f"Warning: Found {len(invalid_econdep)} rows with invalid econdep values")
        print(invalid_econdep)

    # Check for missing values
    missing_fips = df_clean["COUNTY_FIPS"].isnull().sum()
    missing_econdep = df_clean["econdep"].isnull().sum()

    if missing_fips > 0:
        print(f"Warning: {missing_fips} missing COUNTY_FIPS values")
    if missing_econdep > 0:
        print(f"Warning: {missing_econdep} missing econdep values")

    # Remove rows with missing values
    df_clean = df_clean.dropna()

    # Check for duplicate FIPS codes
    duplicates = df_clean[df_clean["COUNTY_FIPS"].duplicated()]
    if len(duplicates) > 0:
        print(f"Warning: Found {len(duplicates)} duplicate COUNTY_FIPS codes")
        print(duplicates)
        # Keep first occurrence
        df_clean = df_clean.drop_duplicates(subset=["COUNTY_FIPS"], keep="first")

    print(f"Cleaned data: {df_clean.shape[0]} counties")

    return df_clean


def generate_summary(df):
    """Generate summary statistics for county typology data"""
    print("\n" + "=" * 80)
    print("County Typology Summary Statistics")
    print("=" * 80)

    print(f"\nTotal counties: {len(df)}")

    # Economic dependency type distribution
    print("\nEconomic Dependency Type Distribution:")
    print("-" * 80)

    econdep_labels = {
        1: "Farming-dependent",
        2: "Mining-dependent",
        3: "Manufacturing-dependent",
        4: "Federal/State government-dependent",
        5: "Services-dependent",
        6: "Nonspecialized",
    }

    econdep_counts = df["econdep"].value_counts().sort_index()

    for econdep_code, label in econdep_labels.items():
        count = econdep_counts.get(econdep_code, 0)
        pct = (count / len(df)) * 100
        print(f"  {econdep_code}. {label:40s}: {count:4d} ({pct:5.1f}%)")

    print("\n" + "=" * 80)

    return econdep_counts


def save_output(df, output_file, summary_file):
    """Save cleaned data and summary"""
    print(f"\nSaving cleaned data to: {output_file}")
    df.to_csv(output_file, index=False)
    print(f"Saved {len(df)} counties")

    # Save summary statistics
    print(f"\nSaving summary statistics to: {summary_file}")

    econdep_labels = {
        1: "Farming-dependent",
        2: "Mining-dependent",
        3: "Manufacturing-dependent",
        4: "Federal/State government-dependent",
        5: "Services-dependent",
        6: "Nonspecialized",
    }

    with open(summary_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write("USDA ERS County Typology (2004) - Summary Statistics\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total counties: {len(df)}\n\n")

        f.write("Economic Dependency Type Distribution:\n")
        f.write("-" * 80 + "\n")

        econdep_counts = df["econdep"].value_counts().sort_index()

        for econdep_code, label in econdep_labels.items():
            count = econdep_counts.get(econdep_code, 0)
            pct = (count / len(df)) * 100
            f.write(f"  {econdep_code}. {label:40s}: {count:4d} ({pct:5.1f}%)\n")

        f.write("\n" + "=" * 80 + "\n\n")

        f.write("Variable Codebook:\n")
        f.write("-" * 80 + "\n")
        f.write("COUNTY_FIPS: 5-digit FIPS code (State + County)\n")
        f.write("econdep: Economic dependency type code (1-6)\n\n")

        f.write("Economic Dependency Type Codes:\n")
        for code, label in econdep_labels.items():
            f.write(f"  {code} = {label}\n")

        f.write("\n" + "=" * 80 + "\n")

    print("Summary statistics saved")


def main():
    """Main cleaning pipeline"""
    print("=" * 80)
    print("USDA ERS County Typology Data Cleaning")
    print("=" * 80)

    # Load configuration
    config = get_config()
    project_root = Path(__file__).resolve().parents[2]

    # File paths
    input_file = os.path.join(
        project_root,
        config["data_sources"]["socioeconomic"]["usda_ers"]["original"],
        "2004CountyTypologyCode.xls",
    )

    output_file = os.path.join(
        project_root,
        config["data_sources"]["socioeconomic"]["usda_ers"]["processed"],
        "County_Typology_2004.csv",
    )

    summary_file = os.path.join(
        project_root,
        config["data_sources"]["socioeconomic"]["usda_ers"]["processed"],
        "County_Typology_2004_Summary.txt",
    )

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Step 1: Load data
    df = load_county_typology_data(input_file)

    # Step 2: Clean data
    df_clean = clean_county_typology_data(df)

    # Step 3: Generate summary
    generate_summary(df_clean)

    # Step 4: Save output
    save_output(df_clean, output_file, summary_file)

    print("\n" + "=" * 80)
    print("Data cleaning completed successfully!")
    print("=" * 80)
    print(f"\nOutput file: {output_file}")
    print(f"Summary file: {summary_file}")


if __name__ == "__main__":
    main()
