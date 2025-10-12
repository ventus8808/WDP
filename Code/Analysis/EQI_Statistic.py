#!/usr/bin/env python3
"""
EQI Statistics Script

This script computes basic statistics for the EQI 2000-2005 dataset and exports
the results to a text file.

Usage:
    python Code/Analysis/EQI_Statistic.py
"""

import pandas as pd
import yaml
from pathlib import Path


def main():
    # Find project root and load config
    script_dir = Path(__file__).resolve().parents[2]  # WDP root
    config_path = script_dir / "config.yaml"

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Get paths from config
    processed_dir = script_dir / config['data_directories']['processed']
    eqi_dir = processed_dir / "EQI"
    eqi_file = eqi_dir / "EQI0005.csv"

    results_dir = script_dir / config['result_directories']['tables']
    output_file = results_dir / "EQI_Statistic.txt"

    # Read the EQI data
    df = pd.read_csv(eqi_file)

    # Compute basic statistics - count of each distinct value for each column
    numeric_cols = df.select_dtypes(include=['number']).columns
    # Exclude COUNTY_FIPS as it has unique values for each county
    cols_to_analyze = numeric_cols.drop('COUNTY_FIPS') if 'COUNTY_FIPS' in numeric_cols else numeric_cols

    # Write to text file
    with open(output_file, 'w') as f:
        f.write("EQI 2000-2005 Basic Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write("Data shape: {}\n".format(df.shape))
        f.write("Number of counties: {}\n\n".format(len(df)))
        f.write("Count of distinct values for each column:\n\n")

        for col in cols_to_analyze:
            f.write(f"{col}:\n")
            counts = df[col].value_counts().sort_index()
            f.write(counts.to_string())
            f.write("\n\n")

    print(f"Statistics exported to {output_file}")


if __name__ == "__main__":
    main()