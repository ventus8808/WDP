#!/usr/bin/env python3
"""
Cluster Result Combiner for BRMS Analysis

This script aggregates Bayesian modeling results from different environmental clusters
into disease-specific CSV files. Each output file contains results for all clusters
for a specific cancer type.

Usage:
    python Code/Analysis/Cluster_result_Combine.py

Input:
    Result/brms_cluster/*.csv - Individual cluster results

Output:
    Result/brms_cluster_combined/ - Directory with combined results
    - C00_C97.csv (all clusters for all cancers)
    - C15_C26.csv
    - etc.
"""

import pandas as pd
import os
import glob
from pathlib import Path

def main():
    # Setup paths
    project_root = Path(__file__).resolve().parents[2]
    input_dir = project_root / "Result" / "brms_cluster"
    output_dir = project_root / "Result" / "brms_cluster_combined"

    # Create output directory
    output_dir.mkdir(exist_ok=True)

    # Find all CSV files
    csv_files = glob.glob(str(input_dir / "*_brms_Cluster*.csv"))
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return

    print(f"Found {len(csv_files)} CSV files to process")

    # Read and combine all data
    all_data = []
    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            all_data.append(df)
            print(f"Loaded: {Path(csv_file).name}")
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")

    if not all_data:
        print("No data loaded")
        return

    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)

    # Sort by ICD_Code, EQI_Period, AAMR_Period, Lag, Model
    combined_df = combined_df.sort_values(['ICD_Code', 'EQI_Period', 'AAMR_Period', 'Lag', 'Model'])

    # Group by ICD_Code and save separate files
    icd_codes = combined_df['ICD_Code'].unique()

    print(f"\nProcessing {len(icd_codes)} ICD codes...")

    for icd_code in sorted(icd_codes):
        # Filter data for this ICD code
        icd_data = combined_df[combined_df['ICD_Code'] == icd_code].copy()

        # Reset index for cleaner output
        icd_data = icd_data.reset_index(drop=True)

        # Save to file
        output_file = output_dir / f"{icd_code}.csv"
        icd_data.to_csv(output_file, index=False)

        # Print summary
        n_scenarios = len(icd_data)
        models = icd_data['Model'].unique()
        print(f"✓ {icd_code}: {n_scenarios} scenarios, {len(models)} model types saved to {output_file.name}")

    print(f"\n✅ All results combined and saved to {output_dir}")
    print(f"Total scenarios processed: {len(combined_df)}")

if __name__ == "__main__":
    main()