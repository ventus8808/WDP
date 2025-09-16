#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP Principal Component Analysis (PCA) Core Module

This script performs the main PCA workflow:
1. Loads the master dataset using the PCA_Data_Loading module.
2. Performs VIF-based variable selection for SVI and Climate domains.
3. Executes PCA on the selected variables.
4. Generates and saves three key outputs:
    a) The final master covariate file with PCA scores.
    b) A comprehensive diagnostics table for all analyses.
    c) The raw data required for generating plots in a separate script.

Author: WDP Analysis Team (with Cascade AI)
Date: 2024-09-04
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor

# Import the data loading function from our new module
from PCA_Data_Loading import load_all_processed_data

warnings.filterwarnings('ignore')

# --- Configuration and Setup ---
# 使用相对于项目根目录的路径
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# --- Core Analysis Functions ---

def select_variables_with_vif(df: pd.DataFrame, variables: list, threshold: float = 10.0) -> tuple[list, pd.DataFrame]:
    """Iteratively select variables by removing those with a VIF above the threshold."""
    print(f"Starting VIF selection with threshold {threshold}...")
    vif_df = df[variables].dropna()
    scaler = StandardScaler()
    vif_df_scaled = pd.DataFrame(scaler.fit_transform(vif_df), columns=vif_df.columns)
    
    selected_vars = list(vif_df_scaled.columns)
    vif_log = []

    while True:
        vif_data = pd.DataFrame()
        vif_data["Variable"] = selected_vars
        vif_data["VIF"] = [variance_inflation_factor(vif_df_scaled[selected_vars].values, i) for i in range(len(selected_vars))]
        
        max_vif = vif_data["VIF"].max()
        if max_vif > threshold:
            max_vif_var = vif_data.sort_values("VIF", ascending=False).iloc[0]["Variable"]
            vif_log.append({"Variable Removed": max_vif_var, "VIF": max_vif, "Status": "Removed"})
            selected_vars.remove(max_vif_var)
            print(f"  - Removing '{max_vif_var}' (VIF: {max_vif:.2f})")
        else:
            print("  All remaining variables are below the VIF threshold.")
            for _, row in vif_data.iterrows():
                vif_log.append({"Variable Removed": row["Variable"], "VIF": row["VIF"], "Status": "Kept"})
            break
            
    print(f"  Final selected variables: {selected_vars}")
    return selected_vars, pd.DataFrame(vif_log)

def run_pca(df: pd.DataFrame, selected_vars: list, analysis_name: str) -> dict:
    """Run PCA and determine components to keep using the Kaiser criterion."""
    print(f"Running PCA for {analysis_name} on {len(selected_vars)} variables...")
    # Create a copy to avoid modifying the original dataframe
    pca_df = df[['COUNTY_FIPS', 'Year'] + selected_vars].copy()
    # Drop rows where ANY of the selected variables for this specific PCA are missing
    pca_df.dropna(subset=selected_vars, inplace=True)  # type: ignore
    
    if pca_df.empty:
        print(f"  Warning: No complete data available for {analysis_name} PCA. Skipping.")
        return None  # type: ignore

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(pca_df[selected_vars])
    
    # First pass to get all eigenvalues for Kaiser criterion
    pca = PCA()
    pca.fit(scaled_data)
    
    eigenvalues = pca.explained_variance_
    n_components_kaiser = np.sum(eigenvalues > 1.0)
    if n_components_kaiser == 0:
        n_components_kaiser = 1 # Always keep at least one component
        print(f"  Kaiser criterion suggests 0 components. Defaulting to 1.")
    else:
        print(f"  Kaiser criterion suggests keeping {n_components_kaiser} components.")
    
    # Refit PCA with the selected number of components
    final_pca = PCA(n_components=n_components_kaiser)
    scores = final_pca.fit_transform(scaled_data)

    # --- SVI Inversion Logic ---
    # If this is the SVI analysis, ensure the first component represents vulnerability.
    # We check the loading of a variable where a higher value means more vulnerable (e.g., Poverty).
    if analysis_name == 'SVI':
        try:
            # Find the index of a key vulnerability indicator
            vulnerability_indicator = 'Poverty_Percent_All_Ages'
            if vulnerability_indicator in selected_vars:
                indicator_index = selected_vars.index(vulnerability_indicator)
                # If the loading for this indicator is negative, flip the component
                if final_pca.components_[0, indicator_index] < 0:
                    print("  Inverting SVI PC1 to align with vulnerability (higher score = more vulnerable).")
                    final_pca.components_[0, :] *= -1
                    scores[:, 0] *= -1
        except (ValueError, IndexError) as e:
            print(f"  Could not perform SVI inversion check: {e}")

    # Create a dataframe with the scores and original identifiers
    score_cols = {f'{analysis_name}_PC{i+1}': scores[:, i] for i in range(n_components_kaiser)}
    scores_df = pca_df[['COUNTY_FIPS', 'Year']].copy()
    scores_df = scores_df.assign(**score_cols)  # type: ignore
    
    return {
        "pca_model": final_pca,
        "scores_df": scores_df,
        "selected_vars": selected_vars,
        "eigenvalues": eigenvalues
    }

# --- Reporting Functions ---

def save_unified_diagnostics_table(reports: dict, filename: Path):
    """Save a unified PCA diagnostics table with standardized format including VIF values."""
    print(f"Creating unified PCA diagnostics table at {filename}...")
    
    unified_data = []
    
    for analysis_name, results in reports.items():
        if results['pca_results'] is None:
            continue
            
        pca_results = results['pca_results']
        pca = pca_results['pca_model']
        selected_vars = pca_results['selected_vars']
        vif_log = results['vif_log']
        
        # Create a dictionary to map variables to their final VIF values
        # For variables that were kept, get their VIF from the log
        # For variables that were removed, get their VIF at removal
        vif_dict = {}
        for _, row in vif_log.iterrows():
            vif_dict[row['Variable Removed']] = row['VIF']
        
        # Get component diagnostics
        eigenvalues = pca.explained_variance_
        explained_variance_ratio = pca.explained_variance_ratio_ * 100
        cumulative_variance = np.cumsum(explained_variance_ratio)
        
        # For each component, add rows for each variable loading
        for comp_idx in range(pca.n_components_):
            component_name = f'PC{comp_idx + 1}'
            eigenvalue = eigenvalues[comp_idx]
            explained_var = explained_variance_ratio[comp_idx]
            cumulative_var = cumulative_variance[comp_idx]
            
            # Get loadings for this component and sort by absolute value
            component_loadings = pca.components_[comp_idx]
            
            # Create list of (variable, loading) tuples and sort by absolute loading value
            var_loading_pairs = [(selected_vars[i], component_loadings[i]) for i in range(len(selected_vars))]
            var_loading_pairs.sort(key=lambda x: abs(x[1]), reverse=True)
            
            for variable_name, loading_value in var_loading_pairs:
                # Get VIF value for this variable (should be from the final kept variables)
                vif_value = vif_dict.get(variable_name, np.nan)
                
                unified_data.append({
                    'PCA_Type': analysis_name,
                    'Component': component_name,
                    'Eigenvalue': round(eigenvalue, 3),
                    'Explained Variance (%)': round(explained_var, 3),
                    'Cumulative Variance (%)': round(cumulative_var, 3),
                    'Variable': variable_name,
                    'Loading': round(loading_value, 3),
                    'VIF': round(vif_value, 3) if not pd.isna(vif_value) else np.nan
                })
    
    # Create DataFrame and save
    unified_df = pd.DataFrame(unified_data)
    unified_df.to_csv(filename, index=False)
    print(f"  Unified diagnostics table saved with {len(unified_df)} rows.")

def save_plot_data(reports: dict, data_dir: Path):
    """Saves the raw data needed for plotting."""
    print(f"Saving all plotting data to {data_dir}...")
    for name, results in reports.items():
        pca_results = results['pca_results']
        pca = pca_results['pca_model']
        
        # Eigenvalues for Scree Plot
        eigen_df = pd.DataFrame({'Component': range(1, len(pca_results['eigenvalues']) + 1), 'Eigenvalue': pca_results['eigenvalues']})
        eigen_df.to_csv(data_dir / f'{name}_scree_plot_data.csv', index=False)
        
        # Loadings and Scores for Biplot/Loading Plot
        loadings_df = pd.DataFrame(pca.components_.T, 
                                   columns=[f'PC{i+1}' for i in range(pca.n_components_)],  # type: ignore
                                   index=pca_results['selected_vars'])
        loadings_df.to_csv(data_dir / f'{name}_loadings_data.csv')
        
        # Save a sample of scores to keep file size manageable
        scores_df = pca_results['scores_df']
        sample_scores = scores_df.sample(n=min(5000, len(scores_df)), random_state=1)
        sample_scores.to_csv(data_dir / f'{name}_scores_sample_data.csv', index=False)
    print("  Plotting data saved.")


# --- Main Execution ---

def main():
    """Main execution function to run the full PCA pipeline."""
    print("\n" + "="*60)
    print("STARTING WDP PCA ANALYSIS PIPELINE")
    print("="*60)

    # --- 1. Setup Paths ---
    print("Setting up output directories...")
    processed_pca_dir = PROJECT_ROOT / "Data/Processed/PCA"
    tables_dir = PROJECT_ROOT / "Result/Tables"
    plot_data_dir = PROJECT_ROOT / "Result/Figure_Original_Data"
    
    processed_pca_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    plot_data_dir.mkdir(parents=True, exist_ok=True)

    # --- 2. Load Data ---
    master_df = load_all_processed_data()

    # --- 3. Define Analysis Tasks ---
    analysis_tasks = {
        'SVI': {
            'vars': ['Poverty_Percent_All_Ages', 'Median_Household_Income', 'Unemployment_Rate', 
                     'Per_Capita_Income', 'Less_Than_High_School_Percent', 'College_Plus_Percent']
        },
        'ENV': {
            'vars': ['tas_mean_annual', 'prcp_sum_annual', 'wind_mean_annual', 'rh_mean_annual', 
                     'swrad_mean_annual', 'potevap_sum_annual', 'tas_mean_DJF', 'tas_mean_MAM', 
                     'tas_mean_JJA', 'tas_mean_SON', 'prcp_sum_DJF', 'prcp_sum_MAM', 
                     'prcp_sum_JJA', 'prcp_sum_SON']
        }
    }

    # --- 4. Run Workflow for Each Task ---
    all_diagnostics = {}
    final_pca_scores = []
    all_raw_vars = []

    for name, task in analysis_tasks.items():
        print(f"\n--- Processing: {name} ---")
        selected_vars, vif_log = select_variables_with_vif(master_df, task['vars'])
        pca_results = run_pca(master_df, selected_vars, name)
        
        all_diagnostics[name] = {'vif_log': vif_log, 'pca_results': pca_results}
        final_pca_scores.append(pca_results['scores_df'])
        all_raw_vars.extend(task['vars'])

    # --- 5. Create and Save Master Covariate File ---
    print("\n--- Creating Final Master Covariate File ---")
    # Start with a clean base of identifiers and key non-PCA variables
    base_cols = ['COUNTY_FIPS', 'Year', 'Total_Population', 'Urbanization_Code']
    final_df = master_df[base_cols].copy()

    # Merge all PCA scores into the base dataframe
    for scores_df in final_pca_scores:
        if scores_df is not None and not scores_df.empty:
            final_df = final_df.merge(scores_df, on=['COUNTY_FIPS', 'Year'], how='left')

    # --- Formatting and Cleanup ---
    # Round PCA columns to 2 decimal places
    pca_cols = [col for col in final_df.columns if '_PC' in col]
    for col in pca_cols:
        final_df[col] = final_df[col].round(2)

    # Convert population and urbanization to integer, handling NaNs
    final_df['Total_Population'] = final_df['Total_Population'].astype('Int64')
    final_df['Urbanization_Code'] = final_df['Urbanization_Code'].astype('Int64')

    output_master_file = processed_pca_dir / "PCA_Master_Covariables.csv"
    final_df.to_csv(output_master_file, index=False)
    print(f"Master covariate file saved to: {output_master_file}")

    # --- 6. Save Diagnostics and Plot Data ---
    save_unified_diagnostics_table(all_diagnostics, tables_dir / "PCA_Diagnose.csv")
    save_plot_data(all_diagnostics, plot_data_dir)

    print("\n" + "="*60)
    print("PCA ANALYSIS PIPELINE COMPLETED SUCCESSFULLY")
    print("="*60)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)


