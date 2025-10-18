"""
Cluster_Map.py

Performs clustering analysis for k=3 to 10 clusters on EQI residuals,
creates a DataFrame with cluster assignments, and generates choropleth maps
for each cluster count.

Data source: EQI standardized data
Output: DataFrame saved to Data/Processed/df_EQI_AAMR/Cluster_Map.csv
Maps: Result/Figures/Cluster_Map/{k}_Cluster_Map.png
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import yaml
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf
import os

def load_config(project_root: Path) -> dict:
    """Load configuration from config.yaml"""
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_and_prepare_data(standard_file_path):
    """Load standardized EQI data and prepare for clustering"""
    print(f"Loading data from: {standard_file_path}")
    df = pd.read_csv(standard_file_path)

    # Select relevant columns
    eqi_columns = ['EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']
    required_cols = ['COUNTY_FIPS', 'RUCC'] + eqi_columns

    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. Expected: {required_cols}")

    df_clean = df[required_cols].dropna()
    # Filter valid RUCC values (1-4)
    df_clean = df_clean[df_clean['RUCC'].isin([1, 2, 3, 4])]
    print(f"Data loaded: {df_clean.shape[0]} counties, {len(eqi_columns)} EQI dimensions")

    return df_clean, eqi_columns

def compute_residuals(df, eqi_columns):
    """Compute residuals after regressing out RUCC effects"""
    print("Computing residuals...")

    residuals_df = df[['COUNTY_FIPS', 'RUCC']].copy()

    for col in eqi_columns:
        formula = f'{col} ~ C(RUCC)'
        model = smf.ols(formula, data=df).fit()
        residuals_df[f'{col}_residual'] = model.resid

    # Standardize residuals
    residual_cols = [f'{col}_residual' for col in eqi_columns]
    scaler = StandardScaler()
    residuals_df[residual_cols] = scaler.fit_transform(residuals_df[residual_cols])

    return residuals_df, residual_cols

def perform_clustering_for_range(X, k_range):
    """Perform K-Means for each k in range and return labels"""
    cluster_results = {}
    for k in k_range:
        print(f"Clustering with k={k}...")
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)
        cluster_results[k] = labels
    return cluster_results

def create_cluster_dataframe(residuals_df, cluster_results):
    """Create DataFrame with cluster assignments for each k"""
    df_clusters = residuals_df[['COUNTY_FIPS']].copy()
    for k, labels in cluster_results.items():
        df_clusters[f'cluster_{k}'] = labels

    # Ensure COUNTY_FIPS is 5-digit string
    df_clusters['COUNTY_FIPS'] = df_clusters['COUNTY_FIPS'].astype(str).str.zfill(5)

    return df_clusters

def save_cluster_dataframe(df_clusters, output_path):
    """Save the cluster DataFrame to CSV"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_clusters.to_csv(output_path, index=False)
    print(f"Cluster DataFrame saved to: {output_path}")

def create_map_for_k(df_clusters, k, shapefile_path, output_dir, color_scheme):
    """Create and save choropleth map for given k"""
    print(f"Creating map for k={k}...")

    # Load shapefile
    counties = gpd.read_file(shapefile_path)
    counties['COUNTY_FIPS'] = counties['STATEFP'] + counties['COUNTYFP']

    # Filter to contiguous US
    contiguous_states = ['01', '04', '05', '06', '08', '09', '10', '11', '12', '13', '16', '17', '18', '19',
                        '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33',
                        '34', '35', '36', '37', '38', '39', '40', '41', '42', '44', '45', '46', '47', '48',
                        '49', '50', '51', '53', '54', '55', '56']
    counties_contiguous = counties[counties['STATEFP'].isin(contiguous_states)].copy()

    # Merge with cluster data
    counties_merged = counties_contiguous.merge(df_clusters[['COUNTY_FIPS', f'cluster_{k}']],
                                               on='COUNTY_FIPS', how='left')

    # Assign colors
    if k == 3:
        # Fixed green gradient for k=3
        color_map = {0: '#ebf0b5', 1: '#a0d292', 2: '#44a05c', 'No Data': '#d3d3d3'}  # Light gray for no data
        cluster_names = {0: 'Disadvantaged', 1: 'Unbalanced', 2: 'Advantageous', 'No Data': 'No Data'}
    else:
        # Gradual color scheme for k>3 using Greens colormap
        cmap = plt.cm.Greens
        colors = [cmap(i / (k - 1)) for i in range(k)]
        color_map = {i: colors[i] for i in range(k)}
        color_map['No Data'] = '#d3d3d3'  # Light gray for no data
        cluster_names = {i: f'Cluster {i}' for i in range(k)}
        cluster_names['No Data'] = 'No Data'

    # Apply colors, handle NaN
    def get_color(cluster):
        if pd.isna(cluster):
            return color_map['No Data']
        return color_map.get(int(cluster), color_map['No Data'])

    counties_merged['cluster_color'] = counties_merged[f'cluster_{k}'].apply(get_color)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    counties_merged.plot(color=counties_merged['cluster_color'], linewidth=0.1, edgecolor='black', ax=ax)

    # State boundaries
    state_boundaries = counties_merged.dissolve(by='STATEFP')
    state_boundaries.boundary.plot(ax=ax, color='black', linewidth=1.2, alpha=0.9)

    ax.set_axis_off()
    ax.set_title(f'Environmental Quality Clusters (k={k})\nContiguous United States Counties', fontsize=18, pad=20)

    # Legend
    legend_elements = [mpatches.Patch(color=color_map[i], label=cluster_names[i]) for i in sorted(cluster_names.keys(), key=lambda x: (isinstance(x, str), x))]
    ax.legend(handles=legend_elements, bbox_to_anchor=(0.02, 0.02), loc='lower left', fontsize=12, frameon=True)

    plt.tight_layout()

    # Save
    output_filename = output_dir / f"{k}_Cluster_Map.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Map saved to: {output_filename}")
    plt.close()

def main():
    # Project setup
    project_root = Path(__file__).resolve().parents[2]
    config = load_config(project_root)

    # Paths
    standard_file = project_root / config["data_directories"]["processed"] / "EQI" / "EQI0005_Standard.csv"
    shapefile_path = project_root / config["data_sources"]["tiger"]["shapefile"]
    output_dir_maps = project_root / "Result" / "Figures" / "Cluster_Map"
    output_dir_maps.mkdir(parents=True, exist_ok=True)
    output_path_df = "/Users/ventus/Repository/WDP/Data/Processed/df_EQI_AAMR/Cluster_Map.csv"

    # Load and prepare data
    df, eqi_columns = load_and_prepare_data(standard_file)

    # Compute residuals
    residuals_df, residual_cols = compute_residuals(df, eqi_columns)

    # Clustering for k=3 to 10
    k_range = range(3, 11)
    cluster_results = perform_clustering_for_range(residuals_df[residual_cols].values, k_range)

    # Create and save DataFrame
    df_clusters = create_cluster_dataframe(residuals_df, cluster_results)
    save_cluster_dataframe(df_clusters, output_path_df)

    # Create maps
    color_scheme = {'name': 'Green Gradient'}  # Placeholder, handled in create_map_for_k
    for k in k_range:
        create_map_for_k(df_clusters, k, shapefile_path, output_dir_maps, color_scheme)

    print("All clustering and mapping completed!")

if __name__ == "__main__":
    main()