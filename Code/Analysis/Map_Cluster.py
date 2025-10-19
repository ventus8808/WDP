"""
Map RUCC and Cluster by County FIPS

Creates a choropleth map of US counties colored by RUCC (Rural-Urban Continuum Code)
and Cluster assignment for the contiguous United States.

Data source: EQI_AAMR_Interval_Clustered.csv
Output: Result/Figures/RUCC_Cluster_Map.png
"""

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import yaml

def load_config(project_root: Path) -> dict:
    """Load configuration from config.yaml"""
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def main():
    # Project setup
    project_root = Path(__file__).resolve().parents[2]
    config = load_config(project_root)

    # Paths
    data_path = project_root / config["eqi_aamr_outputs"]["base_dir"] / "EQI_AAMR_Interval_Clustered.csv"
    shapefile_path = project_root / config["data_sources"]["tiger"]["shapefile"]
    output_dir = project_root / config["result_directories"]["figures"]
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from: {data_path}")
    print(f"Using shapefile: {shapefile_path}")

    # Load data
    df = pd.read_csv(data_path, dtype={'COUNTY_FIPS': str})

    # Take first occurrence of each county (data has multiple scenarios per county)
    df_unique = df.drop_duplicates(subset='COUNTY_FIPS').copy()

    print(f"Loaded {len(df_unique)} unique counties")

    # Load US counties shapefile
    print("Loading US counties shapefile...")
    counties = gpd.read_file(shapefile_path)

    # Convert FIPS to string for merging
    counties['COUNTY_FIPS'] = counties['STATEFP'] + counties['COUNTYFP']

    # Filter to contiguous US (exclude Alaska and Hawaii)
    contiguous_states = ['01', '04', '05', '06', '08', '09', '10', '11', '12', '13', '16', '17', '18', '19',
                        '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33',
                        '34', '35', '36', '37', '38', '39', '40', '41', '42', '44', '45', '46', '47', '48',
                        '49', '50', '51', '53', '54', '55', '56']  # Excludes 02 (AK) and 15 (HI)

    counties_contiguous = counties[counties['STATEFP'].isin(contiguous_states)].copy()

    print(f"Filtered to {len(counties_contiguous)} contiguous US counties")

    # Merge data with shapefile
    counties_merged = counties_contiguous.merge(df_unique[['COUNTY_FIPS', 'RUCC', 'Cluster']],
                                               on='COUNTY_FIPS', how='left')

    print(f"Merged data: {counties_merged['RUCC'].notna().sum()} counties with data")

    # Create categorical color scheme based on RUCC and Cluster
    def create_category(rucc, cluster):
        if pd.isna(rucc) or pd.isna(cluster):
            return 'No Data'
        return f'RUCC{int(rucc)}_C{int(cluster)}'

    counties_merged['Category'] = counties_merged.apply(
        lambda row: create_category(row['RUCC'], row['Cluster']), axis=1
    )

    # Create categorical color scheme based on RUCC and Cluster
    def create_category(rucc, cluster):
        if pd.isna(rucc) or pd.isna(cluster):
            return 'No Data'
        return f'RUCC{int(rucc)}_C{int(cluster)}'

    counties_merged['Category'] = counties_merged.apply(
        lambda row: create_category(row['RUCC'], row['Cluster']), axis=1
    )

    # Use the specified green gradient color scheme
    color_scheme = {
        'name': 'Green Gradient',
        'colors': {
            0: '#ebf0b5',  # Light yellow-green for Disadvantaged
            1: '#a0d292',  # Light green for Unbalanced
            2: '#44a05c'   # Dark green for Advantageous
        }
    }

    cluster_names = {
        0: 'Disadvantaged',
        1: 'Unbalanced',
        2: 'Advantageous'
    }

    print(f"\nCreating map with {color_scheme['name']} color scheme...")

    # Create a single comprehensive map
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))

    # Create color mapping based on cluster (no data -> cluster 2)
    def get_cluster_color(cluster):
        if pd.isna(cluster):
            return color_scheme['colors'][2]  # No data counties -> cluster 2 (Advantageous)
        return color_scheme['colors'].get(int(cluster), color_scheme['colors'][2])

    counties_merged['cluster_color'] = counties_merged['Cluster'].apply(get_cluster_color)

    # Plot all counties with their cluster colors
    counties_merged.plot(
        color=counties_merged['cluster_color'],
        linewidth=0.1,
        edgecolor='black',
        ax=ax
    )

    # Add thicker state boundaries using county data
    # Group counties by state and create state boundaries
    state_boundaries = counties_merged.dissolve(by='STATEFP')

    # Plot state boundaries with thinner lines
    state_boundaries.boundary.plot(
        ax=ax,
        color='black',
        linewidth=1.2,
        alpha=0.9
    )

    # Remove axes
    ax.set_axis_off()

    # Title
    ax.set_title(f'Environmental Quality Clusters - {color_scheme["name"]}\nContiguous United States Counties',
                fontsize=18, pad=20)

    # Create simplified legend
    legend_elements = [
        mpatches.Patch(color=color_scheme['colors'][2], label='Advantageous'),
        mpatches.Patch(color=color_scheme['colors'][1], label='Unbalanced'),
        mpatches.Patch(color=color_scheme['colors'][0], label='Disadvantaged')
    ]

    ax.legend(handles=legend_elements,
             bbox_to_anchor=(0.02, 0.02),
             loc='lower left',
             fontsize=12,
             frameon=True,
             borderaxespad=0)

    plt.tight_layout()

    # Save the map
    output_filename = output_dir / f"Clusters_Map_Green_Gradient.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Map saved to: {output_filename}")
    plt.close()

    # Print cluster statistics
    print("\nCluster Statistics:")
    for cluster_id in [0, 1, 2]:
        cluster_count = (counties_merged['Cluster'] == cluster_id).sum()
        total_count = len(counties_merged)
        percentage = (cluster_count / total_count) * 100
        print(f"Cluster {cluster_id} ({cluster_names[cluster_id]}): {cluster_count} counties ({percentage:.1f}% of total)")

    # Show overall statistics
    print("\n" + "="*50)
    print("OVERALL DATA SUMMARY:")
    print("="*50)
    print(f"Total counties: {len(counties_merged)}")
    print(f"Counties with data: {counties_merged['RUCC'].notna().sum()}")
    print(f"Counties without data: {counties_merged['RUCC'].isna().sum()}")
    print("\nNote: This dataset contains only RUCC codes 1-4 (urban/metro areas)")
    print("RUCC 5-9 (rural areas) are not present in the current data.")

    if counties_merged['RUCC'].notna().sum() > 0:
        print(f"\nRUCC distribution (urban/metro areas only):")
        print(counties_merged['RUCC'].value_counts().sort_index())

        print(f"\nCluster distribution:")
        print(counties_merged['Cluster'].value_counts().sort_index())

if __name__ == "__main__":
    main()