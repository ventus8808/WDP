import os
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# Set matplotlib parameters for consistent styling
# plt.rcParams["font.family"] = "Georgia"  # Set dynamically in loop
plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["mathtext.fontset"] = "stix"

# Colors for EQI Change Categories
# Improved: Blue (Good/Low Risk)
# Worsened: Red (Bad/High Risk)
# Stable: Grey (Neutral)
CHANGE_COLORS = {
    "improved": "#2170b5",  # Blue (EQI Quintile 1)
    "worsened": "#fc9271",  # Red/Orange (EQI Quintile 5)
    "stable": "#cccccc",  # Grey (No Data color)
}


def plot_eqi_delta_map(config):
    """Plot EQI Change Category map"""
    print("Creating EQI Delta map...")

    # Paths
    project_root = Path(__file__).resolve().parents[2]
    shapefile_path = config["data_sources"]["tiger"]["shapefile"]

    # Input data path
    data_file = project_root / "Data/Processed/df_EQI_AAMR/EQI_AAMR_Delta.csv"

    # Output directory
    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load shapefile
    print(f"Loading shapefile from {shapefile_path}...")
    try:
        counties = gpd.read_file(shapefile_path)
    except Exception as e:
        print(f"Error loading shapefile: {e}")
        return

    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]

    # Filter to contiguous US
    contiguous_states = [
        "01",
        "04",
        "05",
        "06",
        "08",
        "09",
        "10",
        "11",
        "12",
        "13",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "23",
        "24",
        "25",
        "26",
        "27",
        "28",
        "29",
        "30",
        "31",
        "32",
        "33",
        "34",
        "35",
        "36",
        "37",
        "38",
        "39",
        "40",
        "41",
        "42",
        "44",
        "45",
        "46",
        "47",
        "48",
        "49",
        "50",
        "51",
        "53",
        "54",
        "55",
        "56",
    ]
    counties_contiguous = counties[counties["STATEFP"].isin(contiguous_states)].copy()

    # Load Data
    if not data_file.exists():
        print(f"Error: Data file not found at {data_file}")
        return

    print(f"Loading data from {data_file}...")
    df = pd.read_csv(data_file)

    # Validate columns
    required_cols = ["COUNTY_FIPS", "EQI_Change_Category"]
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Input CSV must contain columns: {required_cols}")
        return

    # Ensure FIPS is string and padded
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Merge with shapefile
    counties_merged = counties_contiguous.merge(
        df[["COUNTY_FIPS", "EQI_Change_Category"]], on="COUNTY_FIPS", how="left"
    )

    # Assign colors
    def get_color(cat):
        if pd.isna(cat):
            return CHANGE_COLORS["stable"]
        # Normalize to lowercase to match keys
        cat_lower = str(cat).lower()
        return CHANGE_COLORS.get(cat_lower, CHANGE_COLORS["stable"])

    counties_merged["color"] = counties_merged["EQI_Change_Category"].apply(get_color)

    # Prepare state boundaries once
    state_boundaries = counties_merged.dissolve(by="STATEFP")

    fonts = ["Georgia", "Helvetica"]

    for font in fonts:
        plt.rcParams["font.family"] = font

        # Plot
        fig, ax = plt.subplots(1, 1, figsize=(16, 10))

        # Plot counties
        counties_merged.plot(
            color=counties_merged["color"],
            linewidth=0,
            edgecolor="none",
            ax=ax,
            zorder=1,
        )

        # Plot state boundaries
        state_boundaries.boundary.plot(
            ax=ax, color="#000000", linewidth=0.6, alpha=0.6, zorder=2
        )

        ax.set_axis_off()

        # Legend
        legend_elements = [
            mpatches.Patch(color=CHANGE_COLORS["improved"], label="Improved"),
            mpatches.Patch(color=CHANGE_COLORS["stable"], label="Stable"),
            mpatches.Patch(color=CHANGE_COLORS["worsened"], label="Worsened"),
        ]

        ax.legend(
            handles=legend_elements,
            bbox_to_anchor=(0.02, 0.02),
            loc="lower left",
            frameon=True,
            title="EQI Change Status",
        )

        # Save
        output_filename = output_dir / f"Map_EQI_Delta_{font}.png"
        plt.savefig(output_filename, dpi=300, bbox_inches="tight")
        print(f"Map saved to: {output_filename}")
        plt.close()


if __name__ == "__main__":
    # Load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"

    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f)
        plot_eqi_delta_map(config)
    else:
        print(f"Config file not found at {config_path}")
