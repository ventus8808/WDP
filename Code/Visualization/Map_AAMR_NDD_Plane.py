import argparse
import os
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
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


# AAMR 7-color diverging color map based on percentiles
AAMR_COLORS = {
    1: "#3a5dae",  # Low (deep blue) - 0-5%
    2: "#677eba",  # Medium blue - 5-20%
    3: "#93a4cf",  # Light blue - 20-40%
    4: "#e8eaf0",  # Mid (near white) - 40-60%
    5: "#e0bbc0",  # Light red - 60-80%
    6: "#d88a91",  # Medium red - 80-95%
    7: "#be6e73",  # High (darker red) - 95-100%
    "No Data": "#cccccc",
}


def plot_aamr_map(time_period, icd, config):
    """Plot AAMR quintile map for a given time period and ICD code"""
    print(f"Creating AAMR map for time period {time_period}, ICD {icd}...")

    # Paths from config
    shapefile_path = config["data_sources"]["tiger"]["shapefile"]
    aamr_file = "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv"
    output_dir = "Result/Map"
    os.makedirs(output_dir, exist_ok=True)

    # Load shapefile
    counties = gpd.read_file(shapefile_path)
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

    # Load AAMR data
    df_aamr = pd.read_csv(aamr_file)
    df_aamr["COUNTY_FIPS"] = df_aamr["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Filter for the specified ICD
    df_filtered = df_aamr[(df_aamr["Cancer_Type"] == icd)].copy()

    # Filter for the specific time period
    df_period = df_filtered[df_filtered["Time_Period"] == time_period].copy()

    if df_period.empty:
        print(f"No data for time period {time_period}")
        return

    # Compute custom percentiles for AAMR (higher AAMR = higher level)
    percentiles = [0, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0]
    levels, bins = pd.qcut(
        df_period["AAMR"],
        q=percentiles,
        labels=range(1, 8),
        retbins=True,
        duplicates="drop",
    )
    df_period["AAMR_Level"] = levels

    # Get level bins for legend
    percentile_labels = [
        r"$\mathrm{0-5^{th}}$",
        r"$\mathrm{5^{th}-20^{th}}$",
        r"$\mathrm{20^{th}-40^{th}}$",
        r"$\mathrm{40^{th}-60^{th}}$",
        r"$\mathrm{60^{th}-80^{th}}$",
        r"$\mathrm{80^{th}-95^{th}}$",
        r"$\mathrm{95^{th}-100^{th}}$",
    ]
    bin_labels = []
    for i, perc_label in enumerate(percentile_labels):
        lower = f"{bins[i]:.1f}"
        upper = f"{bins[i + 1]:.1f}"
        if i == 0:
            label = f"{perc_label}: <{upper}"
        elif i == 6:
            label = f"{perc_label}: >{lower}"
        else:
            label = f"{perc_label}: {lower}-{upper}"
        bin_labels.append(label)

    # Merge with shapefile
    counties_merged = counties_contiguous.merge(
        df_period[["COUNTY_FIPS", "AAMR_Level"]], on="COUNTY_FIPS", how="left"
    )

    # Assign colors
    counties_merged["color"] = (
        counties_merged["AAMR_Level"]
        .astype(float)
        .map(AAMR_COLORS)
        .fillna(AAMR_COLORS["No Data"])
    )

    # Prepare state boundaries once
    state_boundaries = counties_merged.dissolve(by="STATEFP")

    # Plot for each font
    fonts = ["Georgia", "Helvetica"]

    for font in fonts:
        plt.rcParams["font.family"] = font

        fig, ax = plt.subplots(1, 1, figsize=(16, 10))

        counties_merged.plot(
            color=counties_merged["color"],
            linewidth=0,
            edgecolor="none",
            ax=ax,
            zorder=1,
        )

        state_boundaries.boundary.plot(
            ax=ax, color="#000000", linewidth=0.6, alpha=0.6, zorder=2
        )

        ax.set_axis_off()

        # Legend with ranges
        legend_elements = [
            mpatches.Patch(color=AAMR_COLORS[i + 1], label=bin_labels[i])
            for i in range(7)
        ]
        ax.legend(
            handles=legend_elements,
            bbox_to_anchor=(0.02, 0.02),
            loc="lower left",
            frameon=True,
        )

        # Save
        output_filename = os.path.join(
            output_dir, f"Map_AAMR_{time_period}_{icd}_{font}.png"
        )
        plt.savefig(output_filename, dpi=300, bbox_inches="tight")
        print(f"Map saved to: {output_filename}")
        plt.close()


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Plot Age-Adjusted Mortality Rate (AAMR) maps for neurodegenerative diseases.\n\n"
        "This script generates choropleth maps of AAMR data across US counties, using custom percentiles "
        "and a diverging color scheme. Maps are saved in the Result/Map/ directory.\n\n"
        "Examples:\n"
        "  python Map_AAMR_Climate.py  # Plot for G20_G30_G12.2_F01_F03"
    )
    args = parser.parse_args()

    # Load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Load AAMR data
    aamr_file = "Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv"
    df_aamr = pd.read_csv(aamr_file)

    # Filter for the specific cancer type
    icd = "G20_G30_G12.2_F01_F03"
    df_filtered = df_aamr[(df_aamr["Cancer_Type"] == icd)].copy()

    if df_filtered.empty:
        print(f"No data for Cancer_Type {icd}")
    else:
        unique_time_periods = df_filtered["Time_Period"].unique()

        # Plot for each unique time period
        for time_period in unique_time_periods:
            plot_aamr_map(time_period, icd, config)
