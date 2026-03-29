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

# EQI color map
EQI_COLORS = {
    1: "#2170b5",
    2: "#6baed6",
    3: "#c6dbef",
    4: "#fee0d2",
    5: "#fc9271",
    "No Data": "#cccccc",  # For missing data
}


def plot_eqi_map(period, config):
    """
    Plot EQI distribution map for a given time period
    FOR CONTIGUOUS US ONLY, without projection (plane view).
    """
    print(f"Creating EQI map for period {period} (Plane view)...")

    # Paths from config
    shapefile_path = config["data_sources"]["tiger"]["shapefile"]
    eqi_dir = config["data_sources"]["epa_eqi"]["processed"]
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

    # Load EQI data
    eqi_file = os.path.join(eqi_dir, f"EQI{period}.csv")
    df_eqi = pd.read_csv(eqi_file, usecols=["COUNTY_FIPS", "EQI"])
    df_eqi["COUNTY_FIPS"] = df_eqi["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Merge with shapefile
    counties_merged = counties_contiguous.merge(df_eqi, on="COUNTY_FIPS", how="left")

    # Assign colors
    counties_merged["color"] = (
        counties_merged["EQI"].map(EQI_COLORS).fillna(EQI_COLORS["No Data"])
    )

    # Prepare state boundaries once
    state_boundaries = counties_merged.dissolve(by="STATEFP")

    # Plot for each font
    fonts = ["Georgia", "Helvetica"]

    for font in fonts:
        plt.rcParams["font.family"] = font

        fig, ax = plt.subplots(1, 1, figsize=(16, 10))

        # Plot counties without projection
        counties_merged.plot(
            color=counties_merged["color"], linewidth=0.1, edgecolor="black", ax=ax
        )

        # State boundaries
        state_boundaries.boundary.plot(ax=ax, color="black", linewidth=1.2, alpha=0.9)

        ax.set_axis_off()

        # Legend
        legend_elements = [
            mpatches.Patch(color=EQI_COLORS[i], label=f"EQI {i}") for i in range(1, 6)
        ]
        ax.legend(
            handles=legend_elements,
            bbox_to_anchor=(0.02, 0.02),
            loc="lower left",
            frameon=True,
        )

        # Title
        period_labels = {"0005": "2000-2005", "0610": "2006-2010"}
        title_period = period_labels.get(period, period)
        plt.suptitle(f"EQI Distribution Map ({title_period})", y=0.82)

        # Save
        output_filename = os.path.join(output_dir, f"EQI_{period}_Plane_{font}.png")
        plt.savefig(output_filename, dpi=300, bbox_inches="tight")
        print(f"Map saved to: {output_filename}")
        plt.close()


if __name__ == "__main__":
    # Load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Plot for both periods
    for period in ["0005", "0610"]:
        plot_eqi_map(period, config)
