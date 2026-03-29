import os
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# Set matplotlib parameters for consistent styling
plt.rcParams["font.family"] = "Georgia"
plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 11

# County Typology color map (thematic colors + gray for nonspecialized/missing)
TYPOLOGY_COLORS = {
    1: "#4DAF4A",  # Farming - natural green (vegetation/crops)
    2: "#984EA3",  # Mining - purple (minerals/underground resources)
    3: "#377EB8",  # Manufacturing - reliable blue (industry/technology)
    4: "#E41A1C",  # Government - authoritative red (administration)
    5: "#FF7F00",  # Services - orange (dynamic service sector)
    6: "#CCCCCC",  # Nonspecialized - gray
    "No Data": "#CCCCCC",  # Missing data - gray (same as nonspecialized)
}

TYPOLOGY_LABELS = {
    1: "Farming",
    2: "Mining",
    3: "Manufacturing",
    4: "Government",
    5: "Services",
    6: "Nonspecialized",
}


def plot_county_typology_map(config):
    """
    Plot USDA ERS County Typology 2004 distribution map
    for contiguous US without projection (plane view).
    """
    print("Creating County Typology map (Plane view)...")

    # Paths from config
    shapefile_path = config["data_sources"]["tiger"]["shapefile"]
    typology_file = "Data/Processed/Socioeconomic/County_Typology_2004.csv"
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

    # Load County Typology data
    df_typology = pd.read_csv(typology_file)
    df_typology["COUNTY_FIPS"] = df_typology["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Merge with shapefile
    counties_merged = counties_contiguous.merge(
        df_typology, on="COUNTY_FIPS", how="left"
    )

    # Assign colors
    counties_merged["color"] = (
        counties_merged["econdep"]
        .map(TYPOLOGY_COLORS)
        .fillna(TYPOLOGY_COLORS["No Data"])
    )

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))

    # Plot counties with typology colors (no projection)
    counties_merged.plot(
        color=counties_merged["color"], linewidth=0.1, edgecolor="white", ax=ax
    )

    # State boundaries
    state_boundaries = counties_merged.dissolve(by="STATEFP")
    state_boundaries.boundary.plot(ax=ax, color="gray", linewidth=0.8, alpha=0.7)

    ax.set_axis_off()

    # Legend - arranged in 2 columns for better layout
    legend_elements = [
        mpatches.Patch(color=TYPOLOGY_COLORS[i], label=TYPOLOGY_LABELS[i])
        for i in range(1, 7)
    ]

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        frameon=True,
        ncol=2,
        columnspacing=1.0,
        handletextpad=0.5,
    )

    # Title
    plt.suptitle(
        "USDA ERS County Economic Typology (2004)",
        y=0.82,
        fontsize=18,
        fontweight="normal",
    )

    # Save
    output_filename = os.path.join(output_dir, "County_Typology_2004_Plane.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Map saved to: {output_filename}")
    plt.close()

    # Print summary statistics
    print("\nCounty Type Distribution:")
    type_counts = counties_merged["econdep"].value_counts().sort_index()
    for type_code in range(1, 7):
        count = type_counts.get(type_code, 0)
        pct = (count / len(counties_merged)) * 100
        print(f"  {TYPOLOGY_LABELS[type_code]:25s}: {count:4d} ({pct:5.1f}%)")

    missing = counties_merged["econdep"].isna().sum()
    if missing > 0:
        print(
            f"  {'Missing data':25s}: {missing:4d} ({missing / len(counties_merged) * 100:5.1f}%)"
        )


if __name__ == "__main__":
    # Load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Generate map
    plot_county_typology_map(config)
