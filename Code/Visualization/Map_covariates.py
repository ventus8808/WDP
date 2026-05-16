import os
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# Set matplotlib parameters for consistent styling
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["font.size"] = 14
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 14
plt.rcParams["ytick.labelsize"] = 14
plt.rcParams["legend.fontsize"] = 14
plt.rcParams["mathtext.fontset"] = "stix"

# Color palettes for each covariate
RUCC_COLORS = {
    1: "#08519C",  # Metropolitan urbanized - deep blue
    2: "#3182BD",  # Non-metropolitan urbanized - medium blue
    3: "#6BAED6",  # Less urbanized - light blue
    4: "#BDD7E7",  # Thinly populated - pale blue
    "No Data": "#cccccc",
}

RUCC_LABELS = {
    1: "Metropolitan Urbanized (RUCC1)",
    2: "Non-Metropolitan Urbanized (RUCC2)",
    3: "Less Urbanized (RUCC3)",
    4: "Thinly Populated (RUCC4)",
}

CENSUS_REGION_COLORS = {
    1: "#5B8C5A",  # Northeast - forest green
    2: "#3B7EA1",  # Midwest - steel blue
    3: "#D4A373",  # South - warm tan
    4: "#9B5DE5",  # West - purple
    "No Data": "#cccccc",
}

CENSUS_REGION_LABELS = {
    1: "Northeast",
    2: "Midwest",
    3: "South",
    4: "West",
}

KOPPEN_COLORS = {
    "B": "#E07A5F",  # Dry - terracotta
    "C": "#81B29A",  # Temperate - sage green
    "D": "#3D5A80",  # Continental - slate blue
    "No Data": "#cccccc",
}

KOPPEN_LABELS = {
    "B": "Dry (B)",
    "C": "Temperate (C)",
    "D": "Continental (D)",
}

# County Typology color map (muted, harmonious palette)
TYPOLOGY_COLORS = {
    1: "#7DB892",  # Farming - soft sage green
    2: "#9B8EC2",  # Mining - muted lavender
    3: "#6A9BC3",  # Manufacturing - soft steel blue
    4: "#D4918E",  # Government - dusty rose
    5: "#E8B87D",  # Services - warm sand
    6: "#B8B8B8",  # Nonspecialized - neutral gray
    "No Data": "#cccccc",
}

TYPOLOGY_LABELS = {
    1: "Farming",
    2: "Mining",
    3: "Manufacturing",
    4: "Government",
    5: "Services",
    6: "Nonspecialized",
}

# Contiguous US state FIPS codes
CONTIGUOUS_STATES = [
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


def plot_covariate_map(
    counties_merged, state_boundaries, covariate, colors, labels, output_dir
):
    """Plot a map for a given covariate"""
    print(f"Creating map for {covariate}...")

    fig, ax = plt.subplots(1, 1, figsize=(16, 10))

    counties_merged.plot(
        color=counties_merged["color"],
        linewidth=0,
        edgecolor="none",
        ax=ax,
        zorder=1,
    )

    state_boundaries.boundary.plot(
        ax=ax, color="#000000", linewidth=0.3, alpha=0.6, zorder=2
    )

    ax.set_axis_off()

    # Create legend
    if covariate == "Climate_Zone":
        legend_keys = ["B", "C", "D"]
    elif covariate == "Economic_type":
        legend_keys = [1, 2, 3, 4, 5, 6]
    else:
        legend_keys = [1, 2, 3, 4]

    legend_elements = [
        mpatches.Patch(color=colors[k], label=labels[k]) for k in legend_keys
    ]

    # Use 2 columns for typology legend (6 items)
    ncol = 2 if covariate == "Economic_type" else 1
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        frameon=True,
        ncol=ncol,
    )

    # Save
    output_filename = os.path.join(output_dir, f"Map_{covariate}.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Map saved to: {output_filename}")
    plt.close()


def main():
    # Load config
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Paths
    shapefile_path = config["data_sources"]["tiger"]["shapefile"]
    data_file = project_root / "Data/Processed/df.csv"
    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load shapefile
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]

    # Filter to contiguous US
    counties_contiguous = counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()

    # Load data
    df = pd.read_csv(data_file, low_memory=False)
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Get unique counties with covariates (use first occurrence)
    df_unique = df.drop_duplicates(subset=["COUNTY_FIPS"])[
        ["COUNTY_FIPS", "RUCC", "Census_Region", "Climate_Zone"]
    ].copy()

    # Merge with shapefile
    counties_merged = counties_contiguous.merge(df_unique, on="COUNTY_FIPS", how="left")

    # Prepare state boundaries
    state_boundaries = counties_merged.dissolve(by="STATEFP")

    # Plot RUCC map
    counties_merged["color"] = (
        counties_merged["RUCC"].map(RUCC_COLORS).fillna(RUCC_COLORS["No Data"])
    )
    plot_covariate_map(
        counties_merged, state_boundaries, "RUCC", RUCC_COLORS, RUCC_LABELS, output_dir
    )

    # Plot Census Region map
    counties_merged["color"] = (
        counties_merged["Census_Region"]
        .map(CENSUS_REGION_COLORS)
        .fillna(CENSUS_REGION_COLORS["No Data"])
    )
    plot_covariate_map(
        counties_merged,
        state_boundaries,
        "Census_Region",
        CENSUS_REGION_COLORS,
        CENSUS_REGION_LABELS,
        output_dir,
    )

    # Plot Climate Zone map
    counties_merged["color"] = (
        counties_merged["Climate_Zone"]
        .map(KOPPEN_COLORS)
        .fillna(KOPPEN_COLORS["No Data"])
    )
    plot_covariate_map(
        counties_merged,
        state_boundaries,
        "Climate_Zone",
        KOPPEN_COLORS,
        KOPPEN_LABELS,
        output_dir,
    )

    # Load and plot Economic Type data
    typology_file = project_root / "Data/Processed/Covariate/Economic_type.csv"
    df_typology = pd.read_csv(typology_file, dtype={"COUNTY_FIPS": str})
    df_typology["COUNTY_FIPS"] = df_typology["COUNTY_FIPS"].astype(str).str.zfill(5)
    df_typology["Economic_type"] = df_typology["Economic_type_0005"]

    # Merge typology with shapefile
    counties_typology = counties_contiguous.merge(
        df_typology[["COUNTY_FIPS", "Economic_type"]], on="COUNTY_FIPS", how="left"
    )

    # Prepare state boundaries for typology
    state_boundaries_typology = counties_typology.dissolve(by="STATEFP")

    # Plot Economic Type map
    counties_typology["color"] = (
        counties_typology["Economic_type"]
        .map(TYPOLOGY_COLORS)
        .fillna(TYPOLOGY_COLORS["No Data"])
    )
    plot_covariate_map(
        counties_typology,
        state_boundaries_typology,
        "Economic_type",
        TYPOLOGY_COLORS,
        TYPOLOGY_LABELS,
        output_dir,
    )


if __name__ == "__main__":
    main()
