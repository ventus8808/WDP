import os
from math import pi

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.image import imread

plt.rcParams["font.family"] = "Georgia"
plt.rcParams["font.size"] = 16
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 14
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["legend.fontsize"] = 16

# Unified cluster colors
CLUSTER_COLORS = ["#44a05c", "#5a88c8", "#f49c4a"]


def get_color(cluster, color_map):
    """Get color for cluster"""
    if pd.isna(cluster):
        return color_map["No Data"]
    return color_map.get(int(cluster), color_map["No Data"])


def create_radar_chart(profiles_df, eqi_columns, output_dir, method_name, k):
    """Create radar chart for cluster profiles"""
    print("Creating radar chart...")

    # Prepare data
    categories = [col.replace("EQI_", "") for col in eqi_columns]
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection="polar"))

    angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
    angles += angles[:1]  # Close the loop

    colors = CLUSTER_COLORS
    for _, row in profiles_df.iterrows():
        values = [row[f"{col}_mean"] for col in eqi_columns]
        values += values[:1]  # Close the loop
        cluster_id = int(row["Cluster"])

        ax.plot(
            angles,
            values,
            "o-",
            linewidth=2,
            label=f"Cluster {cluster_id}",
            color=colors[cluster_id],
        )
        ax.fill(angles, values, alpha=0.25, color=colors[cluster_id])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_ylim(-2, 1)

    ax.set_title(f"(A) Radar Chart of {method_name} (k={k})", pad=20, fontsize=24)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.1, -0.1), fontsize=20)
    ax.grid(True)

    radar_path = os.path.join(output_dir, f"{method_name}_{k}_Radar_Chart.png")
    plt.savefig(radar_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Radar chart saved to: {radar_path}")
    return radar_path


def create_box_plot(df, eqi_columns, output_dir, method_name, k):
    """Create box plot for EQI by cluster"""
    print("Creating box plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=eqi_columns,
        var_name="EQI_Dimension",
        value_name="EQI_Value",
    )
    melted_df["EQI_Dimension"] = melted_df["EQI_Dimension"].str.replace("EQI_", "")
    melted_df["Cluster_Label"] = "Cluster " + melted_df["Cluster"].astype(str)

    plt.figure(figsize=(12, 8))
    sns.boxplot(
        data=melted_df,
        x="EQI_Dimension",
        y="EQI_Value",
        hue="Cluster_Label",
        hue_order=["Cluster 0", "Cluster 1", "Cluster 2"],
        palette=CLUSTER_COLORS,
        showfliers=False,
    )
    plt.title(f"(B) Box Plot of {method_name} (k={k})", pad=20, fontsize=24)
    plt.xlabel("EQI Domain")
    plt.ylabel("EQI Value (Standardized)")
    plt.legend(loc="lower left", fontsize=20)
    plt.grid(True, alpha=0.3)

    box_path = os.path.join(output_dir, f"{method_name}_{k}_Box_Plot.png")
    plt.savefig(box_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Box plot saved to: {box_path}")
    return box_path


def create_map_for_k(df_clusters, k, shapefile_path, output_dir, method_name):
    """Create and save choropleth map for given k and method"""
    print(f"Creating map for {method_name} k={k}...")

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

    # Merge with cluster data
    counties_merged = counties_contiguous.merge(
        df_clusters[["COUNTY_FIPS", f"cluster_{k}"]], on="COUNTY_FIPS", how="left"
    )

    # Assign colors
    # Best environment (cluster 0) = #44a05c, worst (cluster k-1) = #ebf0b5, interpolate in between
    best_color = "#44a05c"
    worst_color = "#ebf0b5"

    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(
        "env_gradient", [best_color, worst_color], N=k
    )
    colors = [cmap(i / (k - 1)) for i in range(k)]
    color_map = {i: colors[i] for i in range(k)}
    color_map["No Data"] = color_map[0]  # No data as cluster 0

    cluster_names = {i: f"Cluster {i}" for i in range(k)}
    cluster_names["No Data"] = "No Data"

    counties_merged["cluster_color"] = counties_merged[f"cluster_{k}"].apply(
        lambda x: get_color(x, color_map)
    )

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    counties_merged.plot(
        color=counties_merged["cluster_color"], linewidth=0.1, edgecolor="black", ax=ax
    )

    # State boundaries
    state_boundaries = counties_merged.dissolve(by="STATEFP")
    state_boundaries.boundary.plot(ax=ax, color="black", linewidth=1.2, alpha=0.9)

    ax.set_axis_off()

    # Legend
    legend_elements = [
        mpatches.Patch(color=color_map[i], label=cluster_names[i]) for i in range(k)
    ]
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        frameon=True,
    )

    plt.suptitle(
        f"(C) Contiguous U.S. Counties Map Clustering by {method_name} (k={k})", y=0.82
    )

    # Save
    output_filename = os.path.join(output_dir, f"{method_name}_{k}_Cluster_Map.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Map saved to: {output_filename}")
    plt.close()
    return output_filename


import os

import matplotlib.pyplot as plt
from matplotlib.image import imread


def create_combined_visualization(
    radar_path, box_path, map_path, output_dir, method_name, k
):
    """
    Create combined visualization:
    (A) Radar (top-left, 45% width)
    (B) Box (top-right, 55% width)
    (C) Map (bottom, 55% height)
    Layout proportions:
      - gap = 0.02 between A and B
      - top_h = 0.45 (top section)
      - bottom_h = 0.55 (bottom section)
      - left_w = 0.45 (A width)
      - right_w = 0.55 (B width)
    """
    print(f"Creating combined visualization for {method_name} k={k}...")

    # Load images
    radar_img = imread(radar_path)
    box_img = imread(box_path)
    map_img = imread(map_path)

    # Create figure
    fig = plt.figure(figsize=(20, 14))

    # Layout parameters
    gap = 0.02  # gap between A and B
    top_h = 0.45  # top section height
    bottom_h = 0.55  # bottom section height
    left_w = 0.42
    right_w = 0.58

    # Margins apply ONLY to top panels (A & B)
    margin_left_top = 0.15  # inward margin for A
    margin_right_top = 0.15  # inward margin for B

    # Compute top section usable width
    usable_w_top = 1.0 - margin_left_top - margin_right_top

    # Normalize A and B widths based on their proportions
    total_w = left_w + right_w
    left_frac = (left_w / total_w) * (usable_w_top - gap)
    right_frac = (right_w / total_w) * (usable_w_top - gap)

    # ---- Axes positions ----
    # [left, bottom, width, height]
    axA_pos = [margin_left_top, bottom_h + gap / 2, left_frac, top_h - gap / 2]  # (A)
    axB_pos = [
        margin_left_top + left_frac + gap,
        bottom_h + gap / 2,
        right_frac,
        top_h - gap / 2,
    ]  # (B)
    axC_pos = [0.0, 0.0, 1.0, bottom_h - gap / 2]

    # ---- (A) Radar ----
    axA = fig.add_axes(axA_pos)
    axA.imshow(radar_img)
    axA.axis("off")
    for spine in axA.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(1)

    # ---- (B) Box ----
    axB = fig.add_axes(axB_pos)
    axB.imshow(box_img)
    axB.axis("off")
    for spine in axB.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(1)

    # ---- (C) Map ----
    axC = fig.add_axes(axC_pos)
    axC.imshow(map_img)
    axC.axis("off")
    for spine in axC.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("black")
        spine.set_linewidth(1)

    # Save figure
    combined_path = os.path.join(output_dir, f"{method_name}_{k}_Combined.png")
    plt.savefig(combined_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Combined visualization saved to: {combined_path}")
