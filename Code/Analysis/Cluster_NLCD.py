# Code/Analysis/Cluster_NLCD.py

import os
from math import pi
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# Define cluster type color scheme
CLUSTER_COLORS = {
    "Urban": "#E68785",  # Coral red
    "Water-Sensitive": "#699DCB",  # Sky blue
    "Natural": "#6AAA81",  # Sage green
    "Agricultural": "#EFC085",  # Golden wheat
    "Mixed": "#95a5a6",  # Gray
}


def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_and_prepare_nlcd_data(nlcd_file_path):
    """Load NLCD/JRC data and prepare for clustering"""
    print(f"Loading data from: {nlcd_file_path}")
    df = pd.read_csv(nlcd_file_path)

    print(f"Raw data: {df.shape[0]} rows, {df['COUNTY_FIPS'].nunique()} counties")
    print(f"Year range: {df['Year'].min()} - {df['Year'].max()}")

    # Filter to 2000-2005 period (matching EQI period)
    df_filtered = df[(df["Year"] >= 2000) & (df["Year"] <= 2005)].copy()
    print(f"Filtered to 2000-2005: {df_filtered.shape[0]} rows")

    # Calculate 6-year average for each county
    group_cols = [
        "jrc_permanent_water_km2",
        "nlcd_cropland_km2",
        "nlcd_pasture_km2",
        "nlcd_urban_km2",
        "nlcd_wetland_km2",
        "nlcd_forest_km2",
        "nlcd_shrub_km2",
        "nlcd_grassland_km2",
        "nlcd_barren_km2",
        "total_area_km2",
    ]

    df_avg = df_filtered.groupby("COUNTY_FIPS")[group_cols].mean().reset_index()
    print(f"Averaged data: {df_avg.shape[0]} counties")

    # Feature engineering: Calculate 4 core ratio features
    print("\nCalculating feature ratios...")

    # 1. Agriculture Intensity = (cropland + pasture) / total_area
    df_avg["Agriculture_Intensity"] = (
        df_avg["nlcd_cropland_km2"] + df_avg["nlcd_pasture_km2"]
    ) / df_avg["total_area_km2"]

    # 2. Urban/Developed Intensity = urban / total_area
    df_avg["Urban_Intensity"] = df_avg["nlcd_urban_km2"] / df_avg["total_area_km2"]

    # 3. Hydro Sensitivity = (permanent_water + wetland) / total_area
    df_avg["Hydro_Sensitivity"] = (
        df_avg["jrc_permanent_water_km2"] + df_avg["nlcd_wetland_km2"]
    ) / df_avg["total_area_km2"]

    # 4. Nature Background = (forest + shrub + grassland + barren) / total_area
    df_avg["Nature_Background"] = (
        df_avg["nlcd_forest_km2"]
        + df_avg["nlcd_shrub_km2"]
        + df_avg["nlcd_grassland_km2"]
        + df_avg["nlcd_barren_km2"]
    ) / df_avg["total_area_km2"]

    # Feature columns
    feature_columns = [
        "Agriculture_Intensity",
        "Urban_Intensity",
        "Hydro_Sensitivity",
        "Nature_Background",
    ]

    # Check for missing values
    missing = df_avg[feature_columns].isnull().sum()
    if missing.any():
        print(f"Warning: Missing values detected:\n{missing}")
        df_avg = df_avg.dropna(subset=feature_columns)
        print(f"After dropping NaN: {df_avg.shape[0]} counties")

    # Basic statistics
    print("\nFeature statistics:")
    print(df_avg[feature_columns].describe())

    # Check that ratios sum to approximately 1 (should be close due to complete coverage)
    df_avg["ratio_sum"] = df_avg[feature_columns].sum(axis=1)
    print(f"\nRatio sum statistics (should be ~1.0):")
    print(df_avg["ratio_sum"].describe())

    return df_avg, feature_columns


def standardize_data(df, feature_columns):
    """Standardize feature data for clustering"""
    print("\nStandardizing feature data...")
    scaler = StandardScaler()
    df_standardized = df.copy()
    df_standardized[feature_columns] = scaler.fit_transform(df[feature_columns])

    print("Standardized feature statistics:")
    print(df_standardized[feature_columns].describe())

    return df_standardized, scaler


def evaluate_clusters(X, k_range, output_dir):
    """Evaluate different numbers of clusters using elbow method and silhouette score"""
    print("\nEvaluating optimal number of clusters...")

    wcss = []
    silhouette_scores = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(X)
        wcss.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(X, kmeans.labels_))
        print(
            f"  k={k}: WCSS={kmeans.inertia_:.2f}, Silhouette={silhouette_scores[-1]:.3f}"
        )

    # Elbow method plot
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(k_range, wcss, "bo-", linewidth=2, markersize=8)
    plt.xlabel("Number of Clusters (K)", fontsize=12)
    plt.ylabel("Within-Cluster Sum of Squares (WCSS)", fontsize=12)
    plt.title("Elbow Method for Optimal K", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)

    # Silhouette score plot
    plt.subplot(1, 2, 2)
    plt.plot(k_range, silhouette_scores, "ro-", linewidth=2, markersize=8)
    plt.xlabel("Number of Clusters (K)", fontsize=12)
    plt.ylabel("Silhouette Score", fontsize=12)
    plt.title("Silhouette Analysis for Optimal K", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    eval_plot_path = os.path.join(output_dir, "NLCD_Cluster_Evaluation.png")
    plt.savefig(eval_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Cluster evaluation plots saved to: {eval_plot_path}")

    return wcss, silhouette_scores


def perform_clustering(X, n_clusters):
    """Perform K-Means clustering"""
    print(f"\nPerforming K-Means clustering with {n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    return labels, kmeans


def identify_cluster_type(profile_row, feature_columns):
    """Identify cluster type based on dominant feature"""
    means = {
        col.replace("_mean", ""): profile_row[col]
        for col in [f"{fc}_mean" for fc in feature_columns]
    }
    dominant_feature = max(means, key=means.get)

    # Map feature to cluster type
    type_map = {
        "Agriculture_Intensity": "Agricultural",
        "Urban_Intensity": "Urban",
        "Nature_Background": "Natural",
        "Hydro_Sensitivity": "Water-Sensitive",
    }

    return type_map.get(dominant_feature, "Mixed")


def profile_clusters(df, feature_columns, labels):
    """Create cluster profiles"""
    print("\nCreating cluster profiles...")

    df["Cluster"] = labels
    profiles = []

    for cluster in np.unique(labels):
        cluster_data = df[df["Cluster"] == cluster]
        profile = {"Cluster": cluster, "Count": len(cluster_data)}

        for col in feature_columns:
            profile[f"{col}_mean"] = cluster_data[col].mean()
            profile[f"{col}_std"] = cluster_data[col].std()

        profiles.append(profile)

    profiles_df = pd.DataFrame(profiles)

    # Identify cluster type
    profiles_df["Cluster_Type"] = profiles_df.apply(
        lambda row: identify_cluster_type(row, feature_columns), axis=1
    )

    print("\nCluster profiles:")
    print(profiles_df[["Cluster", "Count", "Cluster_Type"]])

    return df, profiles_df


def sort_clusters_by_intervention(profiles_df, feature_columns):
    """Sort clusters by human intervention level (Natural -> Agricultural -> Urban)"""
    print("\nSorting clusters by human intervention level...")

    # Calculate intervention score: Urban + Agriculture - Nature
    # Higher score = more human intervention
    profiles_df["intervention_score"] = (
        profiles_df["Urban_Intensity_mean"]
        + profiles_df["Agriculture_Intensity_mean"]
        - profiles_df["Nature_Background_mean"]
    )

    # Sort by intervention score (ascending: least to most intervention)
    profiles_df = profiles_df.sort_values("intervention_score").reset_index(drop=True)

    # Create mapping from old cluster labels to new sorted labels
    old_to_new = {old: new for new, old in enumerate(profiles_df["Cluster"])}

    # Reassign cluster numbers
    profiles_df["old_cluster"] = profiles_df["Cluster"]
    profiles_df["Cluster"] = range(len(profiles_df))

    print("Cluster ordering (by intervention level):")
    for idx, row in profiles_df.iterrows():
        print(
            f"  Cluster {row['Cluster']}: {row['Cluster_Type']} "
            f"(Count: {row['Count']}, Intervention: {row['intervention_score']:.3f})"
        )

    return profiles_df, old_to_new


def create_radar_chart(profiles_df, feature_columns, output_dir, k):
    """Create radar chart for cluster profiles"""
    print("\nCreating radar chart...")

    # Prepare data
    categories = [col.replace("_", "\n") for col in feature_columns]
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection="polar"))

    angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
    angles += angles[:1]  # Close the loop

    # Plot each cluster with its type-specific color
    for idx, row in profiles_df.iterrows():
        values = [row[f"{col}_mean"] for col in feature_columns]
        values += values[:1]  # Close the loop

        cluster_type = row["Cluster_Type"]
        color = CLUSTER_COLORS.get(cluster_type, "#95a5a6")
        label = f"Cluster {int(row['Cluster'])}: {cluster_type}"

        ax.plot(angles, values, "o-", linewidth=2.5, label=label, color=color)
        ax.fill(angles, values, alpha=0.25, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=11)
    ax.set_ylim(-2, 2)  # Standardized data range
    ax.set_title(
        "Land Use Cluster Profiles\n(Standardized Features)",
        size=16,
        fontweight="bold",
        pad=20,
    )
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=10)
    ax.grid(True, alpha=0.3)

    radar_path = os.path.join(output_dir, f"{k}_NLCD_Cluster_Radar_Chart.png")
    plt.savefig(radar_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Radar chart saved to: {radar_path}")


def create_heatmap(profiles_df, feature_columns, output_dir, k):
    """Create heatmap for cluster profiles"""
    print("Creating heatmap...")

    # Prepare data for heatmap
    heatmap_data = profiles_df[[f"{col}_mean" for col in feature_columns]].T
    heatmap_data.columns = [
        f"C{i}: {profiles_df.iloc[i]['Cluster_Type']}" for i in range(len(profiles_df))
    ]
    heatmap_data.index = [col.replace("_", "\n") for col in feature_columns]

    plt.figure(figsize=(12, 6))
    im = plt.imshow(heatmap_data, cmap="RdYlGn_r", aspect="auto", vmin=-2, vmax=2)
    plt.colorbar(im, label="Standardized Mean", fraction=0.046, pad=0.04)
    plt.xticks(
        range(len(heatmap_data.columns)), heatmap_data.columns, rotation=45, ha="right"
    )
    plt.yticks(range(len(heatmap_data.index)), heatmap_data.index)
    plt.title(
        "Land Use Cluster Profiles Heatmap\n(Standardized Features)",
        fontweight="bold",
        fontsize=14,
        pad=10,
    )

    # Add value annotations
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            text_color = "white" if abs(heatmap_data.iloc[i, j]) > 1 else "black"
            plt.text(
                j,
                i,
                f"{heatmap_data.iloc[i, j]:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontsize=9,
                fontweight="bold",
            )

    heatmap_path = os.path.join(output_dir, f"{k}_NLCD_Cluster_Heatmap.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to: {heatmap_path}")


def create_swarm_plot(df, feature_columns, output_dir, k):
    """Create swarm plot for features by cluster"""
    print("Creating swarm plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=feature_columns,
        var_name="Feature",
        value_name="Standardized_Value",
    )
    melted_df["Feature"] = melted_df["Feature"].str.replace("_", "\n")

    # Create color palette for clusters based on their types
    cluster_types = df["Cluster"].map(
        lambda c: df[df["Cluster"] == c].iloc[0]["Cluster_Type"]
        if "Cluster_Type" in df.columns
        else "Mixed"
    )
    unique_clusters = sorted(df["Cluster"].unique())
    palette = [
        CLUSTER_COLORS.get(
            df[df["Cluster"] == c].iloc[0].get("Cluster_Type", "Mixed"), "#95a5a6"
        )
        for c in unique_clusters
    ]

    plt.figure(figsize=(14, 8))
    sns.swarmplot(
        data=melted_df,
        x="Feature",
        y="Standardized_Value",
        hue="Cluster",
        palette=palette,
        dodge=True,
        size=3,
        alpha=0.7,
    )
    plt.title(
        "Land Use Features by Cluster\n(Swarm Plot)", fontsize=14, fontweight="bold"
    )
    plt.xlabel("Land Use Feature", fontsize=12)
    plt.ylabel("Standardized Value", fontsize=12)
    plt.xticks(rotation=0)
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()

    swarm_path = os.path.join(output_dir, f"{k}_NLCD_Cluster_Swarm_Plot.png")
    plt.savefig(swarm_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Swarm plot saved to: {swarm_path}")


def create_box_plot(df, feature_columns, output_dir, k):
    """Create box plot for features by cluster"""
    print("Creating box plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=feature_columns,
        var_name="Feature",
        value_name="Standardized_Value",
    )
    melted_df["Feature"] = melted_df["Feature"].str.replace("_", "\n")

    # Create color palette for clusters based on their types
    unique_clusters = sorted(df["Cluster"].unique())
    palette = [
        CLUSTER_COLORS.get(
            df[df["Cluster"] == c].iloc[0].get("Cluster_Type", "Mixed"), "#95a5a6"
        )
        for c in unique_clusters
    ]

    plt.figure(figsize=(14, 8))
    sns.boxplot(
        data=melted_df,
        x="Feature",
        y="Standardized_Value",
        hue="Cluster",
        palette=palette,
        showfliers=False,
    )
    plt.title(
        "Land Use Features by Cluster\n(Box Plot)", fontsize=14, fontweight="bold"
    )
    plt.xlabel("Land Use Feature", fontsize=12)
    plt.ylabel("Standardized Value", fontsize=12)
    plt.xticks(rotation=0)
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()

    box_path = os.path.join(output_dir, f"{k}_NLCD_Cluster_Box_Plot.png")
    plt.savefig(box_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Box plot saved to: {box_path}")


def create_raincloud_plot(df, feature_columns, output_dir, k):
    """Create raincloud plot (violin + swarm) for features by cluster"""
    print("Creating raincloud plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=feature_columns,
        var_name="Feature",
        value_name="Standardized_Value",
    )

    # Create subplots for each feature
    features = melted_df["Feature"].unique()
    n_features = len(features)
    fig, axes = plt.subplots(n_features, 1, figsize=(12, 5 * n_features), sharex=True)
    if n_features == 1:
        axes = [axes]

    for i, feature in enumerate(features):
        ax = axes[i]
        feature_data = melted_df[melted_df["Feature"] == feature]

        # Create color palette for clusters based on their types
        unique_clusters = sorted(df["Cluster"].unique())
        palette = [
            CLUSTER_COLORS.get(
                df[df["Cluster"] == c].iloc[0].get("Cluster_Type", "Mixed"), "#95a5a6"
            )
            for c in unique_clusters
        ]

        # Violin plot
        sns.violinplot(
            data=feature_data,
            x="Cluster",
            y="Standardized_Value",
            ax=ax,
            palette=palette,
            inner=None,
            alpha=0.6,
        )
        # Swarm plot overlay
        sns.swarmplot(
            data=feature_data,
            x="Cluster",
            y="Standardized_Value",
            ax=ax,
            color="black",
            alpha=0.4,
            size=2,
        )

        ax.set_title(f"{feature.replace('_', ' ')}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Cluster", fontsize=11)
        ax.set_ylabel("Standardized Value", fontsize=11)
        ax.grid(True, alpha=0.3, axis="y")

    plt.suptitle(
        "Land Use Features Distribution by Cluster\n(Raincloud Plot)",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()

    raincloud_path = os.path.join(output_dir, f"{k}_NLCD_Cluster_Raincloud_Plot.png")
    plt.savefig(raincloud_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Raincloud plot saved to: {raincloud_path}")


def create_map_for_k(df_clusters, k, profiles_df, shapefile_path, output_dir):
    """Create and save choropleth map for given k"""
    print(f"\nCreating map for k={k}...")

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

    # Map cluster numbers to colors based on cluster type
    color_map = {}
    for idx, row in profiles_df.iterrows():
        cluster_num = int(row["Cluster"])
        cluster_type = row["Cluster_Type"]
        color_map[cluster_num] = CLUSTER_COLORS.get(cluster_type, "#95a5a6")

    # Apply colors, handle NaN
    def get_color(cluster):
        if pd.isna(cluster):
            return "#e0e0e0"  # Light gray for no data
        return color_map.get(int(cluster), "#95a5a6")

    counties_merged["cluster_color"] = counties_merged[f"cluster_{k}"].apply(get_color)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(18, 12))
    counties_merged.plot(
        color=counties_merged["cluster_color"], linewidth=0.1, edgecolor="black", ax=ax
    )

    # State boundaries
    state_boundaries = counties_merged.dissolve(by="STATEFP")
    state_boundaries.boundary.plot(ax=ax, color="black", linewidth=1.5, alpha=0.9)

    ax.set_axis_off()
    ax.set_title(
        f"Land Use Clusters (k={k})\nContiguous United States Counties (2000-2005 Average)",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )

    # Legend
    legend_elements = []
    for idx, row in profiles_df.iterrows():
        cluster_num = int(row["Cluster"])
        cluster_type = row["Cluster_Type"]
        count = int(row["Count"])
        label = f"Cluster {cluster_num}: {cluster_type} (n={count})"
        legend_elements.append(
            mpatches.Patch(color=color_map[cluster_num], label=label)
        )

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        fontsize=11,
        frameon=True,
        fancybox=True,
        shadow=True,
    )

    plt.tight_layout()

    # Save
    output_filename = os.path.join(output_dir, f"{k}_NLCD_Cluster_Map.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Map saved to: {output_filename}")
    plt.close()


def save_feature_data(df, feature_columns, output_dir):
    """Save feature-engineered dataset"""
    print("\nSaving feature-engineered dataset...")

    # Select relevant columns
    output_cols = ["COUNTY_FIPS"] + feature_columns + ["ratio_sum"]
    df_output = df[output_cols].copy()

    # Ensure COUNTY_FIPS is properly formatted
    df_output["COUNTY_FIPS"] = df_output["COUNTY_FIPS"].astype(str).str.zfill(5)

    output_path = os.path.join(output_dir, "LandUse_Clustering_Features.csv")
    df_output.to_csv(output_path, index=False)
    print(f"Feature data saved to: {output_path}")


def save_final_results(df, output_dir):
    """Save final clustered dataset"""
    print("\nSaving final clustered results...")

    # Ensure COUNTY_FIPS is properly formatted
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    final_path = os.path.join(output_dir, "LandUse_Clusters_All_K.csv")
    df.to_csv(final_path, index=False)
    print(f"Final clustered data saved to: {final_path}")


def main():
    """Main clustering analysis pipeline for land use data (k=3 to 10)"""
    print("=" * 80)
    print("Starting Land Use Clustering Analysis (NLCD/JRC Data)")
    print("=" * 80)

    # Configuration
    config = get_config()
    base_path = config["data_directories"]["processed"]
    project_root = Path(__file__).resolve().parents[2]

    output_dir_base = os.path.join(
        project_root, "Result", "Cluster_Visualization_LandUse"
    )
    os.makedirs(output_dir_base, exist_ok=True)

    # File paths
    nlcd_file = os.path.join(base_path, "Environmental", "NLCD_JRC.csv")
    shapefile_path = os.path.join(
        project_root, config["data_sources"]["tiger"]["shapefile"]
    )

    # Step 1: Load and prepare data with feature engineering
    df, feature_columns = load_and_prepare_nlcd_data(nlcd_file)

    # Save feature-engineered dataset
    save_feature_data(df, feature_columns, output_dir_base)

    # Step 2: Standardize data
    df_standardized, scaler = standardize_data(df, feature_columns)

    # Step 3: Evaluate clusters (for reference, but we'll only process k=4)
    k_range = range(3, 11)
    X = df_standardized[feature_columns].values
    wcss, silhouette_scores = evaluate_clusters(X, k_range, output_dir_base)

    # Step 4: Perform clustering for k=4 only (optimal based on silhouette score)
    k = 4
    cluster_results = {}
    profiles_list = []

    if True:  # Keep indentation consistent
        print("\n" + "=" * 80)
        print(f"Processing k={k}")
        print("=" * 80)

        # Perform clustering
        labels, kmeans = perform_clustering(X, k)

        # Profile clusters
        df_clustered = df_standardized.copy()
        df_clustered["Cluster"] = labels
        df_clustered, profiles_df = profile_clusters(
            df_clustered, feature_columns, labels
        )

        # Sort clusters by human intervention level
        profiles_df, old_to_new = sort_clusters_by_intervention(
            profiles_df, feature_columns
        )

        # Reassign cluster labels and add cluster type to dataframe
        df_clustered["Cluster"] = df_clustered["Cluster"].map(old_to_new)

        # Map cluster type to the main dataframe for visualization
        cluster_type_map = profiles_df.set_index("Cluster")["Cluster_Type"].to_dict()
        df_clustered["Cluster_Type"] = df_clustered["Cluster"].map(cluster_type_map)

        labels = df_clustered["Cluster"].values

        # Update profiles_df
        profiles_df = profiles_df.drop(columns=["old_cluster"])
        profiles_path = os.path.join(output_dir_base, f"{k}_NLCD_Cluster_Profiles.csv")
        profiles_df.to_csv(profiles_path, index=False)
        print(f"Cluster profiles saved to: {profiles_path}")

        cluster_results[k] = labels
        profiles_df["k"] = k
        profiles_list.append(profiles_df)

        # Create visualizations for k=4
        create_radar_chart(profiles_df, feature_columns, output_dir_base, k)
        create_heatmap(profiles_df, feature_columns, output_dir_base, k)
        create_swarm_plot(df_clustered, feature_columns, output_dir_base, k)
        create_box_plot(df_clustered, feature_columns, output_dir_base, k)
        create_raincloud_plot(df_clustered, feature_columns, output_dir_base, k)

        # Create map for k=4
        df_clusters = df[["COUNTY_FIPS"]].copy()
        df_clusters[f"cluster_{k}"] = labels
        df_clusters["COUNTY_FIPS"] = df_clusters["COUNTY_FIPS"].astype(str).str.zfill(5)
        create_map_for_k(df_clusters, k, profiles_df, shapefile_path, output_dir_base)

    # Step 5: Save final results
    df_all_clusters = df[["COUNTY_FIPS"]].copy()
    df_all_clusters[f"cluster_{k}"] = cluster_results[k]

    save_final_results(df_all_clusters, output_dir_base)

    # Save combined profiles
    all_profiles = pd.concat(profiles_list, ignore_index=True)
    all_profiles_path = os.path.join(output_dir_base, "All_K_NLCD_Cluster_Profiles.csv")
    all_profiles.to_csv(all_profiles_path, index=False)
    print(f"\nCombined profiles saved to: {all_profiles_path}")

    print("\n" + "=" * 80)
    print("Land Use Clustering Analysis Completed!")
    print("=" * 80)
    print(f"All results saved to: {output_dir_base}")


if __name__ == "__main__":
    main()
