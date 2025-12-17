# Code/Analysis/Cluster.py

import os
from math import pi
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler


def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_and_prepare_data(standard_file_path):
    """Load standardized EQI data and prepare for clustering"""
    print(f"Loading data from: {standard_file_path}")
    df = pd.read_csv(standard_file_path)

    # Select relevant columns
    eqi_columns = ["EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]
    required_cols = ["COUNTY_FIPS"] + eqi_columns

    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. Expected: {required_cols}")

    df_clean = df[required_cols].dropna()
    print(
        f"Data loaded: {df_clean.shape[0]} counties, {len(eqi_columns)} EQI dimensions"
    )

    return df_clean, eqi_columns


def standardize_data(df, eqi_columns):
    """Standardize EQI data for clustering"""
    print("Standardizing EQI data...")
    scaler = StandardScaler()
    df_standardized = df.copy()
    df_standardized[eqi_columns] = scaler.fit_transform(df[eqi_columns])
    return df_standardized, scaler


def evaluate_clusters(X, k_range, output_dir):
    """Evaluate different numbers of clusters using elbow method and silhouette score"""
    print("Evaluating optimal number of clusters...")

    wcss = []
    silhouette_scores = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(X)
        wcss.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(X, kmeans.labels_))

    # Elbow method plot
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(k_range, wcss, "bo-")
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Within-Cluster Sum of Squares (WCSS)")
    plt.title("Elbow Method")
    plt.grid(True)

    # Silhouette score plot
    plt.subplot(1, 2, 2)
    plt.plot(k_range, silhouette_scores, "ro-")
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Silhouette Score")
    plt.title("Silhouette Analysis")
    plt.grid(True)

    plt.tight_layout()
    eval_plot_path = os.path.join(output_dir, "Cluster_Evaluation.png")
    plt.savefig(eval_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Cluster evaluation plots saved to: {eval_plot_path}")

    # Find optimal K
    optimal_k_elbow = k_range[np.argmin(np.diff(wcss))]  # Elbow point
    optimal_k_silhouette = k_range[np.argmax(silhouette_scores)]

    print(f"Optimal K by Elbow method: {optimal_k_elbow}")
    print(f"Optimal K by Silhouette score: {optimal_k_silhouette}")

    # Choose final K (prefer silhouette, but consider interpretability)
    final_k = (
        optimal_k_silhouette if optimal_k_silhouette in range(3, 7) else optimal_k_elbow
    )
    if final_k not in range(3, 7):
        final_k = 4  # Default fallback

    print(f"Selected final K: {final_k}")
    return final_k


def perform_clustering(X, n_clusters):
    """Perform K-Means clustering"""
    print(f"Performing K-Means clustering with {n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    return labels, kmeans


def profile_clusters(df, eqi_columns, labels):
    """Create cluster profiles"""
    print("Creating cluster profiles...")

    df["Cluster"] = labels
    profiles = []

    for cluster in np.unique(labels):
        cluster_data = df[df["Cluster"] == cluster]
        profile = {"Cluster": cluster, "Count": len(cluster_data)}

        for col in eqi_columns:
            profile[f"{col}_mean"] = cluster_data[col].mean()
            profile[f"{col}_std"] = cluster_data[col].std()

        profiles.append(profile)

    profiles_df = pd.DataFrame(profiles)

    return df, profiles_df


def sort_clusters_by_radar_area(profiles_df, eqi_columns):
    """Sort clusters by radar chart area (larger area = worse environment)"""
    # Calculate radar area as sum of standardized means (positive = worse, negative = better)
    profiles_df["radar_area"] = profiles_df[[f"{col}_mean" for col in eqi_columns]].sum(
        axis=1
    )

    # Sort by radar area ascending (smaller sum = better environment = lower cluster number)
    profiles_df = profiles_df.sort_values("radar_area").reset_index(drop=True)

    # Reassign cluster labels from 0 (best) to k-1 (worst)
    old_to_new = {old: new for new, old in enumerate(profiles_df["Cluster"])}
    profiles_df["new_cluster"] = range(len(profiles_df))

    return profiles_df, old_to_new


def create_radar_chart(profiles_df, eqi_columns, output_dir, k):
    """Create radar chart for cluster profiles"""
    print("Creating radar chart...")

    # Prepare data
    categories = [col.replace("EQI_", "") for col in eqi_columns]
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection="polar"))

    angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
    angles += angles[:1]  # Close the loop

    for _, row in profiles_df.iterrows():
        values = [row[f"{col}_mean"] for col in eqi_columns]
        values += values[:1]  # Close the loop

        ax.plot(
            angles, values, "o-", linewidth=2, label=f"Cluster {int(row['Cluster'])}"
        )
        ax.fill(angles, values, alpha=0.25)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_ylim(-2, 2)  # Assuming standardized data
    ax.set_title("Cluster Profiles - Standardized EQI", size=16, fontweight="bold")
    ax.legend(loc="upper right", bbox_to_anchor=(1.2, 1.0))
    ax.grid(True)

    radar_path = os.path.join(output_dir, f"{k}_Cluster_Radar_Chart.png")
    plt.savefig(radar_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Radar chart saved to: {radar_path}")


def create_heatmap(profiles_df, eqi_columns, output_dir, k):
    """Create heatmap for cluster profiles"""
    print("Creating heatmap...")

    # Prepare data for heatmap
    heatmap_data = profiles_df[[f"{col}_mean" for col in eqi_columns]].T
    heatmap_data.columns = [f"Cluster {i}" for i in range(len(profiles_df))]
    heatmap_data.index = [col.replace("EQI_", "") for col in eqi_columns]

    plt.figure(figsize=(10, 6))
    plt.imshow(
        heatmap_data, cmap="RdYlBu", aspect="auto"
    )  # Changed to RdYlBu (not reversed)
    plt.colorbar(label="Standardized EQI")
    plt.xticks(range(len(heatmap_data.columns)), heatmap_data.columns)
    plt.yticks(range(len(heatmap_data.index)), heatmap_data.index)
    plt.title("Cluster Profiles Heatmap", fontweight="bold")

    # Add value annotations
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            plt.text(
                j,
                i,
                f"{heatmap_data.iloc[i, j]:.2f}",
                ha="center",
                va="center",
                color="black",
                fontsize=8,
            )

    heatmap_path = os.path.join(output_dir, f"{k}_Cluster_Heatmap.png")
    plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to: {heatmap_path}")


def create_swarm_plot(df, eqi_columns, output_dir, k):
    """Create swarm plot for EQI by cluster"""
    print("Creating swarm plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=eqi_columns,
        var_name="EQI_Dimension",
        value_name="EQI_Value",
    )
    melted_df["EQI_Dimension"] = melted_df["EQI_Dimension"].str.replace("EQI_", "")

    plt.figure(figsize=(12, 8))
    sns.swarmplot(
        data=melted_df,
        x="EQI_Dimension",
        y="EQI_Value",
        hue="Cluster",
        palette="Set1",
        dodge=True,
    )
    plt.title("Swarm Plot of EQI by Cluster and EQI Dimension")
    plt.xticks(rotation=45)
    plt.legend(title="Cluster")
    plt.grid(True, alpha=0.3)

    swarm_path = os.path.join(output_dir, f"{k}_Cluster_Swarm_Plot.png")
    plt.savefig(swarm_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Swarm plot saved to: {swarm_path}")


def create_box_plot(df, eqi_columns, output_dir, k):
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

    plt.figure(figsize=(12, 8))
    sns.boxplot(
        data=melted_df,
        x="EQI_Dimension",
        y="EQI_Value",
        hue="Cluster",
        palette="Set1",
        showfliers=False,
    )
    plt.title("Box Plot of EQI by Cluster and EQI Dimension")
    plt.xticks(rotation=45)
    plt.legend(title="Cluster")
    plt.grid(True, alpha=0.3)

    box_path = os.path.join(output_dir, f"{k}_Cluster_Box_Plot.png")
    plt.savefig(box_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Box plot saved to: {box_path}")


def create_raincloud_plot(df, eqi_columns, output_dir, k):
    """Create raincloud plot (violin + swarm) for EQI by cluster"""
    print("Creating raincloud plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=eqi_columns,
        var_name="EQI_Dimension",
        value_name="EQI_Value",
    )
    melted_df["EQI_Dimension"] = melted_df["EQI_Dimension"].str.replace("EQI_", "")

    # Create subplots for each EQI dimension
    dimensions = melted_df["EQI_Dimension"].unique()
    n_dims = len(dimensions)
    fig, axes = plt.subplots(n_dims, 1, figsize=(10, 6 * n_dims), sharex=False)
    if n_dims == 1:
        axes = [axes]

    for i, dim in enumerate(dimensions):
        ax = axes[i]
        dim_data = melted_df[melted_df["EQI_Dimension"] == dim]

        # Violin plot
        sns.violinplot(
            data=dim_data, x="Cluster", y="EQI_Value", ax=ax, palette="Set1", inner=None
        )
        # Swarm plot overlay
        sns.swarmplot(
            data=dim_data,
            x="Cluster",
            y="EQI_Value",
            ax=ax,
            color="black",
            alpha=0.6,
            size=3,
        )

        ax.set_title(f"Raincloud Plot: {dim}")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    raincloud_path = os.path.join(output_dir, f"{k}_Cluster_Raincloud_Plot.png")
    plt.savefig(raincloud_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Raincloud plot saved to: {raincloud_path}")


def create_map_for_k(df_clusters, k, shapefile_path, output_dir):
    """Create and save choropleth map for given k"""
    print(f"Creating map for k={k}...")

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

    # Apply colors, handle NaN
    def get_color(cluster):
        if pd.isna(cluster):
            return color_map["No Data"]
        return color_map.get(int(cluster), color_map["No Data"])

    counties_merged["cluster_color"] = counties_merged[f"cluster_{k}"].apply(get_color)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    counties_merged.plot(
        color=counties_merged["cluster_color"], linewidth=0.1, edgecolor="black", ax=ax
    )

    # State boundaries
    state_boundaries = counties_merged.dissolve(by="STATEFP")
    state_boundaries.boundary.plot(ax=ax, color="black", linewidth=1.2, alpha=0.9)

    ax.set_axis_off()
    ax.set_title(
        f"Environmental Quality Clusters (k={k})\nContiguous United States Counties",
        fontsize=18,
        pad=20,
    )

    # Legend
    legend_elements = [
        mpatches.Patch(color=color_map[i], label=cluster_names[i]) for i in range(k)
    ]
    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        fontsize=12,
        frameon=True,
    )

    plt.tight_layout()

    # Save
    output_filename = os.path.join(output_dir, f"{k}_Cluster_Map.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Map saved to: {output_filename}")
    plt.close()


def save_final_results(df, output_dir):
    """Save final clustered dataset"""
    print("Saving final results...")

    # Ensure COUNTY_FIPS is properly formatted as 5-digit string
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    final_path = os.path.join(output_dir, "EQI_Clusters.csv")
    df.to_csv(final_path, index=False)
    print(f"Final clustered data saved to: {final_path}")


def main():
    """Main clustering analysis pipeline for k=3 to 10"""
    print("Starting EQI Clustering Analysis for k=3 to 10...")

    # Configuration
    config = get_config()
    base_path = config["data_directories"]["processed"]
    project_root = Path(__file__).resolve().parents[2]
    output_dir_base = os.path.join(project_root, "Result", "Cluster_Visualization")
    map_output_dir = os.path.join(project_root, "Result", "Cluster_Visualization")
    os.makedirs(output_dir_base, exist_ok=True)
    os.makedirs(map_output_dir, exist_ok=True)

    # File paths
    standard_file = os.path.join(base_path, "EQI", "EQI0005_Standard.csv")
    shapefile_path = os.path.join(
        project_root, config["data_sources"]["tiger"]["shapefile"]
    )

    # Step 1: Load and prepare data
    df, eqi_columns = load_and_prepare_data(standard_file)

    # Step 2: Standardize data
    df_standardized, scaler = standardize_data(df, eqi_columns)

    # Step 3: Perform clustering; allow specifying a single k via CLI
    import argparse

    parser = argparse.ArgumentParser(description="EQI Clustering Analysis")
    parser.add_argument("--k", type=int, help="Run and plot only this k (e.g., 3)")
    args, _ = parser.parse_known_args()
    k_range = [args.k] if args.k is not None else list(range(3, 11))
    cluster_results = {}
    profiles_list = []

    for k in k_range:
        print(f"\nProcessing k={k}...")

        # Perform clustering
        X = df_standardized[eqi_columns].values
        labels, kmeans = perform_clustering(X, k)

        # Profile clusters
        df_clustered = df_standardized.copy()
        df_clustered["Cluster"] = labels
        df_clustered, profiles_df = profile_clusters(df_clustered, eqi_columns, labels)

        # Sort clusters by radar area
        profiles_df, old_to_new = sort_clusters_by_radar_area(profiles_df, eqi_columns)

        # Reassign cluster labels
        df_clustered["Cluster"] = df_clustered["Cluster"].map(old_to_new)
        labels = df_clustered["Cluster"].values

        # Update profiles_df
        profiles_df["Cluster"] = profiles_df["new_cluster"]
        profiles_df = profiles_df.drop(columns=["new_cluster"])
        profiles_path = os.path.join(output_dir_base, f"{k}_Cluster_Profiles.csv")
        profiles_df.to_csv(profiles_path, index=False)

        cluster_results[k] = labels
        profiles_df["k"] = k
        profiles_list.append(profiles_df)

        # Create visualizations for this k
        create_radar_chart(profiles_df, eqi_columns, output_dir_base, k)
        create_heatmap(profiles_df, eqi_columns, output_dir_base, k)
        create_swarm_plot(df_clustered, eqi_columns, output_dir_base, k)
        create_box_plot(df_clustered, eqi_columns, output_dir_base, k)
        create_raincloud_plot(df_clustered, eqi_columns, output_dir_base, k)

        # Create map for this k
        df_clusters = df[["COUNTY_FIPS"]].copy()
        df_clusters[f"cluster_{k}"] = labels
        df_clusters["COUNTY_FIPS"] = df_clusters["COUNTY_FIPS"].astype(str).str.zfill(5)
        create_map_for_k(df_clusters, k, shapefile_path, map_output_dir)

    # Step 4: Save final results for all k
    df_all_clusters = df[["COUNTY_FIPS"]].copy()
    for k in k_range:
        df_all_clusters[f"cluster_{k}"] = cluster_results[k]
    df_all_clusters["COUNTY_FIPS"] = (
        df_all_clusters["COUNTY_FIPS"].astype(str).str.zfill(5)
    )
    final_path = os.path.join(output_dir_base, "EQI_Clusters_All_K.csv")
    df_all_clusters.to_csv(final_path, index=False)
    print(f"Final clustered data for all k saved to: {final_path}")

    print("Clustering analysis completed for all k=3 to 10!")


if __name__ == "__main__":
    main()
