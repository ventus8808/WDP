import pandas as pd
import numpy as np
import os
import sys
import yaml
from sklearn.cluster import KMeans, AgglomerativeClustering, SpectralClustering, Birch
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(os.path.dirname(__file__))
from Cluster_Plot_Function import (
    create_combined_visualization,
    create_radar_chart,
    create_box_plot,
    create_map_for_k,
)


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


def perform_kmeans_clustering(X, n_clusters):
    """Perform K-Means clustering"""
    print(f"Performing K-Means clustering with {n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    return labels


def perform_agglomerative_clustering(X, n_clusters):
    """Perform Agglomerative clustering"""
    print(f"Performing Agglomerative clustering with {n_clusters} clusters...")
    agglo = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")
    labels = agglo.fit_predict(X)
    return labels


def perform_gmm_clustering(X, n_clusters):
    """Perform Gaussian Mixture Model clustering"""
    print(f"Performing GMM clustering with {n_clusters} clusters...")
    gmm = GaussianMixture(
        n_components=n_clusters, random_state=42, covariance_type="full"
    )
    labels = gmm.fit_predict(X)
    return labels


def perform_spectral_clustering(X, n_clusters):
    """Perform Spectral clustering"""
    print(f"Performing Spectral clustering with {n_clusters} clusters...")
    spectral = SpectralClustering(
        n_clusters=n_clusters, random_state=42, affinity="nearest_neighbors"
    )
    labels = spectral.fit_predict(X)
    return labels


def perform_birch_clustering(X, n_clusters):
    """Perform Birch clustering"""
    print(f"Performing Birch clustering with {n_clusters} clusters...")
    birch = Birch(n_clusters=n_clusters)
    labels = birch.fit_predict(X)
    return labels


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


def calculate_radar_area(means):
    """Calculate the geometric area of a radar chart polygon given mean values."""
    n = len(means)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = means * np.cos(angles)
    y = means * np.sin(angles)
    # Append first point to close the polygon
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    # Shoelace formula
    area = 0.5 * np.abs(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))
    return area


def sort_clusters_by_radar_area(profiles_df, eqi_columns):
    """Sort clusters by true radar chart area (larger area = worse environment)"""
    # Calculate EQI domain sum of standardized means (positive = worse, negative = better)
    profiles_df["EQI_domain_SUM"] = profiles_df[
        [f"{col}_mean" for col in eqi_columns]
    ].sum(axis=1)

    # Calculate true radar area using geometric polygon area
    profiles_df["radar_area"] = profiles_df[
        [f"{col}_mean" for col in eqi_columns]
    ].apply(lambda row: calculate_radar_area(row.values), axis=1)

    # Calculate radar area percent relative to maximum possible area (all dimensions = 1)
    max_area = calculate_radar_area([1.0] * len(eqi_columns))
    profiles_df["radar_area_percent"] = (profiles_df["radar_area"] / max_area) * 100

    # Sort by EQI_domain_SUM ascending (smaller sum = better environment = lower cluster number)
    profiles_df = profiles_df.sort_values("EQI_domain_SUM").reset_index(drop=True)

    # Reassign cluster labels from 0 (best) to k-1 (worst)
    old_to_new = {old: new for new, old in enumerate(profiles_df["Cluster"])}
    profiles_df["new_cluster"] = range(len(profiles_df))

    return profiles_df, old_to_new


def main():
    """Compare different clustering methods on EQI data"""
    print("Starting EQI Clustering Method Comparison...")

    # Configuration
    config = get_config()
    base_path = config["data_directories"]["processed"]
    project_root = Path(__file__).resolve().parents[2]
    output_dir = os.path.join(project_root, "Result", "Cluster_Method")
    os.makedirs(output_dir, exist_ok=True)

    # File paths
    standard_file = os.path.join(base_path, "EQI", "EQI0005_Standard.csv")
    shapefile_path = os.path.join(
        project_root, config["data_sources"]["tiger"]["shapefile"]
    )

    # Step 1: Load and prepare data
    df, eqi_columns = load_and_prepare_data(standard_file)

    # Step 2: Standardize data
    df_standardized, scaler = standardize_data(df, eqi_columns)

    # Step 3: Define clustering methods
    methods = {
        "K-means": perform_kmeans_clustering,
        "Agglomerative": perform_agglomerative_clustering,
        "GMM": perform_gmm_clustering,
        "Spectral": perform_spectral_clustering,
        "Birch": perform_birch_clustering,
    }

    k = 3  # Fixed number of clusters for detailed analysis
    k_range = range(3, 11)
    diagnostics = []
    silhouettes_dict = {}

    for method_name, cluster_func in methods.items():
        print(f"\nProcessing {method_name}...")
        silhouettes = []

        for k_val in k_range:
            print(f"  k={k_val}...")

            # Perform clustering
            X = df_standardized[eqi_columns].values
            labels = cluster_func(X, k_val)

            # Compute diagnostics
            silhouette = silhouette_score(X, labels)
            diagnostics.append(
                {"Method": method_name, "k": k_val, "Silhouette_Score": silhouette}
            )
            silhouettes.append(silhouette)

            # Detailed analysis only for k=3
            if k_val == k:
                # Profile clusters
                df_clustered = df_standardized.copy()
                df_clustered["Cluster"] = labels
                df_clustered, profiles_df = profile_clusters(
                    df_clustered, eqi_columns, labels
                )

                # Sort clusters by radar area
                profiles_df, old_to_new = sort_clusters_by_radar_area(
                    profiles_df, eqi_columns
                )

                # Reassign cluster labels
                df_clustered["Cluster"] = df_clustered["Cluster"].map(old_to_new)
                labels = df_clustered["Cluster"].values

                # Update profiles_df
                profiles_df["Cluster"] = profiles_df["new_cluster"]
                profiles_df = profiles_df.drop(columns=["new_cluster"])
                profiles_path = os.path.join(
                    output_dir, f"{method_name}_{k}_Profiles.csv"
                )
                profiles_df.to_csv(profiles_path, index=False)

                # Create visualizations
                radar_path = create_radar_chart(
                    profiles_df, eqi_columns, output_dir, method_name, k
                )
                box_path = create_box_plot(
                    df_clustered, eqi_columns, output_dir, method_name, k
                )

                # Create map
                df_clusters = df[["COUNTY_FIPS"]].copy()
                df_clusters[f"cluster_{k}"] = labels
                df_clusters["COUNTY_FIPS"] = (
                    df_clusters["COUNTY_FIPS"].astype(str).str.zfill(5)
                )
                map_path = create_map_for_k(
                    df_clusters, k, shapefile_path, output_dir, method_name
                )

                # Create combined visualization
                create_combined_visualization(
                    radar_path, box_path, map_path, output_dir, method_name, k
                )

        silhouettes_dict[method_name] = silhouettes
        print(f"Silhouette Scores for {method_name}: {silhouettes}")

    # Create silhouette comparison plot
    plt.figure(figsize=(10, 6))
    for method_name, silhouettes in silhouettes_dict.items():
        plt.plot(list(k_range), silhouettes, marker="o", label=method_name)

    plt.xlabel("Number of Clusters (k)")
    plt.ylabel("Silhouette Score")
    plt.title("Silhouette Scores Comparison Across Clustering Methods (k=3 to 10)")
    plt.legend()
    plt.grid(True)
    silhouette_plot_path = os.path.join(output_dir, "Silhouette_Scores_Comparison.png")
    plt.savefig(silhouette_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Silhouette comparison plot saved to: {silhouette_plot_path}")

    # Save diagnostics
    diagnostics_df = pd.DataFrame(diagnostics)
    diagnostics_path = os.path.join(output_dir, "Clustering_Diagnostics.csv")
    diagnostics_df.to_csv(diagnostics_path, index=False)
    print(f"Diagnostics saved to: {diagnostics_path}")

    print("Clustering method comparison completed!")


if __name__ == "__main__":
    main()
