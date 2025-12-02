# Code/Analysis/Cluster_NLCD_EQI_Analysis.py

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import stats

# Define cluster type color scheme (matching Cluster_NLCD.py)
CLUSTER_COLORS = {
    "Urban": "#E68785",  # Coral red
    "Water-Sensitive": "#699DCB",  # Sky blue
    "Natural": "#6AAA81",  # Sage green
    "Agricultural": "#EFC085",  # Golden wheat
    "Mixed": "#95a5a6",  # Gray
}

# Map cluster numbers to types (from k=4 results)
CLUSTER_TYPE_MAP = {
    0: "Natural",
    1: "Water-Sensitive",
    2: "Agricultural",
    3: "Urban",
}


def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_and_merge_data(cluster_file, eqi_file):
    """Load cluster and EQI data and merge on COUNTY_FIPS"""
    print("Loading data files...")

    # Load cluster assignments
    clusters = pd.read_csv(cluster_file)
    print(f"Clusters: {clusters.shape[0]} counties")

    # Load EQI data
    eqi = pd.read_csv(eqi_file)
    print(f"EQI: {eqi.shape[0]} counties")

    # Ensure COUNTY_FIPS is string and zero-padded
    clusters["COUNTY_FIPS"] = clusters["COUNTY_FIPS"].astype(str).str.zfill(5)
    eqi["COUNTY_FIPS"] = eqi["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Rename cluster column for clarity
    clusters = clusters.rename(columns={"cluster_4": "Cluster"})

    # Merge datasets
    merged = clusters.merge(eqi, on="COUNTY_FIPS", how="inner")
    print(f"Merged: {merged.shape[0]} counties")

    # Add cluster type labels
    merged["Cluster_Type"] = merged["Cluster"].map(CLUSTER_TYPE_MAP)

    print("\nCluster distribution:")
    print(merged["Cluster_Type"].value_counts().sort_index())

    return merged


def calculate_summary_statistics(df):
    """Calculate mean EQI scores by cluster type"""
    print("\nCalculating summary statistics...")

    # EQI domains to analyze
    eqi_domains = ["EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]

    # Group by cluster type and calculate statistics
    summary = []

    for cluster_type in ["Natural", "Water-Sensitive", "Agricultural", "Urban"]:
        cluster_data = df[df["Cluster_Type"] == cluster_type]
        n = len(cluster_data)

        row = {
            "Cluster_Type": cluster_type,
            "N": n,
        }

        for domain in eqi_domains:
            row[f"{domain}_Mean"] = cluster_data[domain].mean()
            row[f"{domain}_SD"] = cluster_data[domain].std()
            row[f"{domain}_Median"] = cluster_data[domain].median()
            row[f"{domain}_Q25"] = cluster_data[domain].quantile(0.25)
            row[f"{domain}_Q75"] = cluster_data[domain].quantile(0.75)

        summary.append(row)

    summary_df = pd.DataFrame(summary)

    return summary_df, eqi_domains


def perform_statistical_tests(df, eqi_domains):
    """Perform ANOVA and post-hoc tests for EQI differences across clusters"""
    print("\nPerforming statistical tests...")

    results = []

    for domain in eqi_domains:
        # Separate data by cluster type
        groups = [
            df[df["Cluster_Type"] == "Natural"][domain].dropna(),
            df[df["Cluster_Type"] == "Water-Sensitive"][domain].dropna(),
            df[df["Cluster_Type"] == "Agricultural"][domain].dropna(),
            df[df["Cluster_Type"] == "Urban"][domain].dropna(),
        ]

        # One-way ANOVA
        f_stat, p_value = stats.f_oneway(*groups)

        # Effect size (eta-squared)
        grand_mean = df[domain].mean()
        ss_between = sum(len(g) * (g.mean() - grand_mean) ** 2 for g in groups)
        ss_total = sum((df[domain] - grand_mean) ** 2)
        eta_squared = ss_between / ss_total if ss_total > 0 else 0

        result = {
            "Domain": domain,
            "F_statistic": f_stat,
            "P_value": p_value,
            "Eta_squared": eta_squared,
            "Significant": "Yes" if p_value < 0.05 else "No",
        }

        results.append(result)

        print(f"{domain}: F={f_stat:.2f}, p={p_value:.4f}, η²={eta_squared:.3f}")

    results_df = pd.DataFrame(results)

    return results_df


def perform_pairwise_tests(df, eqi_domains):
    """Perform pairwise t-tests between cluster types"""
    print("\nPerforming pairwise comparisons...")

    cluster_types = ["Natural", "Water-Sensitive", "Agricultural", "Urban"]
    pairwise_results = []

    for domain in eqi_domains:
        for i, type1 in enumerate(cluster_types):
            for type2 in cluster_types[i + 1 :]:
                group1 = df[df["Cluster_Type"] == type1][domain].dropna()
                group2 = df[df["Cluster_Type"] == type2][domain].dropna()

                # Independent t-test
                t_stat, p_value = stats.ttest_ind(group1, group2)

                # Cohen's d effect size
                pooled_std = np.sqrt(
                    (
                        (len(group1) - 1) * group1.std() ** 2
                        + (len(group2) - 1) * group2.std() ** 2
                    )
                    / (len(group1) + len(group2) - 2)
                )
                cohens_d = (
                    (group1.mean() - group2.mean()) / pooled_std
                    if pooled_std > 0
                    else 0
                )

                pairwise_results.append(
                    {
                        "Domain": domain,
                        "Comparison": f"{type1} vs {type2}",
                        "Mean_Diff": group1.mean() - group2.mean(),
                        "T_statistic": t_stat,
                        "P_value": p_value,
                        "Cohens_d": cohens_d,
                        "Significant": "Yes" if p_value < 0.05 else "No",
                    }
                )

    pairwise_df = pd.DataFrame(pairwise_results)

    return pairwise_df


def create_boxplot_by_cluster(df, eqi_domains, output_dir):
    """Create boxplots showing EQI distribution by cluster type"""
    print("\nCreating boxplot visualization...")

    # Create color palette based on cluster types
    cluster_order = ["Natural", "Water-Sensitive", "Agricultural", "Urban"]
    palette = [CLUSTER_COLORS[ct] for ct in cluster_order]

    # Create subplots for each EQI domain
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, domain in enumerate(eqi_domains):
        ax = axes[idx]

        sns.boxplot(
            data=df,
            x="Cluster_Type",
            y=domain,
            order=cluster_order,
            palette=palette,
            ax=ax,
            showfliers=False,
        )

        # Add individual points
        sns.stripplot(
            data=df,
            x="Cluster_Type",
            y=domain,
            order=cluster_order,
            color="black",
            alpha=0.3,
            size=2,
            ax=ax,
        )

        ax.set_title(domain.replace("_", " "), fontsize=14, fontweight="bold")
        ax.set_xlabel("Land Use Cluster Type", fontsize=11)
        ax.set_ylabel("EQI Score (Standardized)", fontsize=11)
        ax.tick_params(axis="x", rotation=45)
        ax.grid(True, alpha=0.3, axis="y")
        ax.axhline(y=0, color="red", linestyle="--", linewidth=1, alpha=0.5)

    plt.suptitle(
        "Environmental Quality Index by Land Use Cluster Type",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Boxplot.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Boxplot saved to: {output_path}")


def create_heatmap_summary(summary_df, eqi_domains, output_dir):
    """Create heatmap of mean EQI scores by cluster type"""
    print("\nCreating heatmap visualization...")

    # Prepare data for heatmap
    heatmap_data = summary_df.set_index("Cluster_Type")[
        [f"{domain}_Mean" for domain in eqi_domains]
    ]
    heatmap_data.columns = [
        col.replace("_Mean", "").replace("_", " ") for col in heatmap_data.columns
    ]

    # Create heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".3f",
        cmap="RdYlGn_r",
        center=0,
        cbar_kws={"label": "Mean EQI Score"},
        linewidths=1,
        linecolor="white",
    )

    plt.title(
        "Mean EQI Scores by Land Use Cluster Type",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    plt.xlabel("EQI Domain", fontsize=12)
    plt.ylabel("Land Use Cluster Type", fontsize=12)
    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Heatmap.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to: {output_path}")


def create_violin_plot(df, eqi_domains, output_dir):
    """Create violin plots for EQI domains by cluster type"""
    print("\nCreating violin plot visualization...")

    cluster_order = ["Natural", "Water-Sensitive", "Agricultural", "Urban"]
    palette = [CLUSTER_COLORS[ct] for ct in cluster_order]

    # Melt dataframe for easier plotting
    melted = df.melt(
        id_vars=["Cluster_Type"],
        value_vars=eqi_domains,
        var_name="EQI_Domain",
        value_name="EQI_Score",
    )
    melted["EQI_Domain"] = melted["EQI_Domain"].str.replace("_", " ")

    # Create plot
    fig, ax = plt.subplots(figsize=(16, 8))

    sns.violinplot(
        data=melted,
        x="EQI_Domain",
        y="EQI_Score",
        hue="Cluster_Type",
        hue_order=cluster_order,
        palette=palette,
        split=False,
        inner="quartile",
        ax=ax,
    )

    ax.set_title(
        "Distribution of EQI Scores by Land Use Cluster Type",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.set_xlabel("EQI Domain", fontsize=13)
    ax.set_ylabel("EQI Score (Standardized)", fontsize=13)
    ax.legend(title="Land Use Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.axhline(y=0, color="red", linestyle="--", linewidth=1.5, alpha=0.5)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Violin.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Violin plot saved to: {output_path}")


def create_summary_table_figure(summary_df, output_dir):
    """Create a nice table figure showing summary statistics"""
    print("\nCreating summary table figure...")

    # Select key statistics for display
    display_df = summary_df[
        [
            "Cluster_Type",
            "N",
            "EQI_Mean",
            "EQI_Air_Mean",
            "EQI_Water_Mean",
            "EQI_Land_Mean",
            "EQI_Built_Mean",
            "EQI_Social_Mean",
        ]
    ].copy()

    # Round values
    for col in display_df.columns:
        if col not in ["Cluster_Type", "N"]:
            display_df[col] = display_df[col].round(3)

    # Rename columns
    display_df.columns = [
        "Cluster Type",
        "N",
        "Total EQI",
        "Air",
        "Water",
        "Land",
        "Built",
        "Social",
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.axis("tight")
    ax.axis("off")

    table = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        cellLoc="center",
        loc="center",
        colWidths=[0.15, 0.08, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header
    for i in range(len(display_df.columns)):
        cell = table[(0, i)]
        cell.set_facecolor("#4472C4")
        cell.set_text_props(weight="bold", color="white")

    # Color rows by cluster type
    for i, cluster_type in enumerate(display_df["Cluster Type"]):
        cell = table[(i + 1, 0)]
        cell.set_facecolor(CLUSTER_COLORS.get(cluster_type, "#95a5a6"))
        cell.set_text_props(weight="bold", color="white")

    plt.title(
        "Mean EQI Scores by Land Use Cluster Type",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    output_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Summary_Table.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Summary table saved to: {output_path}")


def main():
    """Main analysis pipeline"""
    print("=" * 80)
    print("Land Use Cluster × EQI Analysis")
    print("=" * 80)

    # Configuration
    config = get_config()
    project_root = Path(__file__).resolve().parents[2]
    base_path = config["data_directories"]["processed"]

    # File paths
    cluster_file = os.path.join(
        project_root,
        "Result",
        "Cluster_Visualization_LandUse",
        "LandUse_Clusters_All_K.csv",
    )
    eqi_file = os.path.join(base_path, "EQI", "EQI0005_Standard.csv")
    output_dir = os.path.join(project_root, "Result", "Cluster_Visualization_LandUse")

    # Step 1: Load and merge data
    df = load_and_merge_data(cluster_file, eqi_file)

    # Step 2: Calculate summary statistics
    summary_df, eqi_domains = calculate_summary_statistics(df)

    # Save summary statistics
    summary_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary statistics saved to: {summary_path}")

    # Step 3: Perform ANOVA tests
    anova_results = perform_statistical_tests(df, eqi_domains)

    # Save ANOVA results
    anova_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_ANOVA.csv")
    anova_results.to_csv(anova_path, index=False)
    print(f"ANOVA results saved to: {anova_path}")

    # Step 4: Perform pairwise tests
    pairwise_results = perform_pairwise_tests(df, eqi_domains)

    # Save pairwise results
    pairwise_path = os.path.join(output_dir, "EQI_by_LandUse_Cluster_Pairwise.csv")
    pairwise_results.to_csv(pairwise_path, index=False)
    print(f"Pairwise comparison results saved to: {pairwise_path}")

    # Step 5: Create visualizations
    create_boxplot_by_cluster(df, eqi_domains, output_dir)
    create_heatmap_summary(summary_df, eqi_domains, output_dir)
    create_violin_plot(df, eqi_domains, output_dir)
    create_summary_table_figure(summary_df, output_dir)

    # Step 6: Save merged dataset
    merged_path = os.path.join(output_dir, "LandUse_Cluster_EQI_Merged.csv")
    df.to_csv(merged_path, index=False)
    print(f"\nMerged dataset saved to: {merged_path}")

    print("\n" + "=" * 80)
    print("Analysis completed!")
    print("=" * 80)
    print(f"All results saved to: {output_dir}")


if __name__ == "__main__":
    main()
