# Code/Analysis/CountyTypology_EQI_Analysis.py

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import stats

# Define county typology color scheme (matching Map_CountyTypology.py)
TYPOLOGY_COLORS = {
    "Farming": "#4DAF4A",  # Natural green
    "Mining": "#984EA3",  # Purple
    "Manufacturing": "#377EB8",  # Reliable blue
    "Government": "#E41A1C",  # Authoritative red
    "Services": "#FF7F00",  # Orange
    "Nonspecialized": "#CCCCCC",  # Gray
}

# Map typology codes to names
TYPOLOGY_MAP = {
    1: "Farming",
    2: "Mining",
    3: "Manufacturing",
    4: "Government",
    5: "Services",
    6: "Nonspecialized",
}


def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def load_and_merge_data(typology_file, eqi_file):
    """Load county typology and EQI data and merge on COUNTY_FIPS"""
    print("Loading data files...")

    # Load county typology
    typology = pd.read_csv(typology_file)
    print(f"County Typology: {typology.shape[0]} counties")

    # Load EQI data
    eqi = pd.read_csv(eqi_file)
    print(f"EQI: {eqi.shape[0]} counties")

    # Ensure COUNTY_FIPS is string and zero-padded
    typology["COUNTY_FIPS"] = typology["COUNTY_FIPS"].astype(str).str.zfill(5)
    eqi["COUNTY_FIPS"] = eqi["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Merge datasets
    merged = typology.merge(eqi, on="COUNTY_FIPS", how="inner")
    print(f"Merged: {merged.shape[0]} counties")

    # Add typology type labels
    merged["Typology_Type"] = merged["econdep"].map(TYPOLOGY_MAP)

    print("\nCounty typology distribution:")
    print(merged["Typology_Type"].value_counts().sort_index())

    return merged


def calculate_summary_statistics(df):
    """Calculate mean EQI scores by county typology"""
    print("\nCalculating summary statistics...")

    # EQI domains to analyze
    eqi_domains = ["EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]

    # Group by typology type and calculate statistics
    summary = []

    for typology_type in [
        "Farming",
        "Mining",
        "Manufacturing",
        "Government",
        "Services",
        "Nonspecialized",
    ]:
        typology_data = df[df["Typology_Type"] == typology_type]
        n = len(typology_data)

        row = {
            "Typology_Type": typology_type,
            "N": n,
        }

        for domain in eqi_domains:
            row[f"{domain}_Mean"] = typology_data[domain].mean()
            row[f"{domain}_SD"] = typology_data[domain].std()
            row[f"{domain}_Median"] = typology_data[domain].median()
            row[f"{domain}_Q25"] = typology_data[domain].quantile(0.25)
            row[f"{domain}_Q75"] = typology_data[domain].quantile(0.75)

        summary.append(row)

    summary_df = pd.DataFrame(summary)

    return summary_df, eqi_domains


def perform_statistical_tests(df, eqi_domains):
    """Perform ANOVA and post-hoc tests for EQI differences across county types"""
    print("\nPerforming statistical tests...")

    results = []

    for domain in eqi_domains:
        # Separate data by typology type
        groups = [
            df[df["Typology_Type"] == "Farming"][domain].dropna(),
            df[df["Typology_Type"] == "Mining"][domain].dropna(),
            df[df["Typology_Type"] == "Manufacturing"][domain].dropna(),
            df[df["Typology_Type"] == "Government"][domain].dropna(),
            df[df["Typology_Type"] == "Services"][domain].dropna(),
            df[df["Typology_Type"] == "Nonspecialized"][domain].dropna(),
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
    """Perform pairwise t-tests between county types"""
    print("\nPerforming pairwise comparisons...")

    typology_types = [
        "Farming",
        "Mining",
        "Manufacturing",
        "Government",
        "Services",
        "Nonspecialized",
    ]
    pairwise_results = []

    for domain in eqi_domains:
        for i, type1 in enumerate(typology_types):
            for type2 in typology_types[i + 1 :]:
                group1 = df[df["Typology_Type"] == type1][domain].dropna()
                group2 = df[df["Typology_Type"] == type2][domain].dropna()

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


def create_boxplot_by_typology(df, eqi_domains, output_dir):
    """Create boxplots showing EQI distribution by county typology"""
    print("\nCreating boxplot visualization...")

    # Create color palette based on typology types
    typology_order = [
        "Farming",
        "Mining",
        "Manufacturing",
        "Government",
        "Services",
        "Nonspecialized",
    ]
    palette = [TYPOLOGY_COLORS[ct] for ct in typology_order]

    # Create subplots for each EQI domain
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.flatten()

    for idx, domain in enumerate(eqi_domains):
        ax = axes[idx]

        sns.boxplot(
            data=df,
            x="Typology_Type",
            y=domain,
            order=typology_order,
            palette=palette,
            ax=ax,
            showfliers=False,
        )

        # Add individual points
        sns.stripplot(
            data=df,
            x="Typology_Type",
            y=domain,
            order=typology_order,
            color="black",
            alpha=0.2,
            size=1.5,
            ax=ax,
        )

        ax.set_title(domain.replace("_", " "), fontsize=14, fontweight="bold")
        ax.set_xlabel("County Economic Type", fontsize=11)
        ax.set_ylabel("EQI Score (Standardized)", fontsize=11)
        ax.tick_params(axis="x", rotation=45)
        ax.grid(True, alpha=0.3, axis="y")
        ax.axhline(y=0, color="red", linestyle="--", linewidth=1, alpha=0.5)

    plt.suptitle(
        "Environmental Quality Index by County Economic Typology (USDA ERS 2004)",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_CountyTypology_Boxplot.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Boxplot saved to: {output_path}")


def create_heatmap_summary(summary_df, eqi_domains, output_dir):
    """Create heatmap of mean EQI scores by county typology"""
    print("\nCreating heatmap visualization...")

    # Prepare data for heatmap
    heatmap_data = summary_df.set_index("Typology_Type")[
        [f"{domain}_Mean" for domain in eqi_domains]
    ]
    heatmap_data.columns = [
        col.replace("_Mean", "").replace("_", " ") for col in heatmap_data.columns
    ]

    # Create heatmap
    plt.figure(figsize=(10, 7))
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
        "Mean EQI Scores by County Economic Typology",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    plt.xlabel("EQI Domain", fontsize=12)
    plt.ylabel("County Economic Type", fontsize=12)
    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_CountyTypology_Heatmap.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to: {output_path}")


def create_violin_plot(df, eqi_domains, output_dir):
    """Create violin plots for EQI domains by county typology"""
    print("\nCreating violin plot visualization...")

    typology_order = [
        "Farming",
        "Mining",
        "Manufacturing",
        "Government",
        "Services",
        "Nonspecialized",
    ]
    palette = [TYPOLOGY_COLORS[ct] for ct in typology_order]

    # Melt dataframe for easier plotting
    melted = df.melt(
        id_vars=["Typology_Type"],
        value_vars=eqi_domains,
        var_name="EQI_Domain",
        value_name="EQI_Score",
    )
    melted["EQI_Domain"] = melted["EQI_Domain"].str.replace("_", " ")

    # Create plot
    fig, ax = plt.subplots(figsize=(18, 8))

    sns.violinplot(
        data=melted,
        x="EQI_Domain",
        y="EQI_Score",
        hue="Typology_Type",
        hue_order=typology_order,
        palette=palette,
        split=False,
        inner="quartile",
        ax=ax,
    )

    ax.set_title(
        "Distribution of EQI Scores by County Economic Typology",
        fontsize=16,
        fontweight="bold",
        pad=15,
    )
    ax.set_xlabel("EQI Domain", fontsize=13)
    ax.set_ylabel("EQI Score (Standardized)", fontsize=13)
    ax.legend(
        title="Economic Type", bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10
    )
    ax.axhline(y=0, color="red", linestyle="--", linewidth=1.5, alpha=0.5)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    output_path = os.path.join(output_dir, "EQI_by_CountyTypology_Violin.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Violin plot saved to: {output_path}")


def create_summary_table_figure(summary_df, output_dir):
    """Create a nice table figure showing summary statistics"""
    print("\nCreating summary table figure...")

    # Select key statistics for display
    display_df = summary_df[
        [
            "Typology_Type",
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
        if col not in ["Typology_Type", "N"]:
            display_df[col] = display_df[col].round(3)

    # Rename columns
    display_df.columns = [
        "Economic Type",
        "N",
        "Total EQI",
        "Air",
        "Water",
        "Land",
        "Built",
        "Social",
    ]

    # Create figure
    fig, ax = plt.subplots(figsize=(15, 6))
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

    # Color rows by typology type
    for i, typology_type in enumerate(display_df["Economic Type"]):
        cell = table[(i + 1, 0)]
        cell.set_facecolor(TYPOLOGY_COLORS.get(typology_type, "#95a5a6"))
        cell.set_text_props(weight="bold", color="white")

    plt.title(
        "Mean EQI Scores by County Economic Typology (USDA ERS 2004)",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    output_path = os.path.join(output_dir, "EQI_by_CountyTypology_Summary_Table.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Summary table saved to: {output_path}")


def create_markdown_report(summary_df, anova_results, pairwise_results, output_dir):
    """Generate a markdown report of the analysis"""
    print("\nCreating markdown report...")

    report_lines = [
        "# County Economic Typology × EQI Analysis Report",
        "",
        "## Summary Statistics",
        "",
        "### Mean EQI Scores by County Economic Type",
        "",
    ]

    # Add summary table
    summary_display = summary_df[
        [
            "Typology_Type",
            "N",
            "EQI_Mean",
            "EQI_Air_Mean",
            "EQI_Water_Mean",
            "EQI_Land_Mean",
            "EQI_Built_Mean",
            "EQI_Social_Mean",
        ]
    ].copy()

    report_lines.append(
        "| Economic Type | N | Total EQI | Air | Water | Land | Built | Social |"
    )
    report_lines.append(
        "|---------------|---|-----------|-----|-------|------|-------|--------|"
    )

    for _, row in summary_display.iterrows():
        report_lines.append(
            f"| {row['Typology_Type']} | {row['N']} | "
            f"{row['EQI_Mean']:.3f} | {row['EQI_Air_Mean']:.3f} | "
            f"{row['EQI_Water_Mean']:.3f} | {row['EQI_Land_Mean']:.3f} | "
            f"{row['EQI_Built_Mean']:.3f} | {row['EQI_Social_Mean']:.3f} |"
        )

    report_lines.extend(
        [
            "",
            "## Statistical Tests",
            "",
            "### ANOVA Results",
            "",
            "| Domain | F-statistic | P-value | η² | Significant |",
            "|--------|-------------|---------|-----|-------------|",
        ]
    )

    for _, row in anova_results.iterrows():
        report_lines.append(
            f"| {row['Domain']} | {row['F_statistic']:.2f} | "
            f"{row['P_value']:.6f} | {row['Eta_squared']:.3f} | {row['Significant']} |"
        )

    report_lines.extend(
        [
            "",
            "### Key Findings",
            "",
        ]
    )

    # Add interpretation
    for _, row in anova_results.iterrows():
        if row["Significant"] == "Yes":
            effect_size = (
                "large"
                if row["Eta_squared"] > 0.14
                else "medium"
                if row["Eta_squared"] > 0.06
                else "small"
            )
            report_lines.append(
                f"- **{row['Domain']}**: Significant differences across county types "
                f"(F={row['F_statistic']:.2f}, p<0.001, {effect_size} effect size η²={row['Eta_squared']:.3f})"
            )

    report_lines.extend(
        [
            "",
            "### Pairwise Comparisons (Selected Significant Results)",
            "",
            "| Domain | Comparison | Mean Diff | Cohen's d | P-value |",
            "|--------|------------|-----------|-----------|---------|",
        ]
    )

    # Show only significant results with medium+ effect size
    significant_pairs = (
        pairwise_results[
            (pairwise_results["Significant"] == "Yes")
            & (abs(pairwise_results["Cohens_d"]) > 0.5)
        ]
        .sort_values("Cohens_d", ascending=False)
        .head(20)
    )

    for _, row in significant_pairs.iterrows():
        report_lines.append(
            f"| {row['Domain']} | {row['Comparison']} | "
            f"{row['Mean_Diff']:.3f} | {row['Cohens_d']:.3f} | {row['P_value']:.6f} |"
        )

    report_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis examines the relationship between county economic dependency types ",
            "(as classified by USDA ERS 2004) and environmental quality across multiple domains.",
            "",
            "**Key Observations:**",
            "",
        ]
    )

    # Find types with highest/lowest EQI
    eqi_summary = summary_df.sort_values("EQI_Mean")
    report_lines.extend(
        [
            f"- Counties with the **lowest** overall environmental quality: **{eqi_summary.iloc[0]['Typology_Type']}** (mean EQI = {eqi_summary.iloc[0]['EQI_Mean']:.3f})",
            f"- Counties with the **highest** overall environmental quality: **{eqi_summary.iloc[-1]['Typology_Type']}** (mean EQI = {eqi_summary.iloc[-1]['EQI_Mean']:.3f})",
            "",
        ]
    )

    report_path = os.path.join(output_dir, "CountyTypology_EQI_Analysis_Report.md")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"Markdown report saved to: {report_path}")


def main():
    """Main analysis pipeline"""
    print("=" * 80)
    print("County Economic Typology × EQI Analysis")
    print("=" * 80)

    # Configuration
    config = get_config()
    project_root = Path(__file__).resolve().parents[2]
    base_path = config["data_directories"]["processed"]

    # File paths
    typology_file = os.path.join(base_path, "Socioeconomic", "County_Typology_2004.csv")
    eqi_file = os.path.join(base_path, "EQI", "EQI0005_Standard.csv")
    output_dir = os.path.join(project_root, "Result", "CountyTypology_Analysis")
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Load and merge data
    df = load_and_merge_data(typology_file, eqi_file)

    # Step 2: Calculate summary statistics
    summary_df, eqi_domains = calculate_summary_statistics(df)

    # Save summary statistics
    summary_path = os.path.join(output_dir, "EQI_by_CountyTypology_Summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary statistics saved to: {summary_path}")

    # Step 3: Perform ANOVA tests
    anova_results = perform_statistical_tests(df, eqi_domains)

    # Save ANOVA results
    anova_path = os.path.join(output_dir, "EQI_by_CountyTypology_ANOVA.csv")
    anova_results.to_csv(anova_path, index=False)
    print(f"ANOVA results saved to: {anova_path}")

    # Step 4: Perform pairwise tests
    pairwise_results = perform_pairwise_tests(df, eqi_domains)

    # Save pairwise results
    pairwise_path = os.path.join(output_dir, "EQI_by_CountyTypology_Pairwise.csv")
    pairwise_results.to_csv(pairwise_path, index=False)
    print(f"Pairwise comparison results saved to: {pairwise_path}")

    # Step 5: Create visualizations
    create_boxplot_by_typology(df, eqi_domains, output_dir)
    create_heatmap_summary(summary_df, eqi_domains, output_dir)
    create_violin_plot(df, eqi_domains, output_dir)
    create_summary_table_figure(summary_df, output_dir)

    # Step 6: Create markdown report
    create_markdown_report(summary_df, anova_results, pairwise_results, output_dir)

    # Step 7: Save merged dataset
    merged_path = os.path.join(output_dir, "CountyTypology_EQI_Merged.csv")
    df.to_csv(merged_path, index=False)
    print(f"\nMerged dataset saved to: {merged_path}")

    print("\n" + "=" * 80)
    print("Analysis completed!")
    print("=" * 80)
    print(f"All results saved to: {output_dir}")


if __name__ == "__main__":
    main()
