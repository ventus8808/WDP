import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def load_config(project_root: Path) -> dict:
    """Load config.yaml from project root."""
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_paths(project_root: Path) -> dict:
    """Get input/output directories for Delta Cluster visualization."""
    result_dir = project_root / "Result" / "brms_delta_cluster"
    vis_dir = project_root / "Result" / "brms_delta_cluster_Visualization"
    vis_dir.mkdir(parents=True, exist_ok=True)
    return {"result": result_dir, "vis": vis_dir}


def parse_effect_cell(cell: str) -> Optional[Tuple[float, float, float, bool]]:
    """
    Parse effect string like "-1.54(-4.20, 1.11)***"
    Returns: (mrd, lower_ci, upper_ci, is_significant) or None
    """
    if pd.isna(cell) or not cell or str(cell).strip() == "0.00":
        return None

    s = str(cell).strip()
    # Check for significance markers
    has_stars = bool(re.search(r"\*+$", s))
    # Remove stars
    core = re.sub(r"\*+$", "", s).strip()

    # Parse: MRD(lower, upper) - handle both English and Chinese commas
    pattern = (
        r"^([+-]?\d+\.?\d*)\s*\(\s*([+-]?\d+\.?\d*)\s*[,，]\s*([+-]?\d+\.?\d*)\s*\)$"
    )
    m = re.match(pattern, core)
    if not m:
        return None

    mrd = float(m.group(1))
    lower = float(m.group(2))
    upper = float(m.group(3))

    return mrd, lower, upper, has_stars


def prepare_panel_data(
    df: pd.DataFrame, panel_type: str, cluster_id: Optional[int] = None
) -> pd.DataFrame:
    """
    Extract data for one panel (National or specific Cluster).

    Parameters:
    -----------
    df : DataFrame with all models for one cancer/lag/cluster
    panel_type : 'National' or 'Cluster'
    cluster_id : Cluster ID (0,1,2,... ) if panel_type='Cluster'

    Returns:
    --------
    DataFrame with columns: [Domain, Model, Lag, MRD_Improved, CI_Lower_Imp, CI_Upper_Imp,
                             Sig_Imp, MRD_Worsened, CI_Lower_Wor, CI_Upper_Wor, Sig_Wor]

    Order: 5y EQI + 5y domains + empty row + 10y EQI + 10y domains
    """
    # Filter df based on Model column
    if panel_type == "National":
        df = df[~df["Model"].str.startswith("Cluster_")]
    else:
        df = df[df["Model"].str.startswith(f"Cluster{cluster_id}_")]

    rows = []

    if panel_type == "National":
        # Include base models + EQI: EQI, Air, Water, Land, Built, Social (EQI at the beginning)
        base_models = ["EQI", "Air", "Water", "Land", "Built", "Social"]
        prefix = ""
    else:
        # Cluster stratum: Cluster_{cluster_id}_EQI, etc. (EQI at the beginning)
        base_models = [
            f"Cluster{cluster_id}_EQI",
            f"Cluster{cluster_id}_Air",
            f"Cluster{cluster_id}_Water",
            f"Cluster{cluster_id}_Land",
            f"Cluster{cluster_id}_Built",
            f"Cluster{cluster_id}_Social",
        ]
        prefix = f"Cluster{cluster_id}_"

    # Process each lag first, then each model within the lag
    lag5_eqi_rows = []
    lag5_domain_rows = []
    lag10_eqi_rows = []
    lag10_domain_rows = []

    for lag in [5, 10]:
        for model_name in base_models:
            model_row = df[(df["Model"] == model_name) & (df["Lag"] == lag)]
            if model_row.empty:
                continue

            # Determine domain label
            domain = model_name.replace(prefix, "")

            # Parse Improved and Worsened effects
            improved_str = model_row["MRD_Q_Improved"].iloc[0]
            worsened_str = model_row["MRD_Q_Worsened"].iloc[0]

            improved = parse_effect_cell(improved_str)
            worsened = parse_effect_cell(worsened_str)

            # Only include if at least one effect is non-null
            if improved is None and worsened is None:
                continue

            row_data = {"Domain": domain, "Model": model_name, "Lag": lag}

            if improved:
                row_data["MRD_Improved"] = improved[0]
                row_data["CI_Lower_Imp"] = improved[1]
                row_data["CI_Upper_Imp"] = improved[2]
                row_data["Sig_Imp"] = improved[3]
            else:
                row_data["MRD_Improved"] = np.nan
                row_data["CI_Lower_Imp"] = np.nan
                row_data["CI_Upper_Imp"] = np.nan
                row_data["Sig_Imp"] = False

            if worsened:
                row_data["MRD_Worsened"] = worsened[0]
                row_data["CI_Lower_Wor"] = worsened[1]
                row_data["CI_Upper_Wor"] = worsened[2]
                row_data["Sig_Wor"] = worsened[3]
            else:
                row_data["MRD_Worsened"] = np.nan
                row_data["CI_Lower_Wor"] = np.nan
                row_data["CI_Upper_Wor"] = np.nan
                row_data["Sig_Wor"] = False

            # Separate EQI from domain models
            if domain == "EQI":
                if lag == 5:
                    lag5_eqi_rows.append(row_data)
                else:
                    lag10_eqi_rows.append(row_data)
            else:
                if lag == 5:
                    lag5_domain_rows.append(row_data)
                else:
                    lag10_domain_rows.append(row_data)

    # Combine: 5y EQI + 5y domains + empty row + 10y EQI + 10y domains (EQI at top of each group)
    rows = (
        lag5_eqi_rows
        + lag5_domain_rows
        + [
            {
                "Domain": "",
                "Model": "",
                "Lag": None,
                "MRD_Improved": np.nan,
                "CI_Lower_Imp": np.nan,
                "CI_Upper_Imp": np.nan,
                "Sig_Imp": False,
                "MRD_Worsened": np.nan,
                "CI_Lower_Wor": np.nan,
                "CI_Upper_Wor": np.nan,
                "Sig_Wor": False,
            }
        ]
        + lag10_eqi_rows
        + lag10_domain_rows
    )

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def plot_bidirectional_forest(
    df: pd.DataFrame,
    icd_code: str,
    k: int,
    output_dir: Path,
):
    """
    Create bidirectional forest plot with panels for National and each cluster (0 to k-1).
    Each panel shows single-domain models for both 5-year and 10-year lags.

    Layout:
    - Left side: Improved (negative MRD = beneficial)
    - Right side: Worsened (positive MRD = harmful)
    - Middle: Model labels with lag (e.g., EQI (5y), EQI (10y))
    - Colors: Blue for beneficial direction, Orange for harmful direction
    - Alpha: Solid (1.0) for significant, Semi-transparent (0.4) for non-significant
    """
    # Set global font
    plt.rcParams["font.family"] = "Georgia"
    plt.rcParams["font.size"] = 12

    # Define panel configurations: National + clusters 0 to k-1
    panels = [("National", None, "National Level")] + [
        ("Cluster", i, f"Cluster {i}") for i in range(k)
    ]

    # Create figure with (k+1) horizontal panels
    fig, axes = plt.subplots(
        k + 1, 1, figsize=(18.0, 4.0 * (k + 1)), constrained_layout=True
    )

    # Color scheme: by environmental quality change direction
    color_improved = "#1976D2"  # Blue for quality improved
    color_worsened = "#F57C00"  # Orange for quality worsened

    # Calculate global x_limit based on all panels
    global_all_values = []
    for panel_type, cluster_id, _ in panels:
        panel_df = prepare_panel_data(df, panel_type, cluster_id)
        if not panel_df.empty:
            for col in [
                "MRD_Improved",
                "CI_Lower_Imp",
                "CI_Upper_Imp",
                "MRD_Worsened",
                "CI_Lower_Wor",
                "CI_Upper_Wor",
            ]:
                global_all_values.extend(panel_df[col].dropna().values)

    if global_all_values:
        global_max_abs = max(abs(min(global_all_values)), abs(max(global_all_values)))
        global_x_limit = global_max_abs * 1.2
    else:
        global_x_limit = 10

    for ax_idx, (panel_type, cluster_id, panel_label) in enumerate(panels):
        ax = axes[ax_idx]

        # Prepare data for this panel (includes both lags)
        panel_df = prepare_panel_data(df, panel_type, cluster_id)

        if panel_df.empty:
            ax.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=11,
                color="gray",
            )
            ax.set_title(panel_label, fontsize=12, fontweight="bold", loc="left")
            ax.axis("off")
            continue

        n_models = len(panel_df)
        # Create y_positions for all rows including empty row (creates visual gap)
        y_positions = np.arange(n_models)[::-1]  # Reverse so first model is at top

        # Set y-axis (model labels in the middle)
        ax.set_ylim(-0.5, n_models - 0.5)
        # Create labels dynamically by analyzing panel_df structure
        labels = []
        for _, row in panel_df.iterrows():
            if pd.isna(row["Lag"]):  # Empty row
                labels.append("")
            elif row["Domain"] == "EQI":
                labels.append("EQI")
            else:
                labels.append(row["Domain"])

        # Filter out ticks and labels for empty rows
        filtered_ticks = []
        filtered_labels = []
        for i, (pos, label) in enumerate(zip(y_positions, labels)):
            if not pd.isna(panel_df.iloc[i]["Lag"]):  # Only include non-empty rows
                filtered_ticks.append(pos)
                filtered_labels.append(label)

        ax.set_yticks(filtered_ticks)
        ax.set_yticklabels(filtered_labels, fontsize=12)

        # Use global x_limit for consistent axis range
        x_limit = global_x_limit
        ax.set_xlim(-x_limit, x_limit)

        # Set x-axis ticks to multiples of 5, or 1-5 if small range
        def get_nice_ticks(limit):
            """Generate tick positions."""
            if limit <= 5:
                # Show 1,2,3,4,5
                return [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
            else:
                # Show multiples of 5
                candidates = []
                for i in range(0, int(limit) + 1, 5):
                    if i > 0:
                        candidates.extend([-i, i])
                return sorted(list(set(candidates)))

        x_ticks = get_nice_ticks(x_limit)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{int(x)}" for x in x_ticks], fontsize=12)

        ax.axvline(0, color="black", linestyle="-", linewidth=1.5, zorder=1)
        # Find the position of the empty row for the dashed line
        empty_row_idx = None
        for idx, row in panel_df.iterrows():
            if pd.isna(row["Lag"]):
                empty_row_idx = idx
                break

        # Add dashed line between 5y and 10y (before the empty row)
        if empty_row_idx is not None:
            # Draw dashed line covering 95% of the width
            ax.plot(
                [x_limit * -0.95, x_limit * 0.95],
                [y_positions[empty_row_idx] + 0, y_positions[empty_row_idx] + 0],
                color="gray",
                linestyle="--",
                linewidth=1.0,
                zorder=1,
            )
            # Add colored labels with white background to cover dashed line
            ax.text(
                x_limit * -0.99,
                y_positions[empty_row_idx] + 0,
                "         << Improved Environment",
                fontsize=12,
                fontweight="bold",
                color=color_improved,
                va="center",
                ha="left",
                zorder=10,
                bbox=dict(
                    facecolor="white", edgecolor="none", boxstyle="round,pad=0.2"
                ),
            )
            ax.text(
                x_limit * 0.99,
                y_positions[empty_row_idx] + 0,
                "Worsen Environment >>         ",
                fontsize=12,
                fontweight="bold",
                color=color_worsened,
                va="center",
                ha="right",
                zorder=10,
                bbox=dict(
                    facecolor="white", edgecolor="none", boxstyle="round,pad=0.2"
                ),
            )
        ax.grid(True, axis="x", alpha=0.3, linestyle=":", linewidth=0.5)

        # Plot Improved (left side) and Worsened (right side)
        for idx, row in panel_df.iterrows():
            y_pos = y_positions[idx]

            # Skip plotting for empty row (where Lag is NaN)
            if pd.isna(row["Lag"]):
                continue

            # Plot Improved (left side) - BLUE because environment improved
            if not pd.isna(row["MRD_Improved"]):
                mrd = row["MRD_Improved"]
                lower = row["CI_Lower_Imp"]
                upper = row["CI_Upper_Imp"]
                is_sig = row["Sig_Imp"]

                # Color = BLUE (environment improved)
                color = color_improved

                # Alpha and size by significance
                alpha = 1.0 if is_sig else 0.4
                marker_size = 8 if is_sig else 6
                line_width = 2.0 if is_sig else 1.2

                # Error bar - circle marker
                ax.errorbar(
                    mrd,
                    y_pos,
                    xerr=[[mrd - lower], [upper - mrd]],
                    fmt="o",
                    color=color,
                    markersize=marker_size,
                    alpha=alpha,
                    capsize=3,
                    capthick=line_width,
                    linewidth=line_width,
                    zorder=3,
                )

            # Plot Worsened (right side) - ORANGE because environment worsened
            if not pd.isna(row["MRD_Worsened"]):
                mrd = row["MRD_Worsened"]
                lower = row["CI_Lower_Wor"]
                upper = row["CI_Upper_Wor"]
                is_sig = row["Sig_Wor"]

                # Color = ORANGE (environment worsened)
                color = color_worsened

                # Alpha and size by significance
                alpha = 1.0 if is_sig else 0.4
                marker_size = 8 if is_sig else 6
                line_width = 2.0 if is_sig else 1.2

                # Error bar - square marker
                ax.errorbar(
                    mrd,
                    y_pos,
                    xerr=[[mrd - lower], [upper - mrd]],
                    fmt="s",
                    color=color,
                    markersize=marker_size,
                    alpha=alpha,
                    capsize=3,
                    capthick=line_width,
                    linewidth=line_width,
                    zorder=3,
                )

        # Labels and title
        if ax_idx == len(panels) - 1:  # Bottom panel
            ax.set_xlabel(
                "Mortality Rate Difference (MRD) and 95% Confidence Interval",
                fontsize=12,
                fontweight="bold",
            )

        # Add vertical panel label on the left edge (rotated 90 degrees)
        ax.text(
            -0.07,
            0.5,
            panel_label,
            fontsize=12,
            fontweight="bold",
            rotation=90,
            verticalalignment="center",
            horizontalalignment="center",
            transform=ax.transAxes,
        )

        # Find the ranges for 5y and 10y data to position lag labels
        empty_row_idx = None
        for idx, row in panel_df.iterrows():
            if pd.isna(row["Lag"]):
                empty_row_idx = idx
                break

        if empty_row_idx is not None:
            # 5y data: from 0 to empty_row_idx-1
            if empty_row_idx > 0:
                y_5y_center = (y_positions[0] + y_positions[empty_row_idx - 1]) / 2
            else:
                y_5y_center = y_positions[0]
            # 10y data: from empty_row_idx+1 to end
            if empty_row_idx + 1 < n_models:
                y_10y_center = (y_positions[empty_row_idx + 1] + y_positions[-1]) / 2
            else:
                y_10y_center = y_positions[-1]

            # Add lag period labels on the left side, close to the axis
            ax.text(
                -x_limit + 0.5,
                y_5y_center,
                "Five-Year Lag",
                fontsize=12,
                fontweight="bold",
                rotation=90,
                verticalalignment="center",
                horizontalalignment="center",
                color="black",
            )

            ax.text(
                -x_limit + 0.5,
                y_10y_center,
                "Ten-Year Lag",
                fontsize=12,
                fontweight="bold",
                rotation=90,
                verticalalignment="center",
                horizontalalignment="center",
                color="black",
            )

        # Spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5)

    # Add legend with clearer design
    from matplotlib.lines import Line2D

    # Section 1: Environmental Quality Change Direction (color and marker with error bars)
    legend_quality = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=color_improved,
            markerfacecolor=color_improved,
            markersize=10,
            label="Improved Environment",
            markeredgewidth=1.5,
            markeredgecolor="black",
            alpha=1.0,
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color=color_worsened,
            markerfacecolor=color_worsened,
            markersize=10,
            label="Worsen Environment",
            markeredgewidth=1.5,
            markeredgecolor="black",
            alpha=1.0,
            linewidth=2,
        ),
    ]

    # Section 2: Statistical Significance (filled vs hollow markers with error bars)
    legend_significance = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="#666666",
            markerfacecolor="#666666",
            markersize=12,
            label="Significant",
            markeredgewidth=1.5,
            markeredgecolor="black",
            alpha=1.0,
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="#666666",
            markerfacecolor="#666666",
            markersize=12,
            label="Significant",
            markeredgewidth=1.5,
            markeredgecolor="black",
            alpha=1.0,
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="#666666",
            markerfacecolor="none",
            markersize=12,
            label="Non-significant",
            markeredgewidth=2.0,
            markeredgecolor="#666666",
            alpha=1.0,
            linewidth=2,
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="#666666",
            markerfacecolor="none",
            markersize=12,
            label="Non-significant",
            markeredgewidth=2.0,
            markeredgecolor="#666666",
            alpha=1.0,
            linewidth=2,
        ),
    ]

    # Combine legend elements
    all_legend = legend_quality + legend_significance

    fig.legend(
        handles=all_legend,
        loc="center right",
        bbox_to_anchor=(0.98, 0.52),
        frameon=True,
        fancybox=True,
        shadow=False,
        ncol=1,
        fontsize=12,
        title="Legend",
        title_fontsize=12,
        framealpha=1.0,
        edgecolor="gray",
        facecolor="white",
        borderpad=1.5,
        labelspacing=0.8,
        handletextpad=0.8,
    )

    # Save figure
    output_file = output_dir / f"{icd_code}_Delta_Cluster_k{k}.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"✓ Generated: {output_file.name}")
    return output_file


def list_available_icds(result_dir: Path) -> List[str]:
    """List available ICD codes from delta cluster result files."""
    icds = []
    for p in sorted(result_dir.glob("*_national.csv")):
        icd = p.stem.replace("_national", "")
        if icd not in icds:
            icds.append(icd)
    return icds


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]

    parser = argparse.ArgumentParser(
        description="Generate bidirectional forest plots for Delta Cluster Bayesian analysis (single-domain models with 5y and 10y lags)"
    )
    parser.add_argument("--icd", type=str, help="Specific ICD code (e.g., C00_C97)")
    parser.add_argument(
        "--all", action="store_true", help="Generate for all available ICDs"
    )

    args = parser.parse_args(argv)

    if not args.icd and not args.all:
        print("Error: Please specify --icd <CODE> or --all")
        return 1

    paths = get_paths(project_root)

    # Determine ICD list
    if args.all:
        icds = list_available_icds(paths["result"])
        if not icds:
            print(f"No delta cluster result files found in {paths['result']}")
            return 1
        print(f"Found {len(icds)} cancer types: {', '.join(icds)}")
    else:
        icds = [args.icd]

    # Generate plots
    generated = 0
    for icd in icds:
        for k in [3, 4]:
            # Collect files: national and k-specific
            files = {
                "national": paths["result"] / f"{icd}_national.csv",
                f"k{k}": paths["result"] / f"{icd}_k{k}.csv",
            }

            # Check if files exist
            missing = [key for key, f in files.items() if not f.exists()]
            if missing:
                print(
                    f"⚠ Missing files for {icd} k={k}: {', '.join(missing)}, skipping"
                )
                continue

            # Load and combine data
            dfs = []
            for key, file_path in files.items():
                try:
                    df = pd.read_csv(file_path, dtype={"Lag": int})
                    dfs.append(df)
                except Exception as e:
                    print(f"⚠ Error reading {file_path.name}: {e}")
                    continue

            if not dfs:
                print(f"⚠ No data loaded for {icd} k={k}, skipping")
                continue

            combined_df = pd.concat(dfs, ignore_index=True)

            print(f"\n📊 Processing {icd} k={k}...")

            try:
                plot_bidirectional_forest(combined_df, icd, k, paths["vis"])
                generated += 1
            except Exception as e:
                print(f"  ✗ Error generating plot for {icd} k={k}: {e}")
                import traceback

                traceback.print_exc()

    print(f"\n{'=' * 60}")
    print(f"✓ Successfully generated {generated} plots")
    print(f"📁 Output directory: {paths['vis']}")
    print(f"{'=' * 60}")

    return 0 if generated > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
