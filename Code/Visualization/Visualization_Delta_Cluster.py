import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# Define font families to generate multiple versions
FONT_FAMILIES = [
    {
        "name": "Georgia",
        "family": "serif",
        "fonts": ["Georgia", "DejaVu Serif", "Times New Roman"],
    },
    {
        "name": "Helvetica",
        "family": "sans-serif",
        "fonts": ["Helvetica", "Arial", "DejaVu Sans"],
    },
]


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


def prepare_panel_data_single_lag(
    df: pd.DataFrame, panel_type: str, lag: int, cluster_id: Optional[int] = None
) -> pd.DataFrame:
    """
    Extract data for one panel (National or specific Cluster) for a SINGLE lag period.

    Parameters:
    -----------
    df : DataFrame with all models for one cancer/lag/cluster
    panel_type : 'National' or 'Cluster'
    lag : 5 or 10 (year lag period)
    cluster_id : Cluster ID (0,1,2,... ) if panel_type='Cluster'

    Returns:
    --------
    DataFrame with columns: [Domain, Model, Lag, MRD_Improved, CI_Lower_Imp, CI_Upper_Imp,
                             Sig_Imp, MRD_Worsened, CI_Lower_Wor, CI_Upper_Wor, Sig_Wor]

    Order: EQI + domains (Air, Water, Land, Built, Social)
    """
    # Filter df based on Model column
    if panel_type == "National":
        df = df[~df["Model"].str.startswith("Cluster_")]
    else:
        df = df[df["Model"].str.startswith(f"Cluster{cluster_id}_")]

    rows = []

    if panel_type == "National":
        base_models = ["EQI", "Air", "Water", "Land", "Built", "Social"]
        prefix = ""
    else:
        base_models = [
            f"Cluster{cluster_id}_EQI",
            f"Cluster{cluster_id}_Air",
            f"Cluster{cluster_id}_Water",
            f"Cluster{cluster_id}_Land",
            f"Cluster{cluster_id}_Built",
            f"Cluster{cluster_id}_Social",
        ]
        prefix = f"Cluster{cluster_id}_"

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

        rows.append(row_data)

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def plot_bidirectional_forest(
    df: pd.DataFrame,
    icd_code: str,
    k: int,
    output_dir: Path,
    font_config: dict = None,
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
    if font_config is None:
        font_config = FONT_FAMILIES[0]
    plt.rcParams["font.family"] = font_config["family"]
    plt.rcParams[f"font.{font_config['family']}"] = font_config["fonts"]
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
    color_improved = "#2EAA6F"  # Green for quality improved
    color_worsened = "#FF6F48"  # Orange for quality worsened

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

                # Consistent line width and marker size, alpha and marker style by significance
                alpha = 1.0 if is_sig else 0.4
                marker_size = 7  # Consistent marker size
                line_width = 2.0  # Consistent line width
                # Hollow marker for non-significant
                markerfacecolor = color if is_sig else "none"
                markeredgewidth = 1.5 if is_sig else 2.0

                # Error bar - circle marker
                ax.errorbar(
                    mrd,
                    y_pos,
                    xerr=[[mrd - lower], [upper - mrd]],
                    fmt="o",
                    color=color,
                    markerfacecolor=markerfacecolor,
                    markeredgecolor=color,
                    markeredgewidth=markeredgewidth,
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

                # Consistent line width and marker size, alpha and marker style by significance
                alpha = 1.0 if is_sig else 0.4
                marker_size = 7  # Consistent marker size
                line_width = 2.0  # Consistent line width
                # Hollow marker for non-significant
                markerfacecolor = color if is_sig else "none"
                markeredgewidth = 1.5 if is_sig else 2.0

                # Error bar - square marker
                ax.errorbar(
                    mrd,
                    y_pos,
                    xerr=[[mrd - lower], [upper - mrd]],
                    fmt="s",
                    color=color,
                    markerfacecolor=markerfacecolor,
                    markeredgecolor=color,
                    markeredgewidth=markeredgewidth,
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
                "MRD and 95% CrI",
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

    # Add legend at bottom with 4 items in one row
    from matplotlib.lines import Line2D

    legend_elements = [
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
        Line2D(
            [0],
            [0],
            marker="o",
            color="#666666",
            markerfacecolor="#666666",
            markersize=10,
            label="CrI Excluding 0",
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
            markersize=10,
            label="CrI Including 0",
            markeredgewidth=2.0,
            markeredgecolor="#666666",
            alpha=0.4,
            linewidth=2,
        ),
    ]

    # Place legend at the bottom, 4 items in one row, no title
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.0),
        frameon=True,
        fancybox=False,
        shadow=False,
        ncol=4,
        fontsize=12,
        framealpha=1.0,
        edgecolor="gray",
        facecolor="white",
        borderpad=1.0,
        labelspacing=0.5,
        handletextpad=0.5,
        columnspacing=1.5,
    )

    # Save figure with suffix for font name
    font_suffix = f"_{font_config['name']}" if font_config else ""
    output_file = output_dir / f"{icd_code}_Delta_Cluster_k{k}{font_suffix}.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"✓ Generated: {output_file.name}")
    return output_file


def plot_bidirectional_forest_separate_years(
    df: pd.DataFrame,
    icd_code: str,
    k: int,
    output_dir: Path,
    font_config: dict = None,
):
    """
    Create bidirectional forest plot with SEPARATE columns for 5-year and 10-year lags.
    Layout: (k+1) rows × 2 columns
    - Left column: 5-year lag
    - Right column: 10-year lag
    """
    # Set global font
    if font_config is None:
        font_config = FONT_FAMILIES[0]
    plt.rcParams["font.family"] = font_config["family"]
    plt.rcParams[f"font.{font_config['family']}"] = font_config["fonts"]
    plt.rcParams["font.size"] = 12

    # Define panel configurations: National + clusters 0 to k-1
    panels = [("National", None, "National")] + [
        ("Cluster", i, f"Cluster {i}") for i in range(k)
    ]

    # Color scheme
    color_improved = "#2EAA6F"  # Green
    color_worsened = "#FF6F48"  # Orange

    # Calculate global x_limit based on all panels and both lags
    global_all_values = []
    for panel_type, cluster_id, _ in panels:
        for lag in [5, 10]:
            panel_df = prepare_panel_data_single_lag(df, panel_type, lag, cluster_id)
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

    # Create figure with (k+1) rows × 2 columns
    fig, axes = plt.subplots(
        k + 1, 2, figsize=(14.0, 2.5 * (k + 1)), constrained_layout=False
    )
    # Adjust spacing between subplots
    # To manually adjust horizontal spacing between 5y and 10y plots, change wspace value:
    # wspace=0.05 means 5% of subplot width. Decrease for tighter spacing, increase for wider.
    plt.subplots_adjust(wspace=0.15, hspace=0.3)

    # Y-axis label horizontal offset for 10-year lag plot (negative=left, positive=right)
    # Adjust this value to move y-axis labels left or right
    y_label_pad = -0.055  # Default: -0.02 (slightly to the left)

    # Fixed x-axis limit
    x_limit = 8

    # Helper function for nice ticks
    def get_nice_ticks(limit):
        # For fixed limit of 8, show only -6 to 6 with step 2 (hide -8 and 8)
        return [-6, -4, -2, 0, 2, 4, 6]

    # Plot each panel
    for row_idx, (panel_type, cluster_id, panel_label) in enumerate(panels):
        for col_idx, lag in enumerate([5, 10]):
            ax = axes[row_idx, col_idx]

            # Prepare data for this panel and lag
            panel_df = prepare_panel_data_single_lag(df, panel_type, lag, cluster_id)

            if panel_df.empty:
                ax.text(
                    0.5,
                    0.5,
                    "No data available",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=12,
                    color="gray",
                )
                ax.set_title(
                    f"{panel_label} ({lag}y)",
                    fontsize=12,
                    loc="left",
                )
                ax.axis("off")
                continue

            n_models = len(panel_df)
            y_positions = np.arange(n_models)[::-1]

            # Set y-axis
            ax.set_ylim(-0.5, n_models - 0.5)
            labels = panel_df["Domain"].tolist()

            # For 5-year lag (left column), move y-axis to right and hide labels
            # For 10-year lag (right column), show on left as normal
            if col_idx == 0:  # 5-year lag
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
                ax.set_yticks(y_positions)
                ax.set_yticklabels([""] * n_models)  # Hide labels
            else:  # 10-year lag
                ax.set_yticks(y_positions)
                # Use tick_params to adjust label position (pad parameter)
                # Positive pad moves labels away from axis, negative moves closer
                ax.set_yticklabels(labels, fontsize=12, ha="center")
                # Apply horizontal offset to y-axis labels
                for tick in ax.get_yticklabels():
                    tick.set_horizontalalignment("center")
                    # Get current position and shift horizontally
                    pos = tick.get_position()
                    tick.set_position((pos[0] + y_label_pad, pos[1]))

            # Set x-axis with fixed limit
            ax.set_xlim(-x_limit, x_limit)
            x_ticks = get_nice_ticks(x_limit)
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f"{int(x)}" for x in x_ticks], fontsize=12)

            # Zero line
            ax.axvline(0, color="black", linestyle="-", linewidth=1.0, zorder=1)

            # Add direction labels at top (moved down)
            if row_idx == 0:
                ax.text(
                    -x_limit * 0.75,
                    n_models - 0.3,
                    "<< Improved",
                    fontsize=15,
                    fontweight="bold",
                    color=color_improved,
                    ha="center",
                )
                ax.text(
                    x_limit * 0.75,
                    n_models - 0.3,
                    "Worsened >>",
                    fontsize=15,
                    fontweight="bold",
                    color=color_worsened,
                    ha="center",
                )

            ax.grid(True, axis="x", alpha=0.3, linestyle=":", linewidth=0.5)

            # Plot data points
            for idx, row in panel_df.iterrows():
                y_pos = y_positions[idx]

                # Plot Improved (left side) - BLUE
                if not pd.isna(row["MRD_Improved"]):
                    mrd = row["MRD_Improved"]
                    lower = row["CI_Lower_Imp"]
                    upper = row["CI_Upper_Imp"]
                    is_sig = row["Sig_Imp"]

                    # Consistent line width and marker size, alpha and marker style by significance
                    alpha = 1.0 if is_sig else 0.4
                    marker_size = 7  # Consistent marker size
                    line_width = 2.0  # Consistent line width
                    # Hollow marker for non-significant
                    markerfacecolor = color_improved if is_sig else "none"
                    markeredgewidth = 1.5 if is_sig else 2.0

                    # Clip to axis limits
                    mrd_clipped = np.clip(mrd, -x_limit, x_limit)
                    lower_clipped = np.clip(lower, -x_limit, x_limit)
                    upper_clipped = np.clip(upper, -x_limit, x_limit)

                    # Check if out of range
                    out_of_range = (
                        (mrd < -x_limit) or (upper < -x_limit) or (lower > x_limit)
                    )

                    if out_of_range:
                        # Draw arrow pointing to out-of-range direction
                        if mrd < -x_limit or upper < -x_limit:
                            # Point exceeds left boundary
                            ax.plot(
                                -x_limit,
                                y_pos,
                                marker="<",
                                color=color_improved,
                                markerfacecolor=markerfacecolor,
                                markeredgecolor=color_improved,
                                markeredgewidth=markeredgewidth,
                                markersize=10,
                                alpha=alpha,
                                zorder=4,
                            )
                        elif lower > x_limit:
                            # Point exceeds right boundary (unlikely for improved but handle it)
                            ax.plot(
                                x_limit,
                                y_pos,
                                marker=">",
                                color=color_improved,
                                markerfacecolor=markerfacecolor,
                                markeredgecolor=color_improved,
                                markeredgewidth=markeredgewidth,
                                markersize=10,
                                alpha=alpha,
                                zorder=4,
                            )
                    else:
                        # Normal error bar
                        ax.errorbar(
                            mrd_clipped,
                            y_pos,
                            xerr=[
                                [mrd_clipped - lower_clipped],
                                [upper_clipped - mrd_clipped],
                            ],
                            fmt="o",
                            color=color_improved,
                            markerfacecolor=markerfacecolor,
                            markeredgecolor=color_improved,
                            markeredgewidth=markeredgewidth,
                            markersize=marker_size,
                            alpha=1.0 if is_sig else 0.3,
                            capsize=3,
                            capthick=line_width,
                            linewidth=line_width,
                            zorder=3,
                        )

                # Plot Worsened (right side) - ORANGE
                if not pd.isna(row["MRD_Worsened"]):
                    mrd = row["MRD_Worsened"]
                    lower = row["CI_Lower_Wor"]
                    upper = row["CI_Upper_Wor"]
                    is_sig = row["Sig_Wor"]

                    # Consistent line width and marker size, alpha and marker style by significance
                    alpha = 1.0 if is_sig else 0.4
                    marker_size = 7  # Consistent marker size
                    line_width = 2.0  # Consistent line width
                    # Hollow marker for non-significant
                    markerfacecolor = color_worsened if is_sig else "none"
                    markeredgewidth = 1.5 if is_sig else 2.0

                    # Clip to axis limits
                    mrd_clipped = np.clip(mrd, -x_limit, x_limit)
                    lower_clipped = np.clip(lower, -x_limit, x_limit)
                    upper_clipped = np.clip(upper, -x_limit, x_limit)

                    # Check if out of range
                    out_of_range = (
                        (mrd > x_limit) or (lower > x_limit) or (upper < -x_limit)
                    )

                    if out_of_range:
                        # Draw arrow pointing to out-of-range direction
                        if mrd > x_limit or lower > x_limit:
                            # Point exceeds right boundary
                            ax.plot(
                                x_limit,
                                y_pos,
                                marker=">",
                                color=color_worsened,
                                markerfacecolor=markerfacecolor,
                                markeredgecolor=color_worsened,
                                markeredgewidth=markeredgewidth,
                                markersize=10,
                                alpha=alpha,
                                zorder=4,
                            )
                        elif upper < -x_limit:
                            # Point exceeds left boundary (unlikely for worsened but handle it)
                            ax.plot(
                                -x_limit,
                                y_pos,
                                marker="<",
                                color=color_worsened,
                                markerfacecolor=markerfacecolor,
                                markeredgecolor=color_worsened,
                                markeredgewidth=markeredgewidth,
                                markersize=10,
                                alpha=alpha,
                                zorder=4,
                            )
                    else:
                        # Normal error bar
                        ax.errorbar(
                            mrd_clipped,
                            y_pos,
                            xerr=[
                                [mrd_clipped - lower_clipped],
                                [upper_clipped - mrd_clipped],
                            ],
                            fmt="s",
                            color=color_worsened,
                            markerfacecolor=markerfacecolor,
                            markeredgecolor=color_worsened,
                            markeredgewidth=markeredgewidth,
                            markersize=marker_size,
                            alpha=1.0 if is_sig else 0.3,
                            capsize=3,
                            capthick=line_width,
                            linewidth=line_width,
                            zorder=3,
                        )

            # X-axis label (only on bottom row)
            if row_idx == len(panels) - 1:
                ax.set_xlabel(
                    "MRD and 95% CrI",
                    fontsize=15,
                )

            # Column title (only on top row)
            if row_idx == 0:
                lag_label = "(A) Five-Year Lag" if lag == 5 else "(B) Ten-Year Lag"
                ax.set_title(
                    lag_label,
                    fontsize=15,
                    pad=5,
                )

            # Row label (only on left column)
            if col_idx == 0:
                ax.text(
                    -0.02,
                    0.5,
                    panel_label,
                    fontsize=15,
                    rotation=90,
                    verticalalignment="center",
                    horizontalalignment="center",
                    transform=ax.transAxes,
                )

            # Spines - show all four borders
            ax.spines["top"].set_visible(True)
            ax.spines["right"].set_visible(True)
            ax.spines["left"].set_visible(True)
            ax.spines["bottom"].set_visible(True)
            ax.spines["top"].set_linewidth(0.5)
            ax.spines["right"].set_linewidth(0.5)
            ax.spines["left"].set_linewidth(0.5)
            ax.spines["bottom"].set_linewidth(0.5)

    # Add legend at bottom with 4 items in one row
    from matplotlib.lines import Line2D

    legend_elements = [
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
        Line2D(
            [0],
            [0],
            marker="o",
            color="#666666",
            markerfacecolor="#666666",
            markersize=10,
            label="CrI Excluding 0",
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
            markersize=10,
            label="CrI Including 0",
            markeredgewidth=2.0,
            markeredgecolor="#666666",
            alpha=0.4,
            linewidth=2,
        ),
    ]

    # Place legend at the bottom, 4 items in one row, no title
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.0),
        frameon=True,
        fancybox=False,
        shadow=False,
        ncol=4,
        fontsize=12,
        framealpha=1.0,
        edgecolor="gray",
        facecolor="white",
        borderpad=1.0,
        labelspacing=0.5,
        handletextpad=0.5,
        columnspacing=1.5,
    )

    # Save figure with suffix for font name
    font_suffix = f"_{font_config['name']}" if font_config else ""
    output_file = (
        output_dir / f"{icd_code}_Delta_Cluster_k{k}_SeparateYear{font_suffix}.png"
    )
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

            # Generate plots for each font family
            for font_config in FONT_FAMILIES:
                print(f"  Font: {font_config['name']}")
                try:
                    plot_bidirectional_forest(
                        combined_df, icd, k, paths["vis"], font_config
                    )
                    generated += 1
                except Exception as e:
                    print(
                        f"  ✗ Error generating plot for {icd} k={k} with {font_config['name']}: {e}"
                    )
                    import traceback

                    traceback.print_exc()

                # Also generate separate-year version
                try:
                    plot_bidirectional_forest_separate_years(
                        combined_df, icd, k, paths["vis"], font_config
                    )
                    generated += 1
                except Exception as e:
                    print(
                        f"  ✗ Error generating separate-year plot for {icd} k={k} with {font_config['name']}: {e}"
                    )
                    import traceback

                    traceback.print_exc()

    print(f"\n{'=' * 60}")
    print(f"✓ Successfully generated {generated} plots")
    print(f"📁 Output directory: {paths['vis']}")
    print(f"{'=' * 60}")

    return 0 if generated > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
