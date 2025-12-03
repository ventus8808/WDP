"""
EQI LMM Visualization - County Economic Typology Stratification

- Reads LMM results from Result/brms_Typology_LandUse (paths via config.yaml)
- Builds forest-plot style figures for County Economic Typology stratifications
- Output: {ICD}_{Typology}_{EQI_Period}_{AAMR_Period}_Lag{lag}.png

CLI Usage:
python Code/Visualization/Visualization_Typology.py --icd C00_C97
python Code/Visualization/Visualization_Typology.py --all

Typology Strata:
1. Farming (Farming-dependent counties)
2. Mining (Mining-dependent counties)
3. Manufacturing (Manufacturing-dependent counties)
4. Government (Federal/State government-dependent counties)
5. Services (Services-dependent counties)
6. Nonspecialized (Nonspecialized counties)
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
import pandas as pd
import yaml

matplotlib.use("Agg")  # Use non-interactive backend
import matplotlib as mpl
import matplotlib.pyplot as plt

# Simple scenario order (fixed five scenarios)
SCENARIO_ORDER: List[Tuple[str, str, int]] = [
    ("2000_2005", "2006_2010", 5),
    ("2000_2005", "2011_2015", 10),
    ("2000_2005", "2016_2020", 15),
    ("2006_2010", "2011_2015", 5),
    ("2006_2010", "2016_2020", 10),
]

# Scenario letter mapping for titles
SCENARIO_LABELS = {
    ("2000_2005", "2006_2010", 5): "A",
    ("2000_2005", "2011_2015", 10): "B",
    ("2000_2005", "2016_2020", 15): "C",
    ("2006_2010", "2011_2015", 5): "D",
    ("2006_2010", "2016_2020", 10): "E",
}

# Set global font to Georgia with sensible fallbacks
mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Garamond", "Times New Roman"],
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.titlesize": 15,
    }
)


def load_config(project_root: Path) -> dict:
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_paths(project_root: Path, cfg: dict) -> Dict[str, Path]:
    """Get paths for Typology results and visualization output"""
    result_dir = project_root / "Result" / "brms_Typology_LandUse"
    vis_dir = project_root / "Result" / "brms_Typology_Visualization"

    vis_dir.mkdir(parents=True, exist_ok=True)
    return {"result": result_dir, "vis": vis_dir}


def list_icds_from_results(result_dir: Path) -> List[str]:
    """List ICD codes from Typology result files"""
    icds = []
    for p in sorted(result_dir.glob("*_Typology.csv")):
        stem = p.stem.replace("_Typology", "")
        if stem not in icds:
            icds.append(stem)
    return icds


EFFECT_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?)\s*\(\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)\s*\)\s*$"
)


def parse_effect_cell(cell: str) -> Optional[Tuple[float, float, float, bool, str]]:
    """Parse a cell like "-1.54(-4.20, 1.11)***" -> (mrd, lcl, ucl, is_sig, stars)"""
    if cell is None:
        return None
    s = str(cell).strip()
    if not s or s == "0.00":
        return None
    m_stars = re.search(r"(\*+)$", s)
    stars = m_stars.group(1) if m_stars else ""
    core = s[: -len(stars)] if stars else s
    core = core.strip()
    m = EFFECT_RE.match(core)
    if not m:
        return None
    mrd = float(m.group(1))
    lcl = float(m.group(2))
    ucl = float(m.group(3))
    return mrd, lcl, ucl, bool(stars), stars


def plot_typology_forest(
    df,
    eqi_period,
    aamr_period,
    lag,
    output_dir=None,
    icd_code="C00_C97",
    cfg=None,
):
    """Create forest plot for County Economic Typology stratification"""

    # Define panel configurations
    panel_labels = [
        "Farming",
        "Mining",
        "Manufacturing",
        "Government",
        "Services",
        "Nonspecialized",
    ]

    # Define EQI types
    eqi_types = ["EQI", "Air", "Water", "Land", "Built", "Social"]

    # Define color mapping (grayscale)
    color_map = {
        "Q2": "#B3B3B3",  # Light gray (70%)
        "Q3": "#808080",  # Medium gray (50%)
        "Q4": "#4D4D4D",  # Dark gray (30%)
        "Q5": "#000000",  # Black (0%)
    }

    # Create figure
    fig, axes = plt.subplots(6, 1, figsize=(6, 14))

    # Get scenario label and cancer name
    scenario_key = (eqi_period, aamr_period, lag)
    scenario_label = SCENARIO_LABELS.get(scenario_key, "")

    # Get cancer name from config
    cancer_name = icd_code
    if cfg and "brms_analysis" in cfg and "icd_mapping" in cfg["brms_analysis"]:
        cancer_name = cfg["brms_analysis"]["icd_mapping"].get(icd_code, icd_code)

    # Format title
    title = (
        f"({scenario_label}) {cancer_name} {lag}-year Lag - County Economic Typology"
    )
    fig.suptitle(title, fontsize=14, y=0.95)

    # For each panel (typology stratum)
    for panel_idx, panel_label in enumerate(panel_labels):
        ax = axes[panel_idx]

        # Filter data for current panel
        typology_models = [
            f"Typology_{panel_label}_EQI",
            f"Typology_{panel_label}_EQI_Air",
            f"Typology_{panel_label}_EQI_Water",
            f"Typology_{panel_label}_EQI_Land",
            f"Typology_{panel_label}_EQI_Built",
            f"Typology_{panel_label}_EQI_Social",
        ]

        panel_df = df[df["Model"].isin(typology_models)].copy()

        # Collect all values for y-axis scaling
        panel_values = []

        # Parse data and collect values
        for model_idx, model in enumerate(typology_models):
            model_data = panel_df[panel_df["Model"] == model]
            if model_data.empty:
                continue

            for q_idx, quintile in enumerate(["Q2", "Q3", "Q4", "Q5"]):
                value_str = (
                    model_data[quintile].iloc[0] if not model_data.empty else None
                )
                if (
                    value_str is None
                    or pd.isna(value_str)
                    or str(value_str).strip() == "0.00"
                ):
                    continue

                # Parse value
                pattern = re.compile(r"(-?\d+\.\d+)\s*\((.+?),\s*(.+?)\)([\*†]*)")
                match = pattern.match(str(value_str))
                if not match:
                    continue

                try:
                    mrd = float(match.group(1))
                    lower = float(match.group(2))
                    upper = float(match.group(3))
                    panel_values.extend([mrd, lower, upper])
                except ValueError:
                    continue

        # Parse data and plot
        for model_idx, model in enumerate(typology_models):
            model_data = panel_df[panel_df["Model"] == model]
            if model_data.empty:
                continue

            x_base = model_idx * 0.6

            # Plot each quintile
            for q_idx, quintile in enumerate(["Q2", "Q3", "Q4", "Q5"]):
                value_str = (
                    model_data[quintile].iloc[0] if not model_data.empty else None
                )
                if (
                    value_str is None
                    or pd.isna(value_str)
                    or str(value_str).strip() == "0.00"
                ):
                    continue

                # Parse value
                pattern = re.compile(r"(-?\d+\.\d+)\s*\((.+?),\s*(.+?)\)([\*†]*)")
                match = pattern.match(str(value_str))
                if not match:
                    continue

                try:
                    mrd = float(match.group(1))
                    lower = float(match.group(2))
                    upper = float(match.group(3))
                except ValueError:
                    continue

                # Calculate x position
                x_pos = x_base + (q_idx - 1.5) * 0.065

                # Plot error bar and point
                color = color_map[quintile]
                ax.errorbar(
                    x_pos,
                    mrd,
                    yerr=[[mrd - lower], [upper - mrd]],
                    fmt="o",
                    color=color,
                    markersize=3,
                    capsize=2,
                    capthick=1.2,
                    linewidth=1.2,
                )

        # Set panel style
        ax.set_xlim(-0.3, 3.3)

        # Dynamic y-axis range
        if panel_values:
            max_abs_value = max(abs(v) for v in panel_values)
            dynamic_limit = ((int(max_abs_value) // 5) + 1) * 5
            ax.set_ylim(-dynamic_limit, dynamic_limit)
        else:
            ax.set_ylim(-20, 20)

        ax.axhline(0, color="gray", linestyle="-", linewidth=0.8)
        ax.grid(True, alpha=0.3)

        # Set x-axis labels
        ax.set_xticks([i * 0.6 for i in range(6)])
        ax.set_xticklabels(eqi_types)

        # Set y-axis label
        if panel_idx == 3:  # Middle panel
            ax.set_ylabel(
                "Mortality Rate Difference and 95% Credible Interval",
                fontsize=15,
                labelpad=10,
            )

        # Add black border
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.0)

        # Add panel label in top right corner
        ax.text(
            0.98,
            0.95,
            panel_label,
            fontsize=11,
            verticalalignment="top",
            horizontalalignment="right",
            transform=ax.transAxes,
        )

    plt.tight_layout(rect=(0, 0, 1, 0.96))

    # Save figure
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        output_filename = (
            output_path / f"{icd_code}_Typology_{eqi_period}_{aamr_period}_Lag{lag}.png"
        )
    else:
        output_filename = f"{icd_code}_Typology_{eqi_period}_{aamr_period}_Lag{lag}.png"

    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Typology forest plot saved as {output_filename}")
    plt.close()

    return str(output_filename)


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Generate County Economic Typology stratified forest plots"
    )
    parser.add_argument(
        "--icd", type=str, default=None, help="Specific ICD code, e.g., C00_C97"
    )
    parser.add_argument(
        "--all", action="store_true", help="Generate for all ICDs found in result dir"
    )

    args = parser.parse_args(argv)

    cfg = load_config(project_root)
    paths = get_paths(project_root, cfg)

    if not args.icd and not args.all:
        print("Please specify --icd or --all")
        return 1

    if args.all:
        icds = list_icds_from_results(paths["result"])
    else:
        icds = [args.icd]

    if not icds:
        print(f"No ICD codes found in {paths['result']}")
        return 1

    for icd in icds:
        result_csv = paths["result"] / f"{icd}_Typology.csv"

        if not result_csv.exists():
            print(f"Result CSV not found for {icd}: {result_csv}; skip.")
            continue

        df = pd.read_csv(result_csv, dtype=str)

        # Generate plots for fixed scenarios
        for eqi, aamr, lag in SCENARIO_ORDER:
            try:
                sub = df[
                    (df["EQI_Period"] == eqi)
                    & (df["AAMR_Period"] == aamr)
                    & (df["Lag"].astype(int) == lag)
                ]
            except Exception:
                # If Lag is stored as str
                sub = df[
                    (df["EQI_Period"] == eqi)
                    & (df["AAMR_Period"] == aamr)
                    & (df["Lag"] == str(lag))
                ]

            if sub.empty:
                print(f"[{icd}] no data for scenario {eqi}/{aamr}/Lag{lag}; skipping.")
                continue

            # Generate the plot
            try:
                plot_typology_forest(
                    df=sub,
                    eqi_period=eqi,
                    aamr_period=aamr,
                    lag=lag,
                    output_dir=str(paths["vis"]),
                    icd_code=icd,
                    cfg=cfg,
                )
                print(f"Generated plot: {icd}_Typology_{eqi}_{aamr}_Lag{lag}.png")
            except Exception as e:
                print(f"Error generating plot for {icd} {eqi}/{aamr}/Lag{lag}: {e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
