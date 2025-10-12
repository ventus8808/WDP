"""
EQI Stratified EQI Visualization

- Reads stratified BRMS results, extracts EQI-related models (EQI, RUCC*_EQI)
- Maps stratification groups: Female/Male -> Female/Male, races -> White/Black/Asian/Others
- Builds forest-plot style figures with stratification groups on X-axis instead of EQI types
- Output: C00_C97_EQI_Stratified_{EQI_Period}_{AAMR_Period}_Lag{lag}_brms.png

CLI Usage:
python Code/Analysis/Visualization_Stratified_EQI.py --all
"""
from __future__ import annotations

import re
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec
import matplotlib as mpl

# Simple scenario order (fixed four scenarios)
SCENARIO_ORDER: List[Tuple[str, str, int]] = [
    ("2000_2005", "2006_2010", 5),
    ("2000_2005", "2011_2015", 10),
    ("2006_2010", "2011_2015", 5),
    ("2006_2010", "2016_2020", 10),
]

# Domains to visualize
DOMAINS = ["EQI", "Air", "Water", "Land", "Built", "Social"]

# Stratification group mapping (add Overall)
STRATIFICATION_MAPPING = {
    "Sex_Female": "Female",
    "Sex_Male": "Male",
    "Race_American Indian or Alaska Native": "Others",
    "Race_Asian or Pacific Islander": "Asian",
    "Race_Black or African American": "Black",
    "Race_White": "White",
    "Overall": "Overall"  # For the overall data
}

# X-axis labels for the plot (stratification groups + Overall)
STRAT_GROUP_ORDER = ["Overall", "Male", "Female", "White", "Black", "Asian", "Others"]

# Set global font to Georgia with sensible fallbacks and slightly larger sizes
mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Georgia", "DejaVu Serif", "Times New Roman", "Palatino"],
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 15,
})

# ----------------------------
# Helpers: project paths
# ----------------------------

def load_config(project_root: Path) -> dict:
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_eqi_stratified_paths(project_root: Path) -> Dict[str, Path]:
    stratified_result_dir = project_root / "Result" / "brms_stratified"
    vis_dir = project_root / "Result" / "brms_stratified_EQI_Visualization"
    vis_dir.mkdir(parents=True, exist_ok=True)
    return {"result": stratified_result_dir, "vis": vis_dir}


# ----------------------------
# Parsing utilities (simplified)
# ----------------------------

EFFECT_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*\(\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)\s*\)$")


def parse_effect_cell(cell: str) -> Optional[Tuple[float, float, float, bool, str]]:
    """Parse a cell like "-1.54(-4.20, 1.11)***" -> (mrd, lcl, ucl, is_sig, stars)
    Simplified: rely on the main regex. If it fails, return None.
    """
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


QUINTILE_ORDER = ["Q2", "Q3", "Q4", "Q5"]

# ----------------------------
# Data collection for EQI stratified plot
# ----------------------------

def collect_eqi_data_for_scenario(paths: Dict[str, Path], domain: str, eqi_period: str, aamr_period: str, lag: int, project_root: Path) -> pd.DataFrame:
    """Collect domain data from all stratified CSVs and overall CSV for a specific scenario."""
    all_data = []
    
    # Stratified data
    for csv_path in paths["result"].glob("C00_C97_brms_*.csv"):
        stem = csv_path.stem
        parts = stem.split("_brms_")
        if len(parts) == 2:
            strat_key = parts[1]
            if strat_key in STRATIFICATION_MAPPING:
                df = pd.read_csv(csv_path, dtype=str)
                # Filter for domain models
                domain_df = df[df['Model'].str.contains(domain, na=False)].copy()
                # Filter for the specific scenario
                scenario_df = domain_df[
                    (domain_df["EQI_Period"] == eqi_period) &
                    (domain_df["AAMR_Period"] == aamr_period) &
                    (domain_df["Lag"].astype(int if domain_df["Lag"].str.isdigit().all() else str) == lag)
                ].copy()
                if not scenario_df.empty:
                    scenario_df['Strat_Group'] = STRATIFICATION_MAPPING[strat_key]
                    all_data.append(scenario_df)
    
    # Overall data
    overall_csv = project_root / "Result" / "brms" / "C00_C97_brms.csv"
    if overall_csv.exists():
        df_overall = pd.read_csv(overall_csv, dtype=str)
        domain_df_overall = df_overall[df_overall['Model'].str.contains(domain, na=False)].copy()
        scenario_df_overall = domain_df_overall[
            (domain_df_overall["EQI_Period"] == eqi_period) &
            (domain_df_overall["AAMR_Period"] == aamr_period) &
            (domain_df_overall["Lag"].astype(int if domain_df_overall["Lag"].str.isdigit().all() else str) == lag)
        ].copy()
        if not scenario_df_overall.empty:
            scenario_df_overall['Strat_Group'] = "Overall"
            all_data.append(scenario_df_overall)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df
    return pd.DataFrame()


# ----------------------------
# Main: generate EQI stratified plots
# ----------------------------

def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate domain stratified plots")
    parser.add_argument("--all", action="store_true", help="Generate for all domains and scenarios")

    args = parser.parse_args(argv)

    if not args.all:
        print("Please specify --all")
        return 1

    paths = get_eqi_stratified_paths(project_root)

    for domain in DOMAINS:
        for eqi, aamr, lag in SCENARIO_ORDER:
            df = collect_eqi_data_for_scenario(paths, domain, eqi, aamr, lag, project_root)
            if df.empty:
                print(f"No data for {domain} scenario {eqi}/{aamr}/Lag{lag}; skipping.")
                continue

            out_name = f"C00_C97_{domain}_Stratified_{eqi}_{aamr}_Lag{lag}_brms.png"
            out_path = paths["vis"] / out_name
            title = f"C00_C97 {domain} Stratified {eqi} {aamr} Lag{lag}"

            try:
                plot_eqi_stratified_forest(df, domain, eqi, aamr, lag, str(paths["vis"]), out_name)
                print(f"Generated plot: {out_path}")
            except Exception as e:
                print(f"Error generating plot for {title}: {e}")

    return 0


def plot_eqi_stratified_forest(df: pd.DataFrame, domain: str, eqi_period: str, aamr_period: str, lag: int, output_dir: str, filename: str):
    """
    Create forest plot for domain stratified by groups.
    
    X-axis: Stratification groups + Overall
    Panels: 5 (mostratified + RUCC1-4)
    """
    # Define panel configurations
    panel_labels = [
        'mostratified',
        'metropolitan urban', 
        'nonmetropolitan urban',
        'less urban',
        'thinly populated'
    ]
    
    # Color map for quintiles
    color_map = {
        'Q2': '#B3B3B3',
        'Q3': '#808080',
        'Q4': '#4D4D4D',
        'Q5': '#000000',
    }
    
    # Create figure
    fig, axes = plt.subplots(5, 1, figsize=(14, 12))  # Wider for 7 groups
    fig.suptitle(f'C00_C97 | {domain} Stratified | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}', 
                 fontsize=14, y=0.98)
    
    # For each panel
    for panel_idx, panel_label in enumerate(panel_labels):
        ax = axes[panel_idx]
        
        # Select models for this panel
        if panel_idx == 0:
            if domain == 'EQI':
                panel_models = ['EQI']
            else:
                panel_models = [f'EQI_{domain}']
        else:
            if domain == 'EQI':
                panel_models = [f'RUCC{panel_idx}_EQI']
            else:
                panel_models = [f'RUCC{panel_idx}_EQI_{domain}']
        
        # Filter data
        panel_df = df[df['Model'].isin(panel_models)].copy()
        
        # Collect values for dynamic ylim
        panel_values = []
        
        # Plot for each stratification group
        for group_idx, strat_group in enumerate(STRAT_GROUP_ORDER):
            group_data = panel_df[panel_df['Strat_Group'] == strat_group]
            if group_data.empty:
                continue
            
            x_base = group_idx
            
            # For each quintile
            for q_idx, quintile in enumerate(QUINTILE_ORDER):
                value_str = group_data[quintile].iloc[0] if not group_data.empty else None
                if value_str is None or pd.isna(value_str) or str(value_str).strip() == '0.00':
                    continue
                
                # Parse value
                pattern = re.compile(r'(-?\d+\.\d+)\s*\((.+?),\s*(.+?)\)([\*†]*)')
                match = pattern.match(str(value_str))
                if not match:
                    continue
                
                try:
                    mrd = float(match.group(1))
                    lower = float(match.group(2))
                    upper = float(match.group(3))
                except ValueError:
                    continue
                
                # X position
                x_pos = x_base + (q_idx - 1.5) * 0.08
                
                # Plot
                color = color_map[quintile]
                ax.errorbar(x_pos, mrd, yerr=[[mrd-lower], [upper-mrd]], 
                           fmt='o', color=color, markersize=3, capsize=2, capthick=1.2, linewidth=1.2)
                
                panel_values.extend([mrd, lower, upper])
        
        # Set axes
        ax.set_xlim(-0.5, 6.5)  # For 7 groups
        if panel_values:
            max_abs = max(abs(v) for v in panel_values)
            dynamic_limit = ((int(max_abs) // 5) + 1) * 5
            ax.set_ylim(-dynamic_limit, dynamic_limit)
        else:
            ax.set_ylim(-20, 20)
        
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3)
        
        # X ticks
        ax.set_xticks(range(7))
        ax.set_xticklabels(STRAT_GROUP_ORDER, rotation=0, ha='center')
        
        # Y label
        if panel_idx == 2:
            ax.set_ylabel('Mortality Rate Difference (95% CI)', fontsize=12, labelpad=20)
        
        # Border
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.0)
        
        # Panel label
        ax.text(-0.05, 0.5, panel_label, fontsize=11, fontweight='bold', 
               verticalalignment='center', horizontalalignment='right',
               rotation=90, transform=ax.transAxes)
    
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    
    # Save
    output_path = Path(output_dir) / filename
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return str(output_path)


if __name__ == "__main__":
    import sys
    sys.exit(main())