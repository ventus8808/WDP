"""
EQI BRMS Climate Visualization

- Reads climate BRMS results from Result/brms_Climate (e.g., C00_C97_koppen_major.csv)
- Builds forest-plot style figures with panels per climate category per ICD/scenario.
- Output: {ICD}_{climate_type}_{EQI_Period}_{AAMR_Period}_Lag{lag}_brms.png

CLI Usage:
python Code/Analysis/Visualization_Climate.py --icd C00_C97 --climate koppen_major
python Code/Analysis/Visualization_Climate.py --all --climate koppen_major
"""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import re
import argparse
from typing import Dict, List, Optional, Tuple

import yaml
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import gridspec
import matplotlib as mpl

from Clean.Climate_Mapping import KOPPEN_MAJOR_MAP, CENSUS_REGION_MAP, CENSUS_DIVISION_MAP, DOE_MAJOR_MAP

# Simple scenario order (fixed four scenarios)
SCENARIO_ORDER: List[Tuple[str, str, int]] = [
    ("2000_2005", "2006_2010", 5),
    ("2000_2005", "2011_2015", 10),
    ("2006_2010", "2011_2015", 5),
    ("2006_2010", "2016_2020", 10),
]

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


def get_climate_paths(project_root: Path) -> Dict[str, Path]:
    climate_result_dir = project_root / "Result" / "brms_Climate"
    vis_dir = project_root / "Result" / "brms_Climate_Visualization"
    vis_dir.mkdir(parents=True, exist_ok=True)
    return {"result": climate_result_dir, "vis": vis_dir}


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
# Main: generate panels per climate category
# ----------------------------

def list_icds(result_dir: Path, climate_type: str) -> List[str]:
    """List ICD codes from CSV filenames for a specific climate type"""
    icds = []
    pattern = f"*{climate_type}.csv"
    for p in sorted(result_dir.glob(pattern)):
        icd = p.stem.replace(f"_{climate_type}", "")  # e.g., C00_C97 from C00_C97_koppen_major.csv
        if icd not in icds:
            icds.append(icd)
    return icds


def list_climate_categories(result_dir: Path, icd: str, climate_type: str) -> List[str]:
    """List climate categories from data for a specific ICD and climate type"""
    csv_file = result_dir / f"{icd}_{climate_type}.csv"
    if not csv_file.exists():
        return []

    df = pd.read_csv(csv_file, dtype=str)
    # Extract categories from Model column (format: {climate_type}_{category}_EQI)
    categories = set()
    prefix = f"{climate_type}_"
    for model in df['Model'].unique():
        if model.startswith(prefix):
            parts = model.split('_')
            if len(parts) >= 3 and parts[2] == 'EQI':
                category = parts[1]
                categories.add(category)

    return sorted(list(categories))


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate climate category comparison figures per ICD")
    parser.add_argument("--icd", type=str, default=None, help="Specific ICD code, e.g., C00_C97")
    parser.add_argument("--climate", type=str, required=True, help="Climate type, e.g., koppen_major")
    parser.add_argument("--all", action="store_true", help="Generate for all ICDs found")

    args = parser.parse_args(argv)

    paths = get_climate_paths(project_root)

    if not args.icd and not args.all:
        print("Please specify --icd or --all")
        return 1

    if args.all:
        icds = list_icds(paths["result"], args.climate)
    else:
        icds = [args.icd] if args.icd else []

    for icd in icds:
        result_csv = paths["result"] / f"{icd}_{args.climate}.csv"
        if not result_csv.exists():
            print(f"Result CSV not found for {icd} {args.climate}: {result_csv}; skip.")
            continue

        df = pd.read_csv(result_csv, dtype=str)

        # For each scenario, generate a combined plot showing all climate categories
        for eqi, aamr, lag in SCENARIO_ORDER:
            try:
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"].astype(int) == lag)]
            except Exception:
                # If Lag is stored as str like '5', compare as strings
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"] == str(lag))]
            if sub.empty:
                print(f"[{icd}] no data for scenario {eqi}/{aamr}/Lag{lag}; skipping panel.")
                continue

            out_name = f"{icd}_{args.climate}_{eqi}_{aamr}_Lag{lag}_brms.png"
            out_path = paths["vis"] / out_name
            title = f"{icd} {args.climate} {eqi} {aamr} Lag{lag}"

            # Generate the plot using the climate forest plot function
            try:
                plot_climate_forest(
                    df=sub,
                    eqi_period=eqi,
                    aamr_period=aamr,
                    lag=lag,
                    output_dir=str(paths["vis"]),
                    icd_code=icd,
                    climate_type=args.climate,
                    model_type="brms"
                )
                print(f"Generated plot: {out_path}")
            except Exception as e:
                print(f"Error generating plot for {title}: {e}")

    return 0


def plot_climate_forest(df, eqi_period, aamr_period, lag, output_dir=None, icd_code="C00_C97", climate_type="koppen_major", model_type="brms"):
    """
    创建climate分层的森林图 - 显示所有climate类别的比较

    Parameters:
    -----------
    df : pandas.DataFrame
        包含BRMS climate结果的数据框
    eqi_period : str
        EQI时期，如 '2000_2005'
    aamr_period : str
        AAMR时期，如 '2006_2010'
    lag : int
        滞后年数，如 5
    output_dir : str or Path, optional
        输出目录，默认为当前工作目录
    icd_code : str
        ICD代码，默认为 'C00_C97'
    climate_type : str
        气候类型，如 'koppen_major'

    Returns:
    --------
    str : 生成的图片文件路径
    """
    # Extract unique categories from the data
    categories = set()
    prefix = f"{climate_type}_"
    for model in df['Model'].unique():
        if model.startswith(prefix):
            parts = model.split('_')
            if len(parts) >= 4 and parts[3] == 'EQI':
                category = parts[2]
                categories.add(category)
    categories = sorted(list(categories))

    # Define panel labels based on categories and climate type
    if climate_type == 'koppen_major':
        panel_labels = [KOPPEN_MAJOR_MAP.get(cat, cat) for cat in categories]
    elif climate_type == 'census_region':
        panel_labels = [CENSUS_REGION_MAP.get(int(cat), cat) for cat in categories]
    elif climate_type == 'census_division':
        panel_labels = [CENSUS_DIVISION_MAP.get(int(cat), cat) for cat in categories]
    elif climate_type == 'doe_major':
        panel_labels = [DOE_MAJOR_MAP.get(int(cat), cat) for cat in categories]
    else:
        panel_labels = [f"{cat}" for cat in categories]

    # 定义EQI类型
    eqi_types = ['EQI', 'Air', 'Water', 'Land', 'Built', 'Social']

    # 定义颜色映射
    color_map = {
        'Q2': '#B3B3B3',  # 浅灰（70%）
        'Q3': '#808080',  # 中灰（50%）
        'Q4': '#4D4D4D',  # 深灰（30%）
        'Q5': '#000000',  # 黑（0%）
    }

    # 创建图形 - 面板数量等于类别数量
    num_panels = len(categories)
    fig, axes = plt.subplots(num_panels, 1, figsize=(8, 6 + 2 * num_panels))
    if num_panels == 1:
        axes = [axes]
    fig.suptitle(f'{icd_code} {climate_type} | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}',
                 fontsize=14, y=0.98)

    # 为每个类别面板绘制数据
    for panel_idx, (category, panel_label) in enumerate(zip(categories, panel_labels)):
        ax = axes[panel_idx]

        # 筛选当前类别的数据
        category_df = df[df['Model'].str.startswith(f'{climate_type}_{category}_')].copy()

        # 收集当前面板的所有数值用于设置坐标轴范围
        panel_values = []

        # 为每个EQI类型绘制数据
        for eqi_idx, eqi_type in enumerate(eqi_types):
            model_name = f'{climate_type}_{category}_EQI' if eqi_type == 'EQI' else f'{climate_type}_{category}_EQI_{eqi_type}'

            # 筛选当前EQI类型的数据
            model_data = category_df[category_df['Model'] == model_name]
            if model_data.empty:
                continue

            # 为每个quintile绘制点和误差线
            for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                if value_str is None or pd.isna(value_str) or str(value_str).strip() == '0.00':
                    continue

                # 解析数值
                pattern = re.compile(r'(-?\d+\.\d+)\s*\((.+?),\s*(.+?)\)([\*†]*)')
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

                # 计算x位置（每个quintile在EQI类型内偏移）
                x_pos = eqi_idx + (q_idx - 1.5) * 0.08

                # 绘制误差线和点
                color = color_map[quintile]
                ax.errorbar(x_pos, mrd, yerr=[[mrd-lower], [upper-mrd]],
                           fmt='o', color=color, markersize=3, capsize=2, capthick=1.2, linewidth=1.2)

        # 设置面板样式
        ax.set_xlim(-0.5, 5.5)

        # 为每个子图独立计算坐标轴范围
        if panel_values:
            max_abs_value = max(abs(v) for v in panel_values)
            # 向上取整到5的倍数
            dynamic_limit = ((int(max_abs_value) // 5) + 1) * 5
            ax.set_ylim(-dynamic_limit, dynamic_limit)
        else:
            ax.set_ylim(-20, 20)  # 默认值

        ax.axhline(0, color='gray', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3)

        # 设置x轴标签
        ax.set_xticks(range(6))
        ax.set_xticklabels(eqi_types)

        # 设置y轴标签
        if panel_idx == num_panels // 2:  # 中间面板
            ax.set_ylabel('Mortality Rate Difference (95% CI)', fontsize=12, labelpad=20)

        # 添加黑色边框
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(1.0)

        # 添加竖直的面板标签，放在左边贴边
        ax.text(-0.05, 0.5, panel_label, fontsize=11, fontweight='bold',
               verticalalignment='center', horizontalalignment='right',
               rotation=90, transform=ax.transAxes)

    plt.tight_layout(rect=(0, 0, 1, 0.96))

    # 保存图片
    if output_dir:
        from pathlib import Path
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        output_filename = output_path / f"{icd_code}_{climate_type}_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"
    else:
        output_filename = f"{icd_code}_{climate_type}_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"森林图已保存为 {output_filename}")
    plt.close()

    return str(output_filename)


if __name__ == "__main__":
    import sys
    sys.exit(main())