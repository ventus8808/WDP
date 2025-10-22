"""
EQI BRMS Stratified Visualization

- Reads stratified BRMS results from Result/brms_stratified (e.g., C00_C97_brms_Sex_Female.csv)
- Builds forest-plot style figures for four predefined scenarios per stratification group.
- Output: {ICD}_{Stratification}_{EQI_Period}_{AAMR_Period}_Lag{lag}_brms.png

CLI Usage:
python Code/Analysis/Visualization_Stratified.py --stratification "Sex_Female" --icd C00_C97
python Code/Analysis/Visualization_Stratified.py --all
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


def get_stratified_paths(project_root: Path) -> Dict[str, Path]:
    stratified_result_dir = project_root / "Result" / "brms_stratified"
    vis_dir = project_root / "Result" / "brms_stratified_Visualization"
    combined_dir = project_root / "Result" / "brms_stratified_Visualization_Combined"
    vis_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    return {"result": stratified_result_dir, "vis": vis_dir, "combined": combined_dir}


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
# Main: generate 4 panels per stratification group
# ----------------------------

def list_stratifications(result_dir: Path) -> List[str]:
    """List stratification groups from CSV filenames like C00_C97_brms_{strat}.csv"""
    stratifications = []
    for p in sorted(result_dir.glob("C00_C97_brms_*.csv")):
        stem = p.stem  # e.g., C00_C97_brms_Sex_Female
        parts = stem.split("_brms_")
        if len(parts) == 2:
            strat = parts[1]  # e.g., Sex_Female
            if strat not in stratifications:
                stratifications.append(strat)
    return stratifications


# Define stratification groups for combining
STRATIFICATION_GROUPS = {
    "Sex": ["Female", "Male"],
    "Race": ["American Indian or Alaska Native", "Asian or Pacific Islander", "Black or African American", "White"]
}


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate 4 panel images per stratification group")
    parser.add_argument("--stratification", type=str, default=None, help="Specific stratification, e.g., Sex_Female")
    parser.add_argument("--all", action="store_true", help="Generate for all stratifications found")

    args = parser.parse_args(argv)

    paths = get_stratified_paths(project_root)

    if not args.stratification and not args.all:
        print("Please specify --stratification or --all")
        return 1

    if args.all:
        stratifications = list_stratifications(paths["result"])
    else:
        stratifications = [args.stratification]

    for strat in stratifications:
        result_csv = paths["result"] / f"C00_C97_brms_{strat}.csv"
        if not result_csv.exists():
            print(f"Result CSV not found for {strat}: {result_csv}; skip.")
            continue

        df = pd.read_csv(result_csv, dtype=str)

        # Ensure columns exist and generate panels for fixed scenarios
        for eqi, aamr, lag in SCENARIO_ORDER:
            try:
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"].astype(int) == lag)]
            except Exception:
                # If Lag is stored as str like '5', compare as strings
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"] == str(lag))]
            if sub.empty:
                print(f"[{strat}] no data for scenario {eqi}/{aamr}/Lag{lag}; skipping panel.")
                continue

            out_name = f"C00_C97_{strat}_{eqi}_{aamr}_Lag{lag}_brms.png"
            out_path = paths["vis"] / out_name
            title = f"C00_C97 {strat} {eqi} {aamr} Lag{lag}"

            # Generate the plot using the reference style forest plot function
            try:
                plot_reference_style_forest(
                    df=sub,
                    eqi_period=eqi,
                    aamr_period=aamr,
                    lag=lag,
                    output_dir=str(paths["vis"]),
                    icd_code="C00_C97",
                    stratification=strat,
                    model_type="brms"
                )
                print(f"Generated plot: {out_path}")
            except Exception as e:
                print(f"Error generating plot for {title}: {e}")

    # After generating all individual plots, combine them
    for group_name, group_strats in STRATIFICATION_GROUPS.items():
        for eqi, aamr, lag in SCENARIO_ORDER:
            # Collect image paths for this group and scenario
            image_paths = []
            for sub_strat in group_strats:
                full_strat = f"{group_name}_{sub_strat}"
                img_path = paths["vis"] / f"C00_C97_{full_strat}_{eqi}_{aamr}_Lag{lag}_brms.png"
                if img_path.exists():
                    image_paths.append(img_path)
            if len(image_paths) == len(group_strats):
                # Combine them horizontally
                combined_img = combine_images_horizontally(image_paths)
                combined_name = f"C00_C97_{group_name}_{eqi}_{aamr}_Lag{lag}_brms.png"
                combined_path = paths["combined"] / combined_name
                combined_img.savefig(combined_path, dpi=300, bbox_inches='tight')
                print(f"Combined plot saved: {combined_path}")
                plt.close(combined_img)

    return 0


def plot_reference_style_forest(df, eqi_period, aamr_period, lag, output_dir=None, icd_code="C00_C97", stratification="Overall", model_type="brms"):
    """
    创建与您提供的参考图片完全一致的多面板森林图（分层版本）

    Parameters:
    -----------
    df : pandas.DataFrame
        包含BRMS分层结果的数据框
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
    stratification : str
        分层组，如 'Sex_Female'

    Returns:
    --------
    str : 生成的图片文件路径
    """
    # 定义面板配置 - 按照参考图片显示的顺序
    panel_labels = [
        'National',
        'Metropolitan Urbanized', 
        'Nonmetropolitan Urban',
        'Less Urban',
        'Thinly Populated'
    ]

    # 定义EQI类型
    eqi_types = ['EQI', 'Air', 'Water', 'Land', 'Built', 'Social']

    # 定义颜色映射
    color_map = {
        'Q2': '#B3B3B3',  # 浅灰（70%）
        'Q3': '#808080',  # 中灰（50%）
        'Q4': '#4D4D4D',  # 深灰（30%）
        'Q5': '#000000',  # 黑（0%）
    }

    # 创建图形 - 调整尺寸使其更窄更紧凑
    fig, axes = plt.subplots(5, 1, figsize=(8, 12))
    fig.suptitle(f'{icd_code} | {stratification} | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}',
                 fontsize=14, y=0.98)

    # 为每个面板绘制数据
    for panel_idx, panel_label in enumerate(panel_labels):
        ax = axes[panel_idx]

        # 筛选当前面板的数据
        if panel_idx == 0:  # mostratified - 使用整体模型
            panel_models = ['EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']
        else:  # RUCC模型
            rucc_num = panel_idx
            panel_models = [f'RUCC{rucc_num}_EQI', f'RUCC{rucc_num}_EQI_Air', f'RUCC{rucc_num}_EQI_Water',
                           f'RUCC{rucc_num}_EQI_Land', f'RUCC{rucc_num}_EQI_Built', f'RUCC{rucc_num}_EQI_Social']

        # 筛选数据 - df已经是过滤后的特定场景数据，只需要根据模型筛选
        panel_df = df[df['Model'].isin(panel_models)].copy()

        # 收集当前面板的所有上限值和下限值
        panel_values = []

        # 解析数据并收集数值
        for model_idx, model in enumerate(panel_models):
            model_data = panel_df[panel_df['Model'] == model]
            if model_data.empty:
                continue

            # 为每个quintile收集数值
            for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                # 修复pd.isna的使用问题
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
                    panel_values.append(mrd)
                    panel_values.append(lower)
                    panel_values.append(upper)
                except ValueError:
                    continue

        # 解析数据并绘制
        for model_idx, model in enumerate(panel_models):
            model_data = panel_df[panel_df['Model'] == model]
            if model_data.empty:
                continue

            x_base = model_idx

            # 为每个quintile绘制点和误差线
            for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                # 修复pd.isna的使用问题
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
                except ValueError:
                    continue

                # 计算x位置（每个quintile在EQI类型内偏移）- 使线更加紧凑
                x_pos = x_base + (q_idx - 1.5) * 0.08

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
        if panel_idx == 2:  # 中间面板
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

    # 修复类型错误：将列表改为元组
    plt.tight_layout(rect=(0, 0, 1, 0.96))

    # 保存图片
    if output_dir:
        from pathlib import Path
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        output_filename = output_path / f"{icd_code}_{stratification}_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"
    else:
        output_filename = f"{icd_code}_{stratification}_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"森林图已保存为 {output_filename}")
    plt.close()

    return str(output_filename)


def combine_images_horizontally(image_paths: List[Path]) -> plt.Figure:
    """Combine multiple images horizontally into a single figure."""
    images = [mpimg.imread(str(p)) for p in image_paths]
    n_images = len(images)
    
    # Assume all images have the same height, combine widths
    heights = [img.shape[0] for img in images]
    widths = [img.shape[1] for img in images]
    
    # Use the max height, sum widths
    total_width = sum(widths)
    max_height = max(heights)
    
    # Create a new figure with combined size
    fig, ax = plt.subplots(figsize=(total_width / 100, max_height / 100))  # Scale down for figsize
    ax.axis('off')
    
    current_x = 0
    for img in images:
        # Display each image at the current x position
        ax.imshow(img, extent=[current_x, current_x + img.shape[1], 0, img.shape[0]])
        current_x += img.shape[1]
    
    # Set limits
    ax.set_xlim(0, total_width)
    ax.set_ylim(0, max_height)
    
    return fig


if __name__ == "__main__":
    import sys
    sys.exit(main())