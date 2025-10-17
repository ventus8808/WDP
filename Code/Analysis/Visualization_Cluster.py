"""
EQI BRMS Cluster Visualization

- Reads cluster BRMS results from Result/brms_cluster_combined (e.g., C00_C97.csv)
- Builds forest-plot style figures with 3 panels (one per cluster) per ICD/scenario.
- Output: {ICD}_AllClusters_{EQI_Period}_{AAMR_Period}_Lag{lag}_brms.png

CLI Usage:
python Code/Analysis/Visualization_Cluster.py --icd C00_C97
python Code/Analysis/Visualization_Cluster.py --all
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


def get_cluster_paths(project_root: Path) -> Dict[str, Path]:
    cluster_result_dir = project_root / "Result" / "brms_cluster_combined"
    vis_dir = project_root / "Result" / "brms_cluster_Visualization"
    combined_dir = project_root / "Result" / "brms_cluster_Visualization_Combined"
    vis_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    return {"result": cluster_result_dir, "vis": vis_dir, "combined": combined_dir}


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
# Main: generate 4 panels per cluster
# ----------------------------

def list_icds(result_dir: Path) -> List[str]:
    """List ICD codes from CSV filenames"""
    icds = []
    for p in sorted(result_dir.glob("*.csv")):
        icd = p.stem  # filename without extension
        if icd not in icds:
            icds.append(icd)
    return icds


def list_clusters(result_dir: Path, icd: str) -> List[str]:
    """List cluster IDs from data for a specific ICD"""
    csv_file = result_dir / f"{icd}.csv"
    if not csv_file.exists():
        return []

    df = pd.read_csv(csv_file, dtype=str)
    # Extract cluster IDs from Model column (format: Cluster{N}_...)
    clusters = set()
    for model in df['Model'].unique():
        if model.startswith('Cluster'):
            cluster_id = model.split('_')[0].replace('Cluster', '')
            clusters.add(cluster_id)

    return sorted(list(clusters))


# Define cluster groups for combining
CLUSTER_GROUPS = {
    "All_Clusters": ["0", "1", "2"]
}


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate 3-panel cluster comparison figures per ICD")
    parser.add_argument("--icd", type=str, default=None, help="Specific ICD code, e.g., C00_C97")
    parser.add_argument("--all", action="store_true", help="Generate for all ICDs found")

    args = parser.parse_args(argv)

    paths = get_cluster_paths(project_root)

    if not args.icd and not args.all:
        print("Please specify --icd or --all")
        return 1

    if args.all:
        icds = list_icds(paths["result"])
    else:
        icds = [args.icd] if args.icd else []

    for icd in icds:
        result_csv = paths["result"] / f"{icd}.csv"
        if not result_csv.exists():
            print(f"Result CSV not found for {icd}: {result_csv}; skip.")
            continue

        df = pd.read_csv(result_csv, dtype=str)

        # For each scenario, generate a combined plot showing all clusters
        for eqi, aamr, lag in SCENARIO_ORDER:
            try:
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"].astype(int) == lag)]
            except Exception:
                # If Lag is stored as str like '5', compare as strings
                sub = df[(df["EQI_Period"] == eqi) & (df["AAMR_Period"] == aamr) & (df["Lag"] == str(lag))]
            if sub.empty:
                print(f"[{icd}] no data for scenario {eqi}/{aamr}/Lag{lag}; skipping panel.")
                continue

            out_name = f"{icd}_AllClusters_{eqi}_{aamr}_Lag{lag}_brms.png"
            out_path = paths["vis"] / out_name
            title = f"{icd} All Clusters {eqi} {aamr} Lag{lag}"

            # Generate the plot using the reference style forest plot function
            try:
                plot_reference_style_forest(
                    df=sub,
                    eqi_period=eqi,
                    aamr_period=aamr,
                    lag=lag,
                    output_dir=str(paths["vis"]),
                    icd_code=icd,
                    cluster="all",  # Not used in the function anymore
                    model_type="brms"
                )
                print(f"Generated plot: {out_path}")
            except Exception as e:
                print(f"Error generating plot for {title}: {e}")

    return 0


def plot_reference_style_forest(df, eqi_period, aamr_period, lag, output_dir=None, icd_code="C00_C97", cluster="0", model_type="brms"):
    """
    创建cluster分层的森林图 - 显示所有cluster的比较

    Parameters:
    -----------
    df : pandas.DataFrame
        包含BRMS cluster结果的数据框（单个cluster的数据）
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
    cluster : str
        cluster ID，如 '0'（这里实际上不会用到，因为我们显示所有cluster）

    Returns:
    --------
    str : 生成的图片文件路径
    """
    # 定义面板配置 - 按照cluster分层显示
    panel_labels = [
        'Cluster 0\nDisadvantaged',
        'Cluster 1\nUnbalanced',
        'Cluster 2\nAdvantageous'
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

    # 创建图形 - 3个面板，每个对应一个cluster
    fig, axes = plt.subplots(3, 1, figsize=(8, 10))
    fig.suptitle(f'{icd_code} | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}',
                 fontsize=14, y=0.98)

    # 为每个cluster面板绘制数据
    for panel_idx, panel_label in enumerate(panel_labels):
        ax = axes[panel_idx]
        cluster_id = str(panel_idx)  # 0, 1, 2

        # 筛选当前cluster的数据
        cluster_df = df[df['Model'].str.startswith(f'Cluster{cluster_id}_')].copy()

        # 收集当前面板的所有数值用于设置坐标轴范围
        panel_values = []

        # 为每个EQI类型绘制数据
        for eqi_idx, eqi_type in enumerate(eqi_types):
            model_name = f'Cluster{cluster_id}_EQI' if eqi_type == 'EQI' else f'Cluster{cluster_id}_EQI_{eqi_type}'

            # 筛选当前EQI类型的数据
            model_data = cluster_df[cluster_df['Model'] == model_name]
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
        if panel_idx == 1:  # 中间面板
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
        output_filename = output_path / f"{icd_code}_AllClusters_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"
    else:
        output_filename = f"{icd_code}_AllClusters_{eqi_period}_{aamr_period}_Lag{lag}_{model_type}.png"

    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"森林图已保存为 {output_filename}")
    plt.close()

    return str(output_filename)


if __name__ == "__main__":
    import sys
    sys.exit(main())