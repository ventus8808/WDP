"""
EQI LMM Visualization

- Reads LMM results from Result/EQI_LMM (paths via config.yaml)
- Builds a forest-plot style figure for four predefined scenarios per ICD.
- Output: Result/EQI_LMM_Visulization/{ICD}_{EQI_Period}_{AAMR_Period}_Lag{lag}[_MICE]_panel.png

CLI (new):
python Code/LMM/Visulization.py --icd C00_C97 --model eqi_lmm
python Code/LMM/Visulization.py --all --model mice

If no --model specified, defaults to eqi_lmm.
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


def get_paths(project_root: Path, cfg: dict, model: str = "eqi_lmm") -> Dict[str, Path]:
    # Determine result directory for the chosen model
    if model == "mice":
        result_dir = project_root / cfg["result_directories"].get("eqi_lmm_mi", cfg["result_directories"]["eqi_lmm"])
    elif model == "brms":
        # prefer brms_analysis.results.output_dir if present
        brms_cfg = cfg.get("brms_analysis", {}).get("results", {})
        if brms_cfg.get("output_dir"):
            result_dir = project_root / brms_cfg.get("output_dir")
        else:
            result_dir = project_root / cfg["result_directories"]["eqi_lmm"]
    else:
        result_dir = project_root / cfg["result_directories"]["eqi_lmm"]

    # Use the same visualization directories for consistent naming; file names include model suffix if needed
    vis_dir = project_root / "Result" / "EQI_LMM_Visulization"
    combined_dir = project_root / "Result" / "EQI_LMM_Visulization_Combined"
    vis_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    return {"result": result_dir, "vis": vis_dir, "combined": combined_dir}


def _list_icds_for_brms(result_dir: Path, cfg: dict) -> List[str]:
    # brms filename template likely 'brms_{cancer_type}_Results.csv'
    brms_cfg = cfg.get("brms_analysis", {}).get("results", {})
    template = brms_cfg.get("filename_template", "brms_{cancer_type}_Results.csv")
    # derive glob pattern by replacing {cancer_type} with *
    glob_pat = template.replace("{cancer_type}", "*")
    icds = []
    for p in sorted((result_dir).glob(glob_pat)):
        name = p.name
        # extract cancer_type via regex: text between 'brms_' and '_Results'
        m = re.match(r"brms_(.+?)_Results", name)
        if m:
            icd = m.group(1)
            if icd not in icds:
                icds.append(icd)
    return icds


# ----------------------------
# Parsing utilities (simplified)
# ----------------------------

EFFECT_RE = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*\(\s*([+-]?\d+(?:\.\d+)?)\s*,\s*([+-]?\d+(?:\.\d+)?)\s*\)\s*$")


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
QUINTILE_COLORS = {
    "Q2": "#1b9e77",
    "Q3": "#66a61e",
    "Q4": "#e6ab02",
    "Q5": "#d95f02",
}

# ----------------------------
# Data shaping for plotting (unchanged)
# ----------------------------

# Removed collect_plot_instructions - no longer needed for new grid layout
# Removed unused constants - now using direct labels in draw_panel_grid


def collect_plot_instructions(df: pd.DataFrame) -> List[dict]:
    items: List[dict] = []

    # 1) Overall
    items.append({"type": "heading", "label": "Overall"})
    for domain in DOMAINS_OVERALL:
        sub = df[df["Model"] == domain]
        if sub.empty:
            continue
        items.append({"type": "group_label", "label": domain})
        row = sub.iloc[0]
        for q in QUINTILE_ORDER:
            parsed = parse_effect_cell(row.get(q))
            if parsed is None:
                continue
            mrd, lcl, ucl, is_sig, stars = parsed
            items.append({
                "type": "data_row",
                "group": domain,
                "quintile": q,
                "mrd": mrd,
                "lcl": lcl,
                "ucl": ucl,
                "is_sig": is_sig,
                "stars": stars,
            })

    # 2) RUCC 1..4 groups
    for rucc in [1, 2, 3, 4]:
        items.append({"type": "heading", "label": f"RUCC {rucc}"})
        for domain in DOMAINS_SUFFIX:
            model_name = f"RUCC{rucc}_{domain}"
            sub = df[df["Model"] == model_name]
            if sub.empty:
                continue
            items.append({"type": "group_label", "label": domain})
            row = sub.iloc[0]
            for q in QUINTILE_ORDER:
                parsed = parse_effect_cell(row.get(q))
                if parsed is None:
                    continue
                mrd, lcl, ucl, is_sig, stars = parsed
                items.append({
                    "type": "data_row",
                    "group": domain,
                    "quintile": q,
                    "mrd": mrd,
                    "lcl": lcl,
                    "ucl": ucl,
                    "is_sig": is_sig,
                    "stars": stars,
                })

    return items


# ----------------------------
# Plotting functions (unchanged)
# ----------------------------

# Drawing function removed - no plotting functionality


# ----------------------------
# Main: generate 4 panels per ICD/model
# ----------------------------

def list_icds_from_results(result_dir: Path, cfg: dict, model: str = "eqi_lmm") -> List[str]:
    # For LMM/MICE, expect LMM_{ICD}.csv; exclude _FDR duplicates
    if model == "brms":
        return _list_icds_for_brms(result_dir, cfg)
    icds = []
    for p in sorted(result_dir.glob("LMM_*.csv")):
        stem = p.stem.replace("LMM_", "")
        if stem.endswith("_FDR"):
            stem = stem[:-4]
        if stem not in icds:
            icds.append(stem)
    return icds


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Generate 4 panel images per ICD for a given model")
    parser.add_argument("--icd", type=str, default=None, help="Specific ICD code, e.g., C00_C97")
    parser.add_argument("--all", action="store_true", help="Generate for all ICDs found in result dir")
    parser.add_argument("--model", type=str, default="eqi_lmm", choices=["eqi_lmm", "mice", "brms"], help="Model: eqi_lmm (default), mice, or brms")

    args = parser.parse_args(argv)

    cfg = load_config(project_root)
    paths = get_paths(project_root, cfg, model=args.model)

    if not args.icd and not args.all:
        print("Please specify --icd or --all")
        return 1

    if args.all:
        icds = list_icds_from_results(paths["result"], cfg, model=args.model)
    else:
        icds = [args.icd]

    for icd in icds:
        # determine result csv per model
        if args.model == "brms":
            brms_cfg = cfg.get("brms_analysis", {}).get("results", {})
            template = brms_cfg.get("filename_template", "brms_{cancer_type}_Results.csv")
            fname = template.replace("{cancer_type}", icd)
            result_csv = paths["result"] / fname
        elif args.model == "mice":
            # Try LMM_{ICD}_MICE.csv first, then fallback to LMM_{ICD}.csv
            candidate = paths["result"] / f"LMM_{icd}_MICE.csv"
            if candidate.exists():
                result_csv = candidate
            else:
                result_csv = paths["result"] / f"LMM_{icd}.csv"
        else:
            result_csv = paths["result"] / f"LMM_{icd}.csv"

        if not result_csv.exists():
            print(f"Result CSV not found for {icd}: {result_csv}; skip.")
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
                print(f"[{icd}] no data for scenario {eqi}/{aamr}/Lag{lag}; skipping panel.")
                continue
            suffix = "_MICE" if args.model == "mice" else ""
            out_name = f"{icd}_{eqi}_{aamr}_Lag{lag}{suffix}_panel.png"
            out_path = paths["vis"] / out_name
            title = f"{icd} {eqi} {aamr} Lag{lag}"
            print(f"Skipping plot generation for: {title} (drawing function removed)")

    return 0


def plot_multi_panel_forest(df, eqi_period, aamr_period, lag, title):
    """
    创建与参考图片完全一致的多面板森林图
    """
    # 定义面板配置
    panels = [
        ('mostratified', 'mostratified'),
        ('metropolitan urban', 'metropolitan urban'),
        ('nonmetropolitan urban', 'nonmetropolitan urban'),
        ('less urban', 'less urban'),
        ('thinly populated', 'thinly populated')
    ]
    
    # 定义EQI类型
    eqi_types = ['EQI', 'Air', 'Water', 'Land', 'Built', 'Social']
    
    # 定义颜色映射
    color_map = {
        'Q2': '#2E8B57',  # Sea Green
        'Q3': '#4169E1',  # Royal Blue  
        'Q4': '#FF8C00',  # Dark Orange
        'Q5': '#DC143C',  # Crimson
    }
    
    # 创建图形
    fig, axes = plt.subplots(5, 1, figsize=(12, 16))
    fig.suptitle(f'C00_C97 | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}', 
                 fontsize=14, y=0.98)
    
    # 为每个面板绘制数据
    for panel_idx, (panel_key, panel_label) in enumerate(panels):
        ax = axes[panel_idx]
        
        # 筛选当前面板的数据
        if panel_idx == 0:  # mostratified - 使用整体模型
            panel_models = ['EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']
        else:  # RUCC模型
            rucc_num = panel_idx
            panel_models = [f'RUCC{rucc_num}_EQI', f'RUCC{rucc_num}_EQI_Air', f'RUCC{rucc_num}_EQI_Water', 
                           f'RUCC{rucc_num}_EQI_Land', f'RUCC{rucc_num}_EQI_Built', f'RUCC{rucc_num}_EQI_Social']
        
        # 筛选数据
        panel_df = df[(df['EQI_Period'] == eqi_period) &
                     (df['AAMR_Period'] == aamr_period) &
                     (df['Lag'] == lag) &
                     (df['Model'].isin(panel_models))].copy()
        
        # 解析数据并绘制
        for model_idx, model in enumerate(panel_models):
            model_data = panel_df[panel_df['Model'] == model]
            if model_data.empty:
                continue
                
            x_base = model_idx * 1.0
            
            # 为每个quintile绘制点和误差线
            for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                if pd.isna(value_str) or value_str == '0.00':
                    continue
                    
                # 解析数值
                pattern = re.compile(r'(-?\d+\.\d+)\s*\((.+?),\s*(.+?)\)([\*†]*)')
                match = pattern.match(str(value_str))
                if not match:
                    continue
                    
                mrd = float(match.group(1))
                lower = float(match.group(2))
                upper = float(match.group(3))
                
                # 计算x位置（每个quintile在EQI类型内偏移）
                x_pos = x_base + (q_idx - 1.5) * 0.1
                
                # 绘制误差线和点
                color = color_map[quintile]
                ax.errorbar(x_pos, mrd, yerr=[[mrd-lower], [upper-mrd]], 
                           fmt='o', color=color, markersize=4, capsize=3, capthick=1)
        
        # 设置面板样式
        ax.set_xlim(-0.5, 5.5)
        ax.set_ylim(-20, 20)
        ax.axhline(0, color='gray', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3)
        
        # 设置x轴标签
        ax.set_xticks(range(6))
        ax.set_xticklabels(eqi_types)
        
        # 设置y轴标签
        if panel_idx == 2:  # 中间面板
            ax.set_ylabel('Mortality Rate Difference (95% CI)', fontsize=12)
        
        # 添加面板标签
        ax.text(-0.4, 15, panel_label, fontsize=11, fontweight='bold', 
               verticalalignment='center')
        
        # 移除顶部和右侧边框
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存图片
    output_filename = f"multi_panel_forest_plot_{eqi_period}_{aamr_period}_lag{lag}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"多面板森林图已保存为 {output_filename}")
    plt.show()
    plt.close()


# 示例使用函数
def create_sample_multi_panel_plot():
    """
    创建示例多面板森林图
    """
    try:
        # 从配置文件读取路径
        from pathlib import Path
        import yaml
        
        project_root = Path(__file__).resolve().parents[2]
        config_path = project_root / "config.yaml"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        # 获取结果目录路径
        result_dir = project_root / config['result_directories']['eqi_lmm']
        data_file = result_dir / 'LMM_C00_C97_FDR.csv'
        
        if not data_file.exists():
            print(f"数据文件不存在: {data_file}")
            return
            
        # 读取数据
        df = pd.read_csv(data_file)
        
        # 创建示例图片
        plot_multi_panel_forest(
            df=df,
            eqi_period='2000_2005',
            aamr_period='2006_2010',
            lag=5,
            title='Sample Multi-Panel Forest Plot'
        )
        
    except Exception as e:
        print(f"创建示例图时出错: {e}")


def create_reference_multi_panel_plot():
    """
    创建与参考图片完全一致的多面板森林图
    """
    try:
        # 从配置文件读取路径  
        from pathlib import Path
        import yaml
        
        project_root = Path(__file__).resolve().parents[2]
        config_path = project_root / "config.yaml"
        
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        # 获取结果目录路径
        result_dir = project_root / config['result_directories']['eqi_lmm']
        data_file = result_dir / 'LMM_C00_C97_FDR.csv'
        
        if not data_file.exists():
            print(f"数据文件不存在: {data_file}")
            return
            
        # 读取数据
        df = pd.read_csv(data_file)
        
        # 定义面板配置 - 按照图片显示的顺序
        panel_labels = [
            'mostratified',
            'metropolitan urban', 
            'nonmetropolitan urban',
            'less urban',
            'thinly populated'
        ]
        
        # 定义EQI类型
        eqi_types = ['EQI', 'Air', 'Water', 'Land', 'Built', 'Social']
        
        # 定义颜色映射 - 改为黑白灰度
        color_map = {
    'Q2': '#B3B3B3',  # 稍加深的浅灰（约70%）
    'Q3': '#808080',  # 中灰（50%）
    'Q4': '#4D4D4D',  # 深灰（30%）
    'Q5': '#000000',  # 黑（0%）
}
      
        
        # 设置参数
        eqi_period = '2000_2005'
        aamr_period = '2006_2010'
        lag = 5
        
        # 创建图形 - 调整尺寸使其更窄更紧凑
        fig, axes = plt.subplots(5, 1, figsize=(8, 12))
        fig.suptitle(f'C00_C97 | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}', 
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
            
            # 筛选数据
            panel_df = df[(df['EQI_Period'] == eqi_period) &
                         (df['AAMR_Period'] == aamr_period) &
                         (df['Lag'] == lag) &
                         (df['Model'].isin(panel_models))].copy()
            
            print(f"Panel {panel_idx} ({panel_label}): 找到 {len(panel_df)} 行数据")
            
            # 解析数据并绘制
            for model_idx, model in enumerate(panel_models):
                model_data = panel_df[panel_df['Model'] == model]
                if model_data.empty:
                    continue
                    
                x_base = model_idx
                
                # 为每个quintile绘制点和误差线
                for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                    value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                    if pd.isna(value_str) or value_str == '0.00':
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
            ax.set_ylim(-20, 20)
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
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # 保存图片
        output_filename = f"reference_multi_panel_forest_plot_{eqi_period}_{aamr_period}_lag{lag}.png"
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"参考多面板森林图已保存为 {output_filename}")
        plt.close()
        
    except Exception as e:
        print(f"创建参考图时出错: {e}")
        import traceback
        traceback.print_exc()


def plot_reference_style_forest(df, eqi_period, aamr_period, lag, output_dir=None, icd_code="C00_C97"):
    """
    创建与您提供的参考图片完全一致的多面板森林图
    
    Parameters:
    -----------
    df : pandas.DataFrame
        包含LMM结果的数据框
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
    
    Returns:
    --------
    str : 生成的图片文件路径
    """
    # 定义面板配置 - 按照参考图片显示的顺序
    panel_labels = [
        'mostratified',
        'metropolitan urban', 
        'nonmetropolitan urban',
        'less urban',
        'thinly populated'
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
    fig.suptitle(f'{icd_code} | Lag {lag} years | AAMR {aamr_period} | EQI {eqi_period}', 
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
        
        # 筛选数据
        panel_df = df[(df['EQI_Period'] == eqi_period) &
                     (df['AAMR_Period'] == aamr_period) &
                     (df['Lag'] == lag) &
                     (df['Model'].isin(panel_models))].copy()
        
        # 解析数据并绘制
        for model_idx, model in enumerate(panel_models):
            model_data = panel_df[panel_df['Model'] == model]
            if model_data.empty:
                continue
                
            x_base = model_idx
            
            # 为每个quintile绘制点和误差线
            for q_idx, quintile in enumerate(['Q2', 'Q3', 'Q4', 'Q5']):
                value_str = model_data[quintile].iloc[0] if not model_data.empty else None
                if pd.isna(value_str) or value_str == '0.00':
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
        ax.set_ylim(-20, 20)
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
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存图片
    if output_dir:
        from pathlib import Path
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        output_filename = output_path / f"reference_style_forest_{icd_code}_{eqi_period}_{aamr_period}_lag{lag}.png"
    else:
        output_filename = f"reference_style_forest_{icd_code}_{eqi_period}_{aamr_period}_lag{lag}.png"
    
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"参考样式森林图已保存为 {output_filename}")
    plt.close()
    
    return str(output_filename)


if __name__ == "__main__":
    import sys
    
    # 如果有命令行参数，运行原main函数
    if len(sys.argv) > 1:
        if sys.argv[1] == "--reference":
            # 创建参考图
            create_reference_multi_panel_plot()
        else:
            sys.exit(main())
    else:
        # 否则创建示例多面板图
        create_sample_multi_panel_plot()