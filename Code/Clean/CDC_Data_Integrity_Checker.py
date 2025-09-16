#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==================== 手动指定疾病分组（必须至少填写 1 个） ====================
MANUAL_ICD_GROUPS: str = "C00-C14"
# 例如："C91-C95" 或 "C15-C16,C18"；也可写绝对路径
# =====================================================================

# 输出模式：紧凑摘要（推荐）
COMPACT: bool = True

import sys
from pathlib import Path
import os
import glob
from typing import List, Tuple

try:
    import yaml  # type: ignore
except Exception as exc:  # pragma: no cover
    print(f"[ERROR] 缺少 PyYAML 依赖: {exc}")
    sys.exit(2)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"
if not CONFIG_PATH.exists():
    print(f"[ERROR] 未找到配置文件: {CONFIG_PATH}")
    sys.exit(2)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

cdc_cfg = (cfg.get("data_sources") or {}).get("cdc_wonder") or {}
base_rel = cdc_cfg.get("integrity_base_dir")
if not base_rel:
    print("[ERROR] config.yaml 缺少 data_sources.cdc_wonder.integrity_base_dir")
    sys.exit(2)

BASE_DIR = (PROJECT_ROOT / base_rel).resolve()
if not BASE_DIR.exists():
    print(f"[ERROR] 根目录不存在: {BASE_DIR}")
    sys.exit(2)

if not MANUAL_ICD_GROUPS.strip():
    print("[ERROR] 请指定疾病（在文件顶部设置 MANUAL_ICD_GROUPS，如 'C91-C95' 或 'C15-C16,C18'）")
    sys.exit(2)


def _split_icd_groups(groups: str) -> List[str]:
    return [g.strip() for g in groups.split(',') if g.strip()]


groups = _split_icd_groups(MANUAL_ICD_GROUPS)

selected_dirs: List[Path] = []
for g in groups:
    if g.startswith('/') or g.startswith('~') or '/' in g:
        subdir = Path(os.path.expanduser(g)).resolve()
    else:
        subdir = (BASE_DIR / g).resolve()
    if not subdir.exists():
        print(f"[ERROR] 指定的疾病目录不存在: {subdir}")
        sys.exit(2)
    selected_dirs.append(subdir)

csv_files: List[str] = []
for d in selected_dirs:
    csv_files.extend(glob.glob(os.path.join(str(d), "*.csv")))

if not csv_files:
    print("[ERROR] 未找到任何 CSV 文件")
    sys.exit(2)

print(f"根目录: {BASE_DIR}")
print(f"疾病目录: {[str(d) for d in selected_dirs]}")
print(f"文件数: {len(csv_files)}")

import pandas as pd

key_columns = ['Year', 'County', 'County Code', 'Sex', 'Sex Code', 
               'Race', 'Race Code', 'Ten-Year Age Groups', 'Ten-Year Age Groups Code']
stat_columns = ['Deaths', 'Population', 'Crude Rate Standard Error']

all_years = set()
files_with_missing_keys = 0
files_with_missing_stats = 0

def _year_min_max(df: pd.DataFrame) -> Tuple[str, str]:
    if 'Year' not in df.columns:
        return ("-", "-")
    years = pd.to_numeric(df['Year'], errors='coerce').dropna()
    if years.empty:
        return ("-", "-")
    return (str(int(years.min())), str(int(years.max())))

print("\n检查每个文件 (紧凑模式):")
for file_path in sorted(set(csv_files)):
    filename = os.path.basename(file_path)
    try:
        df = pd.read_csv(file_path)
        df = df.dropna(how='all')

        missing_keys = [col for col in key_columns if col not in df.columns]
        missing_stats = [col for col in stat_columns if col not in df.columns]
        if missing_keys:
            files_with_missing_keys += 1
        if missing_stats:
            files_with_missing_stats += 1

        y_min, y_max = _year_min_max(df)
        if COMPACT:
            print(f"- {filename}: rows={len(df)}, years={y_min}-{y_max}, missing_keys={len(missing_keys)}, missing_stats={len(missing_stats)}")
        else:
            print(f"\n{filename}:")
            print(f"  行数: {len(df)}")
            print(f"  年份范围: {y_min}-{y_max}")
            print(f"  缺少关键列({len(missing_keys)}): {missing_keys}")
            print(f"  缺少统计列({len(missing_stats)}): {missing_stats}")
    except Exception as e:
        print(f"- {filename}: 读取失败 -> {e}")

print("\n" + "="*50)
print("检查总结 (紧凑):")
print("="*50)
if all_years:
    print(f"年份最小-最大: {min(all_years)}-{max(all_years)}")
print(f"文件总数: {len(csv_files)}")
print(f"缺少关键列的文件数: {files_with_missing_keys}")
print(f"缺少统计列的文件数: {files_with_missing_stats}")
