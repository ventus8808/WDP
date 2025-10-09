"""
Visulization_Combine

- Purpose: Combine four panel images per disease (ICD) horizontally into a single figure.
- Inputs: panel PNGs produced by `Code/LMM/Visulization.py` in Result/EQI_LMM_Visulization
  Filename pattern: {ICD}_{EQI_Period}_{AAMR_Period}_Lag{lag}[_MICE]_panel.png
- Outputs: Result/EQI_LMM_Visulization_Combined/{ICD}_combined[_MICE].png

Usage examples (from project root):
- Combine for a single ICD (both non-MICE and MICE if available):
  python Code/LMM/Visulization_Combine.py --icd C00_C97
- Combine for all ICDs (both sets):
  python Code/LMM/Visulization_Combine.py --all
- Only non-MICE (no imputation):
  python Code/LMM/Visulization_Combine.py --icd C00_C97 --non-mice-only
- Only MICE (imputed):
  python Code/LMM/Visulization_Combine.py --icd C00_C97 --mice-only

Notes
- Reads paths from config.yaml for result CSVs to list ICD codes, but image I/O uses the
  established folders Result/EQI_LMM_Visulization and Result/EQI_LMM_Visulization_Combined.
- Requires Pillow (PIL). If missing, install with: pip install pillow
"""
from __future__ import annotations

import re
import sys
import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import yaml

try:
    from PIL import Image  # type: ignore
    _HAS_PIL = True
except Exception:
    Image = None  # type: ignore
    _HAS_PIL = False


SCENARIO_ORDER: List[Tuple[str, str, int]] = [
    ("2000_2005", "2006_2010", 5),
    ("2000_2005", "2011_2015", 10),
    ("2006_2010", "2011_2015", 5),
    ("2006_2010", "2016_2020", 10),
]


def load_config(project_root: Path) -> dict:
    cfg_path = project_root / "config.yaml"
    with cfg_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_paths(project_root: Path, cfg: dict, model: str = "eqi_lmm") -> dict:
    # Choose which result dir to read ICD list from
    if model == "mice":
        result_key = "eqi_lmm_mi" if "eqi_lmm_mi" in cfg.get("result_directories", {}) else "eqi_lmm"
    elif model == "brms":
        result_key = "brms" if "brms" in cfg.get("result_directories", {}) else "eqi_lmm"
    else:
        result_key = "eqi_lmm"

    result_eqi_lmm = project_root / cfg["result_directories"][result_key]

    # Visualization and combined output directories per model
    if model == "brms":
        vis_dir = project_root / "Result" / "brms_Visulization"
        combined_dir = project_root / "Result" / "brms_Visulization_Combined"
    elif model == "mice":
        vis_dir = project_root / "Result" / "EQI_LMM_Visulization"
        combined_dir = project_root / "Result" / "EQI_LMM_Visulization_Combined_MICE"
    else:
        vis_dir = project_root / "Result" / "EQI_LMM_Visulization"
        combined_dir = project_root / "Result" / "EQI_LMM_Visulization_Combined"

    vis_dir.mkdir(parents=True, exist_ok=True)
    combined_dir.mkdir(parents=True, exist_ok=True)
    return {"result": result_eqi_lmm, "vis": vis_dir, "combined": combined_dir}


def list_icds_from_results(result_dir: Path, cfg: dict, model: str = "eqi_lmm") -> List[str]:
    icds: List[str] = []
    if model == "brms":
        # Use brms filename template from config if available
        brms_cfg = cfg.get("brms_analysis", {}).get("results", {})
        template = brms_cfg.get("filename_template", "brms_{cancer_type}_Results.csv")
        glob_pat = template.replace("{cancer_type}", "*")
        for p in sorted(result_dir.glob(glob_pat)):
            m = re.match(r"brms_(.+?)_Results", p.name)
            if m:
                icd = m.group(1)
                if icd not in icds:
                    icds.append(icd)
        return icds

    # Default: look for LMM_*.csv and strip _FDR and _MICE suffixes
    for p in sorted(result_dir.glob("LMM_*.csv")):
        stem = p.stem.replace("LMM_", "")
        if stem.endswith("_FDR"):
            stem = stem[:-4]
        if stem.endswith("_MICE"):
            stem = stem[:-5]
        if stem not in icds:
            icds.append(stem)
    return icds


def parse_panel_filename(icd: str, fname: str) -> Optional[Tuple[str, str, int, bool]]:
    """Return (eqi, aamr, lag, is_mice) from panel filename or None if not match.

    The original implementation used a very strict full-match regex which failed
    for small filename variations. This version is more permissive: it ensures
    the filename starts with the ICD, then searches for two year-range tokens
    ("YYYY_YYYY") and a Lag value. It still returns None if the essential
    pieces aren't found.
    """
    # Quick guard: filename should start with the ICD followed by an underscore
    if not fname.startswith(f"{icd}_"):
        return None

    # Find year-range tokens like '2000_2005' (we expect at least two of them)
    year_ranges = re.findall(r"(\d{4}_\d{4})", fname)
    if len(year_ranges) < 2:
        return None
    eqi = year_ranges[0]
    aamr = year_ranges[1]

    # Find Lag (case-insensitive, allow optional separators)
    lag_m = re.search(r"[Ll]ag[_-]?(\d+)", fname)
    if not lag_m:
        return None
    lag = int(lag_m.group(1))

    # Detect MICE presence (allow variations like '_MICE' before '_panel' or anywhere)
    is_mice = bool(re.search(r"_MICE(?=_panel|\.|$)", fname)) or ("_MICE" in fname)

    return eqi, aamr, lag, is_mice


def sort_scenarios(items: List[Tuple[str, str, int]]) -> List[Tuple[str, str, int]]:
    order_map = {k: i for i, k in enumerate(SCENARIO_ORDER)}
    return sorted(items, key=lambda x: (order_map.get(x, 999), x[0], x[1], x[2]))


def combine_icd(icd: str, vis_dir: Path, out_dir: Path, model: str = "eqi_lmm") -> Optional[Path]:
    # model-aware pattern: MICE panels contain _MICE in filename
    suffix = "_MICE" if model == "mice" else ""
    files = sorted(vis_dir.glob(f"{icd}_*_Lag*{suffix}_panel.png"))
    if not files:
        print(f"[Combine] No panel images for {icd} (model={model}); skip.")
        return None

    parsed = []
    for p in files:
        info = parse_panel_filename(icd, p.name)
        if info is None:
            continue
        eqi, aamr, lag, is_mice = info
        # ensure mice-only for mice model
        if (model == "mice") != is_mice:
            continue
        parsed.append(((eqi, aamr, lag), p))
    if not parsed:
        print(f"[Combine] No matching panel files for {icd} (model={model}); skip.")
        return None

    # Order scenarios by SCENARIO_ORDER
    scenarios = sort_scenarios([k for k, _ in parsed])
    ordered_files: List[Path] = []
    for sc in scenarios:
        for (k, p) in parsed:
            if k == sc and p not in ordered_files:
                ordered_files.append(p)
                break

    if len(ordered_files) < 2:
        print(f"[Combine] Only {len(ordered_files)} panel(s) for {icd} (model={model}); need >=2; skip.")
        return None

    # Use model-specific combined filename to avoid overwriting between models
    if model == "mice":
        out_path = out_dir / f"{icd}_combined_MICE.png"
    elif model == "brms":
        out_path = out_dir / f"{icd}_combined_BRMS.png"
    else:
        out_path = out_dir / f"{icd}_combined.png"

    # Require Pillow for concatenation (simpler, faster, and consistent)
    if not _HAS_PIL:
        print("Pillow is required for image concatenation. Install with: pip install pillow")
        return None

    images = [Image.open(str(p)).convert("RGBA") for p in ordered_files]
    max_h = max(im.height for im in images)
    normalized = []
    for im in images:
        if im.height != max_h:
            new_w = int(im.width * (max_h / im.height))
            resample = getattr(Image, 'LANCZOS', getattr(Image, 'BICUBIC', 2))
            im = im.resize((new_w, max_h), resample)
        normalized.append(im)
    total_w = sum(im.width for im in normalized)
    canvas = Image.new("RGBA", (total_w, max_h), (255, 255, 255, 0))
    x = 0
    for im in normalized:
        canvas.paste(im, (x, 0))
        x += im.width
    rgb = Image.new("RGB", canvas.size, (255, 255, 255))
    rgb.paste(canvas, mask=canvas.split()[3])
    rgb.save(out_path, format="PNG", compress_level=1)
    print(f"[Combine] Saved: {out_path}")
    return out_path


def main(argv: Optional[List[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[2]
    cfg = load_config(project_root)
    parser = argparse.ArgumentParser(description="Combine 4 panel images per ICD into a horizontal strip")
    parser.add_argument("--icd", type=str, default=None, help="Specific ICD code, e.g., C00_C97")
    parser.add_argument("--all", action="store_true", help="Combine for all ICDs found in Result")
    parser.add_argument("--model", type=str, default="eqi_lmm", choices=["eqi_lmm", "mice", "brms"], help="Model: eqi_lmm (default), mice, or brms")

    args = parser.parse_args(argv)

    if not args.icd and not args.all:
        print("Please specify --icd or --all")
        return 1

    paths = get_paths(project_root, cfg, model=args.model)

    # Build ICD list
    if args.all:
        icds = list_icds_from_results(paths["result"], cfg, model=args.model)
    else:
        icds = [args.icd]

    made_any = False
    for icd in icds:
        out = combine_icd(icd, paths["vis"], paths["combined"], model=args.model)
        made_any = made_any or (out is not None)

    return 0 if made_any else 2


if __name__ == "__main__":
    sys.exit(main())
