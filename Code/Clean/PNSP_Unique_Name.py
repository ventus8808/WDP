#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract unique pesticide compound names from all original USGS PNSP files and
write them to a text file 'Unique_Name.txt' under the processed Pesticide folder.

Paths are read ONLY from config.yaml:
- Input:  data_sources.usgs_pnsp.original
- Output: data_sources.usgs_pnsp.processed/Unique_Name.txt
"""

import sys
from pathlib import Path

try:
    import yaml  # type: ignore
except Exception as e:
    print("ERROR: PyYAML is required. pip install pyyaml", file=sys.stderr)
    sys.exit(1)

import pandas as pd


def load_paths() -> tuple[Path, Path]:
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        print(f"ERROR: config.yaml not found at: {config_path}", file=sys.stderr)
        sys.exit(1)
    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    ds = (cfg.get("data_sources") or {}).get("usgs_pnsp") or {}
    original_rel = ds.get("original")
    processed_rel = ds.get("processed")
    if not original_rel or not processed_rel:
        print("ERROR: config.yaml missing data_sources.usgs_pnsp.original/processed", file=sys.stderr)
        sys.exit(1)
    input_dir = (project_root / original_rel).resolve()
    output_dir = (project_root / processed_rel).resolve()
    if not input_dir.exists():
        print(f"ERROR: Input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(1)
    output_dir.mkdir(parents=True, exist_ok=True)
    return input_dir, output_dir


def collect_unique_compounds(input_dir: Path) -> list[str]:
    unique = set()
    encodings = ["utf-8", "latin1", "iso-8859-1", "cp1252", "windows-1252"]
    files = list(sorted(input_dir.glob("*.txt"))) + list(sorted(input_dir.glob("*.csv")))
    for file_path in files:
        df = None
        # Try multiple encodings and delimiter inference first
        for enc in encodings:
            try:
                df = pd.read_csv(file_path, encoding=enc, sep=None, engine="python")
                break
            except UnicodeDecodeError:
                df = None
                continue
            except Exception:
                df = None
                continue
        # Fallback to common delimiters if needed
        if df is None:
            for sep in ["\t", ",", "|"]:
                try:
                    df = pd.read_csv(file_path, encoding="utf-8", sep=sep)
                    break
                except Exception:
                    df = None
                    continue
        if df is None:
            continue
        # Candidate column names often seen in PNSP datasets
        for col in ["Compound", "compound", "COMPOUND", "Compound_Name", "compound_name"]:
            if col in df.columns:
                values = df[col].dropna().astype(str).str.strip()
                unique.update(v for v in values if v)
                break
    return sorted(unique)


def main() -> None:
    input_dir, output_dir = load_paths()
    names = collect_unique_compounds(input_dir)
    output_file = output_dir / "Unique_Name.txt"
    with output_file.open("w", encoding="utf-8") as f:
        for name in names:
            f.write(f"{name}\n")
    print(f"Wrote {len(names)} unique compound names → {output_file}")


if __name__ == "__main__":
    main()


