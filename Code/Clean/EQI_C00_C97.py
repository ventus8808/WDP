#!/usr/bin/env python3
"""
EQI × CDC Triangulation — Clean AAMR tables for C00_C97

Task
1) Remove rows with missing critical fields.
2) For "unreliable/suppressed" cases, adjust intervals using rules referenced from EQI_AAMR_Interval.py:
   - Deaths == 0 → AAMR_Lower = 0, AAMR_Upper = 0
   - Deaths in [1..9] (suppressed-like) and Population > 0 →
       AAMR_Lower = 1 / Population * 100000
       AAMR_Upper = 9 / Population * 100000
   - If CI is missing but point exists → AAMR_Lower = AAMR, AAMR_Upper = AAMR

Inputs (under project root):
- Data/Original/CDC Triangulation/AAMR/2006_2010_C00_C97.csv
- Data/Original/CDC Triangulation/AAMR/2011_2015_C00_C97.csv
- Data/Original/CDC Triangulation/AAMR/2016-2020_C00_C97.csv
(Underscore/hyphen variations tolerated via glob patterns.)

Outputs:
- Same directory as input with filename suffix "_new.csv"
  e.g., 2006_2010_C00_C97_new.csv

Conventions:
- Read project root from this file path; config.yaml is loaded for consistency,
  though paths here use the established AAMR location under Data/Original.
- Keep AAMR point and SE as-is; only adjust intervals per rules above.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

# Load config for consistency with repo convention (not strictly required here)
try:
    with CONFIG_PATH.open("r", encoding="utf-8") as f:
        CFG = yaml.safe_load(f)
except FileNotFoundError:
    CFG = {}

AAMR_DIR = PROJECT_ROOT / "Data/Original/CDC Triangulation/AAMR"

REQUIRED_COLS = [
    "COUNTY_FIPS",
    "Deaths",
    "Population",
    "AAMR",
    "AAMR_SE",
    "AAMR_Lower",
    "AAMR_Upper",
]


def _coerce_and_standardize(df: pd.DataFrame) -> pd.DataFrame:
    # Normalize column names (strip spaces)
    df = df.rename(columns={c: str(c).strip() for c in df.columns})

    # Ensure FIPS is 5-digit string
    if "COUNTY_FIPS" in df.columns:
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.extract(r"(\d+)")[0]
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].fillna("").astype(str).str.zfill(5)

    # Coerce numeric fields
    for col in ["Deaths", "Population", "AAMR", "AAMR_SE", "AAMR_Lower", "AAMR_Upper"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def _drop_missing_rows(df: pd.DataFrame) -> pd.DataFrame:
    # Drop rows with missing critical fields
    need = [c for c in ["COUNTY_FIPS", "Deaths", "Population"] if c in df.columns]
    df = df.dropna(subset=need).copy()

    # Remove invalid Population (<=0)
    df = df[df["Population"] > 0].copy()

    # Remove invalid FIPS
    df = df[df["COUNTY_FIPS"].str.fullmatch(r"\d{5}", na=False)].copy()

    return df


def _adjust_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply interval adjustments mirroring EQI_AAMR_Interval.py logic:
    - Deaths == 0 → [0, 0]
    - Deaths in [1..9] → crude bounds [1/pop, 9/pop]*100000
    - If CI missing but point exists → CI = point
    - If AAMR is unreliable (suppressed-like) or missing, substitute crude rate
    """
    # Prepare working Series
    deaths = pd.to_numeric(df["Deaths"], errors="coerce")
    pop = pd.to_numeric(df["Population"], errors="coerce")
    lower = pd.to_numeric(
        df.get("AAMR_Lower", pd.Series([pd.NA] * len(df))), errors="coerce"
    )
    upper = pd.to_numeric(
        df.get("AAMR_Upper", pd.Series([pd.NA] * len(df))), errors="coerce"
    )
    point = pd.to_numeric(df.get("AAMR", pd.Series([pd.NA] * len(df))), errors="coerce")

    lower = lower.astype("float64")
    upper = upper.astype("float64")

    # Precompute crude rate where feasible
    with pd.option_context("mode.use_inf_as_na", True):
        crude = (deaths / pop) * 100000.0

    # Deaths == 0 → [0, 0]
    zero_mask = deaths.eq(0)
    lower.loc[zero_mask] = 0.0
    upper.loc[zero_mask] = 0.0
    # Also set AAMR to 0 for zero deaths
    point.loc[zero_mask] = 0.0

    # CI missing but point available → CI = point
    no_ci_with_point = lower.isna() & upper.isna() & point.notna()
    lower.loc[no_ci_with_point] = point.loc[no_ci_with_point]
    upper.loc[no_ci_with_point] = point.loc[no_ci_with_point]

    # Suppressed-like: deaths in [1..9] → crude interval and use crude rate as AAMR
    suppressed_mask = deaths.between(1, 9, inclusive="both") & pop.notna() & (pop > 0)
    if suppressed_mask.any():
        lower.loc[suppressed_mask] = (1.0 / pop.loc[suppressed_mask]) * 100000.0
        upper.loc[suppressed_mask] = (9.0 / pop.loc[suppressed_mask]) * 100000.0
        point.loc[suppressed_mask] = crude.loc[suppressed_mask]

    # If AAMR missing but crude available, substitute crude
    point_missing_with_crude = point.isna() & pop.notna() & (pop > 0) & deaths.notna()
    if point_missing_with_crude.any():
        point.loc[point_missing_with_crude] = crude.loc[point_missing_with_crude]

    # Write back (round to 4 decimals to match AAMR precision commonly used here)
    df["AAMR_Lower"] = lower.round(4)
    df["AAMR_Upper"] = upper.round(4)
    df["AAMR"] = point.round(4)

    return df


def _find_target_files(aamr_dir: Path) -> list[Path]:
    """
    Locate target C00_C97 AAMR files allowing underscore/hyphen/space variations.
    """
    patterns = [
        "*2006*2010*_*C00_C97.csv",
        "*2011*2015*_*C00_C97.csv",
        "*2016*2020*_*C00_C97.csv",
        "*2006*2010*C00_C97.csv",
        "*2011*2015*C00_C97.csv",
        "*2016*2020*C00_C97.csv",
    ]
    hits: list[Path] = []
    for pat in patterns:
        hits.extend(sorted(aamr_dir.glob(pat)))

    # Deduplicate by normalized period token + ICD
    def normalize_key(p: Path) -> str:
        base = p.name
        period_match = re.search(r"(2006[-_]?2010|2011[-_]?2015|2016[-_]?2020)", base)
        period = period_match.group(1).replace("_", "-") if period_match else base
        return f"{period}__C00_C97"

    chosen: dict[str, Path] = {}
    for fp in hits:
        if "C00_C97" not in fp.name:
            continue
        key = normalize_key(fp)
        chosen.setdefault(key, fp)

    return list(sorted(chosen.values()))


def process_file(fp: Path) -> pd.DataFrame:
    print(f"  Processing: {fp.name}")
    df = pd.read_csv(fp, dtype={"COUNTY_FIPS": str})

    # Drop any row containing literal 'Missing' (case-insensitive) in any column
    missing_any = df.apply(
        lambda col: col.astype(str).str.contains("missing", case=False, na=False)
    )
    df = df[~missing_any.any(axis=1)].copy()

    # Coerce and clean
    df = _coerce_and_standardize(df)
    missing_before = df.isna().sum().to_dict()

    df = _drop_missing_rows(df)
    df = _adjust_intervals(df)

    # Ensure integer types for Deaths and Population in output
    df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").round(0).astype(int)
    df["Population"] = (
        pd.to_numeric(df["Population"], errors="coerce").round(0).astype(int)
    )

    print(
        f"    ✓ rows={len(df)} | "
        f"missing_before(AAMR_Lower/Upper)={missing_before.get('AAMR_Lower', 0)}/{missing_before.get('AAMR_Upper', 0)}"
    )
    return df


def main():
    print("=" * 70)
    print("EQI × CDC Triangulation — Clean AAMR C00_C97 Tables")
    print("=" * 70)
    print(f"Input directory:  {AAMR_DIR}")

    if not AAMR_DIR.exists():
        print("⚠️ AAMR directory not found. Exiting.")
        sys.exit(1)

    targets = _find_target_files(AAMR_DIR)
    # Filter strictly to C00_C97 just in case
    targets = [fp for fp in targets if "C00_C97" in fp.name]

    # Ensure we only keep the three periods if available
    want_periods = ["2006-2010", "2011-2015", "2016-2020"]
    filtered = []
    for fp in targets:
        m = re.search(r"(2006[-_]?2010|2011[-_]?2015|2016[-_]?2020)", fp.name)
        if not m:
            continue
        period = m.group(1).replace("_", "-")
        if period in want_periods:
            filtered.append(fp)

    if not filtered:
        print("⚠️ No matching C00_C97 AAMR files found.")
        sys.exit(0)

    total_written = 0
    for fp in filtered:
        try:
            out_df = process_file(fp)
            out_name = fp.with_name(fp.stem + "_new.csv")
            out_df.to_csv(out_name, index=False)
            print(f"    💾 Saved cleaned file: {out_name.name}")
            total_written += 1
        except Exception as e:
            print(f"  ⚠️ Failed on {fp.name}: {e}")

    print("=" * 70)
    print(f"Completed. Files written: {total_written}")


if __name__ == "__main__":
    main()
