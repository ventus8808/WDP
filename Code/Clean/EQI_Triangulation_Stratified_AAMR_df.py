#!/usr/bin/env python3
"""
EQI × CDC Triangulation (Stratified) — AAMR long table
Merges stratified triangulation AAMR results with EQI covariates to produce a long CSV.

Output contains rows per county × time_period × (Outcome+Stratum) × Lag_Years,
with columns AAMR, AAMR_SE, AAMR_Lower/AAMR_Upper, plus RUCC/EQI covariates.

Inputs:
- Data/Original/CDC Triangulation/Stratified_AAMR/*.csv
  Filename examples:
    2016-2020_G20_G30_G12.2_F01_F03_White.csv  -> year='2016-2020', icd='G20_G30_G12.2_F01_F03', stratum='White'
    2016-2020_C00_C97_Male.csv                 -> year='2016-2020', icd='C00_C97', stratum='Male'

Output:
- Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Stratifed.csv  (note: 'Stratifed' by request)

Conventions:
- Cancer_Type embeds the stratification: e.g., 'NDD_Male', 'C00_C97_White'
- Valid EQI-AAMR lag mappings:
  * EQI 2000-2005 (code '0005') + 5y lag -> AAMR 2006-2010
  * EQI 2000-2005 (code '0005') + 10y lag -> AAMR 2011-2015
  * EQI 2000-2005 (code '0005') + 15y lag -> AAMR 2016-2020
  * EQI 2006-2010 (code '0610') + 5y lag -> AAMR 2011-2015
  * EQI 2006-2010 (code '0610') + 10y lag -> AAMR 2016-2020
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml

# ---------------- Paths & Config ----------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    CFG = yaml.safe_load(f)

# Input: stratified triangulation AAMR results
TRIANGULATION_AAMR_DIR = (
    PROJECT_ROOT / "Data/Original/CDC Triangulation/Stratified_AAMR"
)

# Output: merged long table (typo 'Stratifed' preserved per request)
OUTPUT_DIR = PROJECT_ROOT / "Data/Processed/df_EQI_AAMR_Triangulation"
OUTPUT_FILE = OUTPUT_DIR / "EQI_AAMR_Stratifed.csv"

# EQI data (from config)
EQI_DIR = PROJECT_ROOT / CFG["data_sources"]["epa_eqi"]["processed"]

# Smoking data
SMOKING_PATH = (
    PROJECT_ROOT
    / CFG["data_directories"]["processed"]
    / "Smoking"
    / "County_Smoking.csv"
)

# Cluster data
CLUSTER_PATH = PROJECT_ROOT / "Result/Cluster_Visualization/EQI_Clusters_All_K.csv"
CLUSTER_K_VALUES = [3, 4, 5]

# Climate zone data
CLIMATE_PATH = (
    PROJECT_ROOT / CFG["eqi_aamr_outputs"]["base_dir"] / "Climate_Zone_Processed.csv"
)

# EQI columns to merge
EQI_COLS = [
    "RUCC",
    "EQI",
    "EQI_Air",
    "EQI_Water",
    "EQI_Land",
    "EQI_Built",
    "EQI_Social",
    "RUCC_EQI",
    "RUCC_EQI_Air",
    "RUCC_EQI_Water",
    "RUCC_EQI_Land",
    "RUCC_EQI_Built",
    "RUCC_EQI_Social",
]

# Cluster columns to merge
CLUSTER_COLS = [f"cluster_{k}" for k in CLUSTER_K_VALUES]

# Climate columns to merge
CLIMATE_COLS = ["census_region", "koppen_major", "doe_major"]


# ---------------- Helpers ----------------


def _parse_stratified_aamr_filename(
    filename: str,
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Parse stratified AAMR filename to (year_range, icd_raw, stratum).
    Expected format: <YYYY-YYYY>_<ICD tokens>_<Stratum>.csv
    """
    name = filename.replace(".csv", "")
    parts = name.split("_")
    if len(parts) < 2:
        return None, None, None

    year = parts[0]
    if not re.fullmatch(r"\d{4}-\d{4}", year):
        return None, None, None

    if len(parts) == 2:
        # Edge case: <YYYY-YYYY>_<Stratum>.csv (unlikely for AAMR outputs)
        return year, None, parts[1]

    # Join middle tokens as ICD block, last token is stratum
    icd = "_".join(parts[1:-1])
    stratum = parts[-1]
    return year, icd, stratum


def _normalize_stratum(s: str) -> str:
    """
    Normalize stratum values into succinct labels used in Cancer_Type.
    """
    if s is None:
        return "Unknown"
    s_clean = s.strip()

    # Map common verbose labels to concise labels
    mapping = {
        "Male": "Male",
        "Female": "Female",
        "White": "White",
        "Black": "Black",
        "Black or African American": "Black",
        "African American": "Black",
        "Asian": "Asian",
        "Asian or Pacific Islander": "Asian",
        "American Indian": "Indian",
        "American Indian or Alaska Native": "Indian",
        "Alaska Native": "Indian",
        "Indian": "Indian",
    }
    if s_clean in mapping:
        return mapping[s_clean]

    # Fallback: simplify to alphanumeric/underscore
    s2 = re.sub(r"[^A-Za-z0-9]+", "_", s_clean)
    s2 = re.sub(r"_+", "_", s2).strip("_")
    return s2 or "Unknown"


def _icd_to_outcome(icd_raw: Optional[str]) -> Optional[str]:
    """
    Map ICD block to a concise outcome label where applicable.

    Known mappings:
    - C00_C97 -> C00_C97 (All-site cancer)
    - G20_G30_G12.2_F01_F03 -> NDD (Neurodegenerative disorders)
    Otherwise: return the ICD block with '-' replaced by '_'
    """
    if not icd_raw:
        return None

    icd_fmt = icd_raw.replace("-", "_")
    mapping = {
        "C00_C97": "C00_C97",
        "G20_G30_G12.2_F01_F03": "NDD",
    }
    return mapping.get(icd_fmt, icd_fmt)


def _load_eqi_by_period() -> dict:
    """Load EQI datasets for both periods (2000-2005, 2006-2010) with column harmonization."""
    d = {}

    # Column name mapping: old -> new
    column_mapping = {
        "EQI_air": "EQI_Air",
        "EQI_water": "EQI_Water",
        "EQI_land": "EQI_Land",
        "EQI_built": "EQI_Built",
        "EQI_Sociodemographic": "EQI_Social",
        "RUCC_EQI_air": "RUCC_EQI_Air",
        "RUCC_EQI_water": "RUCC_EQI_Water",
        "RUCC_EQI_land": "RUCC_EQI_Land",
        "RUCC_EQI_built": "RUCC_EQI_Built",
        "RUCC_EQI_Sociodemographic": "RUCC_EQI_Social",
    }

    for code in ("0005", "0610"):
        fp = EQI_DIR / f"EQI{code}.csv"
        if fp.exists():
            t = pd.read_csv(fp)
            t["COUNTY_FIPS"] = t["COUNTY_FIPS"].astype(str).str.zfill(5)

            # Detect and map column names if needed
            columns_to_rename = {}
            for old_col, new_col in column_mapping.items():
                if old_col in t.columns and new_col not in t.columns:
                    columns_to_rename[old_col] = new_col

            if columns_to_rename:
                t = t.rename(columns=columns_to_rename)
                print(
                    f"  📝 EQI{code}: Mapped columns {list(columns_to_rename.keys())} -> {list(columns_to_rename.values())}"
                )

            d[code] = t
        else:
            print(f"  ⚠️ EQI{code}.csv not found at {fp}")

    return d


def _load_smoking() -> pd.DataFrame | None:
    """Load smoking rate data (SR_Total -> Smoking_Rate)."""
    if SMOKING_PATH.exists():
        df = pd.read_csv(SMOKING_PATH)
        if "COUNTY_FIPS" in df.columns and "SR_Total" in df.columns:
            df = df[["COUNTY_FIPS", "SR_Total"]].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            df = df.rename(columns={"SR_Total": "Smoking_Rate"})
            return df
        else:
            print(f"  ⚠️ Smoking data missing required columns")
    else:
        print(f"  ⚠️ Smoking data not found at {SMOKING_PATH}")
    return None


def _load_clusters() -> pd.DataFrame | None:
    """Load cluster data; select COUNTY_FIPS and cluster columns only."""
    if CLUSTER_PATH.exists():
        df = pd.read_csv(CLUSTER_PATH)
        if "COUNTY_FIPS" in df.columns:
            cols = ["COUNTY_FIPS"] + [c for c in CLUSTER_COLS if c in df.columns]
            df = df[cols].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            return df
        else:
            print(f"  ⚠️ Cluster data missing COUNTY_FIPS column")
    else:
        print(f"  ⚠️ Cluster data not found at {CLUSTER_PATH}")
    return None


def _load_climate() -> pd.DataFrame | None:
    """Load climate zone data; select COUNTY_FIPS and requested climate columns."""
    if CLIMATE_PATH.exists():
        df = pd.read_csv(CLIMATE_PATH, dtype={"COUNTY_FIPS": str})
        if "COUNTY_FIPS" in df.columns:
            cols = ["COUNTY_FIPS"] + [c for c in CLIMATE_COLS if c in df.columns]
            df = df[cols].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            return df
        else:
            print(f"  ⚠️ Climate data missing COUNTY_FIPS column")
    else:
        print(f"  ⚠️ Climate data not found at {CLIMATE_PATH}")
    return None


def _get_valid_lag_combinations(time_period: str) -> list[tuple[int, str, str]]:
    """
    Get valid (lag_years, eqi_code, eqi_period) combinations for a given AAMR time period.

    Returns: list of (lag_years, eqi_code, eqi_period)
    """
    combinations: list[tuple[int, str, str]] = []

    if time_period == "2006-2010":
        combinations.append((5, "0005", "2000-2005"))
    elif time_period == "2011-2015":
        combinations.append((10, "0005", "2000-2005"))
        combinations.append((5, "0610", "2006-2010"))
    elif time_period == "2016-2020":
        combinations.append((15, "0005", "2000-2005"))
        combinations.append((10, "0610", "2006-2010"))

    return combinations


def _process_aamr_file(
    file_path: Path,
    year_range: str,
    icd_raw: str,
    stratum_raw: str,
    eqi_dict: dict,
    smoking_df: pd.DataFrame | None,
    cluster_df: pd.DataFrame | None,
    climate_df: pd.DataFrame | None,
) -> list[pd.DataFrame]:
    """
    Process one stratified AAMR file and create rows for all valid lag combinations.

    Returns: list of DataFrames (one per valid lag combination)
    """
    print(f"  Processing: {file_path.name}")

    # Load AAMR data
    df = pd.read_csv(file_path, dtype={"COUNTY_FIPS": str})
    if "COUNTY_FIPS" not in df.columns:
        # Compatibility with earlier schema
        if "COUNTY" in df.columns:
            df = df.rename(columns={"COUNTY": "COUNTY_FIPS"})
        else:
            raise ValueError("Input AAMR file missing COUNTY_FIPS column")

    # Ensure COUNTY_FIPS is 5-digit string
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Prepare outcome + stratum label
    outcome = _icd_to_outcome(icd_raw)
    stratum_norm = _normalize_stratum(stratum_raw)

    if not outcome:
        print(f"    ⚠️ Could not determine outcome from ICD: {icd_raw}, skipping")
        return []

    cancer_type = f"{outcome}_{stratum_norm}"

    # Get valid lag combinations for this time period
    valid_combinations = _get_valid_lag_combinations(year_range)
    if not valid_combinations:
        print(f"    ⚠️ No valid EQI-AAMR lag combinations for {year_range}")
        return []

    # Ensure required AAMR columns exist
    # Expected columns from stratified AAMR: Deaths, Population, AAMR, AAMR_SE, AAMR_Lower, AAMR_Upper
    for col in ["Deaths", "Population", "AAMR", "AAMR_SE", "AAMR_Lower", "AAMR_Upper"]:
        if col not in df.columns:
            raise ValueError(f"Missing column in {file_path.name}: {col}")

    result_dfs: list[pd.DataFrame] = []

    for lag, eqi_code, eqi_period in valid_combinations:
        # Create base DataFrame
        out = pd.DataFrame(
            {
                "COUNTY_FIPS": df["COUNTY_FIPS"].astype(str),
                "EQI_Period": eqi_period,
                "Time_Period": year_range,
                "Lag_Years": lag,
                "Cancer_Type": cancer_type,
                "Deaths": df["Deaths"],
                "Population": df["Population"],
                "AAMR": df["AAMR"],
                "AAMR_SE": df["AAMR_SE"],
                "AAMR_Lower": df["AAMR_Lower"],
                "AAMR_Upper": df["AAMR_Upper"],
            }
        )

        # Lowercase copies removed; use AAMR_Lower / AAMR_Upper only

        # Merge EQI data
        if eqi_code and eqi_code in eqi_dict:
            eqidf = eqi_dict[eqi_code][["COUNTY_FIPS"] + EQI_COLS].copy()
            out = out.merge(eqidf, on="COUNTY_FIPS", how="left")
        else:
            print(f"    ⚠️ EQI data not available for code {eqi_code}")
            for c in EQI_COLS:
                out[c] = pd.NA

        # Merge Smoking data
        if smoking_df is not None:
            out = out.merge(smoking_df, on="COUNTY_FIPS", how="left")
        else:
            out["Smoking_Rate"] = pd.NA

        # Merge Cluster data
        if cluster_df is not None:
            out = out.merge(cluster_df, on="COUNTY_FIPS", how="left")
        else:
            for c in CLUSTER_COLS:
                out[c] = pd.NA

        # Merge Climate data
        if climate_df is not None:
            out = out.merge(climate_df, on="COUNTY_FIPS", how="left")
        else:
            for c in CLIMATE_COLS:
                out[c] = pd.NA

        result_dfs.append(out)

    print(
        f"    ✓ Created {len(result_dfs)} lag combination(s) with {len(df)} counties each for {cancer_type}"
    )

    return result_dfs


# ---------------- Main ----------------


def main():
    print("=" * 70)
    print("EQI × CDC Triangulation (Stratified) — AAMR Long Table Builder")
    print("=" * 70)

    # Check input directory
    if not TRIANGULATION_AAMR_DIR.exists():
        print(f"\n⚠️ Input directory not found: {TRIANGULATION_AAMR_DIR}")
        print(
            "Please run EQI_Triangulation_AAMR_Stratified.py first to generate stratified AAMR data."
        )
        sys.exit(1)

    # Get all stratified AAMR files
    aamr_files = sorted(TRIANGULATION_AAMR_DIR.glob("*.csv"))
    if not aamr_files:
        print(f"\n⚠️ No stratified AAMR files found in {TRIANGULATION_AAMR_DIR}")
        print(
            "Please run EQI_Triangulation_AAMR_Stratified.py first to generate stratified AAMR data."
        )
        sys.exit(1)

    print(f"\nFound {len(aamr_files)} stratified AAMR files")

    # Load EQI, Smoking, Cluster, and Climate data
    print("\n📊 Loading EQI, Smoking, Cluster, and Climate data...")
    eqi_dict = _load_eqi_by_period()
    smoking_df = _load_smoking()
    cluster_df = _load_clusters()
    climate_df = _load_climate()

    if not eqi_dict:
        print("⚠️ No EQI data loaded. Continuing without EQI covariates.")

    if smoking_df is None:
        print("⚠️ No smoking data loaded. Continuing without smoking rates.")

    if cluster_df is None:
        print("⚠️ No cluster data loaded. Continuing without cluster IDs.")

    if climate_df is None:
        print("⚠️ No climate data loaded. Continuing without climate zones.")

    # Process all AAMR files
    print("\n" + "=" * 70)
    print("Processing Stratified AAMR Files and Merging with EQI")
    print("=" * 70 + "\n")

    all_rows: list[pd.DataFrame] = []

    for file_path in aamr_files:
        year_range, icd_raw, stratum_raw = _parse_stratified_aamr_filename(
            file_path.name
        )

        if not year_range or not stratum_raw or not icd_raw:
            print(f"  ⚠️ Could not parse filename: {file_path.name}, skipping")
            continue

        # Process file and get all lag combinations
        try:
            result_dfs = _process_aamr_file(
                file_path,
                year_range,
                icd_raw,
                stratum_raw,
                eqi_dict,
                smoking_df,
                cluster_df,
                climate_df,
            )
            all_rows.extend(result_dfs)
        except Exception as e:
            print(f"  ❌ Error processing {file_path.name}: {e}")

    # Combine all results
    if not all_rows:
        print("\n⚠️ No rows to write.")
        # Create empty skeleton
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        empty_df = pd.DataFrame(
            columns=[
                "COUNTY_FIPS",
                "EQI_Period",
                "Time_Period",
                "Lag_Years",
                "Cancer_Type",
                "Deaths",
                "Population",
                "AAMR",
                "AAMR_SE",
                "AAMR_Lower",
                "AAMR_Upper",
                "Smoking_Rate",
            ]
            + EQI_COLS
            + CLUSTER_COLS
            + CLIMATE_COLS
        )
        empty_df.to_csv(OUTPUT_FILE, index=False)
        print(f"💾 Wrote empty skeleton to {OUTPUT_FILE}")
        return

    # Concatenate all DataFrames
    final = pd.concat(all_rows, ignore_index=True)

    # Order columns
    first_cols = [
        "COUNTY_FIPS",
        "EQI_Period",
        "Time_Period",
        "Lag_Years",
        "Cancer_Type",
        "Deaths",
        "Population",
        "AAMR",
        "AAMR_SE",
        "AAMR_Lower",
        "AAMR_Upper",
        "Smoking_Rate",
    ]
    ordered = (
        first_cols
        + [c for c in EQI_COLS if c in final.columns]
        + [c for c in CLUSTER_COLS if c in final.columns]
        + [c for c in CLIMATE_COLS if c in final.columns]
    )
    final = final[ordered]

    # Cast RUCC/EQI quintiles to nullable int to avoid 1.0/2.0 formatting
    for c in EQI_COLS:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors="coerce").astype("Int64")

    # Cast cluster columns to nullable int
    for c in CLUSTER_COLS:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors="coerce").astype("Int64")

    # Cast numeric climate columns to nullable int
    numeric_climate_cols = ["census_region", "doe_major"]
    for c in numeric_climate_cols:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors="coerce").astype("Int64")

    # Ensure integer types for Deaths and Population
    if "Deaths" in final.columns:
        final["Deaths"] = (
            pd.to_numeric(final["Deaths"], errors="coerce").round(0).astype("Int64")
        )
    if "Population" in final.columns:
        final["Population"] = (
            pd.to_numeric(final["Population"], errors="coerce").round(0).astype("Int64")
        )

    # Impute Smoking_Rate missing values with mean
    if "Smoking_Rate" in final.columns:
        final["Smoking_Rate"] = pd.to_numeric(final["Smoking_Rate"], errors="coerce")
        sr_mean = final["Smoking_Rate"].mean(skipna=True)
        if pd.notna(sr_mean):
            final["Smoking_Rate"] = final["Smoking_Rate"].fillna(sr_mean)

    # Sort by key columns
    final = final.sort_values(
        ["Time_Period", "Cancer_Type", "Lag_Years", "COUNTY_FIPS"]
    ).reset_index(drop=True)

    # Save output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUTPUT_FILE, index=False)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total rows:        {len(final):,}")
    print(f"Unique counties:   {final['COUNTY_FIPS'].nunique():,}")
    print(f"Time periods:      {sorted(final['Time_Period'].unique())}")
    print(f"Unique Cancer_Type: {final['Cancer_Type'].nunique()}")
    if "Lag_Years" in final.columns:
        print(f"Lag combinations:  {sorted(final['Lag_Years'].unique())}")
    if "Smoking_Rate" in final.columns:
        sr_missing = final["Smoking_Rate"].isna().sum()
        print(f"Smoking_Rate NA:   {sr_missing:,}")
    print(f"\nOutput saved to:   {OUTPUT_FILE}")

    # Show example
    print("\n" + "=" * 70)
    print("Example Output (first 10 rows)")
    print("=" * 70)
    print(final.head(10).to_string(index=False))

    print("\n✓ Stratified EQI × Triangulation AAMR merge completed successfully!")


if __name__ == "__main__":
    main()
