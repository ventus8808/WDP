#!/usr/bin/env python3
"""
CDC SVI × AAMR long table builder

df_SVI.csv: AAMR + SVI trajectory exposure + covariates (period-matched)

The exposure here is the GBTM SVI trajectory category (A/B/C/D), produced by
Code/Analysis/CDC_SVI_GBTM.py and extracted by Code/Clean/CDC_SVI_Extract.py
into Data/Processed/SVI/SVI_Result.csv. Unlike EQI it is a single, static
classification of each county's 2000-2022 trajectory, so there is no EQI-style
lag dimension: every outcome period for a county carries the same SVI value.

Covariate period-matching rules:
  Static    (no period): Census_Region, Climate_Zone, Cluster_EQI, Cluster_NLCD
  Baseline-era          : Smoking_rate, Heavy_Drinking_rate, Physical_Activities_rate,
                          Obesity_rate, Diabetes_Prevalence_rate,
                          Physician_Density_per100k, Forest_Coverage,
                          AQS_Number, Economic_type
                          (measured pre-outcome; use 2006-2010, fall back to 2000-2005)
  Outcome-period matched: Homeownership_rate/tertile, Uninsured_rate
"""

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
with (PROJECT_ROOT / "config.yaml").open(encoding="utf-8") as f:
    CFG = yaml.safe_load(f)

AAMR_DIR  = PROJECT_ROOT / "Data/Original/CDC Triangulation/AAMR"
OUTPUT_DIR = PROJECT_ROOT / "Data/Processed"
OUTPUT_DF  = OUTPUT_DIR / "df_SVI.csv"

SVI_PATH      = PROJECT_ROOT / "Data/Processed/SVI/SVI_Result.csv"
COVARIATE_DIR = PROJECT_ROOT / "Data/Processed/Covariate"
IHME_DIR      = PROJECT_ROOT / "Data/Processed/IHME"
CLUSTER_DIR   = PROJECT_ROOT / CFG["data_directories"]["processed"] / "Cluster"

# Outcome-period suffixes (cover all four AAMR periods)
OUTCOME_SUFFIX = {"2006-2010": "0610", "2011-2015": "1115",
                  "2016-2020": "1620", "2021-2024": "2124"}
# Baseline-era confounders: prefer 2006-2010, fall back to 2000-2005
BASELINE_SUFFIXES = ("0610", "0005")

INT_COLS = ["Deaths", "Population", "Census_Region",
            "Economic_type", "Homeownership_tertile"]


# ─── helpers ──────────────────────────────────────────────────────────────────

def _fips5(s: pd.Series) -> pd.Series:
    return s.astype(str).str.zfill(5)


def _load_csv(path, label="") -> pd.DataFrame | None:
    p = Path(path)
    if p.exists():
        df = pd.read_csv(p, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = _fips5(df["COUNTY_FIPS"])
        return df
    print(f"  Warning: {label or p.name} not found")
    return None


def _pick_period_col(wide_df, col_prefix, suffix):
    """Two-column df (COUNTY_FIPS, col_prefix) for the given suffix, or None."""
    if wide_df is None:
        return None
    col = f"{col_prefix}_{suffix}"
    if col in wide_df.columns:
        return wide_df[["COUNTY_FIPS", col]].rename(columns={col: col_prefix})
    return None


def _pick_baseline_col(wide_df, col_prefix):
    """Baseline-era value: first available of BASELINE_SUFFIXES."""
    for suf in BASELINE_SUFFIXES:
        tmp = _pick_period_col(wide_df, col_prefix, suf)
        if tmp is not None:
            return tmp
    return None


def _cast_int(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


def _parse_filename(filename):
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)
    if len(parts) == 2 and re.match(r"\d{4}-\d{4}", parts[0]):
        return parts[0], parts[1]
    return None, None


# ─── loaders ──────────────────────────────────────────────────────────────────

def _load_static():
    pieces = []
    for fname, col in [("Census_Region.csv", "Census_Region"),
                       ("Climate_Zone.csv", "Climate_Zone")]:
        d = _load_csv(COVARIATE_DIR / fname, col)
        if d is not None and col in d.columns:
            pieces.append(d[["COUNTY_FIPS", col]])
    for fname, col in [("Cluster_EQI.csv", "Cluster_EQI"),
                       ("Cluster_NLCD.csv", "Cluster_NLCD")]:
        d = _load_csv(CLUSTER_DIR / fname, col)
        if d is not None and col in d.columns:
            pieces.append(d[["COUNTY_FIPS", col]])
    if not pieces:
        return None
    static_df = pieces[0]
    for p in pieces[1:]:
        static_df = static_df.merge(p, on="COUNTY_FIPS", how="outer")
    return static_df


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("CDC SVI × AAMR — Long Table Builder")
    print("=" * 70)

    if not AAMR_DIR.exists():
        print(f"\nError: {AAMR_DIR} not found")
        sys.exit(1)
    aamr_files = sorted(AAMR_DIR.glob("*.csv"))
    if not aamr_files:
        print(f"\nError: no AAMR files in {AAMR_DIR}")
        sys.exit(1)
    print(f"\nFound {len(aamr_files)} AAMR files")

    print("\nLoading data...")
    svi = _load_csv(SVI_PATH, "SVI_Result")
    if svi is None or "SVI" not in svi.columns:
        print(f"\nError: SVI exposure not found at {SVI_PATH} "
              f"(run CDC_SVI_Extract.py first)")
        sys.exit(1)
    svi = svi[["COUNTY_FIPS", "SVI"]]

    static_df  = _load_static()
    homeown_df = _load_csv(COVARIATE_DIR / "Homeownership_rate.csv", "Homeownership")
    uninsured_df = _load_csv(COVARIATE_DIR / "Uninsured_rate.csv", "Uninsured_rate")
    baseline_cov = {
        "Smoking_rate":              _load_csv(IHME_DIR / "IHME_Smoking.csv", "Smoking"),
        "Heavy_Drinking_rate":       _load_csv(IHME_DIR / "IHME_Drinking.csv", "Drinking"),
        "Physical_Activities_rate":  _load_csv(IHME_DIR / "IHME_Physical_Activities.csv", "Physical_Activities"),
        "Obesity_rate":              _load_csv(IHME_DIR / "IHME_Obesity.csv", "Obesity"),
        "Diabetes_Prevalence_rate":  _load_csv(IHME_DIR / "IHME_Diabetes_Prevalence.csv", "Diabetes_Prevalence"),
        "Physician_Density_per100k": _load_csv(COVARIATE_DIR / "Physician_Density_per100k.csv", "Physician"),
        "Forest_Coverage":           _load_csv(COVARIATE_DIR / "Forest_Coverage.csv", "Forest"),
        "AQS_Number":                _load_csv(COVARIATE_DIR / "AQS_Number.csv", "AQS_Number"),
        "Economic_type":             _load_csv(COVARIATE_DIR / "Economic_type.csv", "Economic_type"),
    }

    print("\nProcessing files...")
    all_dfs = []
    for fp in aamr_files:
        year_range, icd_code = _parse_filename(fp.name)
        if not year_range or year_range not in OUTCOME_SUFFIX:
            print(f"  Skipping {fp.name} (period not recognised)")
            continue
        df = pd.read_csv(fp, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = _fips5(df["COUNTY_FIPS"])

        row = pd.DataFrame({
            "COUNTY_FIPS": df["COUNTY_FIPS"].astype(str),
            "Time_Period": year_range,
            "Outcome":     icd_code.replace("-", "_"),
            "Deaths":      df.get("Deaths"),
            "Population":  df.get("Population"),
            "AAMR":        df.get("AAMR"),
            "AAMR_Lower":  df.get("AAMR_Lower"),
            "AAMR_Upper":  df.get("AAMR_Upper"),
        })

        # Exposure (static) + static covariates
        row = row.merge(svi, on="COUNTY_FIPS", how="left")
        if static_df is not None:
            row = row.merge(static_df, on="COUNTY_FIPS", how="left")

        # Baseline-era covariates
        for prefix, wide in baseline_cov.items():
            tmp = _pick_baseline_col(wide, prefix)
            if tmp is not None:
                row = row.merge(tmp, on="COUNTY_FIPS", how="left")

        # Outcome-period covariates
        suf = OUTCOME_SUFFIX[year_range]
        if homeown_df is not None:
            rate_col, tert_col = f"Homeownership_rate_{suf}", f"Homeownership_tertile_{suf}"
            hw_cols = [c for c in (rate_col, tert_col) if c in homeown_df.columns]
            if hw_cols:
                row = row.merge(
                    homeown_df[["COUNTY_FIPS"] + hw_cols].rename(columns={
                        rate_col: "Homeownership_rate",
                        tert_col: "Homeownership_tertile"}),
                    on="COUNTY_FIPS", how="left")
        tmp = _pick_period_col(uninsured_df, "Uninsured_rate", suf)
        if tmp is not None:
            row = row.merge(tmp, on="COUNTY_FIPS", how="left")

        all_dfs.append(row)
        print(f"  {fp.name}: {len(df)} counties")

    if not all_dfs:
        print("\nNo rows produced.")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out = pd.concat(all_dfs, ignore_index=True)
    out = _cast_int(out, INT_COLS)
    out = out.sort_values(["Time_Period", "Outcome", "COUNTY_FIPS"]).reset_index(drop=True)
    out.to_csv(OUTPUT_DF, index=False)

    print("\n" + "=" * 70)
    print(f"Output:   {len(out):,} rows  -> {OUTPUT_DF}")
    print(f"Columns:  {list(out.columns)}")
    print(f"Counties: {out['COUNTY_FIPS'].nunique():,}")
    print(f"Periods:  {sorted(out['Time_Period'].unique())}")
    print(f"Outcomes: {out['Outcome'].nunique()}")
    print(f"SVI dist: {out.drop_duplicates('COUNTY_FIPS')['SVI'].value_counts(dropna=False).sort_index().to_dict()}")
    print("\nDone.")


if __name__ == "__main__":
    main()
