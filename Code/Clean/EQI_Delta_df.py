#!/usr/bin/env python3
"""
Generate df_Delta.csv from df.csv by computing period-over-period changes.

Matched pairs (same lag from EQI, consecutive EQI periods):
  Lag=5:  before=(EQI0005, AAMR2006-2010)  after=(EQI0610, AAMR2011-2015)
  Lag=10: before=(EQI0005, AAMR2011-2015)  after=(EQI0610, AAMR2016-2020)
  Lag=15: before=(EQI0005, AAMR2016-2020)  after=(EQI0610, AAMR2021-2024)

Delta AAMR (interval arithmetic for independent intervals):
  delta_lower = AAMR_Lower_after - AAMR_Upper_before
  delta_upper = AAMR_Upper_after - AAMR_Lower_before
  When both periods are exact (Lower==Upper), delta is also exact.

EQI change category convention (higher quintile = worse environment):
  Improved: quintile_after - quintile_before <= -1
  Worsened: quintile_after - quintile_before >= +1
  Stable:   quintile_after - quintile_before == 0
"""

from pathlib import Path
import pandas as pd
import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
with (PROJECT_ROOT / "config.yaml").open(encoding="utf-8") as f:
    CFG = yaml.safe_load(f)

INPUT_DF   = PROJECT_ROOT / "Data/Processed/df.csv"
OUTPUT_DF  = PROJECT_ROOT / "Data/Processed/df_Delta.csv"

SELECTED_OUTCOMES = [
    "I00_I99",
    "J40_J47_J60_J70_J84_D86_C34",
    "K70_K76_C22",
    "N00_N29_C64_C65",
    "X60_X84_Y87.0",
    "G20_G30_G12.2_F01_F03",
    "C00_C97",
]

DELTA_PAIRS = [
    {"lag": 5,
     "eqi_before": "2000-2005", "aamr_before": "2006-2010", "lag_before": 5,
     "eqi_after":  "2006-2010", "aamr_after":  "2011-2015", "lag_after":  5},
    {"lag": 10,
     "eqi_before": "2000-2005", "aamr_before": "2011-2015", "lag_before": 10,
     "eqi_after":  "2006-2010", "aamr_after":  "2016-2020", "lag_after":  10},
    {"lag": 15,
     "eqi_before": "2000-2005", "aamr_before": "2016-2020", "lag_before": 15,
     "eqi_after":  "2006-2010", "aamr_after":  "2021-2024", "lag_after":  15},
]

EQI_DOMAINS = [
    ("EQI",       "EQI_Change_Category"),
    ("EQI_Air",   "Air_Change_Category"),
    ("EQI_Water", "Water_Change_Category"),
    ("EQI_Land",  "Land_Change_Category"),
    ("EQI_Built", "Built_Change_Category"),
    ("EQI_Social","Social_Change_Category"),
]

RUCC_DOMAINS = [
    ("RUCC_EQI",       "RUCC_EQI_Change_Category"),
    ("RUCC_EQI_Air",   "RUCC_EQI_Air_Change_Category"),
    ("RUCC_EQI_Water", "RUCC_EQI_Water_Change_Category"),
    ("RUCC_EQI_Land",  "RUCC_EQI_Land_Change_Category"),
    ("RUCC_EQI_Built", "RUCC_EQI_Built_Change_Category"),
    ("RUCC_EQI_Social","RUCC_EQI_Social_Change_Category"),
]

STATIC_COLS = ["Cluster_EQI", "Cluster_NLCD", "RUCC"]


def assign_change_category(delta: pd.Series) -> pd.Series:
    """
    Assign EQI change category based on quintile delta.
    S = Stable (delta = 0)
    I = Improved (delta <= -1)
    W = Worsened (delta >= 1)
    """
    cat = pd.Series("S", index=delta.index, dtype=object)
    cat[delta <= -1] = "I"
    cat[delta >= 1]  = "W"
    cat[delta.isna()] = np.nan
    return cat


def process_pair(df: pd.DataFrame, pair: dict) -> pd.DataFrame | None:
    lag = pair["lag"]

    before = df[
        (df["EQI_Period"] == pair["eqi_before"]) &
        (df["Time_Period"] == pair["aamr_before"]) &
        (df["Lag_Years"]   == pair["lag_before"])
    ].copy()

    after = df[
        (df["EQI_Period"] == pair["eqi_after"]) &
        (df["Time_Period"] == pair["aamr_after"]) &
        (df["Lag_Years"]   == pair["lag_after"])
    ].copy()

    print(f"  Lag={lag}: before={len(before):,}  after={len(after):,}")
    if before.empty or after.empty:
        return None

    # Columns to carry from each side
    eqi_cols   = [c for c, _ in EQI_DOMAINS   if c in df.columns]
    rucc_cols  = [c for c, _ in RUCC_DOMAINS  if c in df.columns]
    all_eqi    = eqi_cols + rucc_cols
    static_src = [c for c in STATIC_COLS if c in df.columns]

    b_sel = ["COUNTY_FIPS", "Outcome", "AAMR_Lower", "AAMR_Upper",
             "Smoking_rate"] + all_eqi + static_src
    a_sel = ["COUNTY_FIPS", "Outcome", "AAMR_Lower", "AAMR_Upper",
             "Smoking_rate"] + all_eqi

    b_sel = [c for c in b_sel if c in before.columns]
    a_sel = [c for c in a_sel if c in after.columns]

    b = before[b_sel].rename(columns={
        "AAMR_Lower": "AAMR_Lower_b", "AAMR_Upper": "AAMR_Upper_b",
        "Smoking_rate": "Smoking_b",
        **{c: f"{c}_b" for c in all_eqi if c in before.columns}
    })
    a = after[a_sel].rename(columns={
        "AAMR_Lower": "AAMR_Lower_a", "AAMR_Upper": "AAMR_Upper_a",
        "Smoking_rate": "Smoking_a",
        **{c: f"{c}_a" for c in all_eqi if c in after.columns}
    })

    m = b.merge(a, on=["COUNTY_FIPS", "Outcome"], how="inner")
    print(f"         merged={len(m):,}")
    if m.empty:
        return None

    # State from FIPS
    m["State_FIPS"] = m["COUNTY_FIPS"].str[:2]

    # Delta AAMR (interval arithmetic)
    m["delta_AAMR_lower"] = (m["AAMR_Lower_a"] - m["AAMR_Upper_b"]).round(2)
    m["delta_AAMR_upper"] = (m["AAMR_Upper_a"] - m["AAMR_Lower_b"]).round(2)

    # Delta Smoking
    if "Smoking_b" in m.columns and "Smoking_a" in m.columns:
        m["delta_Smoking_Rate"] = (m["Smoking_a"] - m["Smoking_b"]).round(2)
    else:
        m["delta_Smoking_Rate"] = np.nan

    # EQI change categories
    for src_col, cat_col in EQI_DOMAINS:
        b_col, a_col = f"{src_col}_b", f"{src_col}_a"
        if b_col in m.columns and a_col in m.columns:
            m[cat_col] = assign_change_category(
                pd.to_numeric(m[a_col], errors="coerce") -
                pd.to_numeric(m[b_col], errors="coerce")
            )
        else:
            m[cat_col] = np.nan

    # RUCC_EQI change categories
    for src_col, cat_col in RUCC_DOMAINS:
        b_col, a_col = f"{src_col}_b", f"{src_col}_a"
        if b_col in m.columns and a_col in m.columns:
            m[cat_col] = assign_change_category(
                pd.to_numeric(m[a_col], errors="coerce") -
                pd.to_numeric(m[b_col], errors="coerce")
            )
        else:
            m[cat_col] = np.nan

    m["Lag"] = lag

    out_cols = (
        ["COUNTY_FIPS", "State_FIPS", "Outcome", "Lag",
         "delta_AAMR_lower", "delta_AAMR_upper", "delta_Smoking_Rate"]
        + [c for c in static_src if c in m.columns]
        + [cat for _, cat in EQI_DOMAINS]
        + [cat for _, cat in RUCC_DOMAINS if cat in m.columns]
    )
    out_cols = [c for c in out_cols if c in m.columns]

    result = m[out_cols].dropna(subset=["delta_AAMR_lower", "delta_AAMR_upper",
                                        "EQI_Change_Category"])
    print(f"         output={len(result):,} (after dropping key NAs)")
    return result


def main():
    print("=" * 70)
    print("EQI_AAMR_Delta Dataset Generator")
    print("=" * 70)

    if not INPUT_DF.exists():
        raise FileNotFoundError(f"Input not found: {INPUT_DF}")

    print(f"\nLoading {INPUT_DF}...")
    df = pd.read_csv(INPUT_DF, dtype={"COUNTY_FIPS": str})
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
    print(f"Loaded {len(df):,} rows | {df['COUNTY_FIPS'].nunique():,} counties "
          f"| {df['Outcome'].nunique()} outcomes")

    df = df[df["Outcome"].isin(SELECTED_OUTCOMES)]
    print(f"Filtered to {len(SELECTED_OUTCOMES)} selected outcomes: {len(df):,} rows")

    print("\nProcessing delta pairs...")
    parts = []
    for pair in DELTA_PAIRS:
        result = process_pair(df, pair)
        if result is not None:
            parts.append(result)

    if not parts:
        print("\nNo rows produced.")
        return

    out = pd.concat(parts, ignore_index=True)
    out = out.sort_values(["Outcome", "Lag", "COUNTY_FIPS"]).reset_index(drop=True)

    OUTPUT_DF.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_DF, index=False)

    print("\n" + "=" * 70)
    print(f"Output: {len(out):,} rows → {OUTPUT_DF}")
    print(f"Columns ({len(out.columns)}): {list(out.columns)}")
    print(f"Counties:  {out['COUNTY_FIPS'].nunique():,}")
    print(f"Outcomes:  {out['Outcome'].nunique()}")
    print(f"Lags:      {sorted(out['Lag'].unique())}")
    print("\nEQI_Change_Category distribution:")
    print(out["EQI_Change_Category"].value_counts().to_string())
    print("\nDone.")


if __name__ == "__main__":
    main()
