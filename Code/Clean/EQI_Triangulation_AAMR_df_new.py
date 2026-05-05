#!/usr/bin/env python3
"""
EQI × CDC Triangulation — AAMR long table builder

df_Main.csv:        AAMR + EQI + RUCC-stratified EQI  (no other covariates)
df_Sensitivity.csv: AAMR + EQI (no RUCC_EQI) + all covariates (period-matched)

Covariate period-matching rules:
  Static   (no period): Census_Region, Climate_Zone (koppen_major, doe_zone)
  EQI-period matched:   Smoking_rate, Physician_Density_per100k, Forest_Coverage,
                        AQS_Number, Economic_type
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

TRIANGULATION_AAMR_DIR = PROJECT_ROOT / "Data/Original/CDC Triangulation/AAMR"
OUTPUT_DIR  = PROJECT_ROOT / "Data/Processed/df"
OUTPUT_MAIN = OUTPUT_DIR / "df_Main.csv"
OUTPUT_SENS = OUTPUT_DIR / "df_Sensitivity.csv"

EQI_DIR      = PROJECT_ROOT / CFG["data_sources"]["epa_eqi"]["processed"]
COVARIATE_DIR = PROJECT_ROOT / "Data/Processed/Covariate"
IHME_DIR     = PROJECT_ROOT / "Data/Processed/IHME"

EQI_COLS      = ["RUCC", "EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]
RUCC_EQI_COLS = ["RUCC_EQI", "RUCC_EQI_Air", "RUCC_EQI_Water", "RUCC_EQI_Land", "RUCC_EQI_Built", "RUCC_EQI_Social"]

# EQI_Period label → column suffix in covariate files
EQI_SUFFIX = {"2000-2005": "0005", "2006-2010": "0610"}

# Time_Period (outcome) → column suffix for outcome-period covariates
HOMEOWN_SUFFIX   = {"2006-2010": "0610", "2011-2015": "1115", "2016-2020": "1620", "2021-2024": "2124"}
UNINSURED_SUFFIX = {"2006-2010": "0610", "2011-2015": "1115", "2016-2020": "1620", "2021-2024": "2124"}

EQI_COL_RENAMES = {
    "EQI_air": "EQI_Air", "EQI_water": "EQI_Water",
    "EQI_land": "EQI_Land", "EQI_built": "EQI_Built",
    "EQI_Sociodemographic": "EQI_Social",
    "RUCC_EQI_air": "RUCC_EQI_Air", "RUCC_EQI_water": "RUCC_EQI_Water",
    "RUCC_EQI_land": "RUCC_EQI_Land", "RUCC_EQI_built": "RUCC_EQI_Built",
    "RUCC_EQI_Sociodemographic": "RUCC_EQI_Social",
}


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


def _pick_period_col(wide_df: pd.DataFrame | None, col_prefix: str, suffix: str) -> pd.DataFrame | None:
    """Return a two-column df (COUNTY_FIPS, col_prefix) for the given suffix."""
    if wide_df is None:
        return None
    col = f"{col_prefix}_{suffix}"
    if col in wide_df.columns:
        return wide_df[["COUNTY_FIPS", col]].rename(columns={col: col_prefix})
    return None


def _cast_int(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


# ─── loaders ──────────────────────────────────────────────────────────────────

def _load_eqi() -> dict:
    d = {}
    for code in ("0005", "0610"):
        fp = EQI_DIR / f"EQI{code}.csv"
        if fp.exists():
            df = pd.read_csv(fp)
            df["COUNTY_FIPS"] = _fips5(df["COUNTY_FIPS"])
            renames = {k: v for k, v in EQI_COL_RENAMES.items()
                       if k in df.columns and v not in df.columns}
            if renames:
                df = df.rename(columns=renames)
            d[code] = df
        else:
            print(f"  Warning: EQI{code}.csv not found")
    return d


def _load_static() -> pd.DataFrame | None:
    """Merge Census_Region + Climate_Zone into one static lookup."""
    census  = _load_csv(COVARIATE_DIR / "Census_Region.csv",  "Census_Region")
    climate = _load_csv(COVARIATE_DIR / "Climate_Zone.csv",   "Climate_Zone")
    if census is None and climate is None:
        return None
    if census is None:
        return climate[["COUNTY_FIPS", "koppen_major", "doe_zone"]]
    if climate is None:
        return census[["COUNTY_FIPS", "Census_Region"]]
    return census[["COUNTY_FIPS", "Census_Region"]].merge(
        climate[["COUNTY_FIPS", "koppen_major", "doe_zone"]],
        on="COUNTY_FIPS", how="outer"
    )


# ─── filename parser ───────────────────────────────────────────────────────────

def _parse_filename(filename: str) -> tuple[str | None, str | None]:
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)
    if len(parts) == 2 and re.match(r"\d{4}-\d{4}", parts[0]):
        return parts[0], parts[1]
    return None, None


# ─── lag combinations ──────────────────────────────────────────────────────────

def _get_valid_lag_combinations(time_period: str) -> list[tuple[int, str, str]]:
    combos: list[tuple[int, str, str]] = []
    if time_period == "2006-2010":
        combos.append((5,  "0005", "2000-2005"))
    elif time_period == "2011-2015":
        combos += [(10, "0005", "2000-2005"), (5,  "0610", "2006-2010")]
    elif time_period == "2016-2020":
        combos += [(15, "0005", "2000-2005"), (10, "0610", "2006-2010")]
    elif time_period == "2021-2024":
        combos += [(20, "0005", "2000-2005"), (15, "0610", "2006-2010")]
    return combos


# ─── file processor ────────────────────────────────────────────────────────────

def _process_aamr_file(
    file_path: Path,
    year_range: str,
    icd_code: str,
    eqi_dict: dict,
    static_df: pd.DataFrame | None,
    homeown_df: pd.DataFrame | None,
    uninsured_df: pd.DataFrame | None,
    eqi_cov_dfs: dict[str, pd.DataFrame | None],
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """Returns (main_dfs, sens_dfs) — one DataFrame per valid lag combo each."""
    df = pd.read_csv(file_path, dtype={"COUNTY_FIPS": str})
    df["COUNTY_FIPS"] = _fips5(df["COUNTY_FIPS"])
    icd_fmt = icd_code.replace("-", "_")

    valid_combinations = _get_valid_lag_combinations(year_range)
    if not valid_combinations:
        return [], []

    main_dfs: list[pd.DataFrame] = []
    sens_dfs: list[pd.DataFrame] = []

    for lag, eqi_code, eqi_period in valid_combinations:
        # Base columns
        base = pd.DataFrame({
            "COUNTY_FIPS": df["COUNTY_FIPS"].astype(str),
            "EQI_Period":  eqi_period,
            "Time_Period": year_range,
            "Lag_Years":   lag,
            "Outcome":     icd_fmt,
            "Deaths":      df.get("Deaths"),
            "Population":  df.get("Population"),
            "AAMR":        df.get("AAMR"),
            "AAMR_Lower":  df.get("AAMR_Lower"),
            "AAMR_Upper":  df.get("AAMR_Upper"),
        })

        # Merge EQI (all columns: EQI + RUCC_EQI)
        if eqi_code in eqi_dict:
            eqi_src = eqi_dict[eqi_code]
            all_eqi = [c for c in EQI_COLS + RUCC_EQI_COLS if c in eqi_src.columns]
            base = base.merge(eqi_src[["COUNTY_FIPS"] + all_eqi], on="COUNTY_FIPS", how="left")
        else:
            for c in EQI_COLS + RUCC_EQI_COLS:
                base[c] = pd.NA

        # ── Main: base + EQI + RUCC_EQI ──────────────────────────────────────
        main_dfs.append(base.copy())

        # ── Sensitivity: EQI only (no RUCC_EQI) + all covariates ─────────────
        sens = base.drop(columns=[c for c in RUCC_EQI_COLS if c in base.columns]).copy()

        # Static: Census_Region, koppen_major, doe_zone
        if static_df is not None:
            sens = sens.merge(static_df, on="COUNTY_FIPS", how="left")

        # EQI-period covariates
        eqi_suf = EQI_SUFFIX.get(eqi_period, eqi_code)
        for col_prefix, wide_df in eqi_cov_dfs.items():
            tmp = _pick_period_col(wide_df, col_prefix, eqi_suf)
            if tmp is not None:
                sens = sens.merge(tmp, on="COUNTY_FIPS", how="left")

        # Outcome-period covariates: Homeownership
        hw_suf = HOMEOWN_SUFFIX.get(year_range)
        if homeown_df is not None and hw_suf:
            rate_col = f"Homeownership_rate_{hw_suf}"
            tert_col = f"Homeownership_tertile_{hw_suf}"
            hw_cols  = [c for c in [rate_col, tert_col] if c in homeown_df.columns]
            if hw_cols:
                tmp = homeown_df[["COUNTY_FIPS"] + hw_cols].rename(columns={
                    rate_col: "Homeownership_rate",
                    tert_col: "Homeownership_tertile",
                })
                sens = sens.merge(tmp, on="COUNTY_FIPS", how="left")

        # Outcome-period covariates: Uninsured_rate
        ui_suf = UNINSURED_SUFFIX.get(year_range)
        if uninsured_df is not None and ui_suf:
            tmp = _pick_period_col(uninsured_df, "Uninsured_rate", ui_suf)
            if tmp is not None:
                sens = sens.merge(tmp, on="COUNTY_FIPS", how="left")

        sens_dfs.append(sens)

    print(f"  {file_path.name}: {len(valid_combinations)} combo(s), {len(df)} counties")
    return main_dfs, sens_dfs


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("EQI × Triangulation AAMR — Long Table Builder")
    print("=" * 70)

    if not TRIANGULATION_AAMR_DIR.exists():
        print(f"\nError: {TRIANGULATION_AAMR_DIR} not found")
        sys.exit(1)

    aamr_files = sorted(TRIANGULATION_AAMR_DIR.glob("*.csv"))
    if not aamr_files:
        print(f"\nError: No AAMR files found in {TRIANGULATION_AAMR_DIR}")
        sys.exit(1)

    print(f"\nFound {len(aamr_files)} AAMR files")

    # Load all data upfront
    print("\nLoading data...")
    eqi_dict     = _load_eqi()
    static_df    = _load_static()
    homeown_df   = _load_csv(COVARIATE_DIR / "Homeownership_rate.csv",       "Homeownership")
    uninsured_df = _load_csv(COVARIATE_DIR / "Uninsured_rate.csv",            "Uninsured_rate")
    eqi_cov_dfs  = {
        "Smoking_rate":              _load_csv(IHME_DIR      / "IHME_Smoking.csv",              "Smoking"),
        "Heavy_Drinking_rate":       _load_csv(IHME_DIR      / "IHME_Drinking.csv",             "Drinking"),
        "Physical_Activities_rate":  _load_csv(IHME_DIR      / "IHME_Physical_Activities.csv",  "Physical_Activities"),
        "Obesity_rate":              _load_csv(IHME_DIR      / "IHME_Obesity.csv",              "Obesity"),
        "Diabetes_Prevalence_rate":  _load_csv(IHME_DIR      / "IHME_Diabetes_Prevalence.csv",  "Diabetes_Prevalence"),
        "Physician_Density_per100k": _load_csv(COVARIATE_DIR / "Physician_Density_per100k.csv", "Physician"),
        "Forest_Coverage":           _load_csv(COVARIATE_DIR / "Forest_Coverage.csv",           "Forest"),
        "AQS_Number":                _load_csv(COVARIATE_DIR / "AQS_Number.csv",                "AQS_Number"),
        "Economic_type":             _load_csv(COVARIATE_DIR / "Economic_type.csv",             "Economic_type"),
    }

    # Process all AAMR files
    print("\nProcessing files...")
    all_main: list[pd.DataFrame] = []
    all_sens: list[pd.DataFrame] = []

    for fp in aamr_files:
        year_range, icd_code = _parse_filename(fp.name)
        if not year_range:
            print(f"  Skipping {fp.name} (unrecognised filename)")
            continue
        m, s = _process_aamr_file(fp, year_range, icd_code, eqi_dict,
                                   static_df, homeown_df, uninsured_df, eqi_cov_dfs)
        all_main.extend(m)
        all_sens.extend(s)

    if not all_main:
        print("\nNo rows produced.")
        return

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Write main ────────────────────────────────────────────────────────────
    main_df = pd.concat(all_main, ignore_index=True)
    main_df = _cast_int(main_df, EQI_COLS + RUCC_EQI_COLS + ["Deaths", "Population"])
    main_df = main_df.sort_values(["Time_Period", "Outcome", "Lag_Years", "COUNTY_FIPS"]).reset_index(drop=True)
    main_df.to_csv(OUTPUT_MAIN, index=False)

    # ── Write sensitivity ─────────────────────────────────────────────────────
    sens_df = pd.concat(all_sens, ignore_index=True)
    sens_df = _cast_int(sens_df, EQI_COLS + ["Deaths", "Population",
               "Census_Region", "doe_zone", "Economic_type", "Homeownership_tertile"])
    sens_df = sens_df.sort_values(["Time_Period", "Outcome", "Lag_Years", "COUNTY_FIPS"]).reset_index(drop=True)
    sens_df.to_csv(OUTPUT_SENS, index=False)

    print("\n" + "=" * 70)
    print(f"Main:        {len(main_df):,} rows  → {OUTPUT_MAIN}")
    print(f"             cols: {list(main_df.columns)}")
    print(f"Sensitivity: {len(sens_df):,} rows  → {OUTPUT_SENS}")
    print(f"             cols: {list(sens_df.columns)}")
    print(f"\nCounties:    {main_df['COUNTY_FIPS'].nunique():,}")
    print(f"Periods:     {sorted(main_df['Time_Period'].unique())}")
    print(f"Outcomes:    {main_df['Outcome'].nunique()}")
    print(f"Lags:        {sorted(main_df['Lag_Years'].unique())}")
    print("\nDone.")


if __name__ == "__main__":
    main()
