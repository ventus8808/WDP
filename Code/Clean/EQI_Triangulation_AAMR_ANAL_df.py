#!/usr/bin/env python3
"""
ANAL x CDC Triangulation AAMR long table builder

Output:
  Data/Processed/df_ANAL.csv

Design:
  - Aggregate the annual nighttime-light county panel to the exact exposure
    windows in Code/brms_ANAL/ANAL_Research_Plan.md.
  - Merge each exposure window to the matched AAMR outcome window.
  - Keep period-specific ANAL quintiles as the primary exposure, plus tertiles,
    global categories, continuous log-z scores, alternate light metrics, SVI,
    contextual strata, and existing model covariates.

No new data are downloaded.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
with (PROJECT_ROOT / "config.yaml").open(encoding="utf-8") as f:
    CFG = yaml.safe_load(f)

AAMR_DIR = PROJECT_ROOT / "Data/Original/CDC Triangulation/AAMR"
ANAL_PATH = PROJECT_ROOT / "Data/Original/ANAL/ntl_county_panel_2000_2025.csv"
OUTPUT_DIR = PROJECT_ROOT / "Data/Processed"
OUTPUT_DF = OUTPUT_DIR / "df_ANAL.csv"

COVARIATE_DIR = PROJECT_ROOT / "Data/Processed/Covariate"
IHME_DIR = PROJECT_ROOT / "Data/Processed/IHME"
CLUSTER_DIR = PROJECT_ROOT / CFG["data_directories"]["processed"] / "Cluster"
EQI_DIR = PROJECT_ROOT / CFG["data_sources"]["epa_eqi"]["processed"]
SVI_RESULT_PATH = PROJECT_ROOT / "Data/Processed/SVI/SVI_Result.csv"
AIR_SVI_PATH = PROJECT_ROOT / "Data/Processed/SVI/Air_SVI.csv"

ANAL_METRICS = ["popw_mean_rad", "mean_rad", "sol", "lit_area_km2"]

LAG_WINDOWS = [
    {"Lag_Years": 5, "ANAL_Start": 2001, "ANAL_End": 2005, "AAMR_Period": "2006-2010", "Role": "Lag sensitivity"},
    {"Lag_Years": 5, "ANAL_Start": 2006, "ANAL_End": 2010, "AAMR_Period": "2011-2015", "Role": "Lag sensitivity"},
    {"Lag_Years": 5, "ANAL_Start": 2011, "ANAL_End": 2015, "AAMR_Period": "2016-2020", "Role": "Lag sensitivity"},
    {"Lag_Years": 5, "ANAL_Start": 2016, "ANAL_End": 2019, "AAMR_Period": "2021-2024", "Role": "Lag sensitivity"},
    {"Lag_Years": 10, "ANAL_Start": 2001, "ANAL_End": 2005, "AAMR_Period": "2011-2015", "Role": "Primary"},
    {"Lag_Years": 10, "ANAL_Start": 2006, "ANAL_End": 2010, "AAMR_Period": "2016-2020", "Role": "Primary"},
    {"Lag_Years": 10, "ANAL_Start": 2011, "ANAL_End": 2014, "AAMR_Period": "2021-2024", "Role": "Primary, pandemic sensitivity"},
    {"Lag_Years": 15, "ANAL_Start": 2001, "ANAL_End": 2005, "AAMR_Period": "2016-2020", "Role": "Lag sensitivity"},
    {"Lag_Years": 15, "ANAL_Start": 2006, "ANAL_End": 2009, "AAMR_Period": "2021-2024", "Role": "Lag sensitivity"},
]

OUTCOME_SUFFIX = {
    "2006-2010": "0610",
    "2011-2015": "1115",
    "2016-2020": "1620",
    "2021-2024": "2124",
}

COVARIATE_PERIOD_BY_ANAL = {
    "2001-2005": ("0005", "2000-2005"),
    "2006-2010": ("0610", "2006-2010"),
    "2011-2015": ("0610", "2006-2010"),
    "2016-2019": ("0610", "2006-2010"),
    "2011-2014": ("0610", "2006-2010"),
    "2006-2009": ("0610", "2006-2010"),
}

EQI_COL_RENAMES = {
    "EQI_air": "EQI_Air",
    "EQI_water": "EQI_Water",
    "EQI_land": "EQI_Land",
    "EQI_built": "EQI_Built",
    "EQI_Sociodemographic": "EQI_Social",
}

INT_COLS = [
    "Deaths",
    "Population",
    "State_FIPS",
    "ANAL_Quintile",
    "ANAL_Tertile",
    "ANAL_Global_Quintile",
    "ANAL_Global_Tertile",
    "Census_Region",
    "RUCC",
    "Economic_type",
    "Homeownership_tertile",
]


def _fips5(s: pd.Series) -> pd.Series:
    return s.astype(str).str.extract(r"(\d+)", expand=False).str.zfill(5)


def _load_csv(path: Path, label: str = "") -> pd.DataFrame | None:
    if not path.exists():
        print(f"  Warning: {label or path.name} not found")
        return None
    df = pd.read_csv(path, dtype={"COUNTY_FIPS": str, "GEOID": str})
    if "COUNTY_FIPS" in df.columns:
        df["COUNTY_FIPS"] = _fips5(df["COUNTY_FIPS"])
    if "GEOID" in df.columns:
        df["GEOID"] = _fips5(df["GEOID"])
    return df


def _parse_filename(filename: str) -> tuple[str | None, str | None]:
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)
    if len(parts) == 2 and re.match(r"\d{4}-\d{4}$", parts[0]):
        return parts[0], parts[1].replace("-", "_")
    return None, None


def _pick_period_col(wide_df: pd.DataFrame | None, col_prefix: str, suffix: str) -> pd.DataFrame | None:
    if wide_df is None:
        return None
    col = f"{col_prefix}_{suffix}"
    if col not in wide_df.columns:
        return None
    return wide_df[["COUNTY_FIPS", col]].rename(columns={col: col_prefix})


def _cast_int(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    return df


def _rank_category(values: pd.Series, n: int) -> pd.Series:
    """Return 1..n quantile categories, robust to ties and missing values."""
    out = pd.Series(pd.NA, index=values.index, dtype="Int64")
    ok = values.notna()
    if ok.sum() == 0:
        return out
    ranks = values[ok].rank(method="first")
    out.loc[ok] = pd.qcut(ranks, q=n, labels=range(1, n + 1)).astype("Int64")
    return out


def _load_static() -> pd.DataFrame | None:
    pieces: list[pd.DataFrame] = []
    for fname, col in [("Census_Region.csv", "Census_Region"), ("Climate_Zone.csv", "Climate_Zone")]:
        df = _load_csv(COVARIATE_DIR / fname, col)
        if df is not None and col in df.columns:
            pieces.append(df[["COUNTY_FIPS", col]])
    for fname, col in [("Cluster_EQI.csv", "Cluster_EQI"), ("Cluster_NLCD.csv", "Cluster_NLCD")]:
        df = _load_csv(CLUSTER_DIR / fname, col)
        if df is not None and col in df.columns:
            pieces.append(df[["COUNTY_FIPS", col]])
    if not pieces:
        return None
    out = pieces[0]
    for piece in pieces[1:]:
        out = out.merge(piece, on="COUNTY_FIPS", how="outer")
    return out


def _load_rucc() -> dict[str, pd.DataFrame]:
    rucc: dict[str, pd.DataFrame] = {}
    for code, period in [("0005", "2000-2005"), ("0610", "2006-2010")]:
        fp = EQI_DIR / f"EQI{code}.csv"
        df = _load_csv(fp, f"EQI{code}")
        if df is None:
            continue
        df = df.rename(columns={k: v for k, v in EQI_COL_RENAMES.items() if k in df.columns})
        keep = ["COUNTY_FIPS"] + [c for c in ["RUCC"] if c in df.columns]
        rucc[period] = df[keep].copy()
    return rucc


def _load_svi() -> pd.DataFrame | None:
    svi = _load_csv(SVI_RESULT_PATH, "SVI_Result")
    air = _load_csv(AIR_SVI_PATH, "Air_SVI")

    pieces: list[pd.DataFrame] = []
    if svi is not None and "SVI" in svi.columns:
        pieces.append(svi[["COUNTY_FIPS", "SVI"]])
    if air is not None:
        renames = {"SVI": "SVI_Cont", "SVI_2": "SVI_Cont_2", "SVI_3": "SVI_Cont_3"}
        air = air.rename(columns=renames)
        keep = ["COUNTY_FIPS"] + [c for c in ["SVI_Cont", "SVI_Cont_2", "SVI_Cont_3"] if c in air.columns]
        if len(keep) > 1:
            pieces.append(air[keep])
    if not pieces:
        return None
    out = pieces[0]
    for piece in pieces[1:]:
        out = out.merge(piece, on="COUNTY_FIPS", how="outer")
    return out


def _build_anal_windows() -> pd.DataFrame:
    anal = _load_csv(ANAL_PATH, "ANAL panel")
    if anal is None:
        print(f"\nError: {ANAL_PATH} not found")
        sys.exit(1)

    needed = {"GEOID", "year", *ANAL_METRICS}
    missing = needed - set(anal.columns)
    if missing:
        print(f"\nError: ANAL panel missing columns: {sorted(missing)}")
        sys.exit(1)

    anal["COUNTY_FIPS"] = _fips5(anal["GEOID"])
    anal["year"] = pd.to_numeric(anal["year"], errors="coerce").astype("Int64")
    for col in ANAL_METRICS:
        anal[col] = pd.to_numeric(anal[col], errors="coerce")

    county_cols = ["COUNTY_FIPS"]
    if "NAME" in anal.columns:
        county_cols.append("NAME")
    if "STATEFP" in anal.columns:
        anal["State_FIPS"] = pd.to_numeric(anal["STATEFP"], errors="coerce").astype("Int64")
        county_cols.append("State_FIPS")
    county = anal[county_cols].drop_duplicates("COUNTY_FIPS")

    rows: list[pd.DataFrame] = []
    for spec in LAG_WINDOWS:
        start = spec["ANAL_Start"]
        end = spec["ANAL_End"]
        period = f"{start}-{end}"
        sub = anal[(anal["year"] >= start) & (anal["year"] <= end)].copy()
        agg = sub.groupby("COUNTY_FIPS", as_index=False).agg(
            ANAL_years_observed=("year", "nunique"),
            **{col: (col, "mean") for col in ANAL_METRICS},
        )
        agg = county.merge(agg, on="COUNTY_FIPS", how="left")
        agg["ANAL_Period"] = period
        agg["ANAL_Start"] = start
        agg["ANAL_End"] = end
        agg["ANAL_Years_Expected"] = end - start + 1
        agg["AAMR_Period"] = spec["AAMR_Period"]
        agg["Lag_Years"] = spec["Lag_Years"]
        agg["Lag_Role"] = spec["Role"]
        agg["Matched_Pair"] = f"{period}->{spec['AAMR_Period']}"
        rows.append(agg)

    out = pd.concat(rows, ignore_index=True)
    out["ANAL_log1p_popw_mean_rad"] = np.log1p(out["popw_mean_rad"])
    out["ANAL_Quintile"] = out.groupby("ANAL_Period", group_keys=False)["popw_mean_rad"].apply(
        lambda s: _rank_category(s, 5)
    )
    out["ANAL_Tertile"] = out.groupby("ANAL_Period", group_keys=False)["popw_mean_rad"].apply(
        lambda s: _rank_category(s, 3)
    )
    out["ANAL_Global_Quintile"] = _rank_category(out["popw_mean_rad"], 5)
    out["ANAL_Global_Tertile"] = _rank_category(out["popw_mean_rad"], 3)
    mean = out["ANAL_log1p_popw_mean_rad"].mean(skipna=True)
    sd = out["ANAL_log1p_popw_mean_rad"].std(skipna=True)
    out["ANAL_log1p_popw_mean_rad_z"] = (out["ANAL_log1p_popw_mean_rad"] - mean) / sd
    out["ANAL_SVI_Tertile_Label"] = "a" + out["ANAL_Tertile"].astype("string")
    return out


def main() -> None:
    print("=" * 70)
    print("ANAL x CDC Triangulation AAMR - Long Table Builder")
    print("=" * 70)

    if not AAMR_DIR.exists():
        print(f"\nError: {AAMR_DIR} not found")
        sys.exit(1)
    aamr_files = sorted(AAMR_DIR.glob("*.csv"))
    if not aamr_files:
        print(f"\nError: no AAMR files found in {AAMR_DIR}")
        sys.exit(1)

    print(f"\nFound {len(aamr_files)} AAMR files")
    print("\nLoading ANAL windows and covariates...")
    anal_windows = _build_anal_windows()
    static_df = _load_static()
    rucc_by_period = _load_rucc()
    svi_df = _load_svi()
    homeown_df = _load_csv(COVARIATE_DIR / "Homeownership_rate.csv", "Homeownership")
    uninsured_df = _load_csv(COVARIATE_DIR / "Uninsured_rate.csv", "Uninsured_rate")
    covariates = {
        "Smoking_rate": _load_csv(IHME_DIR / "IHME_Smoking.csv", "Smoking"),
        "Heavy_Drinking_rate": _load_csv(IHME_DIR / "IHME_Drinking.csv", "Drinking"),
        "Physical_Activities_rate": _load_csv(IHME_DIR / "IHME_Physical_Activities.csv", "Physical_Activities"),
        "Obesity_rate": _load_csv(IHME_DIR / "IHME_Obesity.csv", "Obesity"),
        "Diabetes_Prevalence_rate": _load_csv(IHME_DIR / "IHME_Diabetes_Prevalence.csv", "Diabetes_Prevalence"),
        "Physician_Density_per100k": _load_csv(COVARIATE_DIR / "Physician_Density_per100k.csv", "Physician"),
        "Forest_Coverage": _load_csv(COVARIATE_DIR / "Forest_Coverage.csv", "Forest"),
        "AQS_Number": _load_csv(COVARIATE_DIR / "AQS_Number.csv", "AQS_Number"),
        "Economic_type": _load_csv(COVARIATE_DIR / "Economic_type.csv", "Economic_type"),
    }

    print("\nProcessing AAMR files...")
    all_dfs: list[pd.DataFrame] = []
    windows_by_aamr = {
        period: df.copy()
        for period, df in anal_windows.groupby("AAMR_Period", sort=False)
    }
    seen_period_outcomes: set[tuple[str, str]] = set()

    for fp in aamr_files:
        aamr_period, outcome = _parse_filename(fp.name)
        if not aamr_period or aamr_period not in windows_by_aamr:
            print(f"  Skipping {fp.name} (period not in ANAL lag design)")
            continue
        period_outcome = (aamr_period, outcome)
        if period_outcome in seen_period_outcomes:
            print(f"  Skipping {fp.name} (duplicate normalized outcome {outcome})")
            continue
        seen_period_outcomes.add(period_outcome)

        aamr = pd.read_csv(fp, dtype={"COUNTY_FIPS": str})
        aamr["COUNTY_FIPS"] = _fips5(aamr["COUNTY_FIPS"])
        base = pd.DataFrame({
            "COUNTY_FIPS": aamr["COUNTY_FIPS"],
            "AAMR_Period": aamr_period,
            "Time_Period": aamr_period,
            "Outcome": outcome,
            "Deaths": aamr.get("Deaths"),
            "Population": aamr.get("Population"),
            "AAMR": aamr.get("AAMR"),
            "AAMR_Lower": aamr.get("AAMR_Lower"),
            "AAMR_Upper": aamr.get("AAMR_Upper"),
        })

        row = windows_by_aamr[aamr_period].merge(base, on=["COUNTY_FIPS", "AAMR_Period"], how="inner")

        if static_df is not None:
            row = row.merge(static_df, on="COUNTY_FIPS", how="left")
        if svi_df is not None:
            row = row.merge(svi_df, on="COUNTY_FIPS", how="left")
            if "SVI" in row.columns:
                row["ANAL_SVI_Tertile"] = row["ANAL_SVI_Tertile_Label"] + "s" + row["SVI"].astype("string")

        row["Covariate_Period"] = pd.NA
        for anal_period, (suffix, cov_period) in COVARIATE_PERIOD_BY_ANAL.items():
            mask = row["ANAL_Period"].eq(anal_period)
            if not mask.any():
                continue
            row.loc[mask, "Covariate_Period"] = cov_period

            rucc = rucc_by_period.get(cov_period)
            if rucc is not None:
                row = row.merge(
                    rucc.add_suffix(f"__{anal_period}").rename(columns={f"COUNTY_FIPS__{anal_period}": "COUNTY_FIPS"}),
                    on="COUNTY_FIPS",
                    how="left",
                )
                for col in ["RUCC"]:
                    period_col = f"{col}__{anal_period}"
                    if period_col in row.columns:
                        if col not in row.columns:
                            row[col] = pd.NA
                        row.loc[mask, col] = row.loc[mask, period_col]
                        row = row.drop(columns=period_col)

            for prefix, wide in covariates.items():
                tmp = _pick_period_col(wide, prefix, suffix)
                if tmp is None:
                    continue
                tmp = tmp.rename(columns={prefix: f"{prefix}__{anal_period}"})
                row = row.merge(tmp, on="COUNTY_FIPS", how="left")
                period_col = f"{prefix}__{anal_period}"
                if prefix not in row.columns:
                    row[prefix] = pd.NA
                row.loc[mask, prefix] = row.loc[mask, period_col]
                row = row.drop(columns=period_col)

        outcome_suffix = OUTCOME_SUFFIX[aamr_period]
        if homeown_df is not None:
            rate_col = f"Homeownership_rate_{outcome_suffix}"
            tert_col = f"Homeownership_tertile_{outcome_suffix}"
            hw_cols = [c for c in [rate_col, tert_col] if c in homeown_df.columns]
            if hw_cols:
                row = row.merge(
                    homeown_df[["COUNTY_FIPS"] + hw_cols].rename(columns={
                        rate_col: "Homeownership_rate",
                        tert_col: "Homeownership_tertile",
                    }),
                    on="COUNTY_FIPS",
                    how="left",
                )
        uninsured = _pick_period_col(uninsured_df, "Uninsured_rate", outcome_suffix)
        if uninsured is not None:
            row = row.merge(uninsured, on="COUNTY_FIPS", how="left")

        all_dfs.append(row)
        print(f"  {fp.name}: {len(row)} rows across {row['Matched_Pair'].nunique()} pair(s)")

    if not all_dfs:
        print("\nNo rows produced.")
        return

    out = pd.concat(all_dfs, ignore_index=True)
    out = out.drop(columns=["ANAL_SVI_Tertile_Label"], errors="ignore")
    out = _cast_int(out, INT_COLS)
    sort_cols = ["AAMR_Period", "Outcome", "Lag_Years", "ANAL_Period", "COUNTY_FIPS"]
    out = out.sort_values(sort_cols).reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_DF, index=False)

    print("\n" + "=" * 70)
    print(f"Output:        {len(out):,} rows -> {OUTPUT_DF}")
    print(f"Columns:       {list(out.columns)}")
    print(f"Counties:      {out['COUNTY_FIPS'].nunique():,}")
    print(f"Outcomes:      {out['Outcome'].nunique():,}")
    print(f"AAMR periods:  {sorted(out['AAMR_Period'].dropna().unique())}")
    print(f"ANAL periods:  {sorted(out['ANAL_Period'].dropna().unique())}")
    print(f"Lags:          {sorted(out['Lag_Years'].dropna().unique())}")
    print(f"Matched pairs: {out['Matched_Pair'].nunique()}")
    print("\nPrimary ANAL quintile counts by exposure period:")
    print(out.drop_duplicates(["COUNTY_FIPS", "ANAL_Period"]).
          groupby(["ANAL_Period", "ANAL_Quintile"]).size().unstack(fill_value=0).to_string())
    if "SVI" in out.columns:
        print("\nSVI class counts:")
        print(out.drop_duplicates("COUNTY_FIPS")["SVI"].value_counts(dropna=False).sort_index().to_string())
    print("\nDone.")


if __name__ == "__main__":
    main()
