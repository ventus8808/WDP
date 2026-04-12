#!/usr/bin/env python3
"""
EQI × CDC Triangulation — AAMR long table
Merges triangulation AAMR results with EQI covariates to produce a long CSV.
Output contains rows per county × time_period × ICD × Lag_Years,
with columns AAMR, AAMR_SE, AAMR_Lower, AAMR_Upper, plus RUCC/EQI covariates.

Input: Data/Original/CDC Triangulation/AAMR/*.csv
Output: Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv
"""

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    CFG = yaml.safe_load(f)

# Input: triangulation AAMR results
TRIANGULATION_AAMR_DIR = PROJECT_ROOT / "Data/Original/CDC Triangulation/AAMR"

# Output: merged long table
OUTPUT_DIR = PROJECT_ROOT / "Data/Processed/df_EQI_AAMR_Triangulation"
OUTPUT_FILE = OUTPUT_DIR / "EQI_AAMR_Cluster_Climate_Typology_LandUse.csv"
SENSITIVITY_FILE = OUTPUT_DIR / "EQI_AAMR_Sensativity.csv"

# EQI data
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

# Typology data (USDA ERS County Economic Typology 2004)
TYPOLOGY_PATH = (
    PROJECT_ROOT
    / CFG["data_directories"]["processed"]
    / "Socioeconomic"
    / "County_Typology_2004.csv"
)

# LandUse cluster data (NLCD/JRC clusters)
LANDUSE_PATH = (
    PROJECT_ROOT / "Result/Cluster_Visualization_LandUse/LandUse_Clusters_All_K.csv"
)

# Sensitivity analysis covariates
INSURANCE_PATH = PROJECT_ROOT / "Data/Processed/Insurance/Insurance.csv"
DOCTOR_PATH = PROJECT_ROOT / "Data/Processed/Doctor/Doctor.csv"
FOREST_PATH = PROJECT_ROOT / "Data/Processed/Environmental/Forest_Coverage.csv"
MONITORING_PATH = PROJECT_ROOT / "Data/Processed/Monitoring/Monitoring.csv"

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

# Typology columns to merge
TYPOLOGY_COLS = ["econdep"]

# LandUse columns to merge (k=4 from NLCD/JRC analysis)
LANDUSE_COLS = ["cluster_4"]


# ---------------- Helpers ----------------


def _parse_filename(filename: str) -> tuple[str | None, str | None]:
    """
    Parse filename to extract year range and ICD code.
    Example: '2016-2020_G30.csv' -> ('2016-2020', 'G30')
    """
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)

    if len(parts) == 2:
        year_range, icd = parts
        # Validate year range format
        if re.match(r"\d{4}-\d{4}", year_range):
            return year_range, icd

    return None, None


def _load_eqi_by_period() -> dict:
    """Load EQI datasets for both periods (2000-2005, 2006-2010)"""
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
    """Load smoking rate data"""
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
    """Load cluster data"""
    if CLUSTER_PATH.exists():
        df = pd.read_csv(CLUSTER_PATH)
        if "COUNTY_FIPS" in df.columns:
            # Select only COUNTY_FIPS and cluster columns
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
    """Load climate zone data"""
    if CLIMATE_PATH.exists():
        df = pd.read_csv(CLIMATE_PATH, dtype={"COUNTY_FIPS": str})
        if "COUNTY_FIPS" in df.columns:
            # Select only COUNTY_FIPS and climate columns
            cols = ["COUNTY_FIPS"] + [c for c in CLIMATE_COLS if c in df.columns]
            df = df[cols].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            return df
        else:
            print(f"  ⚠️ Climate data missing COUNTY_FIPS column")
    else:
        print(f"  ⚠️ Climate data not found at {CLIMATE_PATH}")
    return None


def _load_typology() -> pd.DataFrame | None:
    """Load county economic typology data"""
    if TYPOLOGY_PATH.exists():
        df = pd.read_csv(TYPOLOGY_PATH, dtype={"COUNTY_FIPS": str})
        if "COUNTY_FIPS" in df.columns:
            # Select only COUNTY_FIPS and typology columns
            cols = ["COUNTY_FIPS"] + [c for c in TYPOLOGY_COLS if c in df.columns]
            df = df[cols].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            return df
        else:
            print(f"  ⚠️ Typology data missing COUNTY_FIPS column")
    else:
        print(f"  ⚠️ Typology data not found at {TYPOLOGY_PATH}")
    return None


def _load_landuse() -> pd.DataFrame | None:
    """Load land use cluster data"""
    if LANDUSE_PATH.exists():
        df = pd.read_csv(LANDUSE_PATH, dtype={"COUNTY_FIPS": str})
        if "COUNTY_FIPS" in df.columns:
            # Select only COUNTY_FIPS and landuse columns
            cols = ["COUNTY_FIPS"] + [c for c in LANDUSE_COLS if c in df.columns]
            df = df[cols].copy()
            df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)
            # Rename cluster_4 to landuse_cluster to avoid confusion with EQI cluster_4
            if "cluster_4" in df.columns:
                df = df.rename(columns={"cluster_4": "landuse_cluster"})
            return df
        else:
            print(f"  ⚠️ LandUse data missing COUNTY_FIPS column")
    else:
        print(f"  ⚠️ LandUse data not found at {LANDUSE_PATH}")
    return None


def _load_insurance() -> pd.DataFrame | None:
    """Load county uninsured rate data (SAHIE 2005)"""
    if INSURANCE_PATH.exists():
        df = pd.read_csv(INSURANCE_PATH, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        return df
    print(f"  ⚠️ Insurance data not found at {INSURANCE_PATH}")
    return None


def _load_doctor() -> pd.DataFrame | None:
    """Load county physician density data (AHRF 2000-2005 mean)"""
    if DOCTOR_PATH.exists():
        df = pd.read_csv(DOCTOR_PATH, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        return df
    print(f"  ⚠️ Doctor data not found at {DOCTOR_PATH}")
    return None


def _load_forest() -> pd.DataFrame | None:
    """Load county forest coverage data (NLCD 2000-2005 mean)"""
    if FOREST_PATH.exists():
        df = pd.read_csv(FOREST_PATH, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        return df
    print(f"  ⚠️ Forest coverage data not found at {FOREST_PATH}")
    return None


def _load_monitoring() -> pd.DataFrame | None:
    """Load county EPA monitoring site count data (2000-2006 mean)"""
    if MONITORING_PATH.exists():
        df = pd.read_csv(MONITORING_PATH, dtype={"COUNTY_FIPS": str})
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        return df
    print(f"  ⚠️ Monitoring data not found at {MONITORING_PATH}")
    return None


def _get_valid_lag_combinations(time_period: str) -> list[tuple[int, str, str]]:
    """
    Get valid (lag_years, eqi_code, eqi_period) combinations for a given AAMR time period.

    Valid combinations:
    1. EQI 2000-2005 (code '0005') + 5 years lag -> AAMR 2006-2010
    2. EQI 2000-2005 (code '0005') + 10 years lag -> AAMR 2011-2015
    3. EQI 2000-2005 (code '0005') + 15 years lag -> AAMR 2016-2020
    4. EQI 2006-2010 (code '0610') + 5 years lag -> AAMR 2011-2015
    5. EQI 2006-2010 (code '0610') + 10 years lag -> AAMR 2016-2020

    Returns: list of (lag_years, eqi_code, eqi_period)
    """
    combinations = []

    if time_period == "2006-2010":
        # Only 2000-2005 EQI with 5 years lag
        combinations.append((5, "0005", "2000-2005"))
    elif time_period == "2011-2015":
        # 2000-2005 EQI with 10 years lag AND 2006-2010 EQI with 5 years lag
        combinations.append((10, "0005", "2000-2005"))
        combinations.append((5, "0610", "2006-2010"))
    elif time_period == "2016-2020":
        # 2000-2005 EQI with 15 years lag AND 2006-2010 EQI with 10 years lag
        combinations.append((15, "0005", "2000-2005"))
        combinations.append((10, "0610", "2006-2010"))

    return combinations


def _process_aamr_file(
    file_path: Path,
    year_range: str,
    icd_code: str,
    eqi_dict: dict,
    smoking_df: pd.DataFrame | None,
    cluster_df: pd.DataFrame | None,
    climate_df: pd.DataFrame | None,
    typology_df: pd.DataFrame | None,
    landuse_df: pd.DataFrame | None,
) -> list[pd.DataFrame]:
    """
    Process one AAMR file and create rows for all valid lag combinations.

    Returns: list of DataFrames (one per valid lag combination)
    """
    print(f"  Processing: {file_path.name}")

    # Load AAMR data
    df = pd.read_csv(file_path, dtype={"COUNTY_FIPS": str})

    # Ensure COUNTY_FIPS is 5-digit string
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].astype(str).str.zfill(5)

    # Format ICD code: replace hyphens with underscores
    icd_fmt = icd_code.replace("-", "_")

    # Get valid lag combinations for this time period
    valid_combinations = _get_valid_lag_combinations(year_range)

    if not valid_combinations:
        print(f"    ⚠️ No valid EQI-AAMR lag combinations for {year_range}")
        return []

    result_dfs = []

    for lag, eqi_code, eqi_period in valid_combinations:
        # Create base DataFrame
        out = pd.DataFrame(
            {
                "COUNTY_FIPS": df["COUNTY_FIPS"].astype(str),
                "EQI_Period": eqi_period,
                "Time_Period": year_range,
                "Lag_Years": lag,
                "Cancer_Type": icd_fmt,
                "Deaths": df["Deaths"],
                "Population": df["Population"],
                "AAMR": df["AAMR"],
                "AAMR_SE": df["AAMR_SE"],
                "AAMR_Lower": df["AAMR_Lower"],
                "AAMR_Upper": df["AAMR_Upper"],
            }
        )

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

        # Merge Typology data
        if typology_df is not None:
            out = out.merge(typology_df, on="COUNTY_FIPS", how="left")
        else:
            for c in TYPOLOGY_COLS:
                out[c] = pd.NA

        # Merge LandUse data
        if landuse_df is not None:
            out = out.merge(landuse_df, on="COUNTY_FIPS", how="left")
        else:
            out["landuse_cluster"] = pd.NA

        result_dfs.append(out)

    print(
        f"    ✓ Created {len(result_dfs)} lag combination(s) with {len(df)} counties each"
    )

    return result_dfs


# ---------------- Main ----------------


def main():
    print("=" * 70)
    print("EQI × CDC Triangulation — AAMR Long Table Builder")
    print("=" * 70)

    # Check input directory
    if not TRIANGULATION_AAMR_DIR.exists():
        print(f"\n⚠️ Input directory not found: {TRIANGULATION_AAMR_DIR}")
        print("Please run EQI_Triangulation_AAMR.py first to generate AAMR data.")
        sys.exit(1)

    # Get all AAMR files
    aamr_files = sorted(TRIANGULATION_AAMR_DIR.glob("*.csv"))

    if not aamr_files:
        print(f"\n⚠️ No AAMR files found in {TRIANGULATION_AAMR_DIR}")
        print("Please run EQI_Triangulation_AAMR.py first to generate AAMR data.")
        sys.exit(1)

    print(f"\nFound {len(aamr_files)} AAMR files")

    # Load EQI, Smoking, Cluster, Climate, Typology, LandUse, and sensitivity covariates
    print("\n📊 Loading EQI, Smoking, Cluster, Climate, Typology, LandUse, and sensitivity covariates...")
    eqi_dict = _load_eqi_by_period()
    smoking_df = _load_smoking()
    cluster_df = _load_clusters()
    climate_df = _load_climate()
    typology_df = _load_typology()
    landuse_df = _load_landuse()
    insurance_df = _load_insurance()
    doctor_df = _load_doctor()
    forest_df = _load_forest()
    monitoring_df = _load_monitoring()

    if not eqi_dict:
        print("⚠️ No EQI data loaded. Continuing without EQI covariates.")

    if smoking_df is None:
        print("⚠️ No smoking data loaded. Continuing without smoking rates.")

    if cluster_df is None:
        print("⚠️ No cluster data loaded. Continuing without cluster IDs.")

    if climate_df is None:
        print("⚠️ No climate data loaded. Continuing without climate zones.")

    if typology_df is None:
        print("⚠️ No typology data loaded. Continuing without county typology.")

    if landuse_df is None:
        print("⚠️ No land use data loaded. Continuing without land use clusters.")

    # Process all AAMR files
    print("\n" + "=" * 70)
    print("Processing AAMR Files and Merging with EQI")
    print("=" * 70 + "\n")

    all_rows = []

    for file_path in aamr_files:
        year_range, icd_code = _parse_filename(file_path.name)

        if not year_range or not icd_code:
            print(f"  ⚠️ Could not parse filename: {file_path.name}, skipping")
            continue

        # Process file and get all lag combinations
        result_dfs = _process_aamr_file(
            file_path,
            year_range,
            icd_code,
            eqi_dict,
            smoking_df,
            cluster_df,
            climate_df,
            typology_df,
            landuse_df,
        )

        all_rows.extend(result_dfs)

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
            + TYPOLOGY_COLS
            + ["landuse_cluster"]
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
        + [c for c in TYPOLOGY_COLS if c in final.columns]
        + ["landuse_cluster"]
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
    # Note: koppen_major is a string column, so we don't cast it
    numeric_climate_cols = ["census_region", "doe_major"]
    for c in numeric_climate_cols:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors="coerce").astype("Int64")

    # Cast typology columns to nullable int
    for c in TYPOLOGY_COLS:
        if c in final.columns:
            final[c] = pd.to_numeric(final[c], errors="coerce").astype("Int64")

    # Cast landuse_cluster to nullable int
    if "landuse_cluster" in final.columns:
        final["landuse_cluster"] = pd.to_numeric(
            final["landuse_cluster"], errors="coerce"
        ).astype("Int64")

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
    print(f"Cancer types:      {final['Cancer_Type'].nunique()}")
    print(f"Lag combinations:  {sorted(final['Lag_Years'].unique())}")
    print(f"\nOutput saved to:   {OUTPUT_FILE}")

    # Show example
    print("\n" + "=" * 70)
    print("Example Output (first 10 rows)")
    print("=" * 70)
    print(final.head(10).to_string(index=False))

    # ---- Sensitivity output: EQI 2000-2005 only + healthcare/green space covariates ----
    print("\n" + "=" * 70)
    print("Building Sensitivity Table (EQI 2000-2005 only)")
    print("=" * 70)

    sensitivity = final[final["EQI_Period"] == "2000-2005"].copy()

    for df_cov, label in [
        (insurance_df, "insurance"),
        (doctor_df, "doctor"),
        (forest_df, "forest"),
        (monitoring_df, "monitoring"),
    ]:
        if df_cov is not None:
            sensitivity = sensitivity.merge(df_cov, on="COUNTY_FIPS", how="left")
        else:
            print(f"  ⚠️ {label} data unavailable, column(s) will be missing")

    # Counties with no EPA monitors get 0 (genuinely no monitoring, not missing)
    if "site_number_mean" in sensitivity.columns:
        sensitivity["site_number_mean"] = sensitivity["site_number_mean"].fillna(0)

    sensitivity.to_csv(SENSITIVITY_FILE, index=False)
    print(f"Sensitivity rows:  {len(sensitivity):,}")
    print(f"Output saved to:   {SENSITIVITY_FILE}")

    print("\n✓ EQI × Triangulation AAMR merge completed successfully!")


if __name__ == "__main__":
    main()
