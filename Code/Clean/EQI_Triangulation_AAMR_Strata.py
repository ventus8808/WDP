"""
Stratified AAMR Calculator for RUCC, Climate, Cluster, Typology, and LandUse

This script calculates national and stratum-level Age-Adjusted Mortality Rates (AAMRs)
with confidence intervals using triangulated death data, stratified by RUCC, climate zones,
EQI clusters, county economic typology, and land use clusters.

Stratifications:
- RUCC: Rural-Urban Continuum Codes from EQI data
- Climate: census_region, koppen_major, doe_major
- Cluster: cluster_3, cluster_4, cluster_5
- Typology: econdep (USDA ERS County Economic Typology 2004)
- LandUse: cluster_4 (Land Use Cluster from NLCD/JRC analysis)

Input: Data/Original/CDC Triangulation/Subtracted/*.csv
Output:
    - Result/Tables/Stratified_AAMR_{RUCC,Climate,Cluster,Typology,LandUse}.csv
    - Result/Tables/Top5_Cancer_AAMR_{RUCC,Climate,Cluster,Typology,LandUse}.csv (summary tables)
    - Result/Tables/Top5_NDD_AAMR_{RUCC,Climate,Cluster,Typology,LandUse}.csv (summary tables)
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml
from scipy.stats import gamma

# 2000 US Standard Population weights
AGE_WEIGHTS = {
    "< 1 year": 0.013818,
    "1-4 years": 0.055317,
    "5-14 years": 0.145565,
    "15-24 years": 0.138646,
    "25-34 years": 0.135573,
    "35-44 years": 0.162613,
    "45-54 years": 0.134834,
    "55-64 years": 0.087247,
    "65-74 years": 0.066037,
    "75-84 years": 0.044842,
    "85+ years": 0.015508,
}

AGE_GROUPS = list(AGE_WEIGHTS.keys())

CLUSTER_COLS = ["cluster_3", "cluster_4", "cluster_5"]

CLIMATE_COLS = ["census_region", "koppen_major", "doe_major"]

TYPOLOGY_COLS = ["econdep"]

LANDUSE_COLS = ["cluster_4"]

# NDD ICD code to abbreviation mapping
NDD_ABBR_MAP = {
    "G20_G30_G12.2_F01_F03": "NDD",
    "G30_F01_F03": "Dementia",
    "G30": "AD",
    "G20": "PD",
    "F01": "VD",
    "F03": "UD",
    "G10": "HD",
    "G12.2": "ALS",
}

# All NDD ICD codes
NDD_ALL_CODES = ["G20_G30_G12.2_F01_F03", "G30_F01_F03", "G30", "G20", "F01", "F03", "G10", "G12.2"]

# Constants for Dementia synthesis
_DEMENTIA_CODES = ["G30", "F01", "F03"]
_DEMENTIA_ICD = "G30_F01_F03"
_DEMENTIA_PERIODS = ["2006-2010", "2011-2015", "2016-2020"]


def load_config() -> Tuple[Path, Dict]:
    """Load configuration from config.yaml"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return project_root, config


def apply_stratification_labels(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    """Apply human-readable labels to stratification values based on config mappings."""
    if df.empty:
        return df

    df = df.copy()
    label_mappings = config.get("stratification_labels", {})

    # Apply mappings based on Stratum_Type
    for idx, row in df.iterrows():
        stratum_type = row.get("Stratum_Type", "")
        stratum_value = str(row.get("Stratum_Value", ""))

        # Skip National
        if stratum_value == "National":
            continue

        # Apply appropriate mapping
        if stratum_type == "Typology":
            mapping = label_mappings.get("typology", {})
            if stratum_value in mapping:
                df.at[idx, "Stratum_Label"] = mapping[stratum_value]
            else:
                df.at[idx, "Stratum_Label"] = stratum_value

        elif stratum_type == "LandUse":
            mapping = label_mappings.get("landuse", {})
            if stratum_value in mapping:
                df.at[idx, "Stratum_Label"] = mapping[stratum_value]
            else:
                df.at[idx, "Stratum_Label"] = stratum_value

        elif stratum_type.startswith("Cluster_k"):
            k = stratum_type.replace("Cluster_k", "")
            mapping = label_mappings.get(f"cluster_k{k}", {})
            if stratum_value in mapping:
                df.at[idx, "Stratum_Label"] = mapping[stratum_value]
            else:
                df.at[idx, "Stratum_Label"] = stratum_value
        else:
            # For RUCC and others, keep original value
            df.at[idx, "Stratum_Label"] = stratum_value

    # For National rows, set label to National
    df.loc[df["Stratum_Value"] == "National", "Stratum_Label"] = "National"

    return df


def parse_filename(filename: str) -> Optional[Tuple[str, str]]:
    """Parse filename to extract year range and ICD code."""
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return None


def calculate_age_specific_rate(deaths: float, population: float) -> float:
    """Calculate age-specific death rate per 100,000."""
    if population == 0:
        return 0.0
    return (deaths / population) * 100000


def calculate_aamr_point(agg_df: pd.DataFrame) -> Tuple[float, Dict]:
    """Calculate AAMR point estimate and stats."""
    aamr = 0.0
    total_deaths = 0
    total_population = 0
    for age_group in AGE_GROUPS:
        if age_group in agg_df.index:
            deaths = agg_df.at[age_group, "Deaths"]
            population = agg_df.at[age_group, "Population"]
            mi = calculate_age_specific_rate(deaths, population)
            aamr += mi * AGE_WEIGHTS[age_group]
            total_deaths += int(deaths)
            total_population += int(population)
    stats = {
        "total_deaths": total_deaths,
        "total_population": total_population,
    }
    return aamr, stats


def calculate_aamr_standard_error(agg_df: pd.DataFrame) -> float:
    """Calculate standard error of AAMR using Poisson variance."""
    variance = 0.0
    for age_group in AGE_GROUPS:
        if age_group in agg_df.index:
            weight = AGE_WEIGHTS[age_group]
            deaths = agg_df.at[age_group, "Deaths"]
            population = agg_df.at[age_group, "Population"]
            if population > 0:
                variance += (weight**2) * (deaths / (population**2))
    return (variance**0.5) * 100000


def calculate_aamr_ci(aamr: float, stats: Dict) -> Tuple[float, float]:
    """Calculate 95% CI using Fay-Feuer method."""
    total_deaths = stats["total_deaths"]
    if total_deaths == 0 or aamr == 0:
        return 0.0, 0.0
    p_star = (total_deaths * 100000) / aamr
    if total_deaths == 0:
        d_lower = 0.0
        d_upper = gamma.ppf(0.975, 1)
    else:
        d_lower = gamma.ppf(0.025, total_deaths)
        d_upper = gamma.ppf(0.975, total_deaths + 1)
    aamr_lower = (d_lower / p_star) * 100000
    aamr_upper = (d_upper / p_star) * 100000
    return aamr_lower, aamr_upper


def load_stratification_data(
    project_root: Path, config: Dict
) -> Dict[str, pd.DataFrame]:
    """Load stratification mappings using config paths, only specified columns."""

    # Helper to get path
    def get_path(keys: List[str], default: str) -> Path:
        d = config
        for key in keys[:-1]:
            d = d.get(key, {})
        path_str = d.get(keys[-1], default)
        return project_root / path_str

    # RUCC from EQI
    eqi_path = (
        get_path(["data_directories", "processed"], "Data/Processed")
        / "EQI"
        / "EQI0610.csv"
    )
    rucc_df = pd.DataFrame()
    if eqi_path.exists():
        eqi_df = pd.read_csv(eqi_path, dtype={"COUNTY_FIPS": str})
        eqi_df["COUNTY_FIPS"] = eqi_df["COUNTY_FIPS"].str.zfill(5)
        if "RUCC" in eqi_df.columns:
            rucc_df = eqi_df[["COUNTY_FIPS", "RUCC"]].copy()
    else:
        print(f"Warning: EQI file not found at {eqi_path}")

    # Climate, only specified columns
    climate_path = (
        get_path(["eqi_aamr_outputs", "base_dir"], "Result/Tables")
        / "Climate_Zone_Processed.csv"
    )
    climate_df = pd.DataFrame()
    if climate_path.exists():
        temp_df = pd.read_csv(climate_path, dtype={"COUNTY_FIPS": str})
        temp_df["COUNTY_FIPS"] = temp_df["COUNTY_FIPS"].str.zfill(5)
        avail_cols = [col for col in CLIMATE_COLS if col in temp_df.columns]
        if avail_cols:
            climate_df = temp_df[["COUNTY_FIPS"] + avail_cols].copy()
    else:
        print(f"Warning: Climate file not found at {climate_path}")

    # Cluster, only specified columns
    cluster_path = (
        get_path(
            ["result_directories", "cluster_visualization"],
            "Result/Cluster_Visualization",
        )
        / "EQI_Clusters_All_K.csv"
    )
    cluster_df = pd.DataFrame()
    if cluster_path.exists():
        temp_df = pd.read_csv(cluster_path, dtype={"COUNTY_FIPS": str})
        temp_df["COUNTY_FIPS"] = temp_df["COUNTY_FIPS"].str.zfill(5)
        avail_cols = [col for col in CLUSTER_COLS if col in temp_df.columns]
        if avail_cols:
            cluster_df = temp_df[["COUNTY_FIPS"] + avail_cols].copy()
    else:
        print(f"Warning: Cluster file not found at {cluster_path}")

    # Typology, only specified columns
    typology_path = (
        get_path(["data_directories", "processed"], "Data/Processed")
        / "Socioeconomic"
        / "County_Typology_2004.csv"
    )
    typology_df = pd.DataFrame()
    if typology_path.exists():
        temp_df = pd.read_csv(typology_path, dtype={"COUNTY_FIPS": str})
        temp_df["COUNTY_FIPS"] = temp_df["COUNTY_FIPS"].str.zfill(5)
        avail_cols = [col for col in TYPOLOGY_COLS if col in temp_df.columns]
        if avail_cols:
            typology_df = temp_df[["COUNTY_FIPS"] + avail_cols].copy()
    else:
        print(f"Warning: Typology file not found at {typology_path}")

    # LandUse, only specified columns
    landuse_path = (
        project_root
        / "Result"
        / "Cluster_Visualization_LandUse"
        / "LandUse_Clusters_All_K.csv"
    )
    landuse_df = pd.DataFrame()
    if landuse_path.exists():
        temp_df = pd.read_csv(landuse_path, dtype={"COUNTY_FIPS": str})
        temp_df["COUNTY_FIPS"] = temp_df["COUNTY_FIPS"].str.zfill(5)
        avail_cols = [col for col in LANDUSE_COLS if col in temp_df.columns]
        if avail_cols:
            landuse_df = temp_df[["COUNTY_FIPS"] + avail_cols].copy()
    else:
        print(f"Warning: LandUse file not found at {landuse_path}")

    return {
        "RUCC": rucc_df,
        "Climate": climate_df,
        "Cluster": cluster_df,
        "Typology": typology_df,
        "LandUse": landuse_df,
    }


def process_file(
    file_path: Path, strat_data: Dict[str, pd.DataFrame], config: Dict
) -> Dict[str, List[Dict]]:
    """Process one subtracted file and calculate stratified AAMRs, handling sub-strata."""
    year_range, icd_code = parse_filename(file_path.name) or (None, None)
    if not year_range or not icd_code:
        return {}
    icd_code_norm = icd_code.replace("-", "_")

    df = pd.read_csv(file_path, dtype={"County Code": str})
    df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
    df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
    df["Population"] = pd.to_numeric(df["Population"], errors="coerce").fillna(0)
    df = df.rename(
        columns={"County Code": "COUNTY_FIPS", "Ten-Year Age Groups": "Age_Group"}
    )
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)

    results = {"RUCC": [], "Climate": [], "Cluster": [], "Typology": [], "LandUse": []}

    # Helper to add result
    def add_result(
        res_list: List,
        strat_type: str,
        strat_value: str,
        aamr: float,
        se: float,
        lower: float,
        upper: float,
        stats: Dict,
    ):
        # Normalize stratum value to consistent format
        if strat_value != "National":
            try:
                # Convert to float then format consistently
                val_float = float(strat_value)
                if val_float == int(val_float):
                    strat_value = str(int(val_float))
                else:
                    strat_value = str(val_float)
            except (ValueError, TypeError):
                strat_value = str(strat_value)

        res_list.append(
            {
                "Time_Period": year_range,
                "ICD-10 Code": icd_code,
                "Outcome": config.get("brms_analysis", {})
                .get("icd_mapping", {})
                .get(icd_code_norm, icd_code),
                "Stratum_Type": strat_type,
                "Stratum_Value": strat_value,
                "Total_Deaths": int(stats["total_deaths"]),
                "Total_Population": int(stats["total_population"]),
                "AAMR": round(aamr, 4),
                "AAMR_SE": round(se, 4),
                "AAMR_Lower": round(lower, 4),
                "AAMR_Upper": round(upper, 4),
                "AAMRSE": f"{round(aamr, 2):.2f} ± {round(se, 2):.2f}",
            }
        )

    for strat_key, sub_cols in [
        ("RUCC", ["RUCC"]),
        ("Climate", CLIMATE_COLS),
        ("Cluster", CLUSTER_COLS),
        ("Typology", TYPOLOGY_COLS),
        ("LandUse", LANDUSE_COLS),
    ]:
        if strat_data[strat_key].empty:
            continue
        merge_df = df.merge(strat_data[strat_key], on="COUNTY_FIPS", how="left")
        # Filter sub_cols to only include columns that exist in merged data
        sub_cols = [c for c in sub_cols if c in merge_df.columns]
        if not sub_cols:
            continue
        # Drop rows missing any sub_col
        merge_df = merge_df.dropna(subset=sub_cols)

        # Compute national once per strat_key, but will repeat for each sub_col
        nat_agg = merge_df.groupby("Age_Group")[["Deaths", "Population"]].sum()
        nat_aamr, nat_stats = calculate_aamr_point(nat_agg)
        nat_se = calculate_aamr_standard_error(nat_agg)
        nat_lower, nat_upper = calculate_aamr_ci(nat_aamr, nat_stats)

        for col in sub_cols:
            if strat_key == "RUCC":
                s_type = "RUCC"
            elif strat_key == "Climate":
                s_type = col  # e.g., 'census_region'
            elif strat_key == "Cluster":
                k = col.split("_")[1]
                s_type = f"Cluster_k{k}"
            elif strat_key == "Typology":
                s_type = "Typology"
            else:  # LandUse
                s_type = "LandUse"

            # Add national for this sub-type
            add_result(
                results[strat_key],
                s_type,
                "National",
                nat_aamr,
                nat_se,
                nat_lower,
                nat_upper,
                nat_stats,
            )

            # Stratum-level for this column
            grouped = merge_df.groupby([col, "Age_Group"])[
                ["Deaths", "Population"]
            ].sum()
            strata = sorted(grouped.index.get_level_values(0).unique())
            for stratum in strata:
                stratum_agg = grouped.loc[stratum]
                aamr, stats = calculate_aamr_point(stratum_agg)
                se = calculate_aamr_standard_error(stratum_agg)
                lower, upper = calculate_aamr_ci(aamr, stats)
                add_result(
                    results[strat_key],
                    s_type,
                    str(stratum),
                    aamr,
                    se,
                    lower,
                    upper,
                    stats,
                )

    return results


# ---------------------------------------------------------------------------
# Dementia synthesis helpers
# ---------------------------------------------------------------------------


def _process_clean_df(
    df: pd.DataFrame,
    year_range: str,
    icd_code: str,
    strat_data: Dict[str, pd.DataFrame],
    config: Dict,
) -> Dict[str, List[Dict]]:
    """
    Calculate stratified AAMRs from a pre-cleaned DataFrame.
    Expects columns: COUNTY_FIPS, Age_Group, Deaths (numeric), Population (numeric).
    """
    icd_code_norm = icd_code.replace("-", "_")
    results = {"RUCC": [], "Climate": [], "Cluster": [], "Typology": [], "LandUse": []}

    def add_result(res_list, strat_type, strat_value, aamr, se, lower, upper, stats):
        if strat_value != "National":
            try:
                val_float = float(strat_value)
                strat_value = (
                    str(int(val_float)) if val_float == int(val_float) else str(val_float)
                )
            except (ValueError, TypeError):
                strat_value = str(strat_value)
        res_list.append(
            {
                "Time_Period": year_range,
                "ICD-10 Code": icd_code,
                "Outcome": config.get("brms_analysis", {})
                .get("icd_mapping", {})
                .get(icd_code_norm, icd_code),
                "Stratum_Type": strat_type,
                "Stratum_Value": strat_value,
                "Total_Deaths": int(stats["total_deaths"]),
                "Total_Population": int(stats["total_population"]),
                "AAMR": round(aamr, 4),
                "AAMR_SE": round(se, 4),
                "AAMR_Lower": round(lower, 4),
                "AAMR_Upper": round(upper, 4),
                "AAMRSE": f"{round(aamr, 2):.2f} ± {round(se, 2):.2f}",
            }
        )

    for strat_key, sub_cols in [
        ("RUCC", ["RUCC"]),
        ("Climate", CLIMATE_COLS),
        ("Cluster", CLUSTER_COLS),
        ("Typology", TYPOLOGY_COLS),
        ("LandUse", LANDUSE_COLS),
    ]:
        if strat_data[strat_key].empty:
            continue
        merge_df = df.merge(strat_data[strat_key], on="COUNTY_FIPS", how="left")
        sub_cols = [c for c in sub_cols if c in merge_df.columns]
        if not sub_cols:
            continue
        merge_df = merge_df.dropna(subset=sub_cols)

        nat_agg = merge_df.groupby("Age_Group")[["Deaths", "Population"]].sum()
        nat_aamr, nat_stats = calculate_aamr_point(nat_agg)
        nat_se = calculate_aamr_standard_error(nat_agg)
        nat_lower, nat_upper = calculate_aamr_ci(nat_aamr, nat_stats)

        for col in sub_cols:
            if strat_key == "RUCC":
                s_type = "RUCC"
            elif strat_key == "Climate":
                s_type = col
            elif strat_key == "Cluster":
                k = col.split("_")[1]
                s_type = f"Cluster_k{k}"
            elif strat_key == "Typology":
                s_type = "Typology"
            else:
                s_type = "LandUse"

            add_result(
                results[strat_key], s_type, "National",
                nat_aamr, nat_se, nat_lower, nat_upper, nat_stats,
            )

            grouped = merge_df.groupby([col, "Age_Group"])[["Deaths", "Population"]].sum()
            strata = sorted(grouped.index.get_level_values(0).unique())
            for stratum in strata:
                stratum_agg = grouped.loc[stratum]
                aamr, stats = calculate_aamr_point(stratum_agg)
                se = calculate_aamr_standard_error(stratum_agg)
                lower, upper = calculate_aamr_ci(aamr, stats)
                add_result(
                    results[strat_key], s_type, str(stratum),
                    aamr, se, lower, upper, stats,
                )

    return results


def synthesize_dementia_df(input_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load G30, F01, F03 subtracted files and sum deaths per county × age group.

    Returns {period: cleaned_df} with columns COUNTY_FIPS, Age_Group, Deaths, Population.
    Population is taken from G30 (identical across ICD codes for the same county-age cell).
    """
    result = {}
    for period in _DEMENTIA_PERIODS:
        component_dfs = []
        pop_ref = None
        missing = []

        for code in _DEMENTIA_CODES:
            fp = input_dir / f"{period}_{code}.csv"
            if not fp.exists():
                missing.append(code)
                continue
            t = pd.read_csv(fp, dtype={"County Code": str})
            t = t[t["County Code"].notna() & (t["Population"] != "Missing")].copy()
            t["Deaths"] = pd.to_numeric(t["Deaths"], errors="coerce").fillna(0)
            t["Population"] = pd.to_numeric(t["Population"], errors="coerce").fillna(0)
            t = t.rename(
                columns={"County Code": "COUNTY_FIPS", "Ten-Year Age Groups": "Age_Group"}
            )
            t["COUNTY_FIPS"] = t["COUNTY_FIPS"].str.zfill(5)
            if pop_ref is None:
                pop_ref = t[["COUNTY_FIPS", "Age_Group", "Population"]].copy()
            component_dfs.append(t[["COUNTY_FIPS", "Age_Group", "Deaths"]])

        if missing:
            print(f"  ⚠ Missing dementia component file(s) for {period}: {missing}, skipping")
            continue

        combined = pd.concat(component_dfs, ignore_index=True)
        death_sum = combined.groupby(
            ["COUNTY_FIPS", "Age_Group"], as_index=False
        )["Deaths"].sum()
        merged = death_sum.merge(pop_ref, on=["COUNTY_FIPS", "Age_Group"], how="left")
        result[period] = merged

    return result


def get_top5_cancers(df: pd.DataFrame, time_period: str = "2006-2010") -> List[str]:
    """
    Identify top 5 cancer types based on total deaths in the specified period.

    Args:
        df: Stratified AAMR DataFrame
        time_period: Time period to use for ranking

    Returns:
        List of top 5 cancer ICD codes
    """
    # Filter for the specified time period and National stratum only
    period_data = df[
        (df["Time_Period"] == time_period) & (df["Stratum_Value"] == "National")
    ].copy()

    # Exclude All-site Cancer and NDD-related outcomes
    cancer_exclude = ["C00_C97", "G20_G30_G12.2_F01_F03"]
    ndd_codes = ["F01", "F03", "G10", "G12.2", "G20", "G30"]
    exclude_codes = cancer_exclude + ndd_codes

    cancer_data = period_data[~period_data["ICD-10 Code"].isin(exclude_codes)].copy()

    # Group by ICD code and outcome, sum deaths
    cancer_summary = (
        cancer_data.groupby(["ICD-10 Code", "Outcome"])["Total_Deaths"]
        .sum()
        .reset_index()
    )

    # Sort by deaths and get top 5
    cancer_summary = cancer_summary.sort_values("Total_Deaths", ascending=False)
    top5 = cancer_summary.head(5)["ICD-10 Code"].tolist()

    print(f"\nTop 5 Cancer types by deaths in {time_period}:")
    for _, row in cancer_summary.head(5).iterrows():
        print(
            f"  {row['ICD-10 Code']}: {row['Outcome']} - {row['Total_Deaths']:,} deaths"
        )

    return top5


def get_top5_ndd(df: pd.DataFrame, time_period: str = "2006-2010") -> List[str]:
    """
    Identify top 5 NDD types based on total deaths in the specified period.

    Args:
        df: Stratified AAMR DataFrame
        time_period: Time period to use for ranking

    Returns:
        List of top 5 NDD ICD codes
    """
    # Filter for the specified time period and National stratum only
    period_data = df[
        (df["Time_Period"] == time_period) & (df["Stratum_Value"] == "National")
    ].copy()

    # Include only individual NDD codes (exclude combined NDD)
    ndd_codes = ["F01", "F03", "G10", "G12.2", "G20", "G30"]
    ndd_data = period_data[period_data["ICD-10 Code"].isin(ndd_codes)].copy()

    # Group by ICD code and outcome, sum deaths
    ndd_summary = (
        ndd_data.groupby(["ICD-10 Code", "Outcome"])["Total_Deaths"].sum().reset_index()
    )

    # Sort by deaths and get top 5
    ndd_summary = ndd_summary.sort_values("Total_Deaths", ascending=False)
    top5 = ndd_summary.head(5)["ICD-10 Code"].tolist()

    print(f"\nTop 5 NDD types by deaths in {time_period}:")
    for _, row in ndd_summary.head(5).iterrows():
        print(
            f"  {row['ICD-10 Code']}: {row['Outcome']} - {row['Total_Deaths']:,} deaths"
        )

    return top5


def extract_summary_table(
    df: pd.DataFrame,
    codes: List[str],
    strat_type_filter: str,
    include_national: bool = True,
) -> pd.DataFrame:
    """
    Extract summary table for specified ICD codes from a stratified AAMR DataFrame.

    Args:
        df: Stratified AAMR DataFrame
        codes: List of ICD codes to extract
        strat_type_filter: Stratum type to extract (e.g., 'RUCC', 'koppen_major')
        include_national: Whether to include National row

    Returns:
        Summary DataFrame with Strata, Outcome, and time period columns
    """
    time_periods = ["2006-2010", "2011-2015", "2016-2020"]
    period_labels = ["5", "10", "15"]

    results = []

    # Filter by stratum type
    filtered_df = df[df["Stratum_Type"] == strat_type_filter].copy()

    # Get unique strata values
    if include_national:
        strata_values = ["National"] + sorted(
            [v for v in filtered_df["Stratum_Value"].unique() if v != "National"]
        )
    else:
        strata_values = sorted(
            [v for v in filtered_df["Stratum_Value"].unique() if v != "National"]
        )

    for stratum_val in strata_values:
        for code in codes:
            # Get data for this code and stratum
            subset = filtered_df[
                (filtered_df["ICD-10 Code"] == code)
                & (filtered_df["Stratum_Value"] == stratum_val)
            ].copy()

            if len(subset) == 0:
                continue

            # Create strata label with descriptive names
            if stratum_val == "National":
                strata_label = "National"
            elif strat_type_filter == "RUCC":
                try:
                    rucc_val = int(float(stratum_val))
                    rucc_names = {
                        1: "RUCC: Metropolitan urbanized",
                        2: "RUCC: Non-metropolitan urbanized",
                        3: "RUCC: Less urbanized",
                        4: "RUCC: Thinly populated",
                    }
                    strata_label = rucc_names.get(rucc_val, f"RUCC{rucc_val}")
                except (ValueError, TypeError):
                    strata_label = f"RUCC{stratum_val}"
            elif strat_type_filter == "koppen_major":
                koppen_names = {
                    "B": "Köppen-Geiger Climate Zone: Dry",
                    "C": "Köppen-Geiger Climate Zone: Temperate",
                    "D": "Köppen-Geiger Climate Zone: Continental",
                }
                strata_label = koppen_names.get(stratum_val, f"Climate {stratum_val}")
            elif strat_type_filter == "census_region":
                try:
                    region_val = int(float(stratum_val))
                    region_names = {
                        1: "Census Region: Northeast",
                        2: "Census Region: Midwest",
                        3: "Census Region: South",
                        4: "Census Region: West",
                    }
                    strata_label = region_names.get(
                        region_val, f"Census Region {region_val}"
                    )
                except (ValueError, TypeError):
                    strata_label = f"Census Region {stratum_val}"
            elif strat_type_filter == "doe_major":
                try:
                    doe_val = int(float(stratum_val))
                    doe_names = {
                        2: "DOE Climate Zone: Hot",
                        3: "DOE Climate Zone: Warm",
                        4: "DOE Climate Zone: Mixed",
                        5: "DOE Climate Zone: Cool",
                        6: "DOE Climate Zone: Cold",
                        7: "DOE Climate Zone: Very Cold",
                    }
                    strata_label = doe_names.get(doe_val, f"DOE {doe_val}")
                except (ValueError, TypeError):
                    strata_label = f"DOE {stratum_val}"
            elif strat_type_filter.startswith("Cluster_k"):
                k = strat_type_filter.split("_k")[1]
                try:
                    strata_label = f"Cluster_k{k} {int(float(stratum_val))}"
                except (ValueError, TypeError):
                    strata_label = f"Cluster_k{k} {stratum_val}"
            else:
                strata_label = f"{strat_type_filter} {stratum_val}"

            row_data = {"Strata": strata_label, "Outcome": subset.iloc[0]["Outcome"]}

            # Extract data for each time period
            for period, label in zip(time_periods, period_labels):
                period_data = subset[subset["Time_Period"] == period]

                if len(period_data) > 0:
                    deaths = int(period_data["Total_Deaths"].iloc[0])
                    population = int(period_data["Total_Population"].iloc[0])
                    aamr = period_data["AAMRSE"].iloc[0]

                    # Calculate per 1,000 population
                    per_thousand = (deaths / population) * 1000 if population > 0 else 0

                    row_data[f"{label} Deaths"] = f"{deaths:,} ({per_thousand:.2f}‰)"
                    row_data[f"{label} AAMR"] = aamr
                else:
                    row_data[f"{label} Deaths"] = None
                    row_data[f"{label} AAMR"] = None

            results.append(row_data)

    return pd.DataFrame(results)


def reshape_to_wide_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reshape summary table from long to wide format with outcomes as columns.

    Args:
        df: DataFrame with columns [Strata, Outcome, 5 Deaths, 5 AAMR, ...]

    Returns:
        Wide format DataFrame with strata groups and outcomes as columns
    """
    # Prepare data with strata groups
    result_rows = []

    # Get unique strata in order
    strata_order = df["Strata"].unique()

    # Get unique outcomes (should be consistent across strata)
    outcomes = df["Outcome"].unique()

    for strata in strata_order:
        # Add strata header row
        result_rows.append({"Outcome": strata})

        # Get data for this stratum
        strata_df = df[df["Strata"] == strata]

        # Add each outcome row
        for outcome in outcomes:
            outcome_df = strata_df[strata_df["Outcome"] == outcome]
            if len(outcome_df) > 0:
                row = {"Outcome": outcome}
                row["5 Deaths"] = outcome_df.iloc[0]["5 Deaths"]
                row["5 AAMR"] = outcome_df.iloc[0]["5 AAMR"]
                row["10 Deaths"] = outcome_df.iloc[0]["10 Deaths"]
                row["10 AAMR"] = outcome_df.iloc[0]["10 AAMR"]
                row["15 Deaths"] = outcome_df.iloc[0]["15 Deaths"]
                row["15 AAMR"] = outcome_df.iloc[0]["15 AAMR"]
                result_rows.append(row)

    return pd.DataFrame(result_rows)


def create_summary_tables(all_results: Dict[str, List[Dict]], output_dir: Path):
    """
    Create combined Top 5 Cancer and NDD summary tables with all strata.

    Args:
        all_results: Dictionary with RUCC, Climate, Cluster results
        output_dir: Output directory path
    """
    print("\n" + "=" * 80)
    print("Creating Summary Tables")
    print("=" * 80)

    # Define stratum type mappings
    strat_mappings = {
        "RUCC": ["RUCC"],
        "Climate": ["census_region", "koppen_major", "doe_major"],
        "Cluster": ["Cluster_k5"],  # Using k=5 clusters
    }

    all_cancer_summaries = []
    all_ndd_summaries = []

    for strat_key, strat_types in strat_mappings.items():
        if not all_results[strat_key]:
            continue

        print(f"\nProcessing {strat_key} summaries...")
        df = pd.DataFrame(all_results[strat_key])

        # Get top 5 cancers from 2006-2010 data (only once, using RUCC)
        if strat_key == "RUCC":
            top5_codes = get_top5_cancers(df, "2006-2010")
            cancer_codes = ["C00_C97"] + top5_codes

            top5_ndd_codes = get_top5_ndd(df, "2006-2010")
            ndd_codes = ["G20_G30_G12.2_F01_F03"] + top5_ndd_codes

        # Process each stratum type (RUCC has 1, Climate has 3, Cluster has 1)
        for strat_type in strat_types:
            # Extract Cancer summary for this strata type
            cancer_summary = extract_summary_table(
                df,
                cancer_codes,
                strat_type,
                include_national=(strat_key == "RUCC" and strat_type == "RUCC"),
            )
            # Skip if empty
            if cancer_summary.empty:
                continue
            # Add strata type prefix to non-National rows
            if not (strat_key == "RUCC" and strat_type == "RUCC"):
                # Add National row with strata type prefix
                national_summary = extract_summary_table(
                    df, cancer_codes, strat_type, include_national=True
                )
                if national_summary.empty or "Strata" not in national_summary.columns:
                    continue
                national_only = national_summary[
                    national_summary["Strata"] == "National"
                ].copy()
                # Use specific prefix for climate types
                if strat_key == "Climate":
                    if strat_type == "census_region":
                        prefix = "Census Region National"
                    elif strat_type == "koppen_major":
                        prefix = "Climate National"
                    elif strat_type == "doe_major":
                        prefix = "DOE National"
                    else:
                        prefix = f"{strat_type} National"
                else:
                    prefix = f"{strat_key} National"
                national_only["Strata"] = prefix
                cancer_summary = pd.concat(
                    [national_only, cancer_summary], ignore_index=True
                )

            all_cancer_summaries.append(cancer_summary)

            # Extract NDD summary for this strata type
            ndd_summary = extract_summary_table(
                df,
                ndd_codes,
                strat_type,
                include_national=(strat_key == "RUCC" and strat_type == "RUCC"),
            )
            # Skip if empty
            if ndd_summary.empty:
                continue
            # Add strata type prefix to non-National rows
            if not (strat_key == "RUCC" and strat_type == "RUCC"):
                # Add National row with strata type prefix
                national_summary = extract_summary_table(
                    df, ndd_codes, strat_type, include_national=True
                )
                if national_summary.empty or "Strata" not in national_summary.columns:
                    continue
                national_only = national_summary[
                    national_summary["Strata"] == "National"
                ].copy()
                # Use specific prefix for climate types
                if strat_key == "Climate":
                    if strat_type == "census_region":
                        prefix = "Census Region National"
                    elif strat_type == "koppen_major":
                        prefix = "Climate National"
                    elif strat_type == "doe_major":
                        prefix = "DOE National"
                    else:
                        prefix = f"{strat_type} National"
                else:
                    prefix = f"{strat_key} National"
                national_only["Strata"] = prefix
                ndd_summary = pd.concat([national_only, ndd_summary], ignore_index=True)

            all_ndd_summaries.append(ndd_summary)

    # Combine all cancer summaries
    combined_cancer = pd.concat(all_cancer_summaries, ignore_index=True)

    # Reshape to wide format
    wide_cancer = reshape_to_wide_format(combined_cancer)

    cancer_output = output_dir / "Top5_Cancer_AAMR.csv"
    wide_cancer.to_csv(cancer_output, index=False)
    print(f"\n  Saved combined Cancer summary to {cancer_output}")
    print(f"  Total shape: {wide_cancer.shape}")

    # Combine all NDD summaries
    combined_ndd = pd.concat(all_ndd_summaries, ignore_index=True)

    # Reshape to wide format
    wide_ndd = reshape_to_wide_format(combined_ndd)

    ndd_output = output_dir / "Top5_NDD_AAMR.csv"
    wide_ndd.to_csv(ndd_output, index=False)
    print(f"  Saved combined NDD summary to {ndd_output}")
    print(f"  Total shape: {wide_ndd.shape}")


def create_ndd_strata_table(all_results: Dict[str, List[Dict]], output_dir: Path):
    """
    Create NDD_AAMR_Strata.csv with all NDD outcomes across all stratifications.

    Format:
    Outcome | 2006-2010 Death,n(‰) | AAMR | 2011-2015 Death,n(‰) | AAMR | 2016-2020 Death,n(‰) | AAMR

    Each stratum appears as a header row, followed by NDD outcomes ordered by deaths.
    """
    print("\n" + "=" * 80)
    print("Creating NDD_AAMR_Strata Table")
    print("=" * 80)

    time_periods = ["2006-2010", "2011-2015", "2016-2020"]

    # Define all strata types to process with their display name prefixes
    strata_config = [
        (
            "RUCC",
            "RUCC",
            {
                1: "RUCC: Metropolitan urbanized",
                2: "RUCC: Non-metropolitan urbanized",
                3: "RUCC: Less urbanized",
                4: "RUCC: Thinly populated",
            },
        ),
        (
            "Climate",
            "census_region",
            {
                1: "Census Region: Northeast",
                2: "Census Region: Midwest",
                3: "Census Region: South",
                4: "Census Region: West",
            },
        ),
        (
            "Climate",
            "koppen_major",
            {
                "B": "Köppen-Geiger Climate Zone: Dry",
                "C": "Köppen-Geiger Climate Zone: Temperate",
                "D": "Köppen-Geiger Climate Zone: Continental",
            },
        ),
        (
            "Climate",
            "doe_major",
            {
                2: "DOE Climate Zone: Hot",
                3: "DOE Climate Zone: Warm",
                4: "DOE Climate Zone: Mixed",
                5: "DOE Climate Zone: Cool",
                6: "DOE Climate Zone: Cold",
                7: "DOE Climate Zone: Very Cold",
            },
        ),
        (
            "Typology",
            "Typology",
            {
                1: "County Economic Typology: Farming",
                2: "County Economic Typology: Mining",
                3: "County Economic Typology: Manufacturing",
                4: "County Economic Typology: Government",
                5: "County Economic Typology: Services",
                6: "County Economic Typology: Nonspecialized",
            },
        ),
        ("LandUse", "LandUse", None),
    ]

    result_rows = []

    # Add National results FIRST
    if all_results.get("RUCC"):
        df = pd.DataFrame(all_results["RUCC"])
        national_data = df[
            (df["Stratum_Value"] == "National")
            & (df["ICD-10 Code"].isin(NDD_ALL_CODES))
        ].copy()

        if not national_data.empty:
            result_rows.append({"Outcome": "National"})

            # Calculate total deaths for ordering
            code_deaths = {}
            for code in NDD_ALL_CODES:
                code_data = national_data[national_data["ICD-10 Code"] == code]
                total_deaths = code_data["Total_Deaths"].sum()
                code_deaths[code] = total_deaths

            sorted_codes = sorted(
                code_deaths.keys(), key=lambda x: code_deaths[x], reverse=True
            )

            for code in sorted_codes:
                code_data = national_data[national_data["ICD-10 Code"] == code]
                if code_data.empty:
                    continue

                abbr = NDD_ABBR_MAP.get(code, code)
                row = {"Outcome": abbr}

                for period in time_periods:
                    period_data = code_data[code_data["Time_Period"] == period]
                    if not period_data.empty:
                        deaths = int(period_data["Total_Deaths"].iloc[0])
                        population = int(period_data["Total_Population"].iloc[0])
                        aamr_se = period_data["AAMRSE"].iloc[0]
                        pct = (deaths / population) * 1000 if population > 0 else 0
                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = aamr_se
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""

                result_rows.append(row)

    # Add Sex and Race section
    demographic_rows = process_stratified_ndd_files(project_root)
    if demographic_rows:
        result_rows.append({"Outcome": "Sex and Race"})
        # Add total NDD first (from National)
        if all_results.get("RUCC"):
            df = pd.DataFrame(all_results["RUCC"])
            ndd_national = df[
                (df["Stratum_Value"] == "National")
                & (df["ICD-10 Code"] == "G20_G30_G12.2_F01_F03")
            ].copy()
            if not ndd_national.empty:
                row = {"Outcome": "NDD"}
                for period in time_periods:
                    period_data = ndd_national[ndd_national["Time_Period"] == period]
                    if not period_data.empty:
                        deaths = int(period_data["Total_Deaths"].iloc[0])
                        population = int(period_data["Total_Population"].iloc[0])
                        aamr_se = period_data["AAMRSE"].iloc[0]
                        pct = (deaths / population) * 1000 if population > 0 else 0
                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = aamr_se
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""
                result_rows.append(row)
        # Add demographic rows (already formatted with NDD(Male), etc.)
        result_rows.extend(demographic_rows)

    for strat_key, strat_type, name_map in strata_config:
        if not all_results.get(strat_key):
            continue

        df = pd.DataFrame(all_results[strat_key])

        # Filter for this stratum type and NDD codes only
        filtered = df[
            (df["Stratum_Type"] == strat_type) & (df["ICD-10 Code"].isin(NDD_ALL_CODES))
        ].copy()

        if filtered.empty:
            continue

        # Get unique strata values (excluding National)
        strata_values = [
            v for v in filtered["Stratum_Value"].unique() if v != "National"
        ]

        # Sort strata values
        try:
            strata_values = sorted(strata_values, key=lambda x: float(x))
        except (ValueError, TypeError):
            strata_values = sorted(strata_values)

        for stratum_val in strata_values:
            # Get stratum display name
            if name_map:
                try:
                    key = (
                        int(float(stratum_val))
                        if stratum_val not in ["B", "C", "D"]
                        else stratum_val
                    )
                    stratum_name = name_map.get(key, f"{strat_type}: {stratum_val}")
                except (ValueError, TypeError):
                    stratum_name = name_map.get(
                        stratum_val, f"{strat_type}: {stratum_val}"
                    )
            else:
                # For Cluster, Typology, LandUse without predefined names
                if strat_type.startswith("Cluster_k"):
                    stratum_name = f"EQI {strat_type}: {stratum_val}"
                elif strat_type == "Typology":
                    stratum_name = f"County Economic Typology: {stratum_val}"
                elif strat_type == "LandUse":
                    stratum_name = f"Land Use Cluster: {stratum_val}"
                else:
                    stratum_name = f"{strat_type}: {stratum_val}"

            # Add stratum header row
            result_rows.append({"Outcome": stratum_name})

            # Get data for this stratum
            stratum_data = filtered[filtered["Stratum_Value"] == stratum_val].copy()

            # Calculate total deaths across all periods for ordering
            code_deaths = {}
            for code in NDD_ALL_CODES:
                code_data = stratum_data[stratum_data["ICD-10 Code"] == code]
                total_deaths = code_data["Total_Deaths"].sum()
                code_deaths[code] = total_deaths

            # Sort codes by total deaths (descending)
            sorted_codes = sorted(
                code_deaths.keys(), key=lambda x: code_deaths[x], reverse=True
            )

            # Add rows for each NDD outcome
            for code in sorted_codes:
                code_data = stratum_data[stratum_data["ICD-10 Code"] == code]
                if code_data.empty:
                    continue

                abbr = NDD_ABBR_MAP.get(code, code)
                row = {"Outcome": abbr}

                for period in time_periods:
                    period_data = code_data[code_data["Time_Period"] == period]

                    if not period_data.empty:
                        deaths = int(period_data["Total_Deaths"].iloc[0])
                        population = int(period_data["Total_Population"].iloc[0])
                        aamr_se = period_data["AAMRSE"].iloc[0]

                        # Calculate per mille (deaths / population * 1000)
                        pct = (deaths / population) * 1000 if population > 0 else 0

                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = aamr_se
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""

                result_rows.append(row)

    # Create DataFrame and save
    output_df = pd.DataFrame(result_rows)

    # Reorder columns
    col_order = ["Outcome"]
    for period in time_periods:
        col_order.append(f"{period} Death,n(‰)")
        col_order.append(f"{period} AAMR")

    # Only include columns that exist
    col_order = [c for c in col_order if c in output_df.columns]
    output_df = output_df[col_order]

    output_path = output_dir / "NDD_AAMR_Strata.csv"
    output_df.to_csv(output_path, index=False)
    print(f"\nSaved NDD_AAMR_Strata to {output_path}")
    print(f"Total rows: {len(output_df)}")

    return result_rows


def process_stratified_ndd_files_sex_only(
    project_root: Path, all_results: Dict[str, List[Dict]] = None
) -> List[Dict]:
    """
    Process sex stratified NDD files from Stratified_Subtracted directory.
    Male deaths are calculated as Total - Female to ensure consistency.

    Returns list of result rows with _type field for identification.
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return []

    time_periods = ["2006-2010", "2011-2015", "2016-2020"]

    # Get national totals from all_results
    national_totals = {}
    if all_results and all_results.get("RUCC"):
        df = pd.DataFrame(all_results["RUCC"])
        ndd_national = df[
            (df["Stratum_Value"] == "National")
            & (df["ICD-10 Code"] == "G20_G30_G12.2_F01_F03")
        ].copy()
        for period in time_periods:
            period_data = ndd_national[ndd_national["Time_Period"] == period]
            if not period_data.empty:
                national_totals[period] = {
                    "deaths": int(period_data["Total_Deaths"].iloc[0]),
                    "population": int(period_data["Total_Population"].iloc[0]),
                }

    # First get Female data
    female_data = {}
    for period in time_periods:
        file_path = input_dir / f"{period}_G20_G30_G12.2_F01_F03_Female.csv"
        if file_path.exists():
            df = pd.read_csv(file_path, dtype={"County Code": str})
            df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
            df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
            df["Population"] = pd.to_numeric(df["Population"], errors="coerce").fillna(
                0
            )
            df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})

            agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
            aamr, stats = calculate_aamr_point(agg_df)
            se = calculate_aamr_standard_error(agg_df)

            female_data[period] = {
                "deaths": stats["total_deaths"],
                "population": stats["total_population"],
                "aamr": aamr,
                "se": se,
            }

    result_rows = []

    # Calculate Male as Total - Female
    male_row = {"_type": "male"}
    for period in time_periods:
        if period in national_totals and period in female_data:
            total_deaths = national_totals[period]["deaths"]
            female_deaths = female_data[period]["deaths"]
            male_deaths = total_deaths - female_deaths

            # For population and AAMR, use Male file data
            file_path = input_dir / f"{period}_G20_G30_G12.2_F01_F03_Male.csv"
            if file_path.exists():
                df = pd.read_csv(file_path, dtype={"County Code": str})
                df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
                df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
                df["Population"] = pd.to_numeric(
                    df["Population"], errors="coerce"
                ).fillna(0)
                df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})

                agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
                aamr, stats = calculate_aamr_point(agg_df)
                se = calculate_aamr_standard_error(agg_df)

                population = stats["total_population"]
                pct = (male_deaths / population) * 1000 if population > 0 else 0
                male_row[f"{period} Death,n(‰)"] = f"{male_deaths:,} ({pct:.2f}‰)"
                male_row[f"{period} AAMR"] = f"{aamr:.2f} ± {se:.2f}"

    result_rows.append(male_row)

    # Add Female row
    female_row = {"_type": "female"}
    for period, data in female_data.items():
        deaths = data["deaths"]
        population = data["population"]
        aamr = data["aamr"]
        se = data["se"]

        pct = (deaths / population) * 1000 if population > 0 else 0
        female_row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
        female_row[f"{period} AAMR"] = f"{aamr:.2f} ± {se:.2f}"

    result_rows.append(female_row)

    return result_rows


def process_stratified_ndd_files_race_only(project_root: Path) -> List[Dict]:
    """
    Process race stratified NDD files from Stratified_Subtracted directory.

    Returns list of result rows for NDD by race strata with format:
    NDD (White), NDD (Black), NDD (Asian), NDD (Indian)
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return []

    time_periods = ["2006-2010", "2011-2015", "2016-2020"]

    race_strata = [
        ("White", "NDD (White)"),
        ("Black", "NDD (Black)"),
        ("Asian", "NDD (Asian)"),
        ("Indian", "NDD (Indian)"),
    ]

    result_rows = []

    for stratum_file, display_name in race_strata:
        stratum_data = []
        for period in time_periods:
            file_path = input_dir / f"{period}_G20_G30_G12.2_F01_F03_{stratum_file}.csv"
            if file_path.exists():
                df = pd.read_csv(file_path, dtype={"County Code": str})
                df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
                df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
                df["Population"] = pd.to_numeric(
                    df["Population"], errors="coerce"
                ).fillna(0)
                df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})

                agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
                aamr, stats = calculate_aamr_point(agg_df)
                se = calculate_aamr_standard_error(agg_df)

                stratum_data.append(
                    {
                        "period": period,
                        "deaths": stats["total_deaths"],
                        "population": stats["total_population"],
                        "aamr": aamr,
                        "se": se,
                    }
                )

        if stratum_data:
            row = {"Outcome": display_name}
            for data in stratum_data:
                period = data["period"]
                deaths = data["deaths"]
                population = data["population"]
                aamr = data["aamr"]
                se = data["se"]

                pct = (deaths / population) * 1000 if population > 0 else 0
                row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                row[f"{period} AAMR"] = f"{aamr:.2f} ± {se:.2f}"

            result_rows.append(row)

    return result_rows


def create_ndd_strata_table_with_demographics(
    all_results: Dict[str, List[Dict]], project_root: Path, output_dir: Path
):
    """
    Create NDD_AAMR_Strata.csv with all NDD outcomes including race/sex strata.
    """
    print("\n" + "=" * 80)
    print("Creating NDD_AAMR_Strata Table (with Race/Sex)")
    print("=" * 80)

    time_periods = ["2006-2010", "2011-2015", "2016-2020"]

    # Define all strata types to process with their display name prefixes
    strata_config = [
        (
            "RUCC",
            "RUCC",
            {
                1: "RUCC: Metropolitan urbanized",
                2: "RUCC: Non-metropolitan urbanized",
                3: "RUCC: Less urbanized",
                4: "RUCC: Thinly populated",
            },
        ),
        (
            "Climate",
            "census_region",
            {
                1: "Census Region: Northeast",
                2: "Census Region: Midwest",
                3: "Census Region: South",
                4: "Census Region: West",
            },
        ),
        (
            "Climate",
            "koppen_major",
            {
                "B": "Köppen-Geiger Climate Zone: Dry",
                "C": "Köppen-Geiger Climate Zone: Temperate",
                "D": "Köppen-Geiger Climate Zone: Continental",
            },
        ),
        (
            "Climate",
            "doe_major",
            {
                2: "DOE Climate Zone: Hot",
                3: "DOE Climate Zone: Warm",
                4: "DOE Climate Zone: Mixed",
                5: "DOE Climate Zone: Cool",
                6: "DOE Climate Zone: Cold",
                7: "DOE Climate Zone: Very Cold",
            },
        ),
        (
            "Typology",
            "Typology",
            {
                1: "County Economic Typology: Farming",
                2: "County Economic Typology: Mining",
                3: "County Economic Typology: Manufacturing",
                4: "County Economic Typology: Government",
                5: "County Economic Typology: Services",
                6: "County Economic Typology: Nonspecialized",
            },
        ),
        ("LandUse", "LandUse", None),
    ]

    result_rows = []

    # Add National results FIRST with specific order
    if all_results.get("RUCC"):
        df = pd.DataFrame(all_results["RUCC"])
        national_data = df[
            (df["Stratum_Value"] == "National")
            & (df["ICD-10 Code"].isin(NDD_ALL_CODES))
        ].copy()

        if not national_data.empty:
            result_rows.append({"Outcome": "National"})

            # Get sex stratified data
            sex_rows = process_stratified_ndd_files_sex_only(project_root, all_results)

            # Define specific order: NDD, NDD (Male), NDD (Female), UD, AD, PD, VD, ALS, HD
            national_order = [
                ("G20_G30_G12.2_F01_F03", "NDD"),
                ("sex_male", "NDD (Male)"),
                ("sex_female", "NDD (Female)"),
                ("F03", "UD"),
                ("G30", "AD"),
                ("G20", "PD"),
                ("F01", "VD"),
                ("G12.2", "ALS"),
                ("G10", "HD"),
            ]

            for code, display_name in national_order:
                if code.startswith("sex_"):
                    # Get from sex_rows
                    sex_type = code.replace("sex_", "")
                    for sex_row in sex_rows:
                        if sex_row.get("_type") == sex_type:
                            row = {
                                k: v
                                for k, v in sex_row.items()
                                if not k.startswith("_")
                            }
                            row["Outcome"] = display_name
                            result_rows.append(row)
                            break
                else:
                    code_data = national_data[national_data["ICD-10 Code"] == code]
                    if code_data.empty:
                        continue

                    row = {"Outcome": display_name}

                    for period in time_periods:
                        period_data = code_data[code_data["Time_Period"] == period]
                        if not period_data.empty:
                            deaths = int(period_data["Total_Deaths"].iloc[0])
                            population = int(period_data["Total_Population"].iloc[0])
                            aamr_se = period_data["AAMRSE"].iloc[0]
                            pct = (deaths / population) * 1000 if population > 0 else 0
                            row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                            row[f"{period} AAMR"] = aamr_se
                        else:
                            row[f"{period} Death,n(‰)"] = ""
                            row[f"{period} AAMR"] = ""

                    result_rows.append(row)

    # Add Race section
    race_rows = process_stratified_ndd_files_race_only(project_root)
    if race_rows:
        result_rows.append({"Outcome": "Race"})
        # Add total NDD first (from National)
        if all_results.get("RUCC"):
            df = pd.DataFrame(all_results["RUCC"])
            ndd_national = df[
                (df["Stratum_Value"] == "National")
                & (df["ICD-10 Code"] == "G20_G30_G12.2_F01_F03")
            ].copy()
            if not ndd_national.empty:
                row = {"Outcome": "NDD"}
                for period in time_periods:
                    period_data = ndd_national[ndd_national["Time_Period"] == period]
                    if not period_data.empty:
                        deaths = int(period_data["Total_Deaths"].iloc[0])
                        population = int(period_data["Total_Population"].iloc[0])
                        aamr_se = period_data["AAMRSE"].iloc[0]
                        pct = (deaths / population) * 1000 if population > 0 else 0
                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = aamr_se
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""
                result_rows.append(row)
        # Add race rows
        result_rows.extend(race_rows)

    for strat_key, strat_type, name_map in strata_config:
        if not all_results.get(strat_key):
            continue

        df = pd.DataFrame(all_results[strat_key])

        # Filter for this stratum type and NDD codes only
        filtered = df[
            (df["Stratum_Type"] == strat_type) & (df["ICD-10 Code"].isin(NDD_ALL_CODES))
        ].copy()

        if filtered.empty:
            continue

        # Get unique strata values (excluding National)
        strata_values = [
            v for v in filtered["Stratum_Value"].unique() if v != "National"
        ]

        # Sort strata values
        try:
            strata_values = sorted(strata_values, key=lambda x: float(x))
        except (ValueError, TypeError):
            strata_values = sorted(strata_values)

        for stratum_val in strata_values:
            # Get stratum display name
            if name_map:
                try:
                    key = (
                        int(float(stratum_val))
                        if stratum_val not in ["B", "C", "D"]
                        else stratum_val
                    )
                    stratum_name = name_map.get(key, f"{strat_type}: {stratum_val}")
                except (ValueError, TypeError):
                    stratum_name = name_map.get(
                        stratum_val, f"{strat_type}: {stratum_val}"
                    )
            else:
                # For Cluster, Typology, LandUse without predefined names
                if strat_type.startswith("Cluster_k"):
                    stratum_name = f"EQI {strat_type}: {stratum_val}"
                elif strat_type == "Typology":
                    stratum_name = f"County Economic Typology: {stratum_val}"
                elif strat_type == "LandUse":
                    stratum_name = f"Land Use Cluster: {stratum_val}"
                else:
                    stratum_name = f"{strat_type}: {stratum_val}"

            # Add stratum header row
            result_rows.append({"Outcome": stratum_name})

            # Get data for this stratum
            stratum_data = filtered[filtered["Stratum_Value"] == stratum_val].copy()

            # Calculate total deaths across all periods for ordering
            code_deaths = {}
            for code in NDD_ALL_CODES:
                code_data = stratum_data[stratum_data["ICD-10 Code"] == code]
                total_deaths = code_data["Total_Deaths"].sum()
                code_deaths[code] = total_deaths

            # Sort codes by total deaths (descending)
            sorted_codes = sorted(
                code_deaths.keys(), key=lambda x: code_deaths[x], reverse=True
            )

            # Add rows for each NDD outcome
            for code in sorted_codes:
                code_data = stratum_data[stratum_data["ICD-10 Code"] == code]
                if code_data.empty:
                    continue

                abbr = NDD_ABBR_MAP.get(code, code)
                row = {"Outcome": abbr}

                for period in time_periods:
                    period_data = code_data[code_data["Time_Period"] == period]

                    if not period_data.empty:
                        deaths = int(period_data["Total_Deaths"].iloc[0])
                        population = int(period_data["Total_Population"].iloc[0])
                        aamr_se = period_data["AAMRSE"].iloc[0]

                        # Calculate per mille (deaths / population * 1000)
                        pct = (deaths / population) * 1000 if population > 0 else 0

                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = aamr_se
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""

                result_rows.append(row)

    # Create DataFrame and save
    output_df = pd.DataFrame(result_rows)

    # Reorder columns
    col_order = ["Outcome"]
    for period in time_periods:
        col_order.append(f"{period} Death,n(‰)")
        col_order.append(f"{period} AAMR")

    # Only include columns that exist
    col_order = [c for c in col_order if c in output_df.columns]
    output_df = output_df[col_order]

    output_path = output_dir / "NDD_AAMR_Strata.csv"
    output_df.to_csv(output_path, index=False)
    print(f"\nSaved NDD_AAMR_Strata to {output_path}")
    print(f"Total rows: {len(output_df)}")


def main():
    project_root, config = load_config()
    input_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"
    output_dir = project_root / "Result/Tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    strat_data = load_stratification_data(project_root, config)

    subtracted_files = list(input_dir.glob("*.csv"))
    all_results = {
        "RUCC": [],
        "Climate": [],
        "Cluster": [],
        "Typology": [],
        "LandUse": [],
    }

    for file_path in sorted(subtracted_files):
        file_results = process_file(file_path, strat_data, config)
        for key in all_results:
            if key in file_results:
                all_results[key].extend(file_results[key])

    # Synthesize Dementia (G30+F01+F03) stratified AAMR
    print("\nSynthesizing Dementia (G30+F01+F03) stratified AAMR...")
    dementia_dfs = synthesize_dementia_df(input_dir)
    for period, dem_df in dementia_dfs.items():
        dem_results = _process_clean_df(dem_df, period, _DEMENTIA_ICD, strat_data, config)
        for key in all_results:
            if key in dem_results:
                all_results[key].extend(dem_results[key])

    # Write stratified AAMR outputs
    for strat_type in all_results:
        if all_results[strat_type]:
            df = pd.DataFrame(all_results[strat_type])
            df = df.sort_values(
                ["Time_Period", "ICD-10 Code", "Stratum_Type", "Stratum_Value"]
            )

            # Apply stratification labels
            df = apply_stratification_labels(df, config)

            # Reorder columns to put Stratum_Label after Stratum_Value
            cols = df.columns.tolist()
            if "Stratum_Label" in cols:
                cols.remove("Stratum_Label")
                value_idx = cols.index("Stratum_Value")
                cols.insert(value_idx + 1, "Stratum_Label")
                df = df[cols]

            output_path = output_dir / f"Stratified_AAMR_{strat_type}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved {strat_type} stratified AAMR to {output_path}")

    print("\nStratified AAMR calculation completed.")

    # Create summary tables
    create_summary_tables(all_results, output_dir)

    # Create NDD strata table (with race/sex demographics)
    create_ndd_strata_table_with_demographics(all_results, project_root, output_dir)

    print("\n" + "=" * 80)
    print("All processing completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
