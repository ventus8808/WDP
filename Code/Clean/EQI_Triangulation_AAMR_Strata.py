"""
Stratified AAMR Calculator for RUCC, Climate, and Cluster

This script calculates national and stratum-level Age-Adjusted Mortality Rates (AAMRs)
with confidence intervals using triangulated death data, stratified by RUCC, climate zones,
and EQI clusters.

Stratifications:
- RUCC: Rural-Urban Continuum Codes from EQI data
- Climate: census_region, koppen_major, doe_major
- Cluster: cluster_3, cluster_4, cluster_5

Input: Data/Original/CDC Triangulation/Subtracted/*.csv
Output:
    - Result/Tables/Stratified_AAMR_{RUCC,Climate,Cluster}.csv
    - Result/Tables/Top5_Cancer_AAMR_{RUCC,Climate,Cluster}.csv (summary tables)
    - Result/Tables/Top5_NDD_AAMR_{RUCC,Climate,Cluster}.csv (summary tables)
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


def load_config() -> Tuple[Path, Dict]:
    """Load configuration from config.yaml"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return project_root, config


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

    return {"RUCC": rucc_df, "Climate": climate_df, "Cluster": cluster_df}


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

    results = {"RUCC": [], "Climate": [], "Cluster": []}

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
    ]:
        if strat_data[strat_key].empty:
            continue
        merge_df = df.merge(strat_data[strat_key], on="COUNTY_FIPS", how="left")
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
            else:  # Cluster
                k = col.split("_")[1]
                s_type = f"Cluster_k{k}"

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
            # Add strata type prefix to non-National rows
            if not (strat_key == "RUCC" and strat_type == "RUCC"):
                # Add National row with strata type prefix
                national_summary = extract_summary_table(
                    df, cancer_codes, strat_type, include_national=True
                )
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
            # Add strata type prefix to non-National rows
            if not (strat_key == "RUCC" and strat_type == "RUCC"):
                # Add National row with strata type prefix
                national_summary = extract_summary_table(
                    df, ndd_codes, strat_type, include_national=True
                )
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


def main():
    project_root, config = load_config()
    input_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"
    output_dir = project_root / "Result/Tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    strat_data = load_stratification_data(project_root, config)

    subtracted_files = list(input_dir.glob("*.csv"))
    all_results = {"RUCC": [], "Climate": [], "Cluster": []}

    for file_path in sorted(subtracted_files):
        file_results = process_file(file_path, strat_data, config)
        for key in all_results:
            if key in file_results:
                all_results[key].extend(file_results[key])

    # Write stratified AAMR outputs
    for strat_type in all_results:
        if all_results[strat_type]:
            df = pd.DataFrame(all_results[strat_type])
            df = df.sort_values(
                ["Time_Period", "ICD-10 Code", "Stratum_Type", "Stratum_Value"]
            )
            output_path = output_dir / f"Stratified_AAMR_{strat_type}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved {strat_type} stratified AAMR to {output_path}")

    print("\nStratified AAMR calculation completed.")

    # Create summary tables
    create_summary_tables(all_results, output_dir)

    print("\n" + "=" * 80)
    print("All processing completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
