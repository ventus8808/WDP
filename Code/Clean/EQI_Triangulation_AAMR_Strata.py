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
Output: Result/Tables/Stratified_AAMR_{RUCC,Climate,Cluster}.csv
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

    # Write outputs
    for strat_type in all_results:
        if all_results[strat_type]:
            df = pd.DataFrame(all_results[strat_type])
            df = df.sort_values(
                ["Time_Period", "ICD-10 Code", "Stratum_Type", "Stratum_Value"]
            )
            output_path = output_dir / f"Stratified_AAMR_{strat_type}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved {strat_type} stratified AAMR to {output_path}")

    print("Stratified AAMR calculation completed.")


if __name__ == "__main__":
    main()
