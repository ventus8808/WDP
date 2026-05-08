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

CLIMATE_COLS = ["census_region", "koppen_major"]

TYPOLOGY_COLS = ["econdep"]

# Populated from config in main(); kept as module-level so all functions share them
NDD_ABBR_MAP: Dict[str, str] = {}
NDD_ALL_CODES: List[str] = []
ICD_SHORTNAME_MAP: Dict[str, str] = {}


def _build_shortname_map(config: Dict) -> Dict[str, str]:
    """Build {icd_code: shortname} for every disease entry in config."""
    mapping: Dict[str, str] = {}
    for group_cfg in config.get("diseases", {}).values():
        overall = group_cfg.get("overall", {})
        code = overall.get("icd_code")
        if code:
            mapping[code] = overall.get("shortname", code)
        for subtype in group_cfg.get("subtypes", {}).values():
            sub_code = subtype.get("icd_code")
            if sub_code:
                mapping[sub_code] = subtype.get("shortname", sub_code)
    return mapping


def _build_ndd_maps(config: Dict) -> tuple:
    """Build NDD_ABBR_MAP and NDD_ALL_CODES from diseases config."""
    ndd = config.get("diseases", {}).get("ndd", {})
    abbr_map: Dict[str, str] = {}
    all_codes: List[str] = []

    overall = ndd.get("overall", {})
    if overall.get("icd_code"):
        abbr_map[overall["icd_code"]] = overall.get("shortname", overall["icd_code"])
        all_codes.append(overall["icd_code"])

    for subtype in ndd.get("subtypes", {}).values():
        code = subtype.get("icd_code")
        if code:
            abbr_map.setdefault(code, subtype.get("shortname", code))
            if code not in all_codes:
                all_codes.append(code)

    return abbr_map, all_codes


def load_config() -> Tuple[Path, Dict]:
    """Load configuration from config.yaml"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return project_root, config


def apply_stratification_labels(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    """Add Stratification (category name) and Stratum_Label (sub-label) columns."""
    if df.empty:
        return df

    df = df.copy()

    strat_type_to_map: Dict = {
        strat_type: name_map
        for _, strat_type, name_map in STRATA_DISPLAY_CONFIG
        if name_map is not None
    }

    stratification_col = []
    label_col = []
    for _, row in df.iterrows():
        st = str(row["Stratum_Type"])
        sv = str(row["Stratum_Value"])
        stratification_col.append(STRAT_TYPE_TO_CATEGORY.get(st, st))
        label_col.append(
            "National"
            if sv == "National"
            else _get_stratum_name(st, sv, strat_type_to_map.get(st))
        )

    df["Stratification"] = stratification_col
    df["Stratum_Label"] = label_col
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
        get_path(["data_directories", "processed"], "Data/Processed")
        / "Covariate"
        / "Climate_Zone.csv"
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

    # EQI quintile assignments for both periods
    eqi_dir = project_root / config.get("data_sources", {}).get("epa_eqi", {}).get(
        "processed", "Data/Processed/EQI"
    )
    eqi_0005_df = pd.DataFrame()
    fp = eqi_dir / "EQI0005.csv"
    if fp.exists():
        t = pd.read_csv(fp, dtype={"COUNTY_FIPS": str})
        t["COUNTY_FIPS"] = t["COUNTY_FIPS"].str.zfill(5)
        if "EQI" in t.columns:
            eqi_0005_df = t[["COUNTY_FIPS", "EQI"]].copy()
    else:
        print(f"Warning: EQI0005.csv not found at {fp}")

    return {
        "RUCC": rucc_df,
        "Climate": climate_df,
        "Cluster": cluster_df,
        "Typology": typology_df,
        "EQI_0005": eqi_0005_df,
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

    results = {"RUCC": [], "Climate": [], "Cluster": [], "Typology": [], "EQI": []}

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
                "Outcome_Fullname": config.get("brms_analysis", {})
                .get("icd_mapping", {})
                .get(icd_code_norm, icd_code),
                "Outcome_Shortname": ICD_SHORTNAME_MAP.get(icd_code, icd_code_norm),
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

    # EQI quintile stratification — always EQI 2000-2005 (single exposure period, 3 lags)
    eqi_df = strat_data.get("EQI_0005", pd.DataFrame())

    if not eqi_df.empty:
        eqi_merged = df.merge(eqi_df, on="COUNTY_FIPS", how="left")
        eqi_merged = eqi_merged.dropna(subset=["EQI"])

        nat_agg = eqi_merged.groupby("Age_Group")[["Deaths", "Population"]].sum()
        nat_aamr, nat_stats = calculate_aamr_point(nat_agg)
        nat_se = calculate_aamr_standard_error(nat_agg)
        nat_lower, nat_upper = calculate_aamr_ci(nat_aamr, nat_stats)
        add_result(
            results["EQI"],
            "EQI",
            "National",
            nat_aamr,
            nat_se,
            nat_lower,
            nat_upper,
            nat_stats,
        )

        grouped = eqi_merged.groupby(["EQI", "Age_Group"])[
            ["Deaths", "Population"]
        ].sum()
        for q in sorted(eqi_merged["EQI"].dropna().unique()):
            q_int = int(q)
            try:
                q_agg = grouped.loc[q_int]
            except KeyError:
                continue
            aamr, stats = calculate_aamr_point(q_agg)
            se = calculate_aamr_standard_error(q_agg)
            lower, upper = calculate_aamr_ci(aamr, stats)
            add_result(
                results["EQI"], "EQI", f"Q{q_int}", aamr, se, lower, upper, stats
            )

    return results


def _build_disease_groups(config: Dict) -> Dict:
    """Build ordered display info for each disease group from diseases config."""
    groups: Dict = {}
    for group_key, group_cfg in config.get("diseases", {}).items():
        overall = group_cfg.get("overall", {})
        overall_code = overall.get("icd_code", "")
        shortname = overall.get("shortname", group_key.upper())
        abbr_map: Dict[str, str] = {}
        display_codes: List[str] = []
        hidden_codes: set = set()
        if overall_code:
            abbr_map[overall_code] = shortname
            display_codes.append(overall_code)
            if overall.get("show") is False:
                hidden_codes.add(overall_code)
        for subtype in group_cfg.get("subtypes", {}).values():
            code = subtype.get("icd_code")
            if code and subtype.get("summarize", True):
                abbr_map[code] = subtype.get("shortname", code)
                display_codes.append(code)
                if subtype.get("show") is False:
                    hidden_codes.add(code)
        display_codes_simplified = [c for c in display_codes if c not in hidden_codes]
        groups[group_key] = {
            "shortname": shortname,
            "overall_code": overall_code,
            "display_codes": display_codes,
            "display_codes_simplified": display_codes_simplified,
            "hidden_codes": hidden_codes,
            "abbr_map": abbr_map,
        }
    return groups


# Maps raw Stratum_Type values to human-readable category names (used as "Stratification" column)
STRAT_TYPE_TO_CATEGORY: Dict[str, str] = {
    "RUCC": "RUCC",
    "census_region": "Census Region",
    "koppen_major": "Köppen-Geiger Climate Zone",
    "Typology": "County Economic Typology",
    "EQI": "EQI Quintile",
}

# Strata display order and label maps — values are sub-labels only (no category prefix)
STRATA_DISPLAY_CONFIG = [
    (
        "RUCC",
        "RUCC",
        {
            1: "Metropolitan urbanized",
            2: "Non-metropolitan urbanized",
            3: "Less urbanized",
            4: "Thinly populated",
        },
    ),
    (
        "Climate",
        "census_region",
        {
            1: "Northeast",
            2: "Midwest",
            3: "South",
            4: "West",
        },
    ),
    (
        "Climate",
        "koppen_major",
        {
            "B": "Dry",
            "C": "Temperate",
            "D": "Continental",
        },
    ),
    (
        "Typology",
        "Typology",
        {
            1: "Farming",
            2: "Mining",
            3: "Manufacturing",
            4: "Government",
            5: "Services",
            6: "Nonspecialized",
        },
    ),
    (
        "EQI",
        "EQI",
        {
            "Q1": "Q1 (best)",
            "Q2": "Q2",
            "Q3": "Q3",
            "Q4": "Q4",
            "Q5": "Q5 (worst)",
        },
    ),
]


def _get_stratum_name(
    strat_type: str, stratum_val: str, name_map: Optional[Dict]
) -> str:
    """Convert raw stratum value to display name."""
    if name_map is None:
        return f"{strat_type}: {stratum_val}"
    try:
        key = (
            int(float(stratum_val))
            if stratum_val not in ("B", "C", "D")
            else stratum_val
        )
        return name_map.get(key, f"{strat_type}: {stratum_val}")
    except (ValueError, TypeError):
        return name_map.get(stratum_val, f"{strat_type}: {stratum_val}")


def create_disease_aamr_table(
    group_key: str,
    group_info: Dict,
    all_results: Dict[str, List[Dict]],
    project_root: Path,
    output_dir: Path,
) -> None:
    """
    Create AAMR_{shortname}.csv for one disease group across all stratifications.

    Format: stratum header row followed by disease outcome rows, one block per stratum.
    Within each stratum block, diseases are sorted by their 2006-2010 AAMR descending.
    For the ndd group, a Sex / Race section is inserted after National.

    If any disease has show=false in config, a second _Simplified file is written
    that omits those diseases. The full file always includes all display_codes.
    """
    shortname = group_info["shortname"]
    display_codes = group_info["display_codes"]
    display_codes_simplified = group_info.get("display_codes_simplified", display_codes)
    hidden_codes = group_info.get("hidden_codes", set())
    abbr_map = group_info["abbr_map"]
    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]

    print(f"Creating AAMR_{shortname}.csv")

    def build_rows(codes_to_use: List[str]) -> List[Dict]:
        rows: List[Dict] = []

        def add_stratum_section(header: str, stratum_df: pd.DataFrame) -> None:
            def get_aamr_2006(code: str) -> float:
                mask = (stratum_df["ICD-10 Code"] == code) & (
                    stratum_df["Time_Period"] == "2006-2010"
                )
                period_data = stratum_df[mask]
                return float(period_data["AAMR"].iloc[0]) if not period_data.empty else -1.0

            sorted_codes = sorted(codes_to_use, key=get_aamr_2006, reverse=True)
            rows.append({"Outcome": header})
            for code in sorted_codes:
                code_data = stratum_df[stratum_df["ICD-10 Code"] == code]
                if code_data.empty:
                    continue
                row: Dict = {"Outcome": abbr_map.get(code, code)}
                for period in time_periods:
                    period_row = code_data[code_data["Time_Period"] == period]
                    if not period_row.empty:
                        deaths = int(period_row["Total_Deaths"].iloc[0])
                        population = int(period_row["Total_Population"].iloc[0])
                        pct = (deaths / population) * 1000 if population > 0 else 0
                        row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                        row[f"{period} AAMR"] = period_row["AAMRSE"].iloc[0]
                    else:
                        row[f"{period} Death,n(‰)"] = ""
                        row[f"{period} AAMR"] = ""
                rows.append(row)

        # Determine code sort order once (by 2006-2010 national AAMR) for reuse below
        sorted_codes: List[str] = list(codes_to_use)
        national_df: pd.DataFrame = pd.DataFrame()
        if all_results.get("RUCC"):
            nat_df = pd.DataFrame(all_results["RUCC"])
            national_df = nat_df[
                (nat_df["Stratum_Value"] == "National")
                & (nat_df["ICD-10 Code"].isin(codes_to_use))
            ]
            if not national_df.empty:
                def _get_nat_aamr_2006(code: str) -> float:
                    mask = (national_df["ICD-10 Code"] == code) & (national_df["Time_Period"] == "2006-2010")
                    d = national_df[mask]
                    return float(d["AAMR"].iloc[0]) if not d.empty else -1.0
                sorted_codes = sorted(codes_to_use, key=_get_nat_aamr_2006, reverse=True)

        def _make_period_cells(code_data: pd.DataFrame, row_out: Dict) -> None:
            for period in time_periods:
                period_row = code_data[code_data["Time_Period"] == period]
                if not period_row.empty:
                    deaths = int(period_row["Total_Deaths"].iloc[0])
                    population = int(period_row["Total_Population"].iloc[0])
                    pct = (deaths / population) * 1000 if population > 0 else 0
                    row_out[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
                    row_out[f"{period} AAMR"] = period_row["AAMRSE"].iloc[0]
                else:
                    row_out[f"{period} Death,n(‰)"] = ""
                    row_out[f"{period} AAMR"] = ""

        # National section: disease rows with inline (M)/(F) after each code that has sex data
        if not national_df.empty:
            sex_lookup: Dict[str, pd.DataFrame] = {}
            if all_results.get("Sex"):
                sex_all = pd.DataFrame(all_results["Sex"])
                for sv in ("Male", "Female"):
                    sex_lookup[sv] = sex_all[
                        (sex_all["Stratum_Value"] == sv)
                        & (sex_all["ICD-10 Code"].isin(codes_to_use))
                    ]

            rows.append({"Outcome": "National"})
            for code in sorted_codes:
                code_data = national_df[national_df["ICD-10 Code"] == code]
                if code_data.empty:
                    continue
                row: Dict = {"Outcome": abbr_map.get(code, code)}
                _make_period_cells(code_data, row)
                rows.append(row)

                for sv, suffix in (("Male", "(M)"), ("Female", "(F)")):
                    if sv not in sex_lookup:
                        continue
                    sex_code = sex_lookup[sv][sex_lookup[sv]["ICD-10 Code"] == code]
                    if sex_code.empty:
                        continue
                    sex_row: Dict = {"Outcome": abbr_map.get(code, code) + suffix}
                    _make_period_cells(sex_code, sex_row)
                    rows.append(sex_row)

        # Race section: flat rows named "NDD(White)", "NDD(Black)", "NDD(Others)"
        if all_results.get("Race"):
            race_df = pd.DataFrame(all_results["Race"])
            race_filtered = race_df[race_df["ICD-10 Code"].isin(codes_to_use)]
            if not race_filtered.empty:
                rows.append({"Outcome": "Race"})
                for code in sorted_codes:
                    for race_val, suffix in (("White", "(White)"), ("Black", "(Black)"), ("Others", "(Others)")):
                        race_code = race_filtered[
                            (race_filtered["ICD-10 Code"] == code)
                            & (race_filtered["Stratum_Value"] == race_val)
                        ]
                        if race_code.empty:
                            continue
                        race_row: Dict = {"Outcome": abbr_map.get(code, code) + suffix}
                        _make_period_cells(race_code, race_row)
                        rows.append(race_row)

        # All stratifications
        for strat_key, strat_type, name_map in STRATA_DISPLAY_CONFIG:
            if not all_results.get(strat_key):
                continue
            df = pd.DataFrame(all_results[strat_key])
            filtered = df[
                (df["Stratum_Type"] == strat_type) & (df["ICD-10 Code"].isin(codes_to_use))
            ]
            if filtered.empty:
                continue
            strata_values = [
                v for v in filtered["Stratum_Value"].unique() if v != "National"
            ]
            try:
                strata_values = sorted(strata_values, key=lambda x: float(x))
            except (ValueError, TypeError):
                strata_values = sorted(strata_values)
            for stratum_val in strata_values:
                sub_label = _get_stratum_name(strat_type, str(stratum_val), name_map)
                category = STRAT_TYPE_TO_CATEGORY.get(strat_type, strat_type)
                stratum_name = f"{category} ({sub_label})"
                add_stratum_section(
                    stratum_name, filtered[filtered["Stratum_Value"] == stratum_val]
                )

        return rows

    def write_table(rows: List[Dict], out_path: Path) -> None:
        output_df = pd.DataFrame(rows)
        col_order = ["Outcome"]
        for period in time_periods:
            col_order += [f"{period} Death,n(‰)", f"{period} AAMR"]
        output_df = output_df[[c for c in col_order if c in output_df.columns]]
        output_df.to_csv(out_path, index=False)
        print(f"  Saved {out_path.name}  ({len(output_df)} rows)")

    # Full version (always written)
    write_table(build_rows(display_codes), output_dir / f"AAMR_{shortname}.csv")

    # Simplified version: omit show=false diseases (written only when any are hidden)
    if hidden_codes:
        write_table(
            build_rows(display_codes_simplified),
            output_dir / f"AAMR_{shortname}_Simplified.csv",
        )


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

    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]

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


def process_sex_stratified_for_code(
    project_root: Path, icd_code: str, all_results: Dict[str, List[Dict]] = None
) -> Dict[str, Dict]:
    """
    Process sex stratified files for any ICD code from Stratified_Subtracted directory.
    Male deaths are calculated as Total - Female for consistency.

    Returns dict with 'male' and 'female' keys, each containing a row dict.
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return {}

    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]

    # Get national totals for this ICD code from all_results
    national_totals = {}
    if all_results and all_results.get("RUCC"):
        df = pd.DataFrame(all_results["RUCC"])
        ndd_national = df[
            (df["Stratum_Value"] == "National") & (df["ICD-10 Code"] == icd_code)
        ].copy()
        for period in time_periods:
            period_data = ndd_national[ndd_national["Time_Period"] == period]
            if not period_data.empty:
                national_totals[period] = {
                    "deaths": int(period_data["Total_Deaths"].iloc[0]),
                    "population": int(period_data["Total_Population"].iloc[0]),
                }

    # Get Female data
    female_data = {}
    for period in time_periods:
        file_path = input_dir / f"{period}_{icd_code}_Female.csv"
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

    result = {}

    # Calculate Male as Total - Female
    male_row: Dict = {"_type": "male"}
    for period in time_periods:
        if period in national_totals and period in female_data:
            total_deaths = national_totals[period]["deaths"]
            female_deaths = female_data[period]["deaths"]
            male_deaths = total_deaths - female_deaths

            file_path = input_dir / f"{period}_{icd_code}_Male.csv"
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

    result["male"] = male_row

    # Female row
    female_row: Dict = {"_type": "female"}
    for period, data in female_data.items():
        deaths = data["deaths"]
        population = data["population"]
        aamr = data["aamr"]
        se = data["se"]

        pct = (deaths / population) * 1000 if population > 0 else 0
        female_row[f"{period} Death,n(‰)"] = f"{deaths:,} ({pct:.2f}‰)"
        female_row[f"{period} AAMR"] = f"{aamr:.2f} ± {se:.2f}"

    result["female"] = female_row

    return result


def process_sex_stratified_all_codes(
    project_root: Path,
    all_results: Dict[str, List[Dict]],
    codes: List[str] = None,
) -> List[Dict]:
    """
    Process sex stratified files for all specified ICD codes from Stratified_Subtracted.
    Returns list of result dicts in all_results format with Stratum_Type='Sex'.
    Male deaths are calculated as Total - Female for consistency.
    Defaults to NDD_ALL_CODES when codes is None.
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return []

    codes_to_process = codes if codes is not None else NDD_ALL_CODES
    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]
    results = []

    # Pre-fetch national death totals for all codes (used for Male = Total - Female)
    national_totals: Dict[str, Dict[str, int]] = {}
    if all_results and all_results.get("RUCC"):
        nat_df = pd.DataFrame(all_results["RUCC"])
        nat_df = nat_df[nat_df["Stratum_Value"] == "National"]
        for code in codes_to_process:
            code_df = nat_df[nat_df["ICD-10 Code"] == code]
            national_totals[code] = {}
            for period in time_periods:
                period_df = code_df[code_df["Time_Period"] == period]
                if not period_df.empty:
                    national_totals[code][period] = int(
                        period_df["Total_Deaths"].iloc[0]
                    )

    for code in codes_to_process:
        female_deaths_by_period: Dict[str, int] = {}

        # Female
        for period in time_periods:
            female_path = input_dir / f"{period}_{code}_Female.csv"
            if not female_path.exists():
                continue
            df = pd.read_csv(female_path, dtype={"County Code": str})
            df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
            df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
            df["Population"] = pd.to_numeric(df["Population"], errors="coerce").fillna(
                0
            )
            df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})
            agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
            aamr, stats = calculate_aamr_point(agg_df)
            se = calculate_aamr_standard_error(agg_df)
            lower, upper = calculate_aamr_ci(aamr, stats)
            female_deaths_by_period[period] = stats["total_deaths"]
            results.append(
                {
                    "Time_Period": period,
                    "ICD-10 Code": code,
                    "Outcome": ICD_SHORTNAME_MAP.get(code, code),
                    "Stratum_Type": "Sex",
                    "Stratum_Value": "Female",
                    "Total_Deaths": stats["total_deaths"],
                    "Total_Population": stats["total_population"],
                    "AAMR": round(aamr, 4),
                    "AAMR_SE": round(se, 4),
                    "AAMR_Lower": round(lower, 4),
                    "AAMR_Upper": round(upper, 4),
                    "AAMRSE": f"{round(aamr, 2):.2f} ± {round(se, 2):.2f}",
                }
            )

        # Male (AAMR from Male file; deaths = Total - Female for consistency)
        for period in time_periods:
            male_path = input_dir / f"{period}_{code}_Male.csv"
            if not male_path.exists():
                continue
            df = pd.read_csv(male_path, dtype={"County Code": str})
            df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
            df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
            df["Population"] = pd.to_numeric(df["Population"], errors="coerce").fillna(
                0
            )
            df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})
            agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
            aamr, stats = calculate_aamr_point(agg_df)
            se = calculate_aamr_standard_error(agg_df)
            lower, upper = calculate_aamr_ci(aamr, stats)
            male_deaths = stats["total_deaths"]
            if (
                code in national_totals
                and period in national_totals[code]
                and period in female_deaths_by_period
            ):
                male_deaths = (
                    national_totals[code][period] - female_deaths_by_period[period]
                )
            results.append(
                {
                    "Time_Period": period,
                    "ICD-10 Code": code,
                    "Outcome": ICD_SHORTNAME_MAP.get(code, code),
                    "Stratum_Type": "Sex",
                    "Stratum_Value": "Male",
                    "Total_Deaths": male_deaths,
                    "Total_Population": stats["total_population"],
                    "AAMR": round(aamr, 4),
                    "AAMR_SE": round(se, 4),
                    "AAMR_Lower": round(lower, 4),
                    "AAMR_Upper": round(upper, 4),
                    "AAMRSE": f"{round(aamr, 2):.2f} ± {round(se, 2):.2f}",
                }
            )

    return results


def process_race_stratified_all_codes(
    project_root: Path,
    codes: List[str] = None,
) -> List[Dict]:
    """
    Process race stratified files (White, Black, Others) for all specified ICD codes.
    Returns list of result dicts in all_results format with Stratum_Type='Race'.
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return []

    codes_to_process = codes if codes is not None else NDD_ALL_CODES
    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]
    race_strata = ["White", "Black", "Others"]
    results = []

    for code in codes_to_process:
        for race in race_strata:
            for period in time_periods:
                file_path = input_dir / f"{period}_{code}_{race}.csv"
                if not file_path.exists():
                    continue
                df = pd.read_csv(file_path, dtype={"County Code": str})
                df = df[df["County Code"].notna() & (df["Population"] != "Missing")]
                df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce").fillna(0)
                df["Population"] = pd.to_numeric(df["Population"], errors="coerce").fillna(0)
                df = df.rename(columns={"Ten-Year Age Groups": "Age_Group"})
                agg_df = df.groupby("Age_Group")[["Deaths", "Population"]].sum()
                aamr, stats = calculate_aamr_point(agg_df)
                se = calculate_aamr_standard_error(agg_df)
                lower, upper = calculate_aamr_ci(aamr, stats)
                results.append(
                    {
                        "Time_Period": period,
                        "ICD-10 Code": code,
                        "Outcome": ICD_SHORTNAME_MAP.get(code, code),
                        "Stratum_Type": "Race",
                        "Stratum_Value": race,
                        "Total_Deaths": stats["total_deaths"],
                        "Total_Population": stats["total_population"],
                        "AAMR": round(aamr, 4),
                        "AAMR_SE": round(se, 4),
                        "AAMR_Lower": round(lower, 4),
                        "AAMR_Upper": round(upper, 4),
                        "AAMRSE": f"{round(aamr, 2):.2f} ± {round(se, 2):.2f}",
                    }
                )

    return results


def process_stratified_ndd_files_race_only(project_root: Path) -> List[Dict]:
    """
    Process race stratified NDD files from Stratified_Subtracted directory.

    Returns list of result rows for NDD by race strata with format:
    NDD (White), NDD (Black), NDD (Asian), NDD (Indian)
    """
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    if not input_dir.exists():
        return []

    time_periods = ["2006-2010", "2011-2015", "2016-2020", "2021-2024"]

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


def main():
    project_root, config = load_config()

    global NDD_ABBR_MAP, NDD_ALL_CODES, ICD_SHORTNAME_MAP
    NDD_ABBR_MAP, NDD_ALL_CODES = _build_ndd_maps(config)
    ICD_SHORTNAME_MAP = _build_shortname_map(config)

    disease_groups = _build_disease_groups(config)

    input_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"
    output_dir = project_root / "Result/Table_AAMR"
    output_dir.mkdir(parents=True, exist_ok=True)

    strat_data = load_stratification_data(project_root, config)

    subtracted_files = list(input_dir.glob("*.csv"))
    all_results = {
        "RUCC": [],
        "Climate": [],
        "Cluster": [],
        "Typology": [],
        "EQI": [],
        "Sex": [],
        "Race": [],
    }

    for file_path in sorted(subtracted_files):
        file_results = process_file(file_path, strat_data, config)
        for key in all_results:
            if key in file_results:
                all_results[key].extend(file_results[key])

    # Collect all overall disease codes that have stratified files
    all_strat_codes = [
        info["overall_code"]
        for info in disease_groups.values()
        if info.get("overall_code")
    ]

    # Populate sex and race stratification after all county-level results are collected
    all_results["Sex"] = process_sex_stratified_all_codes(
        project_root, all_results, codes=all_strat_codes
    )
    all_results["Race"] = process_race_stratified_all_codes(
        project_root, codes=all_strat_codes
    )

    # Write raw stratified AAMR data
    # Cluster, Sex, and Race are not written as separate files.
    # Climate is split: census_region → Region, koppen_major → Climate.
    SKIP_STRAT = {"Cluster", "Sex", "Race"}
    for strat_key, results_list in all_results.items():
        if strat_key in SKIP_STRAT or not results_list:
            continue
        df = pd.DataFrame(results_list)
        df = df.sort_values(
            ["Time_Period", "ICD-10 Code", "Stratum_Type", "Stratum_Value"]
        )
        df = apply_stratification_labels(
            df, config
        )  # adds Stratification, Stratum_Label
        df = df.rename(columns={"Stratum_Label": "Strata"})
        df = df.drop(columns=["Stratum_Type", "Stratum_Value"])
        # Canonical column order
        col_order = [
            "Time_Period",
            "ICD-10 Code",
            "Outcome_Fullname",
            "Outcome_Shortname",
            "Stratification",
            "Strata",
            "Total_Deaths",
            "Total_Population",
            "AAMR",
            "AAMR_SE",
            "AAMR_Lower",
            "AAMR_Upper",
            "AAMRSE",
        ]
        df = df[[c for c in col_order if c in df.columns]]

        if strat_key == "Climate":
            region_df = df[df["Stratification"] == "Census Region"].copy()
            climate_df = df[df["Stratification"] == "Köppen-Geiger Climate Zone"].copy()
            if not region_df.empty:
                region_df.to_csv(output_dir / "Stratified_AAMR_Region.csv", index=False)
                print("Saved Region stratified AAMR")
            if not climate_df.empty:
                climate_df.to_csv(
                    output_dir / "Stratified_AAMR_Climate.csv", index=False
                )
                print("Saved Climate stratified AAMR")
        else:
            output_path = output_dir / f"Stratified_AAMR_{strat_key}.csv"
            df.to_csv(output_path, index=False)
            print(f"Saved {strat_key} stratified AAMR to {output_path}")

    print("\nStratified AAMR calculation completed.")

    # Create one summary table per disease group (AAMR_{shortname}.csv)
    for group_key, group_info in disease_groups.items():
        create_disease_aamr_table(
            group_key, group_info, all_results, project_root, output_dir
        )

    print("All processing completed successfully!")


if __name__ == "__main__":
    main()
