"""
CDC Triangulation (Stratified) — Data Integrity Checker, Merger, and Subtraction

This script mirrors EQI_Triangulation.py but supports stratified files, where
the strata (e.g., White, Male, Black) are encoded in the filenames.

For 2021-2024, population data is provided in separate 5-year age group files
(named '{year} ALL_{strat} {range}.csv'). The Sex file contains both Male and
Female populations; Black and White files each contain a single race.
No triangulation subtraction is performed for 2021-2024.

Phase 1: Data Integrity Validation (years, ages, ICD codes, and stratum presence)
Phase 2: Merge age-partitioned files (Age1, Age35, Age65) for each (year, ICD set, stratum)
Phase 3 (2006-2020): Subtract residual ALL from ALL_{ICD} to obtain target-disease deaths
Phase 3b (2021-2024): Attach 5-year population (aggregated to 10-year) to subtracted files

Important filename patterns:
- Death files:      2016-2020 ALL_G20_G30_G12.2_F01_F03_White Age1.csv
- Death files:      2021-2024 ALL_I00_I99_Black Age1.csv
- Population files: 2021-2024 ALL_Black 0-9.csv
- Population files: 2021-2024 ALL_Sex 10-19.csv  (contains Male and Female)
- Population files: 2021-2024 ALL_White 80-89.csv

Directories:
- Input:       Data/Original/CDC Triangulation/Stratified_NotMerged
- Merged:      Data/Original/CDC Triangulation/Stratified_Merged
- Subtracted:  Data/Original/CDC Triangulation/Stratified_Subtracted
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from io import StringIO

import pandas as pd
import yaml

AGE_GROUP_MAPPINGS = {
    "Age1": ["< 1 year", "1-4 years", "5-14 years", "15-24 years", "25-34 years"],
    "Age35": ["35-44 years", "45-54 years", "55-64 years"],
    "Age65": ["65-74 years", "75-84 years", "85+ years"],
}

AGE_GROUP_ORDER = [
    "< 1 year",
    "1-4 years",
    "5-14 years",
    "15-24 years",
    "25-34 years",
    "35-44 years",
    "45-54 years",
    "55-64 years",
    "65-74 years",
    "75-84 years",
    "85+ years",
]

YEAR_RANGE_MAPPINGS = {
    "2006-2010": ["2006", "2007", "2008", "2009", "2010"],
    "2011-2015": ["2011", "2012", "2013", "2014", "2015"],
    "2016-2020": ["2016", "2017", "2018", "2019", "2020"],
    "2021-2024": ["2021", "2022", "2023", "2024"],
}

DIRECT_YEARS = {"2021-2024"}

# Maps 10-year death age groups to their constituent 5-year population groups.
# < 1 year and 1-4 years both use the 0-4 years population (cannot be split further).
FIVE_YEAR_TO_TEN_YEAR = {
    "< 1 year":    ["0-4 years"],
    "1-4 years":   ["0-4 years"],
    "5-14 years":  ["5-9 years", "10-14 years"],
    "15-24 years": ["15-19 years", "20-24 years"],
    "25-34 years": ["25-29 years", "30-34 years"],
    "35-44 years": ["35-39 years", "40-44 years"],
    "45-54 years": ["45-49 years", "50-54 years"],
    "55-64 years": ["55-59 years", "60-64 years"],
    "65-74 years": ["65-69 years", "70-74 years"],
    "75-84 years": ["75-79 years", "80-84 years"],
    "85+ years":   ["85+ years"],
}

# Maps stratum name to its population source file prefix and optional gender filter.
# Male and Female populations are both stored in the "Sex" file.
STRATUM_POP_SOURCE = {
    "Black":  {"file_strat": "Black",  "filter_col": None,     "filter_val": None},
    "White":  {"file_strat": "White",  "filter_col": None,     "filter_val": None},
    "Male":   {"file_strat": "Sex",    "filter_col": "Gender", "filter_val": "Male"},
    "Female": {"file_strat": "Sex",    "filter_col": "Gender", "filter_val": "Female"},
}


def load_config() -> tuple[Path, dict]:
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return project_root, config


def is_population_file_stratified(filename: str) -> bool:
    """
    Detect stratified population-only files: '{year} ALL_{strat} {digit-range}.csv'
    These differ from death files which end with ' AgeN.csv'.
    """
    return bool(re.match(r"^\d{4}-\d{4} ALL_\w+ \d[\d-]*\.csv$", filename))


def parse_filename(filename: str) -> Optional[Dict[str, Optional[str]]]:
    """
    Parse stratified filename to extract year range, ICD code block, stratum, and age partition.

    Expected patterns (examples):
      - '2016-2020 ALL_G20_G30_G12.2_F01_F03_White Age1.csv'
      - '2016-2020 ALL_Indian Age1.csv'
      - '2021-2024 ALL_I00_I99_Black Age35.csv'

    Returns:
      {'year': '2016-2020', 'icd': 'G20_G30_G12.2_F01_F03' or None,
       'stratum': 'White', 'age': 'Age1' | 'Age35' | 'Age65'}
    """
    if not filename.endswith(".csv"):
        return None

    name = filename[:-4]

    if " Age" not in name:
        return None

    base, age = name.rsplit(" ", 1)
    if not re.fullmatch(r"Age\d+", age):
        return None

    if " " not in base:
        return None

    year, rest = base.split(" ", 1)
    if not re.fullmatch(r"\d{4}-\d{4}", year):
        return None

    if not rest.startswith("ALL"):
        # Handle 'All_Cause_of_Death_{stratum}' naming convention from CDC WONDER
        m = re.fullmatch(r"All_Cause_of_Death_(\w+)", rest)
        if m:
            return {"year": year, "icd": None, "stratum": m.group(1), "age": age}
        return None

    suffix = rest[3:]
    if suffix.startswith("_"):
        suffix = suffix[1:]

    if suffix == "":
        icd_block = None
        stratum = None
    else:
        tokens = suffix.split("_")
        if len(tokens) == 1:
            icd_block = None
            stratum = tokens[0]
        else:
            icd_block = "_".join(tokens[:-1]) if len(tokens) > 1 else None
            stratum = tokens[-1]

    return {"year": year, "icd": icd_block, "stratum": stratum, "age": age}


def read_file_metadata(file_path: Path) -> Dict[str, str]:
    with open(file_path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    query_param_idx = -1
    for i in range(len(lines) - 1, max(0, len(lines) - 500), -1):
        if '"Query Parameters:"' in lines[i]:
            query_param_idx = i
            break

    if query_param_idx == -1:
        return {}

    metadata_lines = lines[query_param_idx:]
    metadata_text = " ".join(line.strip().strip('"') for line in metadata_lines)

    metadata: Dict[str, str] = {}

    title_match = re.search(
        r'Title:\s*([^\s][^"]*?)(?=\s+ICD-10 Codes:|\s+Ten-Year Age Groups:|\s+Year/Month:)',
        metadata_text,
    )
    if title_match:
        metadata["title"] = title_match.group(1).strip()

    icd_match = re.search(
        r"ICD-10 Codes:\s*(.*?)(?=\s+Ten-Year Age Groups:|\s+Year/Month:)",
        metadata_text,
        re.DOTALL,
    )
    if icd_match:
        metadata["icd_codes"] = icd_match.group(1).strip()

    age_match = re.search(
        r'Ten-Year Age Groups:\s*([^"]+?)(?=\s+Year/Month:|\s+Group By:)', metadata_text
    )
    if age_match:
        metadata["age_groups"] = age_match.group(1).strip()

    year_match = re.search(
        r'Year/Month:\s*([^"]+?)(?=\s+Group By:|\s+Show)', metadata_text
    )
    if year_match:
        metadata["years"] = year_match.group(1).strip()

    return metadata


def extract_data_rows(
    file_path: Path, keep_metadata: bool = True
) -> Tuple[pd.DataFrame, List[str]]:
    with open(file_path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    metadata_start = -1
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].strip() == '"---"':
            metadata_start = i
            break

    if metadata_start == -1:
        metadata_start = len(lines)

    metadata_rows = (
        lines[metadata_start:] if metadata_start < len(lines) and keep_metadata else []
    )

    csv_data = "".join(lines[:metadata_start])
    df = pd.read_csv(StringIO(csv_data), dtype={"County Code": str})

    columns_to_drop = []
    if "Notes" in df.columns:
        columns_to_drop.append("Notes")
    if "County" in df.columns:
        columns_to_drop.append("County")

    if columns_to_drop:
        df = df.drop(columns=columns_to_drop)

    # Drop footnote/note rows that leaked into the CSV region (County Code is NaN)
    df = df[df["County Code"].notna()].copy()

    return df, metadata_rows


def validate_year_range(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    expected_years = YEAR_RANGE_MAPPINGS.get(filename_info["year"])
    if not expected_years:
        return False, f"Unknown year range: {filename_info['year']}"

    metadata_years = metadata.get("years", "")
    for year in expected_years:
        if year not in metadata_years:
            return False, f"Missing year {year} in metadata. Found: {metadata_years}"

    return True, "OK"


def validate_age_group(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    expected_ages = AGE_GROUP_MAPPINGS.get(filename_info["age"])
    if not expected_ages:
        return False, f"Unknown age group: {filename_info['age']}"

    metadata_ages = metadata.get("age_groups", "")
    for age in expected_ages:
        if age not in metadata_ages:
            return (
                False,
                f"Missing age group '{age}' in metadata. Found: {metadata_ages}",
            )

    return True, "OK"


def validate_icd_codes(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    icd_in_filename = filename_info.get("icd")
    year = filename_info.get("year", "")

    # For 2021-2024, death files are residual (ICD codes excluded from query).
    # The target ICD codes therefore do not appear in the metadata ICD list.
    # Validate using the Title field instead, which always reflects the filename.
    if year in DIRECT_YEARS:
        if not icd_in_filename:
            return True, "OK - 2021-2024 residual ALL file"
        title = metadata.get("title", "")
        if f"ALL_{icd_in_filename}" in title:
            return True, "OK - 2021-2024 residual file confirmed by title"
        return True, f"OK - 2021-2024 file accepted (title: {title!r})"

    if not icd_in_filename:
        icd_text = metadata.get("icd_codes", "")
        title = metadata.get("title", "")
        if "A00-B99" in icd_text or "D50-D89" in icd_text or "All causes" in icd_text:
            return True, "OK - General disease categories"
        if "All_Cause_of_Death" in title:
            return True, "OK - All_Cause file confirmed by title"
        return False, "Expected general disease categories for ALL file"

    icd_codes = icd_in_filename.replace("-", "_").split("_")
    icd_text = metadata.get("icd_codes", "")
    missing_codes = [code for code in icd_codes if code not in icd_text]

    if not missing_codes:
        return True, "OK"

    title = metadata.get("title", "")
    title_alnum = re.sub(r"[^A-Za-z0-9]", "", title)
    icd_alnum = re.sub(r"[^A-Za-z0-9]", "", icd_in_filename)
    if f"ALL{icd_alnum}" in title_alnum:
        return True, "OK - Residual file confirmed by title"

    return False, f"Missing ICD codes in metadata: {', '.join(missing_codes)}"


def validate_stratum(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    stratum = filename_info.get("stratum")
    if not stratum:
        return True, "No stratum in filename"

    haystacks = []
    for key in ("title", "icd_codes", "age_groups", "years"):
        val = metadata.get(key, "")
        if isinstance(val, str) and val:
            haystacks.append(val.lower())

    needle = stratum.lower()

    aliases = {
        "indian": ["american indian", "american indian or alaska native"],
        "white": ["white"],
        "black": ["black", "african american"],
        "male": ["male"],
        "female": ["female"],
        "hispanic": ["hispanic"],
        "asian": ["asian"],
    }
    candidates = [needle] + aliases.get(needle, [])

    for cand in candidates:
        if any(cand in h for h in haystacks):
            return True, "OK"

    return True, f"Stratum '{stratum}' not confirmed in metadata title; proceeding"


def validate_file(file_path: Path) -> Dict:
    filename = file_path.name
    result = {
        "filename": filename,
        "path": str(file_path),
        "parsed": False,
        "year_valid": False,
        "age_valid": False,
        "icd_valid": False,
        "stratum_valid": False,
        "overall_valid": False,
        "errors": [],
    }

    filename_info = parse_filename(filename)
    if not filename_info:
        result["errors"].append(f"Could not parse filename: {filename}")
        return result

    result["parsed"] = True
    result["filename_info"] = filename_info

    try:
        metadata = read_file_metadata(file_path)
        result["metadata"] = metadata
    except Exception as e:
        result["errors"].append(f"Error reading metadata: {e}")
        return result

    year_valid, year_msg = validate_year_range(filename_info, metadata)
    result["year_valid"] = year_valid
    if not year_valid:
        result["errors"].append(f"Year validation: {year_msg}")

    age_valid, age_msg = validate_age_group(filename_info, metadata)
    result["age_valid"] = age_valid
    if not age_valid:
        result["errors"].append(f"Age validation: {age_msg}")

    icd_valid, icd_msg = validate_icd_codes(filename_info, metadata)
    result["icd_valid"] = icd_valid
    if not icd_valid:
        result["errors"].append(f"ICD validation: {icd_msg}")

    stratum_valid, stratum_msg = validate_stratum(filename_info, metadata)
    result["stratum_valid"] = stratum_valid
    if not stratum_valid:
        result["errors"].append(f"Stratum validation: {stratum_msg}")

    result["overall_valid"] = year_valid and age_valid and icd_valid and stratum_valid

    return result


def get_base_filename(filename: str) -> Optional[str]:
    info = parse_filename(filename)
    if not info:
        return None

    year = info["year"]
    icd = info.get("icd")
    stratum = info.get("stratum")

    if stratum and icd:
        return f"{year} ALL_{icd}_{stratum}"
    elif stratum and not icd:
        return f"{year} ALL_{stratum}"
    elif icd and not stratum:
        return f"{year} ALL_{icd}"
    else:
        return f"{year} ALL"


def merge_file_group(file_paths: List[Path], output_path: Path):
    all_dataframes = []
    metadata_rows = None

    print(f"  Merging {len(file_paths)} files...")

    for file_path in sorted(file_paths):
        print(f"    - {file_path.name}")
        df, meta_rows = extract_data_rows(file_path)

        all_dataframes.append(df)

        if metadata_rows is None:
            metadata_rows = meta_rows

    merged_df = pd.concat(all_dataframes, ignore_index=True)

    age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}

    if "Ten-Year Age Groups" in merged_df.columns:
        merged_df["_age_order"] = merged_df["Ten-Year Age Groups"].map(age_order_map)
        merged_df = merged_df.sort_values(
            by=["County Code", "_age_order"], na_position="last"
        )
        merged_df = merged_df.drop(columns=["_age_order"])
    else:
        merged_df = merged_df.sort_values(by=["County Code"])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        merged_df.to_csv(f, index=False)
        if metadata_rows:
            f.writelines(metadata_rows)

    print(f"  ✓ Saved to: {output_path.name} ({len(merged_df)} rows)")


def _safe_int(val, default: int = 0) -> int:
    """Parse a CDC WONDER value that may be NaN, 'Missing', 'Suppressed', or 'Unreliable'."""
    if pd.isna(val):
        return default
    s = str(val).strip()
    if s.lower() in ("missing", "suppressed", "unreliable", "nan", ""):
        return default
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return default


def parse_death_value(value):
    if pd.isna(value):
        return None, False, True

    if isinstance(value, str):
        if "Suppressed" in value:
            return None, True, False
        try:
            return int(float(value)), False, False
        except ValueError:
            return None, False, True

    return int(float(value)), False, False


def calculate_target_deaths(ti_value, ri_value):
    ti, ti_suppressed, ti_missing = parse_death_value(ti_value)
    ri, ri_suppressed, ri_missing = parse_death_value(ri_value)

    if ti_missing:
        return None, "Missing_Ti"
    if ri_missing:
        return None, "Missing_Ri"

    if ti_suppressed and ri_suppressed:
        return 0, "Both_Suppressed"
    if ti_suppressed:
        return 0, "Suppressed_Ti"
    if ri_suppressed:
        return 0, "Suppressed_Ri"

    di = int(max(0, ti - ri))
    return di, "Clean"


def load_population_2021_2024_stratified(
    input_dir: Path, year_range: str, stratum: str
) -> Dict[str, Dict[str, int]]:
    """
    Load stratified 5-year population files for a direct-year range and aggregate to 10-year groups.
    Returns: {county_code: {ten_year_age_group: population}}

    For Male/Female strata, reads the combined Sex file and filters by Gender column.
    Both < 1 year and 1-4 years receive the combined 0-4 years population.
    """
    source = STRATUM_POP_SOURCE.get(stratum)
    if not source:
        print(f"  ⚠ No population source defined for stratum '{stratum}'")
        return {}

    file_strat = source["file_strat"]
    filter_col = source["filter_col"]
    filter_val = source["filter_val"]

    pop_files = sorted(
        f for f in input_dir.glob(f"{year_range} ALL_{file_strat} *.csv")
        if is_population_file_stratified(f.name)
    )

    if not pop_files:
        print(f"  ⚠ No population files found for {year_range} / {stratum}")
        return {}

    county_pop_5yr: Dict[str, Dict[str, int]] = {}

    for pop_file in pop_files:
        with open(pop_file, "r", encoding="utf-8") as f:
            lines = f.readlines()

        metadata_start = len(lines)
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].strip() == '"---"':
                metadata_start = i
                break

        csv_data = "".join(lines[:metadata_start])
        df = pd.read_csv(StringIO(csv_data), dtype={"County Code": str})

        if "Notes" in df.columns:
            df = df[df["Notes"].isna() | (df["Notes"].astype(str) != "Total")]

        if filter_col and filter_val and filter_col in df.columns:
            df = df[df[filter_col].astype(str) == filter_val]

        df = df[df["Age Group"].notna()]

        for _, row in df.iterrows():
            county = str(row["County Code"]).zfill(5)
            age_5yr = str(row["Age Group"])
            pop = row["Population"]

            if pd.isna(pop):
                continue

            if county not in county_pop_5yr:
                county_pop_5yr[county] = {}
            county_pop_5yr[county][age_5yr] = int(float(pop))

    result: Dict[str, Dict[str, int]] = {}
    for county, pop_5yr in county_pop_5yr.items():
        result[county] = {}
        for ten_yr_group, five_yr_groups in FIVE_YEAR_TO_TEN_YEAR.items():
            total = sum(pop_5yr.get(fg, 0) for fg in five_yr_groups)
            result[county][ten_yr_group] = total

    print(f"  ✓ Population loaded ({stratum}): {len(result)} counties from {len(pop_files)} files")
    return result


def subtract_deaths_stratified(merged_dir: Path, output_dir: Path):
    """
    Phase 3 (Stratified): Calculate target disease deaths by subtraction for all year ranges.
    Di = Ti (ALL_{ICD}) - Ri (ALL, no ICD filter)
    Requires {year} ALL_{stratum}.csv in merged_dir for every stratum.
    For 2021-2024, the merged files may lack a Population column; Phase 3b attaches it.
    """
    print("\n" + "=" * 70)
    print("PHASE 3 (Stratified): Calculate Target Disease Deaths (Subtraction)")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    all_files = list(merged_dir.glob("*.csv"))

    group_keys = set()
    all_files_by_group: Dict[Tuple[str, str], Path] = {}
    icd_files_by_group: Dict[Tuple[str, str], Dict[str, Path]] = {}

    for file_path in all_files:
        filename = file_path.name
        info = parse_filename(filename.replace(".csv", " Age1.csv"))
        if not info:
            continue

        year_range = info["year"]

        stratum = info.get("stratum") or ""
        icd_block = info.get("icd")

        key = (year_range, stratum)
        group_keys.add(key)

        if icd_block is None:
            all_files_by_group[key] = file_path
        else:
            icd_files_by_group.setdefault(key, {})
            icd_files_by_group[key][icd_block] = file_path

    print(f"\nFound {len(group_keys)} (year, stratum) groups for subtraction")

    total_processed = 0

    for year_range, stratum in sorted(group_keys):
        print(f"\n{year_range} — {stratum if stratum else '(no stratum)'}:")

        if (year_range, stratum) not in all_files_by_group:
            print(f"  ⚠ Missing ALL file for group, skipping")
            continue

        if (year_range, stratum) not in icd_files_by_group:
            print(f"  ⚠ No ICD files for group, skipping")
            continue

        all_file = all_files_by_group[(year_range, stratum)]
        print(f"  Loading residual (ALL): {all_file.name}")
        df_all, _ = extract_data_rows(all_file, keep_metadata=False)

        for icd_code, icd_file in sorted(
            icd_files_by_group[(year_range, stratum)].items()
        ):
            print(f"    Processing {icd_code}...")

            df_icd, _ = extract_data_rows(icd_file, keep_metadata=False)

            merged = df_all.merge(
                df_icd,
                on=["County Code", "Ten-Year Age Groups"],
                suffixes=("_Total", "_Residual"),
                how="inner",
            )

            results = []
            for _, row in merged.iterrows():
                ti = row["Deaths_Total"]
                ri = row["Deaths_Residual"]

                di, quality_flag = calculate_target_deaths(ti, ri)

                pop_val = row.get("Population_Residual", None)

                results.append(
                    {
                        "County Code": row["County Code"],
                        "Ten-Year Age Groups": row["Ten-Year Age Groups"],
                        "Ten-Year Age Groups Code": row.get(
                            "Ten-Year Age Groups Code_Residual",
                            row.get("Ten-Year Age Groups Code", ""),
                        ),
                        "Deaths": di,
                        "Population": _safe_int(pop_val) if pop_val is not None and not (isinstance(pop_val, float) and pd.isna(pop_val)) else pd.NA,
                        "Quality_Flag": quality_flag,
                    }
                )

            result_df = pd.DataFrame(results)

            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(age_order_map)
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])

            result_df["Deaths"] = result_df["Deaths"].astype("Int64")
            result_df["Population"] = result_df["Population"].astype("Int64")

            output_filename = (
                f"{year_range}_{icd_code}_{stratum}.csv"
                if stratum
                else f"{year_range}_{icd_code}.csv"
            )
            output_path = output_dir / output_filename
            result_df.to_csv(output_path, index=False)

            clean_count = (result_df["Quality_Flag"] == "Clean").sum()
            suppressed_count = result_df["Quality_Flag"].str.contains("Suppressed").sum()
            total_count = len(result_df)

            print(f"      ✓ Saved to: {output_filename}")
            print(f"        {total_count} rows: {clean_count} clean, {suppressed_count} suppressed")

            total_processed += 1

    print(f"\n✓ Stratified subtraction completed: {total_processed} ICD code blocks processed")


def compute_female(
    stratified_subtracted_dir: Path,
    nonstrat_subtracted_dir: Path,
    output_dir: Path,
):
    """
    Compute 'Female' sex stratum for all year ranges: Female = Total - Male.
    Reads Male from stratified_subtracted_dir and Total from nonstrat_subtracted_dir.
    """
    print("\n" + "=" * 70)
    print("PHASE 3d (Stratified): Compute 'Female' Sex Stratum")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    year_icd_male: Dict[str, Dict[str, Path]] = defaultdict(dict)

    for f in sorted(stratified_subtracted_dir.glob("*_Male.csv")):
        m = re.match(r"^(\d{4}-\d{4})_(.+)_Male$", f.stem)
        if m:
            year_icd_male[m.group(1)][m.group(2)] = f

    total_processed = 0

    for year_range, male_icds in sorted(year_icd_male.items()):
        print(f"\n{year_range}: {len(male_icds)} ICD codes")

        for icd, male_file in sorted(male_icds.items()):
            total_file = nonstrat_subtracted_dir / f"{year_range}_{icd}.csv"
            if not total_file.exists():
                print(f"  ⚠ No total file for {icd}, skipping")
                continue

            total_df = pd.read_csv(total_file, dtype={"County Code": str})
            total_df = total_df[total_df["County Code"].notna()].copy()
            male_df = pd.read_csv(male_file, dtype={"County Code": str})
            male_df = male_df[male_df["County Code"].notna()].copy()

            key_cols = ["County Code", "Ten-Year Age Groups"]
            merged = (
                total_df[key_cols + ["Ten-Year Age Groups Code", "Deaths", "Population"]]
                .rename(columns={"Deaths": "Deaths_Total", "Population": "Population_Total"})
                .merge(
                    male_df[key_cols + ["Deaths", "Population"]].rename(
                        columns={"Deaths": "Deaths_Male", "Population": "Population_Male"}
                    ),
                    on=key_cols,
                    how="left",
                )
            )

            results = []
            for _, row in merged.iterrows():
                t_d = _safe_int(row["Deaths_Total"])
                m_d = _safe_int(row["Deaths_Male"])

                t_p = _safe_int(row["Population_Total"])
                m_p = _safe_int(row["Population_Male"])

                results.append({
                    "County Code": str(row["County Code"]).zfill(5),
                    "Ten-Year Age Groups": row["Ten-Year Age Groups"],
                    "Ten-Year Age Groups Code": row.get("Ten-Year Age Groups Code", ""),
                    "Deaths": max(0, t_d - m_d),
                    "Population": max(0, t_p - m_p),
                    "Quality_Flag": "Derived",
                })

            result_df = pd.DataFrame(results)
            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(age_order_map)
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])
            result_df["Deaths"] = result_df["Deaths"].astype("Int64")
            result_df["Population"] = result_df["Population"].astype("Int64")

            output_path = output_dir / f"{year_range}_{icd}_Female.csv"
            result_df.to_csv(output_path, index=False)

            print(f"  ✓ {year_range}_{icd}_Female.csv ({len(result_df)} rows)")
            total_processed += 1

    print(f"\n✓ Female computation: {total_processed} files created")


def compute_others_race(
    stratified_subtracted_dir: Path,
    nonstrat_subtracted_dir: Path,
    output_dir: Path,
):
    """
    Compute 'Others' race category for all year ranges: Others = Total - Black - White.
    Reads from already-computed stratified (Black, White) and non-stratified (Total) subtracted files.
    Requires both scripts to have run first.
    """
    print("\n" + "=" * 70)
    print("PHASE 3c (Stratified): Compute 'Others' Race Category")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    year_icd_black: Dict[str, Dict[str, Path]] = defaultdict(dict)
    year_icd_white: Dict[str, Dict[str, Path]] = defaultdict(dict)

    for f in sorted(stratified_subtracted_dir.glob("*_Black.csv")):
        m = re.match(r"^(\d{4}-\d{4})_(.+)_Black$", f.stem)
        if m:
            year_icd_black[m.group(1)][m.group(2)] = f

    for f in sorted(stratified_subtracted_dir.glob("*_White.csv")):
        m = re.match(r"^(\d{4}-\d{4})_(.+)_White$", f.stem)
        if m:
            year_icd_white[m.group(1)][m.group(2)] = f

    all_years = sorted(set(year_icd_black) | set(year_icd_white))
    total_processed = 0

    for year_range in all_years:
        black_icds = year_icd_black.get(year_range, {})
        white_icds = year_icd_white.get(year_range, {})
        common_icds = sorted(set(black_icds) & set(white_icds))

        if not common_icds:
            continue

        print(f"\n{year_range}: {len(common_icds)} ICD codes")

        for icd in common_icds:
            total_file = nonstrat_subtracted_dir / f"{year_range}_{icd}.csv"
            if not total_file.exists():
                print(f"  ⚠ No total file for {icd}, skipping")
                continue

            total_df = pd.read_csv(total_file, dtype={"County Code": str})
            total_df = total_df[total_df["County Code"].notna()].copy()
            black_df = pd.read_csv(black_icds[icd], dtype={"County Code": str})
            black_df = black_df[black_df["County Code"].notna()].copy()
            white_df = pd.read_csv(white_icds[icd], dtype={"County Code": str})
            white_df = white_df[white_df["County Code"].notna()].copy()

            key_cols = ["County Code", "Ten-Year Age Groups"]
            merged = (
                total_df[key_cols + ["Ten-Year Age Groups Code", "Deaths", "Population"]]
                .rename(columns={"Deaths": "Deaths_Total", "Population": "Population_Total"})
                .merge(
                    black_df[key_cols + ["Deaths", "Population"]].rename(
                        columns={"Deaths": "Deaths_Black", "Population": "Population_Black"}
                    ),
                    on=key_cols,
                    how="left",
                )
                .merge(
                    white_df[key_cols + ["Deaths", "Population"]].rename(
                        columns={"Deaths": "Deaths_White", "Population": "Population_White"}
                    ),
                    on=key_cols,
                    how="left",
                )
            )

            results = []
            for _, row in merged.iterrows():
                t_d = _safe_int(row["Deaths_Total"])
                b_d = _safe_int(row["Deaths_Black"])
                w_d = _safe_int(row["Deaths_White"])

                t_p = _safe_int(row["Population_Total"])
                b_p = _safe_int(row["Population_Black"])
                w_p = _safe_int(row["Population_White"])

                results.append({
                    "County Code": str(row["County Code"]).zfill(5),
                    "Ten-Year Age Groups": row["Ten-Year Age Groups"],
                    "Ten-Year Age Groups Code": row.get("Ten-Year Age Groups Code", ""),
                    "Deaths": max(0, t_d - b_d - w_d),
                    "Population": max(0, t_p - b_p - w_p),
                    "Quality_Flag": "Derived",
                })

            result_df = pd.DataFrame(results)
            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(age_order_map)
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])
            result_df["Deaths"] = result_df["Deaths"].astype("Int64")
            result_df["Population"] = result_df["Population"].astype("Int64")

            output_path = output_dir / f"{year_range}_{icd}_Others.csv"
            result_df.to_csv(output_path, index=False)

            print(f"  ✓ {year_range}_{icd}_Others.csv ({len(result_df)} rows)")
            total_processed += 1

    print(f"\n✓ Others computation: {total_processed} files created")


def attach_population_stratified(subtracted_dir: Path, input_dir: Path):
    """
    Phase 3b (Stratified): Fill in missing Population for 2021-2024 subtracted files.
    For 2021-2024, CDC WONDER death files do not include population; it must be attached
    from the separate 5-year age-group population files in input_dir.
    Only rows with null Population are updated (earlier years already have it from WONDER).
    """
    print("\n" + "=" * 70)
    print("PHASE 3b (Stratified): Attach Population for 2021-2024 Subtracted Files")
    print("=" * 70)

    for year_range in sorted(DIRECT_YEARS):
        subtracted_files = sorted(subtracted_dir.glob(f"{year_range}_*.csv"))
        if not subtracted_files:
            print(f"\n  No subtracted files found for {year_range}, skipping")
            continue

        # Group by stratum to load each population source only once
        stratum_to_files: Dict[str, List[Path]] = defaultdict(list)
        for sf in subtracted_files:
            stem = sf.stem  # e.g. "2021-2024_K74_Male"
            parts = stem.split("_")
            if len(parts) >= 3:
                stratum = parts[-1]
                stratum_to_files[stratum].append(sf)

        for stratum, files in sorted(stratum_to_files.items()):
            if stratum not in STRATUM_POP_SOURCE:
                continue  # Female/Others are derived later; skip here

            print(f"\nLoading population for {year_range} / {stratum}...")
            pop_lookup = load_population_2021_2024_stratified(
                input_dir, year_range, stratum
            )
            if not pop_lookup:
                continue

            for sf in sorted(files):
                df = pd.read_csv(sf, dtype={"County Code": str})
                if "Population" not in df.columns:
                    df["Population"] = None

                # Fill only rows where Population is null
                mask = df["Population"].isna()
                if not mask.any():
                    continue

                def _fill_pop(row):
                    county = str(row["County Code"]).zfill(5)
                    age_grp = row["Ten-Year Age Groups"]
                    return pop_lookup.get(county, {}).get(age_grp, None)

                df.loc[mask, "Population"] = df[mask].apply(_fill_pop, axis=1)
                df["Deaths"] = df["Deaths"].astype("Int64")
                df.to_csv(sf, index=False)

                filled = mask.sum()
                print(f"  ✓ {sf.name}: filled population for {filled} rows")

    print("\n✓ Population attachment completed")


def main():
    print("=" * 70)
    print("CDC Triangulation (Stratified) — Validation, Merge, and Subtraction")
    print("=" * 70)

    project_root, _ = load_config()

    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_NotMerged"
    merged_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Merged"
    subtracted_dir = (
        project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    )
    nonstrat_subtracted_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"

    print(f"\nInput directory:       {input_dir}")
    print(f"Merged directory:      {merged_dir}")
    print(f"Subtracted directory:  {subtracted_dir}")

    all_csv_files = sorted(input_dir.glob("*.csv"))
    csv_files = [f for f in all_csv_files if not is_population_file_stratified(f.name)]
    pop_files = [f for f in all_csv_files if is_population_file_stratified(f.name)]

    print(f"\nFound {len(csv_files)} death files, {len(pop_files)} population files")

    # Phase 1: Validation
    print("\n" + "=" * 70)
    print("PHASE 1 (Stratified): Data Integrity Validation")
    print("=" * 70)

    validation_results = []
    for file_path in csv_files:
        print(f"\nValidating: {file_path.name}")
        result = validate_file(file_path)
        validation_results.append(result)

        if result["overall_valid"]:
            print("  ✓ PASSED")
        else:
            print("  ✗ FAILED")
            for error in result["errors"]:
                print(f"    - {error}")

    print("\n" + "=" * 70)
    print("Validation Summary")
    print("=" * 70)

    total_files = len(validation_results)
    passed_files = sum(1 for r in validation_results if r["overall_valid"])
    failed_files = total_files - passed_files

    print(f"Total files:  {total_files}")
    pct = (passed_files / total_files * 100) if total_files > 0 else 0.0
    print(f"Passed:       {passed_files} ({pct:.1f}%)")
    print(f"Failed:       {failed_files} ({(100 - pct):.1f}%)")

    if failed_files > 0:
        print("\nFiles with errors:")
        for result in validation_results:
            if not result["overall_valid"]:
                print(f"  - {result['filename']}")
                for error in result["errors"]:
                    print(f"      {error}")

    # Phase 2: Merging
    print("\n" + "=" * 70)
    print("PHASE 2 (Stratified): Merge Files by Base Name (year, ICD set, stratum)")
    print("=" * 70)

    file_groups = defaultdict(list)
    for result in validation_results:
        if result["overall_valid"] and result["parsed"]:
            filename = result["filename"]
            file_path = Path(result["path"])
            base_name = get_base_filename(filename)

            if base_name:
                file_groups[base_name].append(file_path)

    print(f"\nFound {len(file_groups)} unique base names to merge")

    for base_name, file_paths in sorted(file_groups.items()):
        print(f"\n{base_name}:")
        output_filename = f"{base_name}.csv"
        output_path = merged_dir / output_filename

        merge_file_group(file_paths, output_path)

    # Phase 3: Subtraction for all year ranges (including 2021-2024)
    # Requires {year} ALL_{stratum}.csv in Stratified_Merged for every stratum and year range.
    subtract_deaths_stratified(merged_dir, subtracted_dir)

    # Phase 3b: Attach population to 2021-2024 subtracted files (WONDER files lack population)
    attach_population_stratified(subtracted_dir, input_dir)

    # Phase 3c: Others = Total - Black - White (all year ranges)
    compute_others_race(subtracted_dir, nonstrat_subtracted_dir, subtracted_dir)

    # Phase 3d: Female = Total - Male (all year ranges)
    compute_female(subtracted_dir, nonstrat_subtracted_dir, subtracted_dir)

    print("\n" + "=" * 70)
    print("SUMMARY (Stratified)")
    print("=" * 70)
    print(f"Validated:      {total_files} death files ({len(pop_files)} population files skipped)")
    print(f"Merged outputs: {len(file_groups)} files")
    print(f"Merged dir:     {merged_dir}")
    print(f"Subtracted dir: {subtracted_dir}")
    print("\n✓ All stratified phases completed successfully!")


if __name__ == "__main__":
    main()
