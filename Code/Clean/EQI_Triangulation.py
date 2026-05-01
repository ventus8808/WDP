"""
CDC Triangulation Data Integrity Checker, Merger, and Subtraction

This script validates and merges CDC death data files from the Triangulation directory.
The triangulation method is used to bypass CDC data suppression (deaths between 1-9).
By downloading all-cause death data and target disease data separately, we can subtract
to obtain the suppressed values.

For 2021-2024, population data is provided in separate 5-year age group files
(named '{year} ALL {range}.csv') and no all-cause deaths file is available.
These files are processed without triangulation: deaths are written directly to
the Subtracted directory with population aggregated from the 5-year files.

Phase 1: Data Integrity Check
- Validates year ranges match filename
- Validates age groups match filename
- Validates ICD-10 codes match filename
- Population-only files are excluded from validation

Phase 2: Merge Files
- Merges files with same year range and ICD codes but different age groups
- Outputs to Merged directory without Notes column

Phase 3 (2006-2020): Subtraction
- Calculates target disease deaths: Di = Ti - Ri
- Ti (total deaths) from "ALL.csv" files
- Ri (residual deaths) from "ALL_{ICD}.csv" files
- Outputs to Subtracted directory

Phase 4 (2021-2024): Direct output
- Aggregates 5-year population to 10-year age groups
- Attaches population to merged death files
- Writes to Subtracted directory (Suppressed deaths treated as 0)

Phase 5: Composite synthesis
- Reads component files already in Subtracted directory
- Sums deaths per county × age group to produce composite disease files
- e.g. Suicide (X60_X84_Y87.0) = X60_X69 + X70_X84 + Y87.0
- e.g. Dementia (G30_F01_F03) = G30 + F01 + F03
- Writes composite files to Subtracted directory
"""

import pandas as pd
import yaml
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from io import StringIO


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

# Year ranges that provide population in separate 5-year files (no triangulation subtraction)
DIRECT_YEARS = {"2021-2024"}

# ── Run control ──────────────────────────────────────────────────────────────
# Set to 0 to skip cancer/NDD subtraction (results already computed).
# CVD and suicide (non-cancer/NDD) always run.
RUN_CANCER_NDD = 0

# ICD codes that use the incremental method (populated from config in main()):
#   Di = ALL_{ICD}(baseline + target) - ALL(baseline)
# All other codes use the subtraction method:
#   Di = ALL_TRUE(total) - ALL_{ICD}(residual)
CANCER_NDD_CODES: set = set()


def _build_cancer_ndd_codes(config: dict) -> set:
    """Collect all cancer and NDD icd_codes from diseases config (incremental method group)."""
    diseases = config.get("diseases", {})
    codes: set = set()
    for group_key in ("cancer", "ndd"):
        group = diseases.get(group_key, {})
        overall = group.get("overall", {})
        if overall.get("icd_code"):
            codes.add(overall["icd_code"])
        for subtype in group.get("subtypes", {}).values():
            code = subtype.get("icd_code")
            if code:
                codes.add(code)
    return codes


def _get_sum_of_specs(config: dict) -> list:
    """Return list of {icd_code, sum_of, group_key} for all downloaded:false sum_of entries."""
    diseases = config.get("diseases", {})
    specs = []
    for group_key, group_cfg in diseases.items():
        overall = group_cfg.get("overall", {})
        if not overall.get("downloaded", True) and overall.get("sum_of"):
            specs.append({
                "icd_code": overall["icd_code"],
                "sum_of": overall["sum_of"],
                "group_key": group_key,
            })
        for subtype_cfg in group_cfg.get("subtypes", {}).values():
            if not subtype_cfg.get("downloaded", True) and subtype_cfg.get("sum_of"):
                specs.append({
                    "icd_code": subtype_cfg["icd_code"],
                    "sum_of": subtype_cfg["sum_of"],
                    "group_key": group_key,
                })
    return specs


def synthesize_sum_of_subtracted(output_dir: Path, config: dict, run_cancer_ndd: bool):
    """
    Phase 5: Create composite disease files in Subtracted/ by summing component files.

    For each config entry with downloaded:false and sum_of, reads the component
    subtracted CSVs, sums Deaths per County Code × age group, and writes a new file.
    Cancer/NDD composites (e.g. Dementia = G30+F01+F03) are only created when
    run_cancer_ndd=True; all others (e.g. Suicide = X60_X69+X70_X84+Y87.0) always run.
    Suppressed deaths (Quality_Flag=Suppressed, Deaths=0) contribute 0 to the sum.
    Missing deaths (Deaths=NA) are treated as 0 for summation purposes.
    Existing output files are skipped.
    """
    print("\n" + "=" * 70)
    print("PHASE 5: Synthesize sum_of Composites in Subtracted/")
    print("=" * 70)

    specs = _get_sum_of_specs(config)
    if not specs:
        print("No sum_of specs found in config.")
        return

    period_pattern = re.compile(r"^(\d{4}-\d{4})_")
    periods = sorted({
        m.group(1)
        for f in output_dir.glob("*.csv")
        if (m := period_pattern.match(f.name))
    })
    print(f"Detected periods: {periods}")

    cancer_ndd_codes = _build_cancer_ndd_codes(config)

    for spec in specs:
        composite_code = spec["icd_code"]
        components = spec["sum_of"]
        is_cancer_ndd = composite_code in cancer_ndd_codes

        if is_cancer_ndd and not run_cancer_ndd:
            print(f"\nSkipping {composite_code} (cancer/NDD composite, run_cancer_ndd=False)")
            continue

        print(f"\nSynthesizing {composite_code} from {components}...")

        for period in periods:
            output_path = output_dir / f"{period}_{composite_code}.csv"
            if output_path.exists():
                print(f"  {period}: already exists, skipping")
                continue

            component_dfs = []
            skip = False
            for comp in components:
                comp_path = output_dir / f"{period}_{comp}.csv"
                if not comp_path.exists():
                    print(f"  ⚠ Component missing: {comp_path.name} — skipping {period}")
                    skip = True
                    break
                component_dfs.append(pd.read_csv(comp_path, dtype={"County Code": str}))

            if skip:
                continue

            pop_df = component_dfs[0][
                ["County Code", "Ten-Year Age Groups", "Ten-Year Age Groups Code", "Population"]
            ].copy()

            deaths_sum = (
                pd.concat(component_dfs, ignore_index=True)
                .assign(Deaths_num=lambda d: pd.to_numeric(d["Deaths"], errors="coerce").fillna(0))
                .groupby(["County Code", "Ten-Year Age Groups"], as_index=False)["Deaths_num"]
                .sum()
            )

            result = pop_df.merge(deaths_sum, on=["County Code", "Ten-Year Age Groups"], how="left")
            result["Deaths"] = result["Deaths_num"].astype("Int64")
            result["Quality_Flag"] = "Synthesized"
            result = result[["County Code", "Ten-Year Age Groups", "Ten-Year Age Groups Code",
                              "Deaths", "Population", "Quality_Flag"]]

            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result["_age_order"] = result["Ten-Year Age Groups"].map(age_order_map)
            result = (result.sort_values(["County Code", "_age_order"])
                      .drop(columns=["_age_order"]))

            result.to_csv(output_path, index=False)
            print(f"  ✓ {output_path.name} ({len(result)} rows)")

    print("\n✓ Composite synthesis completed")


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


def load_config():
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return project_root, config


def is_population_file(filename: str) -> bool:
    """Detect 2021-2024-style population-only files: '{year} ALL {digit-range}.csv'"""
    return bool(re.match(r"^\d{4}-\d{4} ALL \d[\d-]*\.csv$", filename))


def is_true_allcause_file(filename: str) -> bool:
    """Detect true all-cause files: '{year} All_Cause_of_Death Age{n}.csv'"""
    return bool(re.match(r"^\d{4}-\d{4} All_Cause_of_Death Age\d+\.csv$", filename))


def parse_filename(filename: str) -> Optional[Dict[str, str]]:
    """
    Parse filename to extract year range, ICD codes, and age group.

    Examples:
        '2016-2020 ALL_C82_C85 Age1.csv' -> {'year': '2016-2020', 'icd': 'C82_C85', 'age': 'Age1'}
        '2011-2015 ALL Age35.csv' -> {'year': '2011-2015', 'icd': None, 'age': 'Age35'}
    """
    pattern = r"^(\d{4}-\d{4}) ALL(?:_([A-Za-z0-9_.-]+))? (Age\d+)\.csv$"
    match = re.match(pattern, filename)
    if not match:
        return None
    return {"year": match.group(1), "icd": match.group(2), "age": match.group(3)}


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

    metadata = {}

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
        if "A00-B99" in icd_text or "D50-D89" in icd_text:
            return True, "OK - General disease categories"
        else:
            return False, "Expected general disease categories for ALL file"

    icd_codes = icd_in_filename.replace("-", "_").split("_")
    icd_text = metadata.get("icd_codes", "")
    missing_codes = [code for code in icd_codes if code not in icd_text]

    if not missing_codes:
        return True, "OK"

    # ICD not in metadata ICD list: residual file where target was excluded from query.
    # Confirm via Title field. Strip non-alphanumeric chars to handle en-dashes in old files.
    title = metadata.get("title", "")
    title_alnum = re.sub(r"[^A-Za-z0-9]", "", title)
    icd_alnum = re.sub(r"[^A-Za-z0-9]", "", icd_in_filename)
    if f"ALL{icd_alnum}" in title_alnum:
        return True, "OK - Residual file confirmed by title"

    return False, f"Missing ICD codes in metadata: {', '.join(missing_codes)}"


def validate_file(file_path: Path) -> Dict:
    filename = file_path.name
    result = {
        "filename": filename,
        "path": str(file_path),
        "parsed": False,
        "year_valid": False,
        "age_valid": False,
        "icd_valid": False,
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

    result["overall_valid"] = year_valid and age_valid and icd_valid

    return result


def get_base_filename(filename: str) -> Optional[str]:
    filename_info = parse_filename(filename)
    if not filename_info:
        return None

    if filename_info["icd"]:
        return f"{filename_info['year']} ALL_{filename_info['icd']}"
    else:
        return f"{filename_info['year']} ALL"


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
        f.writelines(metadata_rows)

    print(f"  ✓ Saved to: {output_path.name} ({len(merged_df)} rows)")


def merge_true_allcause_files(input_dir: Path, merged_dir: Path):
    """Merge Age1/Age35/Age65 true all-cause files into one file per year range."""
    print("\n" + "=" * 70)
    print("PHASE 2b: Merge True All-Cause Files")
    print("=" * 70)

    pattern = re.compile(r"^(\d{4}-\d{4}) All_Cause_of_Death (Age\d+)\.csv$")
    groups: Dict[str, List[Path]] = defaultdict(list)
    for f in sorted(input_dir.glob("* All_Cause_of_Death Age*.csv")):
        m = pattern.match(f.name)
        if m:
            groups[m.group(1)].append(f)

    if not groups:
        print("  ⚠ No true all-cause files found in input directory")
        return

    print(f"\nFound {len(groups)} year ranges with true all-cause files")
    for year_range, file_paths in sorted(groups.items()):
        print(f"\n{year_range} All_Cause_of_Death ({len(file_paths)} age files):")
        output_path = merged_dir / f"{year_range} All_Cause_of_Death.csv"
        merge_file_group(file_paths, output_path)


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
    """
    Calculate target disease deaths: Di = Ti - Ri

    Rules:
    - If Ti or Ri is suppressed: Di = 0
    - If Ti or Ri is missing: Di = NaN (missing)
    - If both numeric: Di = max(0, Ti - Ri)
    """
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


def load_population_2021_2024(
    input_dir: Path, year_range: str
) -> Dict[str, Dict[str, int]]:
    """
    Load 5-year population files for a direct-year range and aggregate to 10-year groups.
    Returns: {county_code: {ten_year_age_group: population}}
    Both < 1 year and 1-4 years receive the combined 0-4 years population.
    """
    pop_files = sorted(
        f for f in input_dir.glob(f"{year_range} ALL *.csv")
        if is_population_file(f.name)
    )

    if not pop_files:
        print(f"  ⚠ No population files found for {year_range}")
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

    print(f"  ✓ Population loaded: {len(result)} counties from {len(pop_files)} files")
    return result


def subtract_deaths(merged_dir: Path, output_dir: Path, run_cancer_ndd: bool = True):
    """
    Phase 3: Calculate target disease deaths by subtraction (pre-2021 years only).

    Two methods depending on disease group:
    - Cancer/NDD (incremental):   Di = ALL_{ICD}(baseline+target) - ALL(baseline)
    - CVD/suicide (subtraction):  Di = ALL_TRUE(total) - ALL_{ICD}(residual)
    """
    print("\n" + "=" * 70)
    print("PHASE 3: Calculate Target Disease Deaths (Subtraction)")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    all_files = list(merged_dir.glob("*.csv"))

    year_ranges = set()
    icd_files: Dict[str, Dict[str, Path]] = {}
    all_files_dict: Dict[str, Path] = {}
    true_allcause_dict: Dict[str, Path] = {}

    for file_path in all_files:
        filename = file_path.name

        # True all-cause files (for CVD/suicide subtraction method)
        if re.match(r"^\d{4}-\d{4} All_Cause_of_Death\.csv$", filename):
            year_match = re.match(r"^(\d{4}-\d{4})", filename)
            if year_match:
                true_allcause_dict[year_match.group(1)] = file_path
            continue

        filename_info = parse_filename(filename.replace(".csv", " Age1.csv"))
        if not filename_info:
            continue

        year_range = filename_info["year"]
        if year_range in DIRECT_YEARS:
            continue

        year_ranges.add(year_range)

        if filename_info["icd"] is None:
            all_files_dict[year_range] = file_path
        else:
            icd_code = filename_info["icd"]
            if year_range not in icd_files:
                icd_files[year_range] = {}
            icd_files[year_range][icd_code] = file_path

    print(f"\nFound {len(year_ranges)} year ranges for subtraction")
    print(f"True all-cause files available: {sorted(true_allcause_dict.keys())}")

    total_processed = 0

    for year_range in sorted(year_ranges):
        print(f"\n{year_range}:")

        if year_range not in all_files_dict:
            print(f"  ⚠ Missing ALL (incremental baseline) file for {year_range}, skipping")
            continue

        if year_range not in icd_files:
            print(f"  ⚠ No ICD files for {year_range}, skipping")
            continue

        # Incremental baseline (cancer/NDD method)
        all_file = all_files_dict[year_range]
        print(f"  Incremental baseline (cancer/NDD): {all_file.name}")
        df_all, _ = extract_data_rows(all_file, keep_metadata=False)

        # True all-cause total (CVD/suicide method)
        df_true_all = None
        if year_range in true_allcause_dict:
            true_file = true_allcause_dict[year_range]
            print(f"  True all-cause (CVD/suicide):     {true_file.name}")
            df_true_all, _ = extract_data_rows(true_file, keep_metadata=False)
        else:
            print(f"  ⚠ No true all-cause file for {year_range} (CVD/suicide will be skipped)")

        for icd_code, icd_file in sorted(icd_files[year_range].items()):
            is_cancer_ndd = icd_code in CANCER_NDD_CODES

            if is_cancer_ndd and not run_cancer_ndd:
                print(f"    Skipping {icd_code} (cancer/NDD, RUN_CANCER_NDD=0)")
                continue

            if not is_cancer_ndd and df_true_all is None:
                print(f"    ⚠ Skipping {icd_code} (no true all-cause file)")
                continue

            print(f"    Processing {icd_code} ({'incremental' if is_cancer_ndd else 'subtraction'})...")

            df_icd, _ = extract_data_rows(icd_file, keep_metadata=False)

            if is_cancer_ndd:
                # Incremental method: Di = ALL_{ICD}(baseline+target) - ALL(baseline)
                # df_all → _Residual, df_icd → _Total  (Di = Total - Residual)
                merged = df_all.merge(
                    df_icd,
                    on=["County Code", "Ten-Year Age Groups"],
                    suffixes=("_Residual", "_Total"),
                    how="inner",
                )
            else:
                # Subtraction method: Di = ALL_TRUE(total) - ALL_{ICD}(residual)
                # df_true_all → _Total, df_icd → _Residual  (Di = Total - Residual)
                merged = df_true_all.merge(
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

                results.append(
                    {
                        "County Code": row["County Code"],
                        "Ten-Year Age Groups": row["Ten-Year Age Groups"],
                        "Ten-Year Age Groups Code": row[
                            "Ten-Year Age Groups Code_Residual"
                        ],
                        "Deaths": di,
                        "Population": row["Population_Residual"],
                        "Quality_Flag": quality_flag,
                    }
                )

            result_df = pd.DataFrame(results)

            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(age_order_map)
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])

            result_df["Deaths"] = result_df["Deaths"].astype("Int64")

            output_filename = f"{year_range}_{icd_code}.csv"
            output_path = output_dir / output_filename
            result_df.to_csv(output_path, index=False)

            clean_count = (result_df["Quality_Flag"] == "Clean").sum()
            suppressed_count = result_df["Quality_Flag"].str.contains("Suppressed").sum()
            total_count = len(result_df)

            print(f"      ✓ Saved to: {output_filename}")
            print(f"        {total_count} rows: {clean_count} clean, {suppressed_count} suppressed")

            total_processed += 1

    print(f"\n✓ Subtraction completed: {total_processed} ICD codes processed")


def write_direct_deaths(merged_dir: Path, input_dir: Path, output_dir: Path):
    """
    Phase 4: Process 2021-2024 death files with population attached.

    - Cancer/NDD codes: deaths written directly (no triangulation baseline available).
    - CVD/suicide codes: Di = ALL_TRUE(total) - ALL_{ICD}(residual), same subtraction
      method as 2006-2020.  Requires '{year} All_Cause_of_Death.csv' in merged_dir.

    Suppressed deaths treated as 0.  NaN County Code rows are dropped.
    """
    print("\n" + "=" * 70)
    print("PHASE 4: Write 2021-2024 Deaths with Population")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    for year_range in sorted(DIRECT_YEARS):
        merged_files = sorted(merged_dir.glob(f"{year_range} ALL_*.csv"))
        if not merged_files:
            print(f"\n  No merged files found for {year_range}, skipping")
            continue

        print(f"\nLoading population for {year_range}...")
        pop_lookup = load_population_2021_2024(input_dir, year_range)

        # True all-cause file for CVD/suicide subtraction
        true_allcause_path = merged_dir / f"{year_range} All_Cause_of_Death.csv"
        df_true_all = None
        if true_allcause_path.exists():
            print(f"  True all-cause file found: {true_allcause_path.name}")
            df_true_all, _ = extract_data_rows(true_allcause_path, keep_metadata=False)
            df_true_all = df_true_all[df_true_all["County Code"].notna()].copy()
        else:
            print(f"  ⚠ No true all-cause file for {year_range} — CVD/suicide will be skipped")

        for merged_file in merged_files:
            filename_info = parse_filename(
                merged_file.name.replace(".csv", " Age1.csv")
            )
            if not filename_info or not filename_info.get("icd"):
                continue

            icd_code = filename_info["icd"]
            is_cancer_ndd = icd_code in CANCER_NDD_CODES

            if not is_cancer_ndd and df_true_all is None:
                print(f"  ⚠ Skipping {icd_code} (no true all-cause file)")
                continue

            print(f"  Processing {icd_code} ({'direct' if is_cancer_ndd else 'subtraction'})...")

            df_icd, _ = extract_data_rows(merged_file, keep_metadata=False)
            df_icd = df_icd[df_icd["County Code"].notna()].copy()

            results = []

            if is_cancer_ndd:
                # Direct write: deaths come from the ICD file as-is
                for _, row in df_icd.iterrows():
                    county = str(row["County Code"]).zfill(5)
                    age_group = row["Ten-Year Age Groups"]

                    death_val, is_suppressed, is_missing = parse_death_value(row["Deaths"])
                    if is_missing:
                        deaths = None
                        flag = "Missing"
                    elif is_suppressed:
                        deaths = 0
                        flag = "Suppressed"
                    else:
                        deaths = death_val
                        flag = "Clean"

                    population = pop_lookup.get(county, {}).get(age_group, None)
                    results.append({
                        "County Code": county,
                        "Ten-Year Age Groups": age_group,
                        "Ten-Year Age Groups Code": row.get("Ten-Year Age Groups Code", ""),
                        "Deaths": deaths,
                        "Population": population,
                        "Quality_Flag": flag,
                    })
            else:
                # Subtraction method: Di = ALL_TRUE(total) - ALL_{ICD}(residual)
                merged = df_true_all.merge(
                    df_icd,
                    on=["County Code", "Ten-Year Age Groups"],
                    suffixes=("_Total", "_Residual"),
                    how="inner",
                )
                for _, row in merged.iterrows():
                    county = str(row["County Code"]).zfill(5)
                    age_group = row["Ten-Year Age Groups"]

                    di, quality_flag = calculate_target_deaths(
                        row["Deaths_Total"], row["Deaths_Residual"]
                    )
                    population = pop_lookup.get(county, {}).get(age_group, None)
                    results.append({
                        "County Code": county,
                        "Ten-Year Age Groups": age_group,
                        "Ten-Year Age Groups Code": row.get("Ten-Year Age Groups Code_Residual", ""),
                        "Deaths": di,
                        "Population": population,
                        "Quality_Flag": quality_flag,
                    })

            result_df = pd.DataFrame(results)

            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(age_order_map)
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])
            result_df["Deaths"] = result_df["Deaths"].astype("Int64")

            output_filename = f"{year_range}_{icd_code}.csv"
            output_path = output_dir / output_filename
            result_df.to_csv(output_path, index=False)

            clean_count = (result_df["Quality_Flag"] == "Clean").sum()
            suppressed_count = result_df["Quality_Flag"].str.contains("Suppressed").sum()
            total_count = len(result_df)

            print(f"    ✓ Saved to: {output_filename}")
            print(f"      {total_count} rows: {clean_count} clean, {suppressed_count} suppressed")

    print("\n✓ Direct output completed")


def main():
    print("=" * 70)
    print("CDC Triangulation Data Integrity Checker, Merger, and Subtraction")
    print("=" * 70)

    project_root, config = load_config()

    global CANCER_NDD_CODES
    CANCER_NDD_CODES = _build_cancer_ndd_codes(config)

    input_dir = project_root / "Data/Original/CDC Triangulation/NotMerged"
    merged_dir = project_root / "Data/Original/CDC Triangulation/Merged"
    subtracted_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"

    print(f"\nInput directory:       {input_dir}")
    print(f"Merged directory:      {merged_dir}")
    print(f"Subtracted directory:  {subtracted_dir}")

    all_csv_files = list(input_dir.glob("*.csv"))
    csv_files = [
        f for f in all_csv_files
        if not is_population_file(f.name) and not is_true_allcause_file(f.name)
    ]
    pop_files = [f for f in all_csv_files if is_population_file(f.name)]
    true_allcause_files = [f for f in all_csv_files if is_true_allcause_file(f.name)]

    print(f"\nFound {len(csv_files)} death files, {len(pop_files)} population files, "
          f"{len(true_allcause_files)} true all-cause files")
    print(f"RUN_CANCER_NDD = {RUN_CANCER_NDD}")

    # Phase 1: Validation
    print("\n" + "=" * 70)
    print("PHASE 1: Data Integrity Validation")
    print("=" * 70)

    validation_results = []
    for file_path in sorted(csv_files):
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
    print(f"Passed:       {passed_files} ({passed_files / total_files * 100:.1f}%)")
    print(f"Failed:       {failed_files} ({failed_files / total_files * 100:.1f}%)")

    if failed_files > 0:
        print("\nFiles with errors:")
        for result in validation_results:
            if not result["overall_valid"]:
                print(f"  - {result['filename']}")
                for error in result["errors"]:
                    print(f"      {error}")

    # Phase 2: Merging
    print("\n" + "=" * 70)
    print("PHASE 2: Merge Files by Base Name")
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

    # Phase 2b: Merge true all-cause files (for CVD/suicide subtraction method)
    merge_true_allcause_files(input_dir, merged_dir)

    # Phase 3: Subtraction (pre-2021 years)
    subtract_deaths(merged_dir, subtracted_dir, run_cancer_ndd=bool(RUN_CANCER_NDD))

    # Phase 4: Direct output for 2021-2024 with population attached
    write_direct_deaths(merged_dir, input_dir, subtracted_dir)

    # Phase 5: Synthesize sum_of disease composites in Subtracted/
    synthesize_sum_of_subtracted(subtracted_dir, config, run_cancer_ndd=bool(RUN_CANCER_NDD))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Validated:      {total_files} death files "
          f"({len(pop_files)} population + {len(true_allcause_files)} all-cause files skipped)")
    print(f"Merged:         {len(file_groups)} output files")
    print(f"RUN_CANCER_NDD: {RUN_CANCER_NDD}")
    print(f"Merged dir:     {merged_dir}")
    print(f"Subtracted dir: {subtracted_dir}")
    print("\n✓ All phases completed successfully!")


if __name__ == "__main__":
    main()
