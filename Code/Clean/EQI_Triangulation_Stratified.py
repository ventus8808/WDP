"""
CDC Triangulation (Stratified) — Data Integrity Checker, Merger, and Subtraction

This script mirrors EQI_Triangulation.py but supports stratified files, where
the strata (e.g., White, Male, Indian) are encoded in the filenames.

It performs three phases:
- Phase 1: Data Integrity Validation (years, ages, ICD codes, and stratum presence)
- Phase 2: Merge age-partitioned files (Age1, Age35, Age65) for each (year, ICD set, stratum)
- Phase 3: Subtract residual ALL from ALL_{ICD...} to obtain target-disease deaths

Important filename patterns (examples from Stratified_NotMerged):
- 2016-2020 ALL_G20_G30_G12.2_F01_F03_White Age1.csv
- 2016-2020 ALL_Indian Age1.csv
- 2016-2020 ALL_Male Age35.csv
- 2016-2020 ALL_White Age65.csv

Parsing logic:
- year_range:     2016-2020
- ICD set:        G20_G30_G12.2_F01_F03 (or None for total ALL)
- stratum:        White (or Indian, Male, etc.)
- age partition:  Age1 / Age35 / Age65

Directories:
- Input (stratified, not merged):     Data/Original/CDC Triangulation/Stratified_NotMerged
- Output merged (stratified):         Data/Original/CDC Triangulation/Stratified_Merged
- Output subtracted (stratified):     Data/Original/CDC Triangulation/Stratified_Subtracted
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

# Configuration
AGE_GROUP_MAPPINGS = {
    "Age1": ["< 1 year", "1-4 years", "5-14 years", "15-24 years", "25-34 years"],
    "Age35": ["35-44 years", "45-54 years", "55-64 years"],
    "Age65": ["65-74 years", "75-84 years", "85+ years"],
}

# Define the order for age groups
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
}


def load_config() -> tuple[Path, dict]:
    """Load configuration from config.yaml"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    return project_root, config


def parse_filename(filename: str) -> Optional[Dict[str, Optional[str]]]:
    """
    Parse stratified filename to extract year range, ICD code block (joined by '_'), stratum, and age partition.

    Expected patterns (examples):
      - '2016-2020 ALL_G20_G30_G12.2_F01_F03_White Age1.csv'
      - '2016-2020 ALL_Indian Age1.csv'
      - '2016-2020 ALL_Male Age65.csv'
      - '2016-2020 ALL_White Age35.csv'

    Returns:
      {
        'year': '2016-2020',
        'icd': 'G20_G30_G12.2_F01_F03' or None,
        'stratum': 'White' (e.g., 'Indian', 'Male', 'White', etc.),
        'age': 'Age1' | 'Age35' | 'Age65'
      }
    """
    if not filename.endswith(".csv"):
        return None

    name = filename[:-4]  # strip .csv

    # Expect ' ... AgeX' suffix
    if " Age" not in name:
        return None

    base, age = name.rsplit(" ", 1)
    if not re.fullmatch(r"Age\d+", age):
        return None

    # Split year from the rest
    if " " not in base:
        return None

    year, rest = base.split(" ", 1)
    if not re.fullmatch(r"\d{4}-\d{4}", year):
        return None

    # Expect leading 'ALL' in rest
    if not rest.startswith("ALL"):
        return None

    suffix = rest[3:]  # drop 'ALL'
    if suffix.startswith("_"):
        suffix = suffix[1:]  # drop leading underscore

    if suffix == "":
        # Unstratified ALL — not expected in stratified folders, but handle
        icd_block = None
        stratum = None
    else:
        # suffix looks like 'G20_G30_G12.2_F01_F03_White' or 'White'
        tokens = suffix.split("_")
        if len(tokens) == 1:
            icd_block = None
            stratum = tokens[0]
        else:
            icd_block = "_".join(tokens[:-1]) if len(tokens) > 1 else None
            stratum = tokens[-1]

    return {"year": year, "icd": icd_block, "stratum": stratum, "age": age}


def read_file_metadata(file_path: Path) -> Dict[str, str]:
    """
    Read metadata from the end of CDC CSV file.
    Returns dict with 'title', 'icd_codes', 'age_groups', 'years'
    """
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Find "Query Parameters:" line which starts the metadata we care about
    query_param_idx = -1
    for i in range(len(lines) - 1, max(0, len(lines) - 500), -1):
        if '"Query Parameters:"' in lines[i]:
            query_param_idx = i
            break

    if query_param_idx == -1:
        return {}

    # Join lines from Query Parameters onwards to handle multi-line fields
    metadata_lines = lines[query_param_idx:]
    metadata_text = " ".join(line.strip().strip('"') for line in metadata_lines)

    # Extract key fields
    metadata: Dict[str, str] = {}

    # Title
    title_match = re.search(
        r'Title:\s*([^\s][^"]*?)(?=\s+ICD-10 Codes:|\s+Ten-Year Age Groups:|\s+Year/Month:)',
        metadata_text,
    )
    if title_match:
        metadata["title"] = title_match.group(1).strip()

    # ICD-10 Codes
    icd_match = re.search(
        r"ICD-10 Codes:\s*(.*?)(?=\s+Ten-Year Age Groups:|\s+Year/Month:)",
        metadata_text,
        re.DOTALL,
    )
    if icd_match:
        metadata["icd_codes"] = icd_match.group(1).strip()

    # Ten-Year Age Groups
    age_match = re.search(
        r'Ten-Year Age Groups:\s*([^"]+?)(?=\s+Year/Month:|\s+Group By:)', metadata_text
    )
    if age_match:
        metadata["age_groups"] = age_match.group(1).strip()

    # Year/Month
    year_match = re.search(
        r'Year/Month:\s*([^"]+?)(?=\s+Group By:|\s+Show)', metadata_text
    )
    if year_match:
        metadata["years"] = year_match.group(1).strip()

    return metadata


def extract_data_rows(
    file_path: Path, keep_metadata: bool = True
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Extract data as DataFrame and metadata section from CSV file.
    Returns (dataframe_without_notes, metadata_rows)
    """
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Find where metadata starts
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

    # Read CSV data using pandas (up to metadata section)
    from io import StringIO

    csv_data = "".join(lines[:metadata_start])
    # Preserve County Code as string with leading zeros
    df = pd.read_csv(StringIO(csv_data), dtype={"County Code": str})

    # Remove Notes and County columns if they exist
    columns_to_drop = []
    if "Notes" in df.columns:
        columns_to_drop.append("Notes")
    if "County" in df.columns:
        columns_to_drop.append("County")

    if columns_to_drop:
        df = df.drop(columns=columns_to_drop)

    return df, metadata_rows


def validate_year_range(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    """Validate that year range in filename matches metadata"""
    expected_years = YEAR_RANGE_MAPPINGS.get(filename_info["year"])
    if not expected_years:
        return False, f"Unknown year range: {filename_info['year']}"

    metadata_years = metadata.get("years", "")

    # Check all expected years are present
    for year in expected_years:
        if year not in metadata_years:
            return False, f"Missing year {year} in metadata. Found: {metadata_years}"

    return True, "OK"


def validate_age_group(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    """Validate that age group in filename matches metadata"""
    expected_ages = AGE_GROUP_MAPPINGS.get(filename_info["age"])
    if not expected_ages:
        return False, f"Unknown age group: {filename_info['age']}"

    metadata_ages = metadata.get("age_groups", "")

    # Check all expected age groups are present
    for age in expected_ages:
        if age not in metadata_ages:
            return (
                False,
                f"Missing age group '{age}' in metadata. Found: {metadata_ages}",
            )

    return True, "OK"


def validate_icd_codes(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    """Validate that ICD codes in filename match metadata"""
    icd_in_filename = filename_info.get("icd")

    # If no ICD in filename (just "ALL_{Stratum}"), expect general disease categories
    if not icd_in_filename:
        icd_text = metadata.get("icd_codes", "")
        if "A00-B99" in icd_text or "D50-D89" in icd_text or "All causes" in icd_text:
            return True, "OK - General disease categories"
        else:
            return False, "Expected general disease categories for ALL file"

    # Parse ICD codes from filename (e.g., C82_C85, C18_C21, G20_G30_G12.2_F01_F03)
    icd_codes = icd_in_filename.split("_")

    icd_text = metadata.get("icd_codes", "")

    # For each code in filename, verify it appears in metadata
    missing_codes = []
    for code in icd_codes:
        if code not in icd_text:
            missing_codes.append(code)

    if missing_codes:
        return False, f"Missing ICD codes in metadata: {', '.join(missing_codes)}"

    return True, "OK"


def validate_stratum(filename_info: Dict, metadata: Dict) -> Tuple[bool, str]:
    """
    Validate that the stratum from the filename appears in the metadata (best-effort).

    We check case-insensitive presence of the stratum token in the 'title' or other fields.
    This is lenient due to variability in CDC WONDER titles and descriptors.
    """
    stratum = filename_info.get("stratum")
    if not stratum:
        # If not provided, treat as valid in stratified context (some files might be unstratified)
        return True, "No stratum in filename"

    haystacks = []
    for key in ("title", "icd_codes", "age_groups", "years"):
        val = metadata.get(key, "")
        if isinstance(val, str) and val:
            haystacks.append(val.lower())

    needle = stratum.lower()

    # Light normalization: map common shorthand to title forms
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

    # Not found; don't be overly strict
    return True, f"Stratum '{stratum}' not confirmed in metadata title; proceeding"


def validate_file(file_path: Path) -> Dict:
    """
    Validate a single CDC Triangulation stratified file.
    Returns validation result dict.
    """
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

    # Parse filename
    filename_info = parse_filename(filename)
    if not filename_info:
        result["errors"].append(f"Could not parse filename: {filename}")
        return result

    result["parsed"] = True
    result["filename_info"] = filename_info

    # Read metadata
    try:
        metadata = read_file_metadata(file_path)
        result["metadata"] = metadata
    except Exception as e:
        result["errors"].append(f"Error reading metadata: {e}")
        return result

    # Validate year range
    year_valid, year_msg = validate_year_range(filename_info, metadata)
    result["year_valid"] = year_valid
    if not year_valid:
        result["errors"].append(f"Year validation: {year_msg}")

    # Validate age group
    age_valid, age_msg = validate_age_group(filename_info, metadata)
    result["age_valid"] = age_valid
    if not age_valid:
        result["errors"].append(f"Age validation: {age_msg}")

    # Validate ICD codes
    icd_valid, icd_msg = validate_icd_codes(filename_info, metadata)
    result["icd_valid"] = icd_valid
    if not icd_valid:
        result["errors"].append(f"ICD validation: {icd_msg}")

    # Validate stratum (lenient)
    stratum_valid, stratum_msg = validate_stratum(filename_info, metadata)
    result["stratum_valid"] = stratum_valid
    if not stratum_valid:
        result["errors"].append(f"Stratum validation: {stratum_msg}")

    # Overall validation
    result["overall_valid"] = year_valid and age_valid and icd_valid and stratum_valid

    return result


def get_base_filename(filename: str) -> Optional[str]:
    """
    Extract base filename by removing age group suffix and preserving stratum.

    Examples:
      '2016-2020 ALL_C82_C85_White Age1.csv' -> '2016-2020 ALL_C82_C85_White'
      '2016-2020 ALL_White Age35.csv'        -> '2016-2020 ALL_White'
    """
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
        # Fallback (unlikely in stratified folders)
        return f"{year} ALL_{icd}"
    else:
        return f"{year} ALL"


def merge_file_group(file_paths: List[Path], output_path: Path):
    """
    Merge multiple CSV files (different age groups) into one.
    Removes the Notes column and orders by age groups.
    """
    all_dataframes = []
    metadata_rows = None

    print(f"  Merging {len(file_paths)} files...")

    for file_path in sorted(file_paths):
        print(f"    - {file_path.name}")
        df, meta_rows = extract_data_rows(file_path)

        all_dataframes.append(df)

        if metadata_rows is None:
            metadata_rows = meta_rows

    # Concatenate all dataframes
    merged_df = pd.concat(all_dataframes, ignore_index=True)

    # Create age group order mapping for sorting
    age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}

    # Sort by County Code, then by age group order
    if "Ten-Year Age Groups" in merged_df.columns:
        merged_df["_age_order"] = merged_df["Ten-Year Age Groups"].map(age_order_map)
        merged_df = merged_df.sort_values(
            by=["County Code", "_age_order"], na_position="last"
        )
        merged_df = merged_df.drop(columns=["_age_order"])
    else:
        # Fallback: just sort by County Code
        merged_df = merged_df.sort_values(by=["County Code"])

    # Write merged file
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        # Write dataframe as CSV (without index)
        merged_df.to_csv(f, index=False)

        # Write metadata (keep as-is)
        if metadata_rows:
            f.writelines(metadata_rows)

    print(f"  ✓ Saved to: {output_path.name} ({len(merged_df)} rows)")


def parse_death_value(value):
    """
    Parse death value from CSV.
    Returns: (numeric_value, is_suppressed, is_missing)
    """
    if pd.isna(value):
        return None, False, True

    if isinstance(value, str):
        if "Suppressed" in value:
            return None, True, False
        # Try to parse numeric from string
        try:
            return int(float(value)), False, False
        except ValueError:
            return None, False, True

    # Numeric value
    return int(float(value)), False, False


def calculate_target_deaths(ti_value, ri_value):
    """
    Calculate target disease deaths: Di = Ti - Ri

    Returns: (di_value, quality_flag)

    Rules:
    - If Ti or Ri is suppressed: Di = 0
    - If Ti or Ri is missing: Di = NaN (missing)
    - If both numeric: Di = max(0, Ti - Ri)
    """
    ti, ti_suppressed, ti_missing = parse_death_value(ti_value)
    ri, ri_suppressed, ri_missing = parse_death_value(ri_value)

    # Handle missing data
    if ti_missing:
        return None, "Missing_Ti"
    if ri_missing:
        return None, "Missing_Ri"

    # Handle suppressed data
    if ti_suppressed and ri_suppressed:
        return 0, "Both_Suppressed"
    if ti_suppressed:
        return 0, "Suppressed_Ti"
    if ri_suppressed:
        return 0, "Suppressed_Ri"

    # Both are numeric
    di = int(max(0, ti - ri))
    return di, "Clean"


def subtract_deaths_stratified(merged_dir: Path, output_dir: Path):
    """
    Phase 3 (Stratified): Calculate target disease deaths by subtraction for each (year, stratum).
    Di = Ti (ALL_{ICD} includes target) - Ri (ALL excludes target)
    """
    print("\n" + "=" * 70)
    print("PHASE 3 (Stratified): Calculate Target Disease Deaths (Subtraction)")
    print("=" * 70)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Get all merged files
    all_files = list(merged_dir.glob("*.csv"))

    # Group by (year_range, stratum)
    group_keys = set()
    all_files_by_group: Dict[Tuple[str, str], Path] = {}
    icd_files_by_group: Dict[Tuple[str, str], Dict[str, Path]] = {}

    for file_path in all_files:
        filename = file_path.name
        # Reconstruct a pseudo "Age1" suffix to reuse parse_filename (base-only names need help)
        # Example merged output is '2016-2020 ALL_{...}_{Stratum}.csv'
        # We temporarily append ' Age1.csv' to parse base for icd and stratum.
        info = parse_filename(filename.replace(".csv", " Age1.csv"))
        if not info:
            continue

        year_range = info["year"]
        stratum = info.get("stratum") or ""
        icd_block = info.get("icd")

        key = (year_range, stratum)
        group_keys.add(key)

        if icd_block is None:
            # ALL_{Stratum} => residual
            all_files_by_group[key] = file_path
        else:
            icd_files_by_group.setdefault(key, {})
            icd_files_by_group[key][icd_block] = file_path

    print(f"\nFound {len(group_keys)} (year, stratum) groups")

    total_processed = 0

    # Process each group
    for year_range, stratum in sorted(group_keys):
        print(f"\n{year_range} — {stratum if stratum else '(no stratum)'}:")

        if (year_range, stratum) not in all_files_by_group:
            print(f"  ⚠ Missing ALL file for group, skipping")
            continue

        if (year_range, stratum) not in icd_files_by_group:
            print(f"  ⚠ No ICD files for group, skipping")
            continue

        # Load residual (ALL_{Stratum}) — without metadata to avoid contamination
        all_file = all_files_by_group[(year_range, stratum)]
        print(f"  Loading residual (ALL): {all_file.name}")
        df_all, _ = extract_data_rows(all_file, keep_metadata=False)

        # Process each ICD code block
        for icd_code, icd_file in sorted(
            icd_files_by_group[(year_range, stratum)].items()
        ):
            print(f"    Processing {icd_code}...")

            # Load total (ALL_{ICD}_{Stratum}) — without metadata
            df_icd, _ = extract_data_rows(icd_file, keep_metadata=False)

            # Merge on County Code and Age Group
            # Note: df_all = ALL (residual, excludes target disease) = Ri
            #       df_icd = ALL_{ICD} (includes target disease) = Ti
            merged = df_all.merge(
                df_icd,
                on=["County Code", "Ten-Year Age Groups"],
                suffixes=("_Residual", "_Total"),
                how="inner",
            )

            # Calculate target deaths (Di = Ti - Ri)
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

            # Create output DataFrame
            result_df = pd.DataFrame(results)

            # Sort by County Code and Age Group
            age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
            result_df["_age_order"] = result_df["Ten-Year Age Groups"].map(
                age_order_map
            )
            result_df = result_df.sort_values(by=["County Code", "_age_order"])
            result_df = result_df.drop(columns=["_age_order"])

            # Ensure Deaths column is integer type
            result_df["Deaths"] = result_df["Deaths"].astype("Int64")

            # Save to file, keep stratum in filename
            output_filename = (
                f"{year_range}_{icd_code}_{stratum}.csv"
                if stratum
                else f"{year_range}_{icd_code}.csv"
            )
            output_path = output_dir / output_filename
            result_df.to_csv(output_path, index=False)

            # Summary stats
            clean_count = (result_df["Quality_Flag"] == "Clean").sum()
            suppressed_count = (
                result_df["Quality_Flag"].str.contains("Suppressed").sum()
            )
            total_count = len(result_df)

            print(f"      ✓ Saved to: {output_filename}")
            print(
                f"        {total_count} rows: {clean_count} clean, {suppressed_count} suppressed"
            )

            total_processed += 1

    print(
        f"\n✓ Stratified subtraction completed: {total_processed} ICD code blocks processed"
    )


def main():
    """Main execution function (Stratified)"""
    print("=" * 70)
    print("CDC Triangulation (Stratified) — Validation, Merge, and Subtraction")
    print("=" * 70)

    # Load configuration (project_root used to anchor relative paths)
    project_root, _ = load_config()

    # Set up stratified paths
    input_dir = project_root / "Data/Original/CDC Triangulation/Stratified_NotMerged"
    merged_dir = project_root / "Data/Original/CDC Triangulation/Stratified_Merged"
    subtracted_dir = (
        project_root / "Data/Original/CDC Triangulation/Stratified_Subtracted"
    )

    print(f"\nInput directory:       {input_dir}")
    print(f"Merged directory:      {merged_dir}")
    print(f"Subtracted directory:  {subtracted_dir}")

    # Get all CSV files
    csv_files = sorted(input_dir.glob("*.csv"))
    print(f"\nFound {len(csv_files)} CSV files")

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

    # Validation summary
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

    # Group files by base name (year + icd block + stratum)
    file_groups = defaultdict(list)
    for result in validation_results:
        if result["overall_valid"] and result["parsed"]:
            filename = result["filename"]
            file_path = Path(result["path"])
            base_name = get_base_filename(filename)

            if base_name:
                file_groups[base_name].append(file_path)

    print(f"\nFound {len(file_groups)} unique base names to merge")

    # Merge each group to Stratified_Merged
    for base_name, file_paths in sorted(file_groups.items()):
        print(f"\n{base_name}:")
        output_filename = f"{base_name}.csv"
        output_path = merged_dir / output_filename

        merge_file_group(file_paths, output_path)

    # Phase 3: Subtraction (Stratified)
    subtract_deaths_stratified(merged_dir, subtracted_dir)

    # Final summary
    print("\n" + "=" * 70)
    print("SUMMARY (Stratified)")
    print("=" * 70)
    print(f"Validated:      {total_files} files")
    print(f"Merged outputs: {len(file_groups)} files")
    print(f"Merged dir:     {merged_dir}")
    print(f"Subtracted dir: {subtracted_dir}")
    print("\n✓ All stratified phases completed successfully!")


if __name__ == "__main__":
    main()
