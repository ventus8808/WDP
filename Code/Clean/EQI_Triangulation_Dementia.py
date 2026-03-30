"""
CDC Dementia Composite Death Calculator

This script combines three individual ICD-10 code death files (G30, F01, F03)
from the Subtracted directory to produce a composite Dementia death count.

Source files (already produced by EQI_Triangulation.py):
    {year}_G30.csv   - Alzheimer's disease
    {year}_F01.csv   - Vascular dementia
    {year}_F03.csv   - Unspecified dementia

Method: Direct summation of Deaths across G30 + F01 + F03 per County Code and Age Group.

Output: {year}_G30_F01_F03.csv in the Subtracted directory (same format as other subtracted files)
"""

import pandas as pd
import yaml
from pathlib import Path


YEAR_RANGES = ["2006-2010", "2011-2015", "2016-2020"]
DEMENTIA_CODES = ["G30", "F01", "F03"]
OUTPUT_CODE = "G30_F01_F03"

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


def load_config():
    """Load configuration from config.yaml"""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    return project_root, config


def combine_quality_flags(flags: pd.Series) -> str:
    """
    Combine quality flags from multiple ICD components into one.
    Returns "Clean" only if all are "Clean", otherwise summarises the issues.
    """
    unique_flags = set(flags.dropna())
    if unique_flags == {"Clean"}:
        return "Clean"
    non_clean = [f for f in unique_flags if f != "Clean"]
    return "Mixed_" + "_".join(sorted(non_clean))


def combine_deaths(deaths: pd.Series) -> int | None:
    """
    Sum death counts across ICD codes.
    Treats NaN (missing) as NaN (propagates), otherwise sums integers.
    """
    if deaths.isna().any():
        return None
    return int(deaths.sum())


def build_dementia_csv(subtracted_dir: Path, year_range: str) -> pd.DataFrame:
    """
    Load G30, F01, F03 subtracted files for a given year range and sum Deaths.
    Returns a DataFrame in the same column format as other subtracted files.
    """
    merge_key = ["County Code", "Ten-Year Age Groups", "Ten-Year Age Groups Code"]
    frames = []
    for code in DEMENTIA_CODES:
        file_path = subtracted_dir / f"{year_range}_{code}.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"Required file not found: {file_path}")

        df = pd.read_csv(file_path, dtype={"County Code": str})
        # Drop trailing NaN rows (CDC files have empty footer rows)
        df = df[df["County Code"].notna()].reset_index(drop=True)
        df = df.rename(columns={
            "Deaths": f"Deaths_{code}",
            "Quality_Flag": f"Quality_Flag_{code}",
        })
        frames.append(df)

    # Merge all three on County Code + Age Group; keep Population from G30
    merged = frames[0]
    for df in frames[1:]:
        merged = merged.merge(
            df[[c for c in df.columns if c != "Population"]],
            on=merge_key,
            how="inner",
        )

    # Sum deaths and combine quality flags row by row
    death_cols = [f"Deaths_{c}" for c in DEMENTIA_CODES]
    flag_cols = [f"Quality_Flag_{c}" for c in DEMENTIA_CODES]

    merged["Deaths"] = merged[death_cols].apply(
        lambda row: combine_deaths(row), axis=1
    )
    merged["Quality_Flag"] = merged[flag_cols].apply(
        lambda row: combine_quality_flags(row), axis=1
    )

    result = merged[[
        "County Code",
        "Ten-Year Age Groups",
        "Ten-Year Age Groups Code",
        "Deaths",
        "Population",
        "Quality_Flag",
    ]].copy()

    # Sort by County Code then Age Group order
    age_order_map = {age: idx for idx, age in enumerate(AGE_GROUP_ORDER)}
    result["_age_order"] = result["Ten-Year Age Groups"].map(age_order_map)
    result = result.sort_values(by=["County Code", "_age_order"])
    result = result.drop(columns=["_age_order"])

    # Ensure Deaths is nullable integer
    result["Deaths"] = result["Deaths"].astype("Int64")

    return result


def main():
    """Main execution function"""
    print("=" * 70)
    print("CDC Dementia Composite Death Calculator (G30 + F01 + F03)")
    print("=" * 70)

    project_root, _ = load_config()
    subtracted_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"

    print(f"\nSubtracted directory: {subtracted_dir}")
    print(f"ICD codes to combine: {' + '.join(DEMENTIA_CODES)}")
    print(f"Output code:          {OUTPUT_CODE}")

    for year_range in YEAR_RANGES:
        print(f"\n{year_range}:")

        try:
            result_df = build_dementia_csv(subtracted_dir, year_range)
        except FileNotFoundError as e:
            print(f"  ✗ Skipped: {e}")
            continue

        output_filename = f"{year_range}_{OUTPUT_CODE}.csv"
        output_path = subtracted_dir / output_filename
        result_df.to_csv(output_path, index=False)

        clean_count = (result_df["Quality_Flag"] == "Clean").sum()
        mixed_count = result_df["Quality_Flag"].str.startswith("Mixed_").sum()
        total_count = len(result_df)

        print(f"  ✓ Saved: {output_filename}")
        print(f"    {total_count} rows: {clean_count} clean, {mixed_count} mixed")

    print("\n✓ Dementia composite files generated successfully!")


if __name__ == "__main__":
    main()
