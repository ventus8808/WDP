"""
CDC Triangulation AAMR Calculator

This script calculates Age-Adjusted Mortality Rate (AAMR) and confidence intervals
using the triangulation method results.

Uses proper statistical methods:
- Poisson distribution for variance/standard error
- Gamma distribution for confidence intervals (Fay-Feuer method)
  This produces asymmetric CIs appropriate for small death counts

Input: Subtracted death counts from Data/Original/CDC Triangulation/Subtracted/
Output: AAMR with CI to Data/Original/CDC Triangulation/AAMR/
"""

from __future__ import annotations

import pandas as pd
import yaml
from pathlib import Path
from scipy.stats import gamma
from typing import Any


# ── Run control ──────────────────────────────────────────────────────────────
# Set a disease group to 0 to skip its AAMR calculation. Any other value (or absent key) runs it.
# All groups run by default. Uncomment lines below to skip:
SKIP_GROUPS: dict = {
    # "cancer": 0,
    # "ndd": 0,
}

# Maps icd_code → disease group key; populated from config in main().
ICD_TO_GROUP: dict = {}


def _build_icd_to_group(config: dict) -> dict:
    """Build {icd_code: group_key} from diseases config for group-based skip checks."""
    mapping: dict = {}
    for group_key, group_cfg in config.get("diseases", {}).items():
        overall = group_cfg.get("overall", {})
        if overall.get("icd_code"):
            mapping[overall["icd_code"]] = group_key
        for subtype in group_cfg.get("subtypes", {}).values():
            code = subtype.get("icd_code")
            if code:
                mapping[code] = group_key
    return mapping


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

# Age groups in order
AGE_GROUPS = [
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


def parse_filename(filename: str) -> tuple[str | None, str | None]:
    """
    Parse filename to extract year range and ICD code.

    Example: '2016-2020_G30.csv' -> ('2016-2020', 'G30')
    """
    name = filename.replace(".csv", "")
    parts = name.split("_", 1)

    if len(parts) == 2:
        return parts[0], parts[1]
    return None, None


def calculate_age_specific_rate(deaths: float, population: float) -> float:
    """
    Calculate age-specific death rate per 100,000 population.

    mi = (Di / Pi) × 100,000
    """
    if pd.isna(deaths) or pd.isna(population) or population == 0:
        return 0.0

    return (float(deaths) / float(population)) * 100000


def calculate_aamr_point(df: pd.DataFrame) -> tuple[float, dict[str, Any]]:
    """
    Calculate AAMR point estimate.

    AAMR = Σ(mi × wi)

    Returns: (aamr, stats_dict)
    """
    aamr = 0.0
    total_deaths = 0
    total_population = 0
    age_rates = {}

    for age_group in AGE_GROUPS:
        weight = AGE_WEIGHTS[age_group]

        # Get data for this age group
        age_data = df[df["Ten-Year Age Groups"] == age_group]

        if len(age_data) == 0:
            # Missing age group
            age_rates[age_group] = 0.0
            continue

        deaths = age_data["Deaths"].iloc[0]
        population = age_data["Population"].iloc[0]

        # Handle missing deaths (keep as 0 for AAMR calculation)
        if pd.isna(deaths):
            deaths = 0.0

        if pd.isna(population):
            population = 0.0

        # Calculate age-specific rate
        mi = calculate_age_specific_rate(deaths, population)
        age_rates[age_group] = mi

        # Add to AAMR
        aamr += mi * weight

        # Track totals
        total_deaths += deaths
        total_population += population

    stats = {
        "aamr": aamr,
        "total_deaths": total_deaths,
        "total_population": total_population,
        "age_rates": age_rates,
    }

    return aamr, stats


def calculate_aamr_standard_error(df: pd.DataFrame) -> float:
    """
    Calculate standard error of AAMR using Poisson variance.

    SE(AAMR) = sqrt(Σ(w_i² × D_i / P_i²)) × 100,000

    Returns: standard_error
    """
    variance = 0.0

    for age_group in AGE_GROUPS:
        weight = AGE_WEIGHTS[age_group]

        age_data = df[df["Ten-Year Age Groups"] == age_group]
        if len(age_data) == 0:
            continue

        deaths = age_data["Deaths"].iloc[0]
        population = age_data["Population"].iloc[0]

        if pd.isna(deaths) or pd.isna(population) or population == 0:
            continue

        # Variance component: w_i² × D_i / P_i²
        variance += (weight**2) * (deaths / (population**2))

    # Standard error
    se = (variance**0.5) * 100000

    return se


def calculate_aamr_ci(
    df: pd.DataFrame, aamr: float, stats: dict[str, Any]
) -> tuple[float, float]:
    """
    Calculate 95% confidence interval for AAMR using Fay-Feuer method.

    This method uses Gamma distribution for small death counts (Poisson-based).
    Steps:
    1. Calculate effective population: P* = D_total × 100,000 / AAMR
    2. Get CI for total death count using Gamma distribution
    3. Convert death count CI back to AAMR CI

    This produces asymmetric CIs appropriate for small counts.

    Returns: (aamr_lower, aamr_upper)
    """
    total_deaths = stats["total_deaths"]

    # Handle edge cases
    if total_deaths == 0:
        return 0.0, 0.0

    if aamr == 0:
        return 0.0, 0.0

    # Step B: Calculate effective population (P*)
    # P* = D_total × 100,000 / AAMR
    # This treats the age-adjusted rate as if it were a crude rate
    p_star = (total_deaths * 100000) / aamr

    # Step C: Calculate CI for death count using Gamma distribution
    # For Poisson parameter λ, the 95% CI is:
    # Lower: gamma.ppf(0.025, total_deaths) / 1
    # Upper: gamma.ppf(0.975, total_deaths + 1) / 1

    if total_deaths == 0:
        d_lower = 0.0
        d_upper = gamma.ppf(0.975, 1)  # ~3.69
    else:
        # Lower bound: 2.5th percentile of Gamma(D, 1)
        d_lower = gamma.ppf(0.025, total_deaths)
        # Upper bound: 97.5th percentile of Gamma(D+1, 1)
        d_upper = gamma.ppf(0.975, total_deaths + 1)

    # Step D: Convert death count CI to AAMR CI
    aamr_lower = (d_lower / p_star) * 100000
    aamr_upper = (d_upper / p_star) * 100000

    return aamr_lower, aamr_upper


def process_subtracted_file(file_path: Path) -> pd.DataFrame:
    """
    Process one subtracted file and calculate AAMR for each county.

    Returns DataFrame with: COUNTY_FIPS, Deaths, Population, AAMR, AAMR_Lower, AAMR_Upper
    """
    print(f"  Processing: {file_path.name}")

    # Load data
    df = pd.read_csv(file_path, dtype={"County Code": str})

    # Filter out rows with missing County Code (NaN rows at end of file)
    df = df[df["County Code"].notna()].copy()

    # Filter out rows with 'Missing' population values
    df = df[df["Population"] != "Missing"].copy()

    # Convert numeric columns to proper types
    df["Deaths"] = pd.to_numeric(df["Deaths"], errors="coerce")
    df["Population"] = pd.to_numeric(df["Population"], errors="coerce")

    # Group by county
    counties = df["County Code"].unique()

    results = []

    for county_code in counties:
        county_data = df[df["County Code"] == county_code].copy()

        # Calculate AAMR
        aamr, stats = calculate_aamr_point(county_data)

        # Calculate confidence interval using Fay-Feuer method
        aamr_lower, aamr_upper = calculate_aamr_ci(county_data, aamr, stats)

        # Calculate standard error (for reference)
        se = calculate_aamr_standard_error(county_data)

        results.append(
            {
                "COUNTY_FIPS": county_code,
                "Deaths": int(stats["total_deaths"]),
                "Population": int(stats["total_population"]),
                "AAMR": round(aamr, 4),
                "AAMR_SE": round(se, 4),
                "AAMR_Lower": round(aamr_lower, 4),
                "AAMR_Upper": round(aamr_upper, 4),
            }
        )

    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values("COUNTY_FIPS")

    print(f"    ✓ Processed {len(result_df)} counties")
    print(
        f"      AAMR range: {result_df['AAMR'].min():.2f} - {result_df['AAMR'].max():.2f}"
    )

    return result_df


# ---------------------------------------------------------------------------
# Synthesis of sum_of diseases (config-driven)
# ---------------------------------------------------------------------------


def _compute_aamr_from_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-county AAMR from a cleaned DataFrame.
    Expects columns: 'County Code', 'Ten-Year Age Groups', 'Deaths', 'Population' (numeric).
    """
    counties = df["County Code"].unique()
    results = []
    for county_code in counties:
        county_data = df[df["County Code"] == county_code].copy()
        aamr, stats = calculate_aamr_point(county_data)
        aamr_lower, aamr_upper = calculate_aamr_ci(county_data, aamr, stats)
        se = calculate_aamr_standard_error(county_data)
        results.append(
            {
                "COUNTY_FIPS": county_code,
                "Deaths": int(stats["total_deaths"]),
                "Population": int(stats["total_population"]),
                "AAMR": round(aamr, 4),
                "AAMR_SE": round(se, 4),
                "AAMR_Lower": round(aamr_lower, 4),
                "AAMR_Upper": round(aamr_upper, 4),
            }
        )
    return pd.DataFrame(results).sort_values("COUNTY_FIPS")


def main():
    """Main execution function"""
    print("=" * 70)
    print("CDC Triangulation AAMR Calculator")
    print("=" * 70)

    # Load configuration
    project_root, config = load_config()

    global ICD_TO_GROUP
    ICD_TO_GROUP = _build_icd_to_group(config)

    # Set up paths
    input_dir = project_root / "Data/Original/CDC Triangulation/Subtracted"
    output_dir = project_root / "Data/Original/CDC Triangulation/AAMR"

    print(f"\nInput directory:  {input_dir}")
    print(f"Output directory: {output_dir}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get all subtracted files
    subtracted_files = list(input_dir.glob("*.csv"))

    if len(subtracted_files) == 0:
        print("\n⚠ No subtracted files found!")
        print("Please run EQI_Triangulation.py first to generate subtracted data.")
        return

    print(f"\nFound {len(subtracted_files)} subtracted files")
    print(f"SKIP_GROUPS = {SKIP_GROUPS}")

    # Process each file
    print("\n" + "=" * 70)
    print("Calculating AAMR with Confidence Intervals")
    print("=" * 70 + "\n")

    processed_count = 0
    skipped_count = 0

    for file_path in sorted(subtracted_files):
        year_range, icd_code = parse_filename(file_path.name)

        if not year_range or not icd_code:
            print(f"  ⚠ Could not parse filename: {file_path.name}, skipping")
            continue

        group_key = ICD_TO_GROUP.get(icd_code, "")
        if SKIP_GROUPS.get(group_key) == 0:
            skipped_count += 1
            continue

        # Process file
        result_df = process_subtracted_file(file_path)

        # Save output with same filename
        output_filename = file_path.name
        output_path = output_dir / output_filename
        result_df.to_csv(output_path, index=False)

        print(f"    ✓ Saved to: {output_filename}\n")

        processed_count += 1

    # Final summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Processed:     {processed_count} files")
    print(f"Skipped:       {skipped_count} files (SKIP_GROUPS={SKIP_GROUPS})")
    print(f"Output dir:    {output_dir}")
    print("\n✓ AAMR calculation completed successfully!")

    # Print example output
    if processed_count > 0:
        example_file = sorted(output_dir.glob("*.csv"))[0]
        print(f"\nExample output from {example_file.name}:")
        example_df = pd.read_csv(example_file)
        print(example_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
