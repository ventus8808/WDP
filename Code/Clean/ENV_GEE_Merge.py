#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merges and processes GEE data (NLCD Land Cover and JRC Surface Water).

This script reads the raw GEE CSV files, merges them, creates a complete
county-year panel, and applies linear interpolation to fill missing yearly data for
NLCD land cover categories.
"""

import sys
from pathlib import Path
from typing import List, Optional

import pandas as pd
import yaml


def load_paths() -> tuple[Path, Path]:
    """Loads input and output paths from config.yaml."""
    project_root = Path(__file__).resolve().parents[2]
    config_path = project_root / "config.yaml"
    if not config_path.exists():
        sys.exit(f"ERROR: Configuration file not found: {config_path}")

    with config_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    try:
        gee_config = cfg["data_sources"]["socioeconomic"]["gee"]
        input_rel = gee_config["original"]
    except (KeyError, TypeError):
        sys.exit("ERROR: config.yaml is missing the required path for data_sources.socioeconomic.gee")

    input_dir = (project_root / input_rel).resolve()
    output_path = (project_root / "Data/Processed/Environmental/NLCD_JRC.csv").resolve()

    if not input_dir.exists():
        sys.exit(f"ERROR: Input directory not found: {input_dir}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    return input_dir, output_path


def find_latest_csv(directory: Path, prefix: str) -> Optional[Path]:
    """Finds the most recently modified CSV file with a given prefix."""
    candidates = sorted(directory.glob(f"{prefix}*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def read_and_concat_csvs(paths: List[Path], required_cols: set) -> pd.DataFrame:
    """Reads multiple CSVs, validates columns, and concatenates them."""
    frames = []
    for path in paths:
        df = pd.read_csv(path, dtype={"GEOID": str})
        if not required_cols.issubset(df.columns):
            raise ValueError(f"File {path.name} is missing required columns: {required_cols - set(df.columns)}")
        df["Year"] = pd.to_numeric(df["Year"], errors="coerce").astype("Int64")
        frames.append(df[list(required_cols)])

    if not frames:
        raise FileNotFoundError(f"No valid data files found.")

    combined_df = pd.concat(frames, ignore_index=True)
    return combined_df.drop_duplicates(subset=["GEOID", "Year"], keep="last").reset_index(drop=True)


def merge_gee_data(county_path: Path, jrc_paths: List[Path], nlcd_paths: List[Path]) -> pd.DataFrame:
    """Merges county base, JRC water, and NLCD land cover data."""
    county_df = pd.read_csv(county_path, dtype={"GEOID": str})[["GEOID", "total_area_km2"]]

    jrc_cols = {"GEOID", "Year", "jrc_permanent_water_km2", "jrc_seasonal_water_km2"}
    jrc_df = read_and_concat_csvs(jrc_paths, jrc_cols)

    nlcd_cols = {
        "GEOID", "Year", "nlcd_forest_km2", "nlcd_water_km2", "nlcd_urban_km2",
        "nlcd_agriculture_km2", "nlcd_cropland_km2", "nlcd_pasture_km2", "nlcd_wetland_km2",
        "nlcd_wetland_woody_km2", "nlcd_wetland_herb_km2", "nlcd_shrub_km2",
        "nlcd_grassland_km2", "nlcd_barren_km2",
    }
    nlcd_df = read_and_concat_csvs(nlcd_paths, nlcd_cols)

    merged = pd.merge(jrc_df, nlcd_df, on=["GEOID", "Year"], how="outer")
    merged = pd.merge(merged, county_df, on="GEOID", how="left")
    return merged.sort_values(by=["GEOID", "Year"]).reset_index(drop=True)


def create_full_panel(df: pd.DataFrame, year_min: int, year_max: int) -> pd.DataFrame:
    """Creates a complete county-year panel for the specified date range."""
    geoids = df["GEOID"].unique()
    years = range(year_min, year_max + 1)
    grid = pd.MultiIndex.from_product([geoids, years], names=["GEOID", "Year"]).to_frame(index=False)

    # Preserve total_area_km2 by merging it separately
    area_df = df[["GEOID", "total_area_km2"]].drop_duplicates(subset=["GEOID"])

    panel = pd.merge(grid, df.drop(columns=["total_area_km2"]), on=["GEOID", "Year"], how="left")
    panel = pd.merge(panel, area_df, on="GEOID", how="left")
    return panel.sort_values(["GEOID", "Year"]).reset_index(drop=True)


def interpolate_panel(df: pd.DataFrame) -> pd.DataFrame:
    """Applies linear interpolation to fill missing data within each county group."""
    # All columns except identifiers are interpolated
    value_cols = [col for col in df.columns if col not in ["GEOID", "Year"]]

    def interpolate_group(g):
        g[value_cols] = g[value_cols].interpolate(method="linear", limit_direction="both")
        g[value_cols] = g[value_cols].ffill().bfill() # Fill any remaining gaps at boundaries
        return g

    interpolated = df.groupby("GEOID", group_keys=False).apply(interpolate_group)
    interpolated[value_cols] = interpolated[value_cols].round(4).clip(lower=0)
    return interpolated


def main():
    """Main execution function."""
    input_dir, output_path = load_paths()

    try:
        county_path = find_latest_csv(input_dir, "county_base")
        jrc_paths = sorted(input_dir.glob("jrc_water_*.csv"))
        nlcd_paths = sorted(input_dir.glob("nlcd_landuse_*.csv"))

        if not county_path or not jrc_paths or not nlcd_paths:
            raise FileNotFoundError("Required input files (county_base, jrc_water, nlcd_landuse) not found.")

        print(f"Processing {len(jrc_paths)} JRC files and {len(nlcd_paths)} NLCD files...")

        # Merge raw data
        merged_data = merge_gee_data(county_path, jrc_paths, nlcd_paths)

        # Create a complete panel from 1999-2020
        panel_data = create_full_panel(merged_data, 1999, 2020)

        # Interpolate missing values
        final_data = interpolate_panel(panel_data)

        # Rename GEOID to COUNTY_FIPS for consistency with other datasets
        final_data = final_data.rename(columns={"GEOID": "COUNTY_FIPS"})

        # Save the final interpolated data
        final_data.to_csv(output_path, index=False)
        print(f"\nSuccessfully merged and interpolated data.")
        print(f"Output saved to: {output_path}")
        print(f"Final data shape: {final_data.shape}")

    except (FileNotFoundError, ValueError) as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
