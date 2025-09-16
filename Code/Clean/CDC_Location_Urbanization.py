#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDC Location & Urbanization cleaner (two-table output)

Paths are loaded from config.yaml (data_sources.cdc_wonder.*). If config is
missing or unreadable, fall back to relative paths via get_data_dir().

- Input (from config): data_sources.cdc_wonder.location_urbanization_original
- Outputs (from config):
  - data_sources.cdc_wonder.location_output_file  -> Location.csv
  - data_sources.cdc_wonder.urbanization_output_file -> Urbanization.csv

Tables produced:
  - Location.csv        (county static: COUNTY_FIPS, County, HHS_Region, Census_Region, Census_Division)
  - Urbanization.csv    (county×year: COUNTY_FIPS, Year, County, Urbanization_Code, Urbanization_Type)
"""

import sys
from pathlib import Path
import pandas as pd

# Optional: read paths from YAML
try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover
    yaml = None  # fallback below

# fallback helper from previous config module
sys.path.append(str(Path(__file__).parent.parent))
from config import get_data_dir  # noqa: E402

# ---------- resolve paths from config.yaml (preferred) ----------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

SRC_DIR: Path
LOCATION_OUT: Path
URBAN_OUT: Path

_loaded_from_yaml = False
if yaml is not None and CONFIG_PATH.exists():
    try:
        with CONFIG_PATH.open("r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        ds = (cfg or {}).get("data_sources", {})
        cdc = ds.get("cdc_wonder", {})
        src_rel = cdc.get("location_urbanization_original")
        loc_rel = cdc.get("location_output_file")
        urb_rel = cdc.get("urbanization_output_file")
        if src_rel and loc_rel and urb_rel:
            SRC_DIR = PROJECT_ROOT / src_rel
            LOCATION_OUT = PROJECT_ROOT / loc_rel
            URBAN_OUT = PROJECT_ROOT / urb_rel
            _loaded_from_yaml = True
    except Exception:
        _loaded_from_yaml = False

# ---------- fallback to previous relative resolution ----------
if not _loaded_from_yaml:
    # keep the literal space to match on-disk path
    SRC_DIR = get_data_dir("original") / " CDC WONDER" / "Location and Urbanization"
    LOCATION_OUT = get_data_dir("processed") / "CDC" / "Location.csv"
    URBAN_OUT = get_data_dir("processed") / "CDC" / "Urbanization.csv"

URBAN_PREFIX = "Location_County_Urbanization"
HHS_FILE = "Location_HHS_State.csv"
CENSUS_FILE = "Location_Region_Division_State.csv"


def _standardize_fips(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace('.0', '', regex=False).str.zfill(5)


def load_urbanization_panel() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for fp in sorted(SRC_DIR.glob(f"{URBAN_PREFIX}*.csv")):
        try:
            df = pd.read_csv(fp)
        except Exception as exc:
            print(f"[URB] read failed: {fp.name}: {exc}")
            continue
        needed = {'Year','County','County Code','2013 Urbanization','2013 Urbanization Code'}
        if not needed.issubset(df.columns):
            print(f"[URB] missing columns, skip: {fp.name}")
            continue
        sub = df[['Year','County','County Code','2013 Urbanization','2013 Urbanization Code']].copy()
        sub['COUNTY_FIPS'] = _standardize_fips(sub['County Code'])
        sub = sub.rename(columns={
            '2013 Urbanization': 'Urbanization_Type',
            '2013 Urbanization Code': 'Urbanization_Code'
        })
        sub['Year'] = pd.to_numeric(sub['Year'], errors='coerce').astype('Int64')
        sub = sub[['COUNTY_FIPS','Year','County','Urbanization_Code','Urbanization_Type']]
        frames.append(sub)
    if not frames:
        return pd.DataFrame(columns=['COUNTY_FIPS','Year','County','Urbanization_Code','Urbanization_Type'])
    urb = pd.concat(frames, ignore_index=True)
    urb = urb.dropna(subset=['Year'])
    urb = (urb.sort_values(['COUNTY_FIPS','Year'])
              .drop_duplicates(subset=['COUNTY_FIPS','Year'], keep='first')
              .reset_index(drop=True))
    return urb


def load_location_static() -> pd.DataFrame:
    # HHS
    hhs = pd.DataFrame(columns=['COUNTY_FIPS','County','HHS_Region'])
    fp = SRC_DIR / HHS_FILE
    if fp.exists():
        try:
            df = pd.read_csv(fp)
            need = {'County','County Code','HHS Region'}
            if need.issubset(df.columns):
                sub = df[['County','County Code','HHS Region']].copy()
                sub['COUNTY_FIPS'] = _standardize_fips(sub['County Code'])
                sub = sub.rename(columns={'HHS Region': 'HHS_Region'})
                hhs = sub[['COUNTY_FIPS','County','HHS_Region']]
        except Exception as exc:
            print(f"[LOC] HHS read failed: {exc}")

    # Census Region/Division
    cen = pd.DataFrame(columns=['COUNTY_FIPS','County','Census_Region','Census_Division'])
    fp = SRC_DIR / CENSUS_FILE
    if fp.exists():
        try:
            df = pd.read_csv(fp)
            need = {'County','County Code','Census Region','Census Division'}
            if need.issubset(df.columns):
                sub = df[['County','County Code','Census Region','Census Division']].copy()
                sub['COUNTY_FIPS'] = _standardize_fips(sub['County Code'])
                sub = sub.rename(columns={
                    'Census Region': 'Census_Region',
                    'Census Division': 'Census_Division',
                })
                cen = sub[['COUNTY_FIPS','County','Census_Region','Census_Division']]
        except Exception as exc:
            print(f"[LOC] Census read failed: {exc}")

    if hhs.empty and cen.empty:
        return pd.DataFrame(columns=['COUNTY_FIPS','County','HHS_Region','Census_Region','Census_Division'])

    merged = pd.merge(hhs, cen, on='COUNTY_FIPS', how='outer', suffixes=("_hhs","_cen"))
    merged['County'] = merged.get('County_hhs').fillna(merged.get('County_cen'))
    merged = merged.drop(columns=[c for c in merged.columns if c.startswith('County_') and c != 'County'])
    keep = ['COUNTY_FIPS','County','HHS_Region','Census_Region','Census_Division']
    merged = merged[[c for c in keep if c in merged.columns]]
    merged = (merged.sort_values(['COUNTY_FIPS'])
                    .drop_duplicates(subset=['COUNTY_FIPS'], keep='first')
                    .reset_index(drop=True))
    return merged


def main() -> None:
    print(f"Source: {SRC_DIR}")
    # Ensure output directories
    LOCATION_OUT.parent.mkdir(parents=True, exist_ok=True)
    URBAN_OUT.parent.mkdir(parents=True, exist_ok=True)

    loc = load_location_static()
    urb = load_urbanization_panel()

    if not loc.empty:
        loc.to_csv(LOCATION_OUT, index=False)
        print(f"[SAVE] {LOCATION_OUT} ({len(loc)})")
    else:
        print("[WARN] no location data")

    if not urb.empty:
        urb.to_csv(URBAN_OUT, index=False)
        print(f"[SAVE] {URBAN_OUT} ({len(urb)})")
    else:
        print("[WARN] no urbanization data")


if __name__ == '__main__':
    main()
