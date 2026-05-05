#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDC WONDER Census Region extractor
Source: Location_Region_Division_State.csv (county-level static mapping)
Output: Census_Region.csv — COUNTY_FIPS + Census_Region (1-4)
"""

import yaml
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

with open(PROJECT_ROOT / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

cdc = cfg["data_sources"]["cdc_wonder"]
SRC = PROJECT_ROOT / cdc["census_region_original"] / cdc["census_region_source_file"]
OUT = PROJECT_ROOT / cdc["census_region_output"]

df = pd.read_csv(SRC, dtype=str)
df = df[df["County Code"].notna() & (df["County Code"].str.strip() != "")]
df["COUNTY_FIPS"] = df["County Code"].str.strip().str.zfill(5)
df["Census_Region"] = df["Census Region Code"].str.extract(r"CENS-R(\d)").astype("Int64")

result = (
    df[["COUNTY_FIPS", "Census_Region"]]
    .drop_duplicates("COUNTY_FIPS")
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

OUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT, index=False)
print(f"Saved {len(result):,} counties → {OUT}")
print(result["Census_Region"].value_counts().sort_index())
