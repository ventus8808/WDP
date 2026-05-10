#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climate zone extractor
Source:
  - koeppengeigerUScounty.txt  → Köppen-Geiger code (dominant class per county)
Output: Climate_Zone.csv — COUNTY_FIPS, Climate_Zone
"""

import yaml
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

with open(PROJECT_ROOT / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

cz = cfg["data_sources"]["climate_zone"]
DATA_IN = PROJECT_ROOT / cz["original"]
OUT     = PROJECT_ROOT / cz["output_file"]

# Köppen-Geiger: pick dominant class (highest PROP) per county
koppen = pd.read_csv(DATA_IN / cz["koppen_file"], sep="\t", dtype={"FIPS": str})
print(f"Loaded Köppen data: {len(koppen):,} records")

koppen = koppen[["FIPS", "CLS", "PROP"]].rename(columns={"FIPS": "COUNTY_FIPS", "CLS": "koppen_code"})
koppen["COUNTY_FIPS"] = koppen["COUNTY_FIPS"].str.zfill(5)
koppen["Climate_Zone"] = koppen["koppen_code"].str[0]

print("\nBefore deduplication (by zone):")
print(koppen["Climate_Zone"].value_counts().sort_index())

koppen = (
    koppen.sort_values("PROP", ascending=False)
    .drop_duplicates("COUNTY_FIPS", keep="first")
)

print("\nAfter deduplication (before merging A→B, E→D):")
print(koppen["Climate_Zone"].value_counts().sort_index())

# Merge A into B, E into D
koppen["Climate_Zone"] = koppen["Climate_Zone"].replace({"A": "B", "E": "D"})

print("\nAfter merging A→B and E→D:")
print(koppen["Climate_Zone"].value_counts().sort_index())
print(f"\nTotal unique counties: {len(koppen):,}")

result = (
    koppen[["COUNTY_FIPS", "Climate_Zone"]]
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

OUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT, index=False)
print(f"Saved {len(result):,} counties → {OUT}")
print(result.describe(include="all"))
