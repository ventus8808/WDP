#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climate zone extractor
Sources:
  - koeppengeigerUScounty.txt  → Köppen-Geiger code (dominant class per county)
  - climate_zones.csv          → DOE/IECC climate zone (1-8)
Output: Climate_Zone.csv — COUNTY_FIPS, koppen_code, koppen_major, doe_zone
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
koppen = koppen[["FIPS", "CLS", "PROP"]].rename(columns={"FIPS": "COUNTY_FIPS", "CLS": "koppen_code"})
koppen["COUNTY_FIPS"] = koppen["COUNTY_FIPS"].str.zfill(5)
koppen = (
    koppen.sort_values("PROP", ascending=False)
    .drop_duplicates("COUNTY_FIPS", keep="first")
    [["COUNTY_FIPS", "koppen_code"]]
)
koppen["koppen_major"] = koppen["koppen_code"].str[0]

# DOE/IECC climate zone
doe = pd.read_csv(DATA_IN / cz["doe_file"], dtype={"State FIPS": str, "County FIPS": str})
doe["COUNTY_FIPS"] = doe["State FIPS"].str.zfill(2) + doe["County FIPS"].str.zfill(3)
doe = doe[["COUNTY_FIPS", "IECC Climate Zone"]].rename(columns={"IECC Climate Zone": "doe_zone"})
doe["doe_zone"] = pd.to_numeric(doe["doe_zone"], errors="coerce").astype("Int64")

result = (
    koppen.merge(doe, on="COUNTY_FIPS", how="outer")
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

OUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT, index=False)
print(f"Saved {len(result):,} counties → {OUT}")
print(result.describe(include="all"))
