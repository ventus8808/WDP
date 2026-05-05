#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
USDA ERS County Typology cleaner
Source: 2004CountyTypologyCode.xls (econdep 1-6)
Output: Economic_type.csv — COUNTY_FIPS + Economic_type_0005 + Economic_type_0610
Note: same 2004 classification used for both periods (structural, not annual).
"""

import yaml
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

with open(PROJECT_ROOT / "config.yaml", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

ers = cfg["data_sources"]["socioeconomic"]["usda_ers"]
SRC = PROJECT_ROOT / ers["original"] / ers["typology_source"]
OUT = PROJECT_ROOT / ers["typology_output"]

df = pd.read_excel(SRC, dtype={"FIPSTXT": str})
df["COUNTY_FIPS"] = df["FIPSTXT"].str.zfill(5)
df["econdep"] = pd.to_numeric(df["econdep"], errors="coerce").astype("Int64")
df = df[df["econdep"].notna() & df["COUNTY_FIPS"].str.len().eq(5)].copy()
df = df.drop_duplicates("COUNTY_FIPS", keep="first")

result = (
    df[["COUNTY_FIPS", "econdep"]]
    .rename(columns={"econdep": "Economic_type_0005"})
    .assign(Economic_type_0610=lambda x: x["Economic_type_0005"])
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

OUT_PATH = PROJECT_ROOT / OUT if not str(OUT).startswith("/") else Path(str(OUT))
OUT_PATH = PROJECT_ROOT / ers["typology_output"]
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT_PATH, index=False)
print(f"Saved {len(result):,} counties → {OUT_PATH}")
print(result["Economic_type_0005"].value_counts().sort_index())
