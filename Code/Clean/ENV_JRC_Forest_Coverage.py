#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from pathlib import Path

BASE    = Path(__file__).resolve().parents[2]
INPUT   = BASE / "Data/Processed/Environmental/NLCD_JRC.csv"
OUTPUT  = BASE / "Data/Processed/Covariate/Forest_Coverage.csv"
EQI_DIR = BASE / "Data/Processed/EQI"

PERIODS = {
    "0005": (2000, 2005),
    "0610": (2006, 2010),
}

df = pd.read_csv(INPUT, dtype={"COUNTY_FIPS": str})
df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
df = df[df["nlcd_forest_km2"].notna() & df["total_area_km2"].notna()]
df = df[df["total_area_km2"] > 0]
df["forest_pct"] = df["nlcd_forest_km2"] / df["total_area_km2"] * 100

parts = []
for suffix, (yr_start, yr_end) in PERIODS.items():
    part = (
        df[df["Year"].between(yr_start, yr_end)]
        .groupby("COUNTY_FIPS")["forest_pct"]
        .mean()
        .round(1)
        .reset_index()
        .rename(columns={"forest_pct": f"Forest_Coverage_{suffix}"})
    )
    col = f"Forest_Coverage_{suffix}"
    ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
    ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
    part = ref.merge(part, on="COUNTY_FIPS", how="left")
    n_missing = part[col].isna().sum()
    part[col] = part[col].fillna(part[col].mean()).round(1)
    parts.append(part)
    print(f"{suffix}: {len(part)} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

result = parts[0]
for part in parts[1:]:
    result = result.merge(part, on="COUNTY_FIPS", how="outer")

result = result.sort_values("COUNTY_FIPS").reset_index(drop=True)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, index=False)
print(f"\nDone: {len(result)} counties → {OUTPUT}")
print(result.describe())
