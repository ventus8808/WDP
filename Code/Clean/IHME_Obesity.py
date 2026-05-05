#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 肥胖率数据清洗脚本
数据源：IHME USA BMI County Race Ethnicity 2000-2019 (OBESITY_BOTH)
输出：IHME_Obesity.csv — 县级肥胖率 0005/0610 均值
"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/IHME Obesity/OBESITY_BOTH"
DATA_OUT = BASE / "Data/Processed/IHME"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

PERIODS = {
    "0005": set(range(2000, 2006)),
    "0610": set(range(2006, 2011)),
}
ALL_YEARS = set().union(*PERIODS.values())


def pad_fips(series: pd.Series) -> pd.Series:
    return series.dropna().astype(float).astype(int).astype(str).str.zfill(5)


def main():
    files = sorted(DATA_IN.glob("*.CSV"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {DATA_IN}")

    chunks = []
    for f in files:
        df = pd.read_csv(f, dtype={"fips": str})
        df = df[df["year"].isin(ALL_YEARS)]
        if df.empty:
            continue
        chunks.append(df)

    raw = pd.concat(chunks, ignore_index=True)
    raw["fips"] = pad_fips(raw["fips"])

    base = raw[
        (raw["race_name"] == "Total")
        & (raw["age_name"] == "20 plus")   # "All Ages" does not exist in this file
        & (raw["fips"].str.len() == 5)
    ].copy()

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            base[base["year"].isin(years)]
            .groupby("fips")["val"]
            .mean()
            .mul(100)
            .round(2)
            .reset_index()
            .rename(columns={"fips": "COUNTY_FIPS", "val": f"Obesity_rate_{suffix}"})
        )
        col = f"Obesity_rate_{suffix}"
        ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    county = parts[0].merge(parts[1], on="COUNTY_FIPS", how="outer").sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "IHME_Obesity.csv"
    county.to_csv(out_path, index=False)
    print(f"Saved {len(county):,} counties → {out_path}")
    print(county.describe())


if __name__ == "__main__":
    main()
