#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 体力活动数据清洗脚本
数据源：IHME USA Obesity Physical Activity 2001-2011
注意：无 Both 合并行，取 Male + Female 均值作为两性合并估计。
输出：IHME_Physical_Activities.csv — 县级体力活动率 0005/0610 均值
"""

import re
import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/IHME Physical Activities"
DATA_OUT = BASE / "Data/Processed/IHME"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

PA_FILE = DATA_IN / "IHME_USA_OBESITY_PHYSICAL_ACTIVITY_2001_2011.csv"

PERIODS = {
    "0005": set(range(2001, 2006)),  # data starts 2001
    "0610": set(range(2006, 2011)),
}
ALL_YEARS = set().union(*PERIODS.values())


def main():
    df = pd.read_csv(PA_FILE, dtype={"fips": str})

    df = df.dropna(subset=["fips"])
    df["fips"] = df["fips"].astype(float).astype(int).astype(str).str.zfill(5)
    df = df[(df["fips"].str.len() == 5) & (df["Outcome"] == "Any PA")].copy()

    year_cols = {
        col: int(m.group(1))
        for col in df.columns
        if (m := re.match(r"^Prevalence (\d{4}) \(%\)$", str(col)))
        and int(m.group(1)) in ALL_YEARS
    }
    if not year_cols:
        raise ValueError(f"No prevalence columns found for years {ALL_YEARS}")

    long = df[["fips", "Sex"] + list(year_cols)].melt(
        id_vars=["fips", "Sex"], var_name="year_col", value_name="pa_rate"
    )
    long["year"] = long["year_col"].map(year_cols)

    # average Male + Female per county-year first
    by_year = (
        long.dropna(subset=["pa_rate"])
        .groupby(["fips", "year"])["pa_rate"]
        .mean()
        .reset_index()
    )

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            by_year[by_year["year"].isin(years)]
            .groupby("fips")["pa_rate"]
            .mean()
            .round(2)
            .reset_index()
            .rename(columns={"fips": "COUNTY_FIPS", "pa_rate": f"Physical_Activities_rate_{suffix}"})
        )
        col = f"Physical_Activities_rate_{suffix}"
        ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    county = parts[0].merge(parts[1], on="COUNTY_FIPS", how="outer").sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "IHME_Physical_Activities.csv"
    county.to_csv(out_path, index=False)
    print(f"Saved {len(county):,} counties → {out_path}")
    print(county.describe())


if __name__ == "__main__":
    main()
