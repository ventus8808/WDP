#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 吸烟率数据清洗脚本
数据源：IHME US County Total and Daily Smoking Prevalence 1996-2012
FIPS 通过 IHME Diabetes 文件中的 (State, Location)→FIPS 映射表获取（直接匹配）。
输出：IHME_Smoking.csv — 县级吸烟率 0005/0610 均值（%，保留2位小数）
"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/IHME Smoking"
DIAB_IN  = BASE / "Data/Original/IHME Diabetes Prevalence"
DATA_OUT = BASE / "Data/Processed/IHME"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

SMOKE_FILE = DATA_IN / "IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012.csv"
DIAB_FILE  = DIAB_IN / "IHME_USA_COUNTY_DIABETES_PREVALENCE_1999_2012_NATIONAL_Y2016M08D23.XLSX"

PERIODS = {
    "0005": set(range(2000, 2006)),
    "0610": set(range(2006, 2011)),
}
ALL_YEARS = set().union(*PERIODS.values())

STATE_FIPS = {
    "01": "Alabama", "02": "Alaska", "04": "Arizona", "05": "Arkansas", "06": "California",
    "08": "Colorado", "09": "Connecticut", "10": "Delaware", "11": "District of Columbia",
    "12": "Florida", "13": "Georgia", "15": "Hawaii", "16": "Idaho", "17": "Illinois",
    "18": "Indiana", "19": "Iowa", "20": "Kansas", "21": "Kentucky", "22": "Louisiana",
    "23": "Maine", "24": "Maryland", "25": "Massachusetts", "26": "Michigan", "27": "Minnesota",
    "28": "Mississippi", "29": "Missouri", "30": "Montana", "31": "Nebraska", "32": "Nevada",
    "33": "New Hampshire", "34": "New Jersey", "35": "New Mexico", "36": "New York",
    "37": "North Carolina", "38": "North Dakota", "39": "Ohio", "40": "Oklahoma", "41": "Oregon",
    "42": "Pennsylvania", "44": "Rhode Island", "45": "South Carolina", "46": "South Dakota",
    "47": "Tennessee", "48": "Texas", "49": "Utah", "50": "Vermont", "51": "Virginia",
    "53": "Washington", "54": "West Virginia", "55": "Wisconsin", "56": "Wyoming",
}


def build_fips_lookup() -> pd.DataFrame:
    df = pd.read_excel(DIAB_FILE, sheet_name="Diagnosed", header=1)
    df = df[["Location", "FIPS"]].dropna(subset=["FIPS"])
    df["FIPS"] = df["FIPS"].astype(float).astype(int).astype(str).str.zfill(5)
    df = df[df["FIPS"].str.len() == 5].copy()
    df["state"] = df["FIPS"].str[:2].map(STATE_FIPS)
    return df[["state", "Location", "FIPS"]].rename(columns={"Location": "county", "FIPS": "COUNTY_FIPS"})


def main():
    fips_ref = build_fips_lookup()

    df = pd.read_csv(SMOKE_FILE)

    base = (
        df[
            df["county"].notna()
            & (df["county"] != "")
            & (df["sex"] == "Both")
            & (df["year"].isin(ALL_YEARS))
        ]
        .merge(fips_ref, on=["state", "county"], how="left")
    )

    unmatched = base["COUNTY_FIPS"].isna().sum()
    if unmatched:
        print(f"Warning: {unmatched} rows unmatched")

    base = base.dropna(subset=["COUNTY_FIPS"])

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            base[base["year"].isin(years)]
            .groupby("COUNTY_FIPS")["total_mean"]
            .mean()
            .round(2)
            .reset_index()
            .rename(columns={"total_mean": f"Smoking_rate_{suffix}"})
        )
        col = f"Smoking_rate_{suffix}"
        ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    result = parts[0].merge(parts[1], on="COUNTY_FIPS", how="outer").sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "IHME_Smoking.csv"
    result.to_csv(out_path, index=False)
    print(f"Saved {len(result):,} counties → {out_path}")
    print(result.describe())


if __name__ == "__main__":
    main()
