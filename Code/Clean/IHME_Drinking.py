#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 饮酒率数据清洗脚本
数据源：IHME USA County Alcohol Use Prevalence 2002-2012 (Heavy sheet)
注意：Heavy drinking 数据仅覆盖 2005-2012，故 _0005 均值实为 2005 单年值。
FIPS 通过 IHME Diabetes 文件中的 (State, Location)→FIPS 映射表获取（直接匹配）。
输出：IHME_Drinking.csv — 县级重度饮酒率 0005/0610 均值（%，保留2位小数）
"""

import re
import pandas as pd
from pathlib import Path

BASE     = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/IHME Drinking"
DIAB_IN  = BASE / "Data/Original/IHME Diabetes Prevalence"
DATA_OUT = BASE / "Data/Processed/IHME"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

DRINK_FILE = DATA_IN / "IHME_USA_COUNTY_ALCOHOL_USE_PREVALENCE_2002_2012_NATIONAL_Y2015M04D23.XLSX"
DIAB_FILE  = DIAB_IN / "IHME_USA_COUNTY_DIABETES_PREVALENCE_1999_2012_NATIONAL_Y2016M08D23.XLSX"

PERIODS = {
    "0005": {2005},
    "0610": set(range(2006, 2011)),
}

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
    df["State"] = df["FIPS"].str[:2].map(STATE_FIPS)
    return df[["State", "Location", "FIPS"]].rename(columns={"FIPS": "COUNTY_FIPS"})


def main():
    fips_ref = build_fips_lookup()
    raw = pd.read_excel(DRINK_FILE, sheet_name="Heavy")

    all_target = set().union(*PERIODS.values())
    year_cols = {
        col: int(m.group(1))
        for col in raw.columns
        if (m := re.match(r"^(\d{4}) Both Sexes$", str(col))) and int(m.group(1)) in all_target
    }

    df = raw[raw["Location"] != raw["State"]].copy()
    long = df[["State", "Location"] + list(year_cols)].melt(
        id_vars=["State", "Location"], var_name="year_col", value_name="rate"
    )
    long["year"] = long["year_col"].map(year_cols)
    long = long.merge(fips_ref, on=["State", "Location"], how="left")

    unmatched = long["COUNTY_FIPS"].isna().sum()
    if unmatched:
        print(f"Warning: {unmatched} rows unmatched (multi-county merged areas)")

    long = long.dropna(subset=["COUNTY_FIPS"])

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            long[long["year"].isin(years)]
            .groupby("COUNTY_FIPS")["rate"]
            .mean()
            .round(2)
            .reset_index()
            .rename(columns={"rate": f"Heavy_Drinking_rate_{suffix}"})
        )
        col = f"Heavy_Drinking_rate_{suffix}"
        ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    county = parts[0].merge(parts[1], on="COUNTY_FIPS", how="outer").sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "IHME_Drinking.csv"
    county.to_csv(out_path, index=False)
    print(f"Saved {len(county):,} counties → {out_path}")
    print(county.describe())


if __name__ == "__main__":
    main()
