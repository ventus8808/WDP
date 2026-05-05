#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAHIE 无保险率数据清洗脚本
数据源及格式：
  2000     — Excel (.xls)，All Ages Percent uninsured（col 7，含国家/州行需过滤）
  2005     — 空格分隔 txt，"Datalines Follow:" 标记，pctelig at parts[12]
  2006/07  — 管道分隔 txt，PCTUI at parts[14]（iprcat=0 时等价于 PCTELIG）
  2008-23  — CSV，header 自动定位，PCTELIG，geocat=50/agecat=0/racecat=0/sexcat=0/iprcat=0
输出：Insurance.csv — 县级无保险率五期均值（%，保留2位小数）
"""

import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/SAHIE"
DATA_OUT = BASE / "Data/Processed/Covariate"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

PERIODS = {
    "0005": [2000, 2005],
    "0610": [2006, 2007, 2008, 2009, 2010],
    "1115": list(range(2011, 2016)),
    "1620": list(range(2016, 2021)),
    "2124": [2021, 2022, 2023],
}


def read_year_2000() -> pd.DataFrame:
    df = pd.read_excel(DATA_IN / "sahie2000.xls", header=4, engine="xlrd")
    # col layout (0-indexed): 0=postal, 1=state_fips, 2=county_fips, 3=name,
    #   4=NI_allages, 5=NU_allages, 6=CI, 7=Pct_uninsured_allages, ...
    df.columns = [str(c) for c in df.columns]
    state_col  = df.columns[1]
    county_col = df.columns[2]
    rate_col   = df.columns[7]
    df = df[[state_col, county_col, rate_col]].copy()
    df.columns = ["state_fips", "county_fips", "rate"]
    df["rate"] = pd.to_numeric(df["rate"], errors="coerce")
    df = df.dropna(subset=["state_fips", "county_fips", "rate"])
    df = df[pd.to_numeric(df["county_fips"], errors="coerce") > 0].copy()
    df["COUNTY_FIPS"] = (
        df["state_fips"].astype(float).astype(int).astype(str).str.zfill(2)
        + df["county_fips"].astype(float).astype(int).astype(str).str.zfill(3)
    )
    df["year"] = 2000
    return df[["COUNTY_FIPS", "year", "rate"]]


def read_year_2005() -> pd.DataFrame:
    rows = []
    in_data = False
    with open(DATA_IN / "sahie-2005.txt") as f:
        for line in f:
            if "Datalines Follow:" in line:
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.strip().split()
            if len(parts) < 13:
                continue
            try:
                sf, cf, geocat, agecat, racecat, sexcat, iprcat = (int(parts[i]) for i in range(7))
                rate = float(parts[12])
            except ValueError:
                continue
            if geocat == 50 and agecat == 0 and racecat == 0 and sexcat == 0 and iprcat == 0:
                rows.append({"COUNTY_FIPS": f"{sf:02d}{cf:03d}", "year": 2005, "rate": round(rate, 2)})
    return pd.DataFrame(rows)


def read_year_pipe(year: int) -> pd.DataFrame:
    """2006/2007 pipe-delimited: PCTUI at index 14 (iprcat=0 → equivalent to PCTELIG)"""
    rows = []
    with open(DATA_IN / f"sahie-{year}.txt") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 15:
                continue
            try:
                sf = int(parts[1])
                cf = int(parts[2])
                geocat, agecat, racecat, sexcat, iprcat = (int(parts[i]) for i in (3, 4, 5, 6, 7))
                rate = float(parts[14])
            except ValueError:
                continue
            if geocat == 50 and agecat == 0 and racecat == 0 and sexcat == 0 and iprcat == 0:
                rows.append({"COUNTY_FIPS": f"{sf:02d}{cf:03d}", "year": year, "rate": round(rate, 2)})
    return pd.DataFrame(rows)


def read_year_csv(year: int) -> pd.DataFrame:
    """2008-2023 CSV: auto-detect header line, use PCTELIG"""
    filepath = DATA_IN / f"sahie_{year}.csv"
    with open(filepath) as f:
        lines = f.readlines()
    header_idx = next(i for i, l in enumerate(lines) if l.strip().startswith("year,"))
    df = pd.read_csv(filepath, skiprows=header_idx, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(lambda s: s.str.strip() if s.dtype == object else s)
    df["geocat"]  = pd.to_numeric(df["geocat"],  errors="coerce")
    df["agecat"]  = pd.to_numeric(df["agecat"],  errors="coerce")
    df["racecat"] = pd.to_numeric(df["racecat"], errors="coerce")
    df["sexcat"]  = pd.to_numeric(df["sexcat"],  errors="coerce")
    df["iprcat"]  = pd.to_numeric(df["iprcat"],  errors="coerce")
    df = df[
        (df["geocat"] == 50) & (df["agecat"] == 0) &
        (df["racecat"] == 0) & (df["sexcat"] == 0) & (df["iprcat"] == 0)
    ].copy()
    df["COUNTY_FIPS"] = (
        df["statefips"].astype(int).astype(str).str.zfill(2)
        + df["countyfips"].astype(int).astype(str).str.zfill(3)
    )
    df["rate"] = pd.to_numeric(df["PCTELIG"], errors="coerce").round(2)
    df["year"] = year
    return df[["COUNTY_FIPS", "year", "rate"]].dropna(subset=["rate"])


def main():
    all_dfs = []
    all_dfs.append(read_year_2000())
    all_dfs.append(read_year_2005())
    for yr in [2006, 2007]:
        all_dfs.append(read_year_pipe(yr))
    for yr in range(2008, 2024):
        all_dfs.append(read_year_csv(yr))

    long = pd.concat(all_dfs, ignore_index=True)

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            long[long["year"].isin(years)]
            .groupby("COUNTY_FIPS")["rate"]
            .mean()
            .round(2)
            .reset_index()
            .rename(columns={"rate": f"Uninsured_rate_{suffix}"})
        )
        col = f"Uninsured_rate_{suffix}"
        eqi_key = "0005" if suffix == "0005" else "0610"
        ref = pd.read_csv(EQI_DIR / f"EQI{eqi_key}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    result = parts[0]
    for part in parts[1:]:
        result = result.merge(part, on="COUNTY_FIPS", how="outer")

    result = result.sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "Uninsured_rate.csv"
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result):,} counties → {out_path}")
    print(result.describe())


if __name__ == "__main__":
    main()
