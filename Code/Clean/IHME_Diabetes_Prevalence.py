#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
IHME 糖尿病患病率数据清洗脚本
数据源：IHME USA County Diabetes Prevalence 1999-2012 (Diagnosed sheet)
输出：IHME_Diabetes_Prevalence.csv — 县级已诊断糖尿病率 0005/0610 均值
"""

import re
import pandas as pd
from pathlib import Path

BASE     = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/IHME Diabetes Prevalence"
DATA_OUT = BASE / "Data/Processed/IHME"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

DIAB_FILE = DATA_IN / "IHME_USA_COUNTY_DIABETES_PREVALENCE_1999_2012_NATIONAL_Y2016M08D23.XLSX"

PERIODS = {
    "0005": set(range(2000, 2006)),
    "0610": set(range(2006, 2011)),
}
ALL_YEARS = set().union(*PERIODS.values())


def main():
    df = pd.read_excel(DIAB_FILE, sheet_name="Diagnosed", header=1)
    df = df.dropna(subset=["FIPS"])
    df["FIPS"] = df["FIPS"].astype(float).astype(int).astype(str).str.zfill(5)
    df = df[df["FIPS"].str.len() == 5].copy()

    year_cols = {
        col: int(m.group(1))
        for col in df.columns
        if (m := re.match(r"^Prevalence, (\d{4}), Both Sexes$", str(col)))
        and int(m.group(1)) in ALL_YEARS
    }
    if not year_cols:
        raise ValueError(f"No prevalence columns found for years {ALL_YEARS}")

    long = df[["FIPS"] + list(year_cols)].melt(
        id_vars="FIPS", var_name="year_col", value_name="rate"
    )
    long["year"] = long["year_col"].map(year_cols)

    parts = []
    for suffix, years in PERIODS.items():
        part = (
            long[long["year"].isin(years)]
            .dropna(subset=["rate"])
            .groupby("FIPS")["rate"]
            .mean()
            .round(2)
            .reset_index()
            .rename(columns={"FIPS": "COUNTY_FIPS", "rate": f"Diabetes_Prevalence_rate_{suffix}"})
        )
        col = f"Diabetes_Prevalence_rate_{suffix}"
        ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
        ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
        part = ref.merge(part, on="COUNTY_FIPS", how="left")
        n_missing = part[col].isna().sum()
        part[col] = part[col].fillna(part[col].mean()).round(2)
        parts.append(part)
        print(f"{suffix}: {len(part):,} counties ({n_missing} imputed), mean={part[col].mean():.1f}%")

    county = parts[0].merge(parts[1], on="COUNTY_FIPS", how="outer").sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "IHME_Diabetes_Prevalence.csv"
    county.to_csv(out_path, index=False)
    print(f"Saved {len(county):,} counties → {out_path}")
    print(county.describe())


if __name__ == "__main__":
    main()
