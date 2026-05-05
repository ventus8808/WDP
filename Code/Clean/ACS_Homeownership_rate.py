#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Census Bureau ACS 住房自有率数据清洗脚本
数据源：NHGIS ACS 5年期估计 - B25003 Housing Tenure（住房产权）
变量：自有住房率 = 自有占用住房单元 / 全部占用住房单元
三分位在 EQI 参考县填补缺失值后独立计算（低/中/高）
输出：Homeownership_rate.csv — 县级自有住房率及三分位，四个结局期
"""

from pathlib import Path
import pandas as pd

BASE     = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/Census Bureau ACS"
DATA_OUT = BASE / "Data/Processed/Covariate"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

PERIOD_FILES = [
    ("nhgis0002_ds176_20105_county_E.csv", "0610", "JRKE001",  "JRKE002"),
    ("nhgis0002_ds215_20155_county_E.csv", "1115", "ADP0E001", "ADP0E002"),
    ("nhgis0002_ds249_20205_county_E.csv", "1620", "AMUFE001", "AMUFE002"),
    ("nhgis0002_ds272_20245_county_E.csv", "2124", "AUUEE001", "AUUEE002"),
]

EQI_REF = pd.read_csv(EQI_DIR / "EQI0610.csv", usecols=["COUNTY_FIPS"], dtype=str)
EQI_REF["COUNTY_FIPS"] = EQI_REF["COUNTY_FIPS"].str.zfill(5)


def read_period(filename: str, total_col: str, owner_col: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_IN / filename, dtype={"STATEA": str, "COUNTYA": str})
    df["COUNTY_FIPS"] = df["STATEA"].str.zfill(2) + df["COUNTYA"].str.zfill(3)
    df = df[df["COUNTY_FIPS"].str.len() == 5].copy()
    df["total"] = pd.to_numeric(df[total_col], errors="coerce")
    df["owned"] = pd.to_numeric(df[owner_col], errors="coerce")
    df = df[(df["total"] > 0) & df["owned"].notna()].copy()
    df["rate"] = (df["owned"] / df["total"] * 100).round(2)
    return df[["COUNTY_FIPS", "rate"]]


def assign_tertile(series: pd.Series) -> pd.Series:
    return pd.qcut(series.rank(method="first"), q=3, labels=[1, 2, 3]).astype("Int64")


def main():
    parts = []
    for filename, suffix, total_col, owner_col in PERIOD_FILES:
        df = read_period(filename, total_col, owner_col)

        # Align to EQI reference; impute missing with column mean
        df = EQI_REF.merge(df, on="COUNTY_FIPS", how="left")
        n_missing = df["rate"].isna().sum()
        df["rate"] = df["rate"].fillna(df["rate"].mean()).round(2)

        # Tertile computed after imputation
        df[f"Homeownership_tertile_{suffix}"] = assign_tertile(df["rate"])
        df = df.rename(columns={"rate": f"Homeownership_rate_{suffix}"})

        parts.append(df)
        print(
            f"{suffix}: {len(df):,} counties ({n_missing} imputed), "
            f"mean={df[f'Homeownership_rate_{suffix}'].mean():.1f}%, "
            f"tertile={df[f'Homeownership_tertile_{suffix}'].value_counts().sort_index().to_dict()}"
        )

    result = parts[0]
    for part in parts[1:]:
        result = result.merge(part, on="COUNTY_FIPS", how="outer")

    result = result.sort_values("COUNTY_FIPS").reset_index(drop=True)

    out_path = DATA_OUT / "Homeownership_rate.csv"
    result.to_csv(out_path, index=False)
    print(f"\nSaved {len(result):,} counties → {out_path}")
    print(result.describe())


if __name__ == "__main__":
    main()
