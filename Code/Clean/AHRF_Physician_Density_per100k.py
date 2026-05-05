#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AHRF 医师密度数据清洗脚本
数据源：AHRF2006.txt (2000-2005), AHRF2008.txt (2006-2007)
注意：AHRF2008 仅含 AMA MD 数据至 2007（AMA 数据滞后一年），故 _0610 基于 2006-2007 均值。
输出：Physician_Density_per100k.csv — 县级医师密度（每10万人）两期均值
"""

import pandas as pd
from pathlib import Path

BASE     = Path(__file__).resolve().parents[2]
DATA_IN  = BASE / "Data/Original/AHRF"
DATA_OUT = BASE / "Data/Processed/Covariate"
EQI_DIR  = BASE / "Data/Processed/EQI"
DATA_OUT.mkdir(parents=True, exist_ok=True)

MD_MISSING  = 99999
POP_MISSING = 99999999

FILE_0005 = DATA_IN / "AHRF2006.txt"
SLICES_0005 = [
    ((839, 844),  (19374, 19382)),   # 2005
    ((844, 849),  (19382, 19390)),   # 2004
    ((849, 854),  (19390, 19398)),   # 2003
    ((854, 859),  (19398, 19406)),   # 2002
    ((859, 864),  (19406, 19414)),   # 2001
    ((864, 869),  (19414, 19422)),   # 2000
]
MIN_LEN_0005 = 19423

FILE_0610 = DATA_IN / "AHRF2008.txt"
SLICES_0610 = [
    ((884, 889),  (18459, 18467)),   # 2007
    ((889, 894),  (18467, 18475)),   # 2006
]
MIN_LEN_0610 = 18476


def read_ahrf(filepath: Path, slices: list, min_len: int) -> pd.DataFrame:
    rows = []
    with open(filepath, "r", encoding="latin-1") as f:
        for line in f:
            if len(line) < min_len:
                continue
            fips = line[1:6].strip()
            if not fips.isdigit():
                continue
            densities = []
            for md_sl, pop_sl in slices:
                md_raw  = line[md_sl[0]:md_sl[1]].strip()
                pop_raw = line[pop_sl[0]:pop_sl[1]].strip()
                try:
                    md  = int(md_raw)
                    pop = int(pop_raw)
                except ValueError:
                    continue
                if md == MD_MISSING or pop == POP_MISSING or pop == 0:
                    continue
                densities.append(md / pop * 100_000)
            if densities:
                rows.append({"COUNTY_FIPS": fips, "density": round(sum(densities) / len(densities), 1)})
    return pd.DataFrame(rows)


def _eqi_ref(suffix: str) -> pd.DataFrame:
    ref = pd.read_csv(EQI_DIR / f"EQI{suffix}.csv", usecols=["COUNTY_FIPS"], dtype=str)
    ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)
    return ref


def main():
    df_0005 = read_ahrf(FILE_0005, SLICES_0005, MIN_LEN_0005).rename(columns={"density": "Physician_Density_per100k_0005"})
    df_0610 = read_ahrf(FILE_0610, SLICES_0610, MIN_LEN_0610).rename(columns={"density": "Physician_Density_per100k_0610"})

    for df, suffix in [(df_0005, "0005"), (df_0610, "0610")]:
        ref = _eqi_ref(suffix)
        col = f"Physician_Density_per100k_{suffix}"
        df_aligned = ref.merge(df, on="COUNTY_FIPS", how="left")
        n_missing = df_aligned[col].isna().sum()
        df_aligned[col] = df_aligned[col].fillna(df_aligned[col].mean()).round(1)
        print(f"{suffix}: {len(df_aligned):,} counties ({n_missing} imputed), mean={df_aligned[col].mean():.1f}")
        if suffix == "0005":
            df_0005 = df_aligned
        else:
            df_0610 = df_aligned

    result = (
        df_0005.merge(df_0610, on="COUNTY_FIPS", how="outer")
        .sort_values("COUNTY_FIPS")
        .reset_index(drop=True)
    )

    out_path = DATA_OUT / "Physician_Density_per100k.csv"
    result.to_csv(out_path, index=False)
    print(f"Saved {len(result):,} counties → {out_path}")
    print(result.describe())


if __name__ == "__main__":
    main()
