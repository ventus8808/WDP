import pandas as pd
from pathlib import Path

SVI_DIR = Path("Data/Original/CDC SVI")
OUT_PATH = Path("Data/Processed/SVI/SVI.csv")

# Config: (year, fips_col, svi_col, skiprows)
YEAR_CONFIG = [
    (2000, "STCOFIPS", "USTP",        [1]),   # row 1 is alternate header
    (2010, "FIPS",     "R_PL_THEMES", None),
    (2014, "FIPS",     "RPL_THEMES",  None),
    (2016, "FIPS",     "RPL_THEMES",  None),
    (2018, "FIPS",     "RPL_THEMES",  None),
    (2020, "FIPS",     "RPL_THEMES",  None),
    (2022, "FIPS",     "RPL_THEMES",  None),
]

frames = []
for year, fips_col, svi_col, skiprows in YEAR_CONFIG:
    df = pd.read_csv(
        SVI_DIR / f"SVI_{year}_US_county.csv",
        encoding="utf-8-sig",
        skiprows=skiprows,
        usecols=[fips_col, svi_col],
        dtype={fips_col: str},
    )
    df[fips_col] = df[fips_col].str.zfill(5)
    df[svi_col] = pd.to_numeric(df[svi_col], errors="coerce")
    df.loc[df[svi_col] == -999, svi_col] = pd.NA
    df = df.rename(columns={fips_col: "COUNTY_FIPS", svi_col: f"{year}_SVI_RPL"})
    frames.append(df)

result = frames[0]
for df in frames[1:]:
    result = result.merge(df, on="COUNTY_FIPS", how="outer")

result = result.sort_values("COUNTY_FIPS").reset_index(drop=True)
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT_PATH, index=False)
print(f"Saved {len(result)} counties → {OUT_PATH}")
print(result.head())
