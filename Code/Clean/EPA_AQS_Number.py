import pandas as pd
from pathlib import Path

BASE      = Path(__file__).resolve().parents[2]
INPUT_DIR = BASE / "Data/Original/EPA AIR/annual_conc_by_monitor"
EQI_DIR   = BASE / "Data/Processed/EQI"
OUTPUT    = BASE / "Data/Processed/Covariate/AQS_Number.csv"

VALID_STATES = {str(i).zfill(2) for i in range(1, 57)}

PERIODS = {
    "0005": (range(2000, 2006), EQI_DIR / "EQI0005.csv"),
    "0610": (range(2006, 2011), EQI_DIR / "EQI0610.csv"),
}


def count_sites_for_years(years) -> pd.DataFrame:
    yearly = []
    for year in years:
        fp = INPUT_DIR / f"annual_conc_by_monitor_{year}.csv"
        df = pd.read_csv(fp, dtype=str, usecols=["State Code", "County Code", "Site Num"])
        df = df[df["State Code"].isin(VALID_STATES)]
        df["COUNTY_FIPS"] = df["State Code"].str.zfill(2) + df["County Code"].str.zfill(3)
        sites = (
            df.groupby("COUNTY_FIPS")["Site Num"]
            .nunique()
            .rename(f"sites_{year}")
        )
        yearly.append(sites)
        print(f"  {year}: {len(sites)} counties with monitors")
    return pd.concat(yearly, axis=1)


parts = []
for suffix, (years, eqi_file) in PERIODS.items():
    print(f"\nPeriod {suffix}:")
    combined = count_sites_for_years(years)
    aqs = (
        combined.mean(axis=1)
        .round(1)
        .rename(f"AQS_Number_{suffix}")
        .reset_index()
    )

    ref = pd.read_csv(eqi_file, usecols=["COUNTY_FIPS"], dtype=str)
    ref["COUNTY_FIPS"] = ref["COUNTY_FIPS"].str.zfill(5)

    part = ref.merge(aqs, on="COUNTY_FIPS", how="left")
    part[f"AQS_Number_{suffix}"] = part[f"AQS_Number_{suffix}"].fillna(0)

    parts.append(part)
    n_zero = (part[f"AQS_Number_{suffix}"] == 0).sum()
    print(f"  → {len(part)} counties ({n_zero} with no monitors → set to 0)")

result = parts[0]
for part in parts[1:]:
    result = result.merge(part, on="COUNTY_FIPS", how="outer")

result = result.sort_values("COUNTY_FIPS").reset_index(drop=True)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
result.to_csv(OUTPUT, index=False)
print(f"\nDone: {len(result)} counties → {OUTPUT}")
print(result.describe())
