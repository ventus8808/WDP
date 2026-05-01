import pandas as pd
from pathlib import Path

INPUT_DIR = Path("Data/Original/EPA AIR/annual_conc_by_monitor")
OUTPUT    = Path("Data/Processed/Monitoring/Monitoring.csv")

# Valid US state FIPS codes (01-56, excluding territories 60,66,69,72,78)
VALID_STATES = {str(i).zfill(2) for i in range(1, 57)}

yearly_counts = []

for year in range(2000, 2007):
    fp = INPUT_DIR / f"annual_conc_by_monitor_{year}.csv"
    df = pd.read_csv(fp, dtype=str, usecols=["State Code", "County Code", "Site Num"])

    # Filter to US states only
    df = df[df["State Code"].isin(VALID_STATES)]

    # Build 5-digit COUNTY_FIPS
    df["COUNTY_FIPS"] = df["State Code"].str.zfill(2) + df["County Code"].str.zfill(3)

    # Count unique sites per county this year
    sites = (
        df.groupby("COUNTY_FIPS")["Site Num"]
        .nunique()
        .rename(f"sites_{year}")
    )
    yearly_counts.append(sites)
    print(f"{year}: {len(sites)} counties")

# Combine all years and compute mean
combined = pd.concat(yearly_counts, axis=1)
combined["site_number_mean"] = combined.mean(axis=1).round(1)

out = (
    combined[["site_number_mean"]]
    .reset_index()
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUTPUT, index=False)
print(f"\nDone: {len(out)} counties → {OUTPUT}")
print(out.head())

# Print mean site number by RUCC stratum
eqi = pd.read_csv(
    "Data/Original/EPA EQI/00_05_EQI.csv",
    usecols=["stfips", "cat_rucc"],
    dtype={"stfips": str, "cat_rucc": str},
)
eqi["COUNTY_FIPS"] = eqi["stfips"].str.zfill(5)
merged = out.merge(eqi[["COUNTY_FIPS", "cat_rucc"]], on="COUNTY_FIPS", how="inner")
print("\nMean site number by RUCC stratum:")
for rucc in sorted(merged["cat_rucc"].dropna().unique()):
    mean_val = merged.loc[merged["cat_rucc"] == rucc, "site_number_mean"].mean()
    n = merged["cat_rucc"].eq(rucc).sum()
    print(f"  RUCC {rucc}: {mean_val:.2f} (n={n} counties)")
