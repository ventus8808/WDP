import pandas as pd

INPUT  = "Data/Processed/Environmental/NLCD_JRC.csv"
OUTPUT = "Data/Processed/Environmental/Forest_Coverage.csv"

df = pd.read_csv(INPUT, dtype={"COUNTY_FIPS": str})

# Restrict to EQI exposure window 2000-2005
df = df[df["Year"].between(2000, 2005)]

# Skip rows where either field is missing or total_area is zero
df = df[df["nlcd_forest_km2"].notna() & df["total_area_km2"].notna()]
df = df[df["total_area_km2"] > 0]

# Annual forest coverage rate per county-year
df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
df["forest_pct"] = df["nlcd_forest_km2"] / df["total_area_km2"] * 100

# Mean across years per county
out = (
    df.groupby("COUNTY_FIPS")["forest_pct"]
    .mean()
    .round(1)
    .reset_index()
    .rename(columns={"forest_pct": "forest_coverage"})
    .sort_values("COUNTY_FIPS")
    .reset_index(drop=True)
)

out.to_csv(OUTPUT, index=False)
print(f"Done: {len(out)} counties → {OUTPUT}")
print(out.head())
