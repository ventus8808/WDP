"""
Extract the GBTM trajectory classification into a compact lookup table.

Input
    Data/Processed/SVI/SVI_GBTM.csv   (COUNTY_FIPS, GBTM_Class = "Class A".."Class D", ...)

Output
    Data/Processed/SVI/SVI_Result.csv (COUNTY_FIPS, SVI = A/B/C/D)
"""

from pathlib import Path

import pandas as pd

IN_PATH = Path("Data/Processed/SVI/SVI_GBTM.csv")
OUT_PATH = Path("Data/Processed/SVI/SVI_Result.csv")

df = pd.read_csv(IN_PATH, dtype={"COUNTY_FIPS": str})
df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)

out = pd.DataFrame({
    "COUNTY_FIPS": df["COUNTY_FIPS"],
    "SVI": df["GBTM_Class"].str.split().str[-1],   # "Class A" -> "A"
})

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_PATH, index=False)
print(f"Saved {len(out)} rows -> {OUT_PATH}")
print(out["SVI"].value_counts(dropna=False).sort_index())
