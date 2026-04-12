import os

import pandas as pd

INPUT = "Data/Original/AHRF/AHRF2006.txt"
OUTPUT = "Data/Processed/Doctor/Doctor.csv"

# 1-indexed positions from SAS layout → Python slices [start-1 : start-1+width]
# MD counts (width=5): f1121505..f1121500
# Pop estimates (width=8): f1198405..f1198401, f0453000 (Census 2000)
YEARS = [2005, 2004, 2003, 2002, 2001, 2000]
MD_SLICES = [(839, 844), (844, 849), (849, 854), (854, 859), (859, 864), (864, 869)]
POP_SLICES = [
    (19374, 19382),
    (19382, 19390),
    (19390, 19398),
    (19398, 19406),
    (19406, 19414),
    (19414, 19422),
]

MD_MISSING = 99999
POP_MISSING = 99999999

rows = []
with open(INPUT, "r", encoding="latin-1") as f:
    for line in f:
        if len(line) < 19423:
            continue
        fips = line[1:6].strip()
        if not fips.isdigit():
            continue

        densities = []
        for md_sl, pop_sl in zip(MD_SLICES, POP_SLICES):
            md_raw = line[md_sl[0] : md_sl[1]].strip()
            pop_raw = line[pop_sl[0] : pop_sl[1]].strip()

            try:
                md = int(md_raw)
                pop = int(pop_raw)
            except ValueError:
                continue

            if md == MD_MISSING or pop == POP_MISSING:
                continue
            if pop == 0:
                continue

            densities.append(md / pop * 100_000)

        if densities:
            rows.append(
                {
                    "COUNTY_FIPS": fips,
                    "physician_density_per100k": round(
                        sum(densities) / len(densities), 1
                    ),
                }
            )

df = pd.DataFrame(rows).sort_values("COUNTY_FIPS").reset_index(drop=True)

os.makedirs("Data/Processed/Doctor", exist_ok=True)
df.to_csv(OUTPUT, index=False)
print(f"Done: {len(df)} counties → {OUTPUT}")
print(df.head())
