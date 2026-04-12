import pandas as pd

INPUT = "Data/Original/SAHIE/sahie-2005.txt"
OUTPUT = "Data/Processed/Insurance/Insurance.csv"

rows = []
in_data = False

with open(INPUT, "r") as f:
    for line in f:
        if "Datalines Follow:" in line:
            in_data = True
            continue
        if not in_data:
            continue
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 13:
            continue
        statefips, countyfips, geocat, agecat, racecat, sexcat, iprcat = (
            int(parts[0]), int(parts[1]), int(parts[2]),
            int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6])
        )
        if geocat == 50 and agecat == 0 and racecat == 0 and sexcat == 0 and iprcat == 0:
            pctelig = float(parts[12])
            fips = f"{statefips:02d}{countyfips:03d}"
            rows.append({"COUNTY_FIPS": fips, "uninsured_rate_2005": round(pctelig, 2)})

df = pd.DataFrame(rows).sort_values("COUNTY_FIPS").reset_index(drop=True)
df.to_csv(OUTPUT, index=False)
print(f"Done: {len(df)} counties → {OUTPUT}")
print(df.head())
