"""
Build a county-level joint EQI-Air x SVI exposure table.

Air domain index comes from the EPA EQI 2000-2005 file (continuous
air_EQI_22July2013, no RUCC stratification). SVI is summarised as each county's
mean percentile across the 2000-2022 trajectory (Data/Processed/SVI/SVI.csv).

For each exposure we add:
  *_2  binary  classification (median split -> 1 = low, 2 = high)
  *_3  tertile classification (1 = low, 2 = mid, 3 = high)
and two cross-classified joint variables:
  EQI_Air_SVI_2  2x2 joint, e.g. a1s1 / a1s2 / a2s1 / a2s2
  EQI_Air_SVI_3  3x3 joint, e.g. a1s1 ... a3s3

Output: Data/Processed/SVI/Air_SVI.csv
    COUNTY_FIPS, EQI_Air, EQI_Air_2, EQI_Air_3,
    SVI, SVI_2, SVI_3, EQI_Air_SVI_2, EQI_Air_SVI_3
"""

import os

import numpy as np
import pandas as pd
import yaml

HERE = os.path.dirname(__file__)
PROJECT_ROOT = os.path.dirname(os.path.dirname(HERE))
OUTPUT_PATH = os.path.join(PROJECT_ROOT, "Data/Processed/SVI/Air_SVI.csv")
SVI_PATH = os.path.join(PROJECT_ROOT, "Data/Processed/SVI/SVI.csv")


def get_config():
    with open(os.path.join(PROJECT_ROOT, "config.yaml")) as f:
        return yaml.safe_load(f)


def _fips5(s):
    return s.astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(5)


def _binary(x):
    """Median split -> 1 (low) / 2 (high), as Int64 with NaN preserved."""
    q = pd.qcut(x, 2, labels=[1, 2], duplicates="drop")
    return q.astype("Int64")


def _tertile(x):
    """Tertiles -> 1 (low) / 2 (mid) / 3 (high), as Int64 with NaN preserved."""
    q = pd.qcut(x, 3, labels=[1, 2, 3], duplicates="drop")
    return q.astype("Int64")


def _joint(air, svi):
    """Cross-classify two category Series into 'a{air}s{svi}' (NaN if either NA)."""
    out = pd.Series(pd.NA, index=air.index, dtype="object")
    ok = air.notna() & svi.notna()
    out.loc[ok] = ("a" + air[ok].astype(int).astype(str)
                   + "s" + svi[ok].astype(int).astype(str))
    return out


def main():
    cfg = get_config()
    eqi_src = os.path.join(cfg["data_directories"]["original"],
                           "EPA EQI", "00_05_EQI.csv")

    # ── EQI Air (continuous, no RUCC) ──────────────────────────────────────
    eqi = pd.read_csv(eqi_src)
    air = pd.DataFrame({
        "COUNTY_FIPS": _fips5(eqi["stfips"]),
        "EQI_Air": eqi["air_EQI_22July2013"],
    })

    # ── SVI (continuous = mean percentile across the 2000-2022 trajectory) ──
    svi_raw = pd.read_csv(SVI_PATH, dtype={"COUNTY_FIPS": str})
    svi_raw["COUNTY_FIPS"] = svi_raw["COUNTY_FIPS"].str.zfill(5)
    rpl_cols = [c for c in svi_raw.columns if c.endswith("_SVI_RPL")]
    svi = pd.DataFrame({
        "COUNTY_FIPS": svi_raw["COUNTY_FIPS"],
        "SVI": svi_raw[rpl_cols].mean(axis=1, skipna=True).round(4),
    })

    # ── Merge and classify ─────────────────────────────────────────────────
    df = air.merge(svi, on="COUNTY_FIPS", how="inner")
    df["EQI_Air_2"] = _binary(df["EQI_Air"])
    df["EQI_Air_3"] = _tertile(df["EQI_Air"])
    df["SVI_2"] = _binary(df["SVI"])
    df["SVI_3"] = _tertile(df["SVI"])
    df["EQI_Air_SVI_2"] = _joint(df["EQI_Air_2"], df["SVI_2"])
    df["EQI_Air_SVI_3"] = _joint(df["EQI_Air_3"], df["SVI_3"])

    df = df[["COUNTY_FIPS", "EQI_Air", "EQI_Air_2", "EQI_Air_3",
             "SVI", "SVI_2", "SVI_3", "EQI_Air_SVI_2", "EQI_Air_SVI_3"]]

    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print(f"Saved {len(df)} counties -> {OUTPUT_PATH}")
    print("EQI_Air range:", round(df['EQI_Air'].min(), 3), "to", round(df['EQI_Air'].max(), 3))
    print("SVI range:", round(df['SVI'].min(), 3), "to", round(df['SVI'].max(), 3))
    print("2x2 joint:", df["EQI_Air_SVI_2"].value_counts(dropna=False).sort_index().to_dict())
    print("3x3 joint:", df["EQI_Air_SVI_3"].value_counts(dropna=False).sort_index().to_dict())


if __name__ == "__main__":
    main()
