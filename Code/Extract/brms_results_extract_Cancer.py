"""
Extract MRD + MRR results for top-5 cancer types, combined into one wide table.

Output: Result/Tables/Cancer_MRR_MRD_disease_stratified.csv

Row structure (one EQI_Period × cancer × lag × metric):
  ICD_Code | Outcome | EQI_Period | Lag | Effect | Q1 | Q2 | Q3 | Q4 | Q5

Order per cancer: 5 MRD, 10 MRD, 15 MRD, then 5 MRR, 10 MRR, 15 MRR
  MRD rows  — Q1 = 0.0 (ref), Q2–Q5 = "mean(lower,upper)"  (asterisks stripped)
  MRR rows  — Q1 = 1.0 (ref), Q2–Q5 = "mean(lower,upper)"  (2 decimal places)

Sources: Result/brms_MRR_lag/*_main.csv  (MRD)
         Result/brms_MRR_lag/*_MRR.csv   (MRR)
"""

import os
import re

import pandas as pd
import yaml

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

config_path = os.path.join(base_dir, "config.yaml")
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

icd_mapping = config["brms_analysis"]["icd_mapping"]

input_dir = os.path.join(base_dir, "Result", "brms_MRR_lag")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

TOP5_CANCER = ["C00_C97", "C34", "C18_C21", "C50", "C25", "C61"]
CANCER_ORDER = {icd: i for i, icd in enumerate(TOP5_CANCER)}
QUINTILES = ["Q1", "Q2", "Q3", "Q4", "Q5"]


def strip_stars(cell):
    """Remove trailing asterisks/ns from a formatted cell string."""
    if pd.isna(cell):
        return cell
    return re.sub(r"[\*]+\s*$|ns\s*$", "", str(cell).strip())


# ---------------------------------------------------------------------------
# MRD — read *_main.csv, strip asterisks, reshape to output rows
# ---------------------------------------------------------------------------
mrd_rows = []
for icd in TOP5_CANCER:
    fpath = os.path.join(input_dir, f"{icd}_main.csv")
    if not os.path.exists(fpath):
        print(f"MRD file not found: {fpath}")
        continue
    raw = pd.read_csv(fpath)
    if "Model" in raw.columns:
        raw["Model"] = raw["Model"].astype(str).str.strip('"')
    raw = raw[raw["Model"] == "EQI"] if "Model" in raw.columns else raw

    for _, r in raw.iterrows():
        row = {
            "ICD_Code":   icd,
            "Outcome":    icd_mapping.get(icd, icd),
            "EQI_Period": r["EQI_Period"],
            "Lag":        int(r["Lag"]),
            "Effect":     "MRD",
            "_lag_num":   int(r["Lag"]),
            "_metric":    0,  # MRD sorts before MRR
        }
        for q in QUINTILES:
            row[q] = strip_stars(r.get(q)) if q in raw.columns else None
        mrd_rows.append(row)

# ---------------------------------------------------------------------------
# MRR — read *_MRR.csv, pivot quintiles to wide, format mean(lower,upper)
# ---------------------------------------------------------------------------
mrr_rows_wide = {}  # key: (icd, eqi_period, lag) → {Q2: "x(l,u)", ...}
for icd in TOP5_CANCER:
    fpath = os.path.join(input_dir, f"{icd}_MRR.csv")
    if not os.path.exists(fpath):
        print(f"MRR file not found: {fpath}")
        continue
    raw = pd.read_csv(fpath)
    for _, r in raw.iterrows():
        key = (icd, r["EQI_Period"], int(r["Lag"]))
        if key not in mrr_rows_wide:
            mrr_rows_wide[key] = {}
        q = r["Quintile"]
        mrr_rows_wide[key][q] = (
            f"{r['MRR_mean']:.2f}({r['MRR_lower']:.2f},{r['MRR_upper']:.2f})"
        )

mrr_rows = []
for (icd, eqi_period, lag), q_vals in mrr_rows_wide.items():
    row = {
        "ICD_Code":   icd,
        "Outcome":    icd_mapping.get(icd, icd),
        "EQI_Period": eqi_period,
        "Lag":        lag,
        "Effect":     "MRR",
        "_lag_num":   lag,
        "_metric":    1,  # MRR sorts after MRD
        "Q1":         "1.0",
    }
    for q in ["Q2", "Q3", "Q4", "Q5"]:
        row[q] = q_vals.get(q)
    mrr_rows.append(row)

# ---------------------------------------------------------------------------
# Combine, sort, save
# ---------------------------------------------------------------------------
combined = pd.concat(
    [pd.DataFrame(mrd_rows), pd.DataFrame(mrr_rows)], ignore_index=True
)

combined["_order"] = combined["ICD_Code"].map(CANCER_ORDER).fillna(99)
combined = (
    combined.sort_values(["_order", "EQI_Period", "_metric", "_lag_num"])
    .drop(columns=["_order", "_lag_num", "_metric"])
    .reset_index(drop=True)
)

combined = combined[["ICD_Code", "Outcome", "EQI_Period", "Effect", "Lag"] + QUINTILES]

# Insert a blank row between each disease group
blank = pd.DataFrame([[""] * len(combined.columns)], columns=combined.columns)
groups = [combined[combined["ICD_Code"] == icd] for icd in TOP5_CANCER if icd in combined["ICD_Code"].values]
combined_out = pd.concat(
    [g for pair in zip(groups, [blank] * len(groups)) for g in pair][:-1],
    ignore_index=True,
)

output_path = os.path.join(output_dir, "Cancer_MRR_MRD_disease_stratified.csv")
combined_out.to_csv(output_path, index=False)
print(f"Saved Cancer_MRR_MRD_disease_stratified.csv ({len(combined)} rows)")
print(combined.head(12).to_string(index=False))
