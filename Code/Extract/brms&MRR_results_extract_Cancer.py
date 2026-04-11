"""
Extract MRD + MRR results for top-5 cancer types, combined into one wide table.

Output: Result/Tables/Cancer_MRR_MRD_disease_stratified.csv

Row structure (one EQI_Period × cancer × metric × lag):
  ICD_Code | Outcome | EQI_Period | Effect | Lag | Q1 | Q2 | Q3 | Q4 | Q5

Order per cancer: MRD rows (lag 5→10→15) then MRR rows (lag 5→10→15)
  MRD rows — Q1 = 0.0 (ref), Q2–Q5 = "mean(lower,upper)"  (asterisks stripped)
  MRR rows — Q1 = 1.0 (ref), Q2–Q5 = "mean(lower,upper)"  (2 decimal places)

MRD source : Result/brms/*_main.csv       (Model == EQI)
MRR source : Result/brms_MRR_lag/*_MRR.csv
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

mrd_dir = os.path.join(base_dir, "Result", "brms")
mrr_dir = os.path.join(base_dir, "Result", "brms_MRR_lag")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

TOP5_CANCER = ["C00_C97", "C34", "C18_C21", "C50", "C25", "C61"]
CANCER_ORDER = {icd: i for i, icd in enumerate(TOP5_CANCER)}
EQI_PERIODS = ["2000_2005"]
QUINTILES = ["Q1", "Q2", "Q3", "Q4", "Q5"]


def strip_stars(cell):
    if pd.isna(cell):
        return cell
    return re.sub(r"[\*]+\s*$|ns\s*$", "", str(cell).strip())


# ---------------------------------------------------------------------------
# MRD — Result/brms/*_main.csv, Model == EQI
# ---------------------------------------------------------------------------
mrd_rows = []
for icd in TOP5_CANCER:
    fpath = os.path.join(mrd_dir, f"{icd}_main.csv")
    if not os.path.exists(fpath):
        print(f"MRD file not found: {fpath}")
        continue
    raw = pd.read_csv(fpath)
    if "Model" in raw.columns:
        raw["Model"] = raw["Model"].astype(str).str.strip('"')
    raw = raw[(raw["Model"] == "EQI") & (raw["EQI_Period"].isin(EQI_PERIODS))]

    for _, r in raw.iterrows():
        row = {
            "ICD_Code":   icd,
            "Outcome":    icd_mapping.get(icd, icd),
            "EQI_Period": r["EQI_Period"],
            "Effect":     "MRD",
            "Lag":        int(r["Lag"]),
            "_metric":    0,
            "_order":     CANCER_ORDER[icd],
        }
        for q in QUINTILES:
            row[q] = strip_stars(r.get(q)) if q in raw.columns else None
        mrd_rows.append(row)

# ---------------------------------------------------------------------------
# MRR — Result/brms_MRR_lag/*_MRR.csv, pivot quintiles to wide
# ---------------------------------------------------------------------------
# The CSV may contain appended runs; keep the last occurrence of each
# (ICD_Code, EQI_Period, Lag, Quintile) to get the most recent results.
# Q1 rows for Lag 10 and Lag 15 are dropped — Lag 5 Q1 is the universal
# reference (shown as "1.0") so inter-lag Q1 values are not displayed.
mrr_wide = {}  # key: (icd, eqi_period, lag) → {Q2: "x(l,u)", ...}
for icd in TOP5_CANCER:
    fpath = os.path.join(mrr_dir, f"{icd}_MRR.csv")
    if not os.path.exists(fpath):
        print(f"MRR file not found: {fpath}")
        continue
    raw = pd.read_csv(fpath)
    raw = raw[raw["EQI_Period"].isin(EQI_PERIODS)].copy()
    # Keep only the last row for each unique key (handles appended duplicates)
    raw = raw.drop_duplicates(
        subset=["ICD_Code", "EQI_Period", "Lag", "Quintile"], keep="last"
    )
    # Drop Q1 for Lag 10 and Lag 15 — not shown in output table
    raw = raw[~((raw["Quintile"] == "Q1") & (raw["Lag"] != 5))]
    for _, r in raw.iterrows():
        key = (icd, r["EQI_Period"], int(r["Lag"]))
        if key not in mrr_wide:
            mrr_wide[key] = {}
        mrr_wide[key][r["Quintile"]] = (
            f"{r['MRR_mean']:.2f}({r['MRR_lower']:.2f},{r['MRR_upper']:.2f})"
        )

mrr_rows = []
for (icd, eqi_period, lag), q_vals in mrr_wide.items():
    row = {
        "ICD_Code":   icd,
        "Outcome":    icd_mapping.get(icd, icd),
        "EQI_Period": eqi_period,
        "Effect":     "MRR",
        "Lag":        lag,
        "_metric":    1,
        "_order":     CANCER_ORDER[icd],
        "Q1":         "1.0",  # Lag 5 Q1 is the universal reference
    }
    for q in ["Q2", "Q3", "Q4", "Q5"]:
        row[q] = q_vals.get(q)
    mrr_rows.append(row)

# ---------------------------------------------------------------------------
# Combine, sort, insert blank rows between diseases, save
# ---------------------------------------------------------------------------
combined = pd.concat(
    [pd.DataFrame(mrd_rows), pd.DataFrame(mrr_rows)], ignore_index=True
)

combined = (
    combined.sort_values(["_order", "EQI_Period", "_metric", "Lag"])
    .drop(columns=["_order", "_metric"])
    .reset_index(drop=True)
)

combined["Outcome"] = combined["Outcome"] + " (" + combined["Effect"] + ")"
combined = combined[["Outcome", "Lag"] + QUINTILES]

blank = pd.DataFrame([[""] * len(combined.columns)], columns=combined.columns)
cancer_names = [icd_mapping.get(icd, icd) for icd in TOP5_CANCER]
groups = [combined[combined["Outcome"].str.startswith(name)] for name in cancer_names
          if any(combined["Outcome"].str.startswith(name))]
combined_out = pd.concat(
    [g for pair in zip(groups, [blank] * len(groups)) for g in pair][:-1],
    ignore_index=True,
)

output_path = os.path.join(output_dir, "Cancer_MRR_MRD_disease_stratified.csv")
combined_out.to_csv(output_path, index=False)
print(f"Saved Cancer_MRR_MRD_disease_stratified.csv ({len(combined)} rows)")
print(combined.head(12).to_string(index=False))
