"""
Extract lag comparison test results from brms_MRR_lag/*_lag_test.csv.

Each file contains pairwise lag comparisons:
    ICD_Code, comparison, diff_mean, diff_lower, diff_upper, P_a_gt_b

Output: Result/Tables/lag_test.csv (long format, sorted by disease then comparison)
"""

import os

import pandas as pd
import yaml

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

config_path = os.path.join(base_dir, "config.yaml")
with open(config_path, "r") as f:
    config = yaml.safe_load(f)

icd_mapping = config["brms_analysis"]["icd_mapping"]
icd_mapping["G30_F01_F03"] = "Dementia"

abbr_mapping = {
    "G20_G30_G12.2_F01_F03": "NDD",
    "G30_F01_F03": "Dementia",
    "G30": "AD",
    "G20": "PD",
    "F01": "VD",
    "G12.2": "ALS",
    "G10": "HD",
}

input_dir = os.path.join(base_dir, "Result", "brms_MRR_lag")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

SKIP_CODES = {"F03"}

DISEASE_ORDER = [
    "G20_G30_G12.2_F01_F03",
    "G30_F01_F03",
    "G30",
    "G20",
    "F01",
    "G12.2",
    "G10",
]

COMPARISON_ORDER = ["lag5_vs_lag10", "lag10_vs_lag15", "lag15_vs_lag5"]

frames = []
for fname in sorted(f for f in os.listdir(input_dir) if f.endswith("_lag_test.csv")):
    icd_code = fname.replace("_lag_test.csv", "")
    if icd_code in SKIP_CODES:
        continue
    df = pd.read_csv(os.path.join(input_dir, fname))
    frames.append(df)

if not frames:
    print("No *_lag_test.csv files found.")
    raise SystemExit(1)

combined = pd.concat(frames, ignore_index=True)
combined["Disease"] = combined["ICD_Code"].map(lambda x: abbr_mapping.get(x, x))

# Sort by disease order then comparison order
combined["_sd"] = (
    combined["ICD_Code"].map({d: i for i, d in enumerate(DISEASE_ORDER)}).fillna(99)
)
combined["_sc"] = (
    combined["comparison"]
    .map({c: i for i, c in enumerate(COMPARISON_ORDER)})
    .fillna(99)
)
combined = (
    combined.sort_values(["_sd", "_sc"])
    .drop(columns=["_sd", "_sc"])
    .reset_index(drop=True)
)

combined["Comparison"] = combined["comparison"].str.replace("_", " ", regex=False)
combined["Difference (95% CI)"] = combined.apply(
    lambda r: f"{r['diff_mean']:.2f}({r['diff_lower']:.2f},{r['diff_upper']:.2f})",
    axis=1,
)
combined["P"] = combined["P_a_gt_b"].apply(
    lambda p: "< 0.05" if p < 0.05 else f"{p:.2f}"
)

combined = combined[["Disease", "Comparison", "Difference (95% CI)", "P"]]

output_path = os.path.join(output_dir, "lag_test.csv")
combined.to_csv(output_path, index=False)
print(f"Saved lag_test.csv ({len(combined)} rows)")
print(combined.to_string(index=False))
