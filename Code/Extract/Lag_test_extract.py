"""
Extract lag comparison test results from brms_MRR_lag/*_lag_test.csv.

Each file contains pairwise lag comparisons:
    ICD_Code, comparison, diff_mean, diff_lower, diff_upper, P_a_gt_b

P_a_gt_b = posterior probability that lag_a MRD > lag_b MRD.
  - P >= 0.95 → lag_a significantly greater than lag_b  (*)
  - P <= 0.05 → lag_b significantly greater than lag_a  (*)
  - 0.05 < P < 0.95 → no significant difference

Output (wide format — one row per disease):
  Result/Tables/NDD_lag_test.csv
  Result/Tables/Cancer_lag_test.csv

Columns: Disease | lag5 vs. lag10 | lag10 vs. lag15 | lag15 vs. lag5 | Direction
  * appended to the MRD cell when the comparison is significant.
  Direction is derived from all three pairwise significance results.
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

NDD_ABBR = {
    "G20_G30_G12.2_F01_F03": "NDD",
    "G30_F01_F03": "Dementia",
    "G30": "AD",
    "G20": "PD",
    "F01": "VD",
    "G12.2": "ALS",
    "G10": "HD",
}

CANCER_ABBR = {
    "C00_C97": "All Cancer",
    "C34":     "Lung",
    "C18_C21": "Colorectal",
    "C50":     "Breast",
    "C25":     "Pancreas",
    "C61":     "Prostate",
}

input_dir = os.path.join(base_dir, "Result", "brms_MRR_lag")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

SKIP_CODES = {"F03"}

NDD_ORDER = [
    "G20_G30_G12.2_F01_F03",
    "G30_F01_F03",
    "G30",
    "G20",
    "F01",
    "G12.2",
    "G10",
]

CANCER_ORDER = [
    "C00_C97",
    "C34",
    "C18_C21",
    "C50",
    "C25",
    "C61",
]

COMPARISON_ORDER = ["lag5_vs_lag10", "lag10_vs_lag15", "lag15_vs_lag5"]
COMPARISON_LABELS = {
    "lag5_vs_lag10":  "lag5 vs. lag10",
    "lag10_vs_lag15": "lag10 vs. lag15",
    "lag15_vs_lag5":  "lag15 vs. lag5",
}


def to_float_p(p):
    try:
        return float(p)
    except (ValueError, TypeError):
        return 0.0001  # e.g. "p<0.0001"


def sig_sign(p):
    """+1 if a>b significant, -1 if b>a significant, 0 if no diff."""
    if p >= 0.95:
        return 1
    if p <= 0.05:
        return -1
    return 0


def fmt_cell(diff_mean, diff_lower, diff_upper, p_num):
    """Format MRD difference; append * when significant."""
    cell = f"{diff_mean:.2f} ({diff_lower:.2f}, {diff_upper:.2f})"
    if p_num >= 0.95 or p_num <= 0.05:
        cell += "*"
    return cell


def determine_direction(s5v10, s10v15, s15v5):
    """
    Derive a plain-language ordering from three pairwise significance signs.

    s5v10  : +1 = lag5 > lag10 sig,  -1 = lag10 > lag5 sig,  0 = NS
    s10v15 : +1 = lag10 > lag15 sig, -1 = lag15 > lag10 sig, 0 = NS
    s15v5  : +1 = lag15 > lag5 sig,  -1 = lag5 > lag15 sig,  0 = NS
    """
    if s5v10 == 0 and s10v15 == 0 and s15v5 == 0:
        return "No significant difference"

    # --- Two-comparison chains (s5v10 & s10v15) ---
    if s5v10 == 1 and s10v15 == 1:
        return "lag5 > lag10 > lag15"
    if s5v10 == -1 and s10v15 == -1:
        return "lag15 > lag10 > lag5"
    if s5v10 == 1 and s10v15 == -1:          # lag5 > lag10, lag15 > lag10 → trough at lag10
        if s15v5 == 1:
            return "lag10 < lag5 < lag15"
        if s15v5 == -1:
            return "lag10 < lag15 < lag5"
        return "lag5 > lag10, lag15 > lag10"
    if s5v10 == -1 and s10v15 == 1:          # lag10 > lag5, lag10 > lag15 → peak at lag10
        if s15v5 == 1:
            return "lag5 < lag15 < lag10"
        if s15v5 == -1:
            return "lag15 < lag5 < lag10"
        return "lag10 > lag5, lag10 > lag15"

    # --- Two-comparison chains (s10v15 & s15v5), s5v10 == 0 ---
    if s10v15 == 1 and s15v5 == 1:
        return "lag10 > lag15 > lag5"
    if s10v15 == -1 and s15v5 == -1:
        return "lag5 > lag15 > lag10"
    if s10v15 == -1 and s15v5 == 1:          # lag15 > lag10, lag15 > lag5 → peak at lag15
        return "lag15 > lag10, lag15 > lag5"
    if s10v15 == 1 and s15v5 == -1:          # lag10 > lag15, lag5 > lag15 → trough at lag15
        return "lag10 > lag15, lag5 > lag15"

    # --- Two-comparison chains (s5v10 & s15v5), s10v15 == 0 ---
    if s5v10 == 1 and s15v5 == 1:
        return "lag15 > lag5 > lag10"
    if s5v10 == -1 and s15v5 == -1:
        return "lag10 > lag5 > lag15"
    if s5v10 == 1 and s15v5 == -1:           # lag5 > lag10, lag5 > lag15 → peak at lag5
        return "lag5 > lag10, lag5 > lag15"
    if s5v10 == -1 and s15v5 == 1:           # lag10 > lag5, lag15 > lag5 → trough at lag5
        return "lag5 < lag10, lag5 < lag15"

    # --- Single significant comparison ---
    if s5v10 == 1:   return "lag5 > lag10"
    if s5v10 == -1:  return "lag10 > lag5"
    if s10v15 == 1:  return "lag10 > lag15"
    if s10v15 == -1: return "lag15 > lag10"
    if s15v5 == 1:   return "lag15 > lag5"
    if s15v5 == -1:  return "lag5 > lag15"

    return "No significant difference"


def format_and_save(frames, abbr_map, disease_order, output_path, label):
    if not frames:
        print(f"No files found for {label}, skipping.")
        return

    combined = pd.concat(frames, ignore_index=True)
    combined["_p_num"] = combined["P_a_gt_b"].map(to_float_p)

    rows = []
    for icd_code in disease_order:
        sub = combined[combined["ICD_Code"] == icd_code]
        if sub.empty:
            continue

        row = {"Disease": abbr_map.get(icd_code, icd_code)}
        p_signs = {}

        for comp in COMPARISON_ORDER:
            comp_row = sub[sub["comparison"] == comp]
            col = COMPARISON_LABELS[comp]
            if comp_row.empty:
                row[col] = ""
                p_signs[comp] = 0
            else:
                r = comp_row.iloc[0]
                p_num = r["_p_num"]
                p_signs[comp] = p_num
                row[col] = fmt_cell(r["diff_mean"], r["diff_lower"], r["diff_upper"], p_num)

        row["Direction"] = determine_direction(
            sig_sign(p_signs.get("lag5_vs_lag10", 0)),
            sig_sign(p_signs.get("lag10_vs_lag15", 0)),
            sig_sign(p_signs.get("lag15_vs_lag5", 0)),
        )
        rows.append(row)

    result = pd.DataFrame(rows, columns=["Disease"] + list(COMPARISON_LABELS.values()) + ["Direction"])
    result.to_csv(output_path, index=False)
    fname = os.path.basename(output_path)
    print(f"Saved {fname} ({len(result)} rows)")
    print(result.to_string(index=False))


ndd_frames = []
cancer_frames = []

for fname in sorted(f for f in os.listdir(input_dir) if f.endswith("_lag_test.csv")):
    icd_code = fname.replace("_lag_test.csv", "")
    if icd_code in SKIP_CODES:
        continue
    df = pd.read_csv(os.path.join(input_dir, fname))
    if icd_code in CANCER_ORDER:
        cancer_frames.append(df)
    elif not icd_code.startswith("C"):
        ndd_frames.append(df)

format_and_save(
    ndd_frames,
    NDD_ABBR,
    NDD_ORDER,
    os.path.join(output_dir, "NDD_lag_test.csv"),
    "NDD",
)

format_and_save(
    cancer_frames,
    CANCER_ABBR,
    CANCER_ORDER,
    os.path.join(output_dir, "Cancer_lag_test.csv"),
    "Cancer",
)
