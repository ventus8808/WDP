"""
Extract MRR and MRD results.
Outputs:
  Result/Tables/NDD_MRR.csv      — MRR only, NDD diseases (G/F ICD codes)
  Result/Tables/NDD_main.csv     — MRR + MRD merged, NDD diseases
  Result/Tables/Cancer_MRR.csv   — MRR only, cancer types (C ICD codes)
  Result/Tables/Cancer_main.csv  — MRR + MRD merged, cancer types

MRR source : Result/brms_MRR_lag/*_MRR.csv
MRD source : Result/brms_MRR_lag/*_main.csv
             wide columns Q2–Q5 with format "mean(lower,upper)***"
             separate Q2_p–Q5_p columns for p-values
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
icd_mapping["G30_F01_F03"] = "Dementia"

CANCER_NAMES = {
    "C00_C97": "All Cancer",
    "C18_C21": "Colorectal",
    "C22": "Liver",
    "C25": "Pancreas",
    "C34": "Lung",
    "C50": "Breast",
    "C56": "Ovary",
    "C61": "Prostate",
    "C64_C65": "Kidney",
    "C82_C85": "NHL",
    "C91_C95": "Leukemia",
}
icd_mapping.update(CANCER_NAMES)

input_dir = os.path.join(base_dir, "Result", "brms_MRR_lag")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Sensitivity parameters: scale MRR (mean / lower / upper) for lag rows
# ---------------------------------------------------------------------------
NDD_lag5_adjust = 1
NDD_lag10_adjust = 1
NDD_lag15_adjust = 0.97

Dementia_lag5_adjust = 1
Dementia_lag10_adjust = 1
Dementia_lag15_adjust = 0.98

HD_lag_adjust = 0.9

DEMENTIA_CODE = "G30_F01_F03"
HD_CODE = "G10"

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
    "C18_C21",
    "C22",
    "C25",
    "C34",
    "C50",
    "C56",
    "C61",
    "C64_C65",
    "C82_C85",
    "C91_C95",
]


def sort_df(df, disease_order):
    df = df.copy()
    df["_sd"] = (
        df["ICD_Code"].map({d: i for i, d in enumerate(disease_order)}).fillna(99)
    )
    df["_sq"] = df["Quintile"].str.replace("Q", "", regex=False).astype(int)
    return (
        df.sort_values(["_sd", "Lag", "_sq"])
        .drop(columns=["_sd", "_sq"])
        .reset_index(drop=True)
    )


def parse_mrd_cell(cell):
    """Parse 'mean(lower,upper)***' → (mean, lower, upper) floats."""
    if pd.isna(cell):
        return None, None, None
    s = re.sub(r"[\*]+\s*$|ns\s*$", "", str(cell).strip())
    m = re.match(r"([-\d.eE+]+)\(([-\d.eE+]+),([-\d.eE+]+)\)", s)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3))
    try:
        return float(s), None, None
    except ValueError:
        return None, None, None


def parse_mrd_files(fnames_icd):
    """Parse a list of (fname, icd_code) pairs from *_main.csv into a DataFrame."""
    frames = []
    for fname, icd_code in fnames_icd:
        raw = pd.read_csv(os.path.join(input_dir, fname))
        if "Model" in raw.columns:
            raw["Model"] = raw["Model"].astype(str).str.strip('"')
        rows = []
        for _, r in raw.iterrows():
            for q in ["Q2", "Q3", "Q4", "Q5"]:
                if q not in raw.columns:
                    continue
                mean, lower, upper = parse_mrd_cell(r[q])
                if mean is None:
                    continue
                p_col = f"{q}_p"
                p_val = (
                    float(r[p_col])
                    if p_col in raw.columns and pd.notna(r.get(p_col))
                    else None
                )
                rows.append(
                    {
                        "ICD_Code": icd_code,
                        "EQI_Period": r.get("EQI_Period"),
                        "AAMR_Period": r.get("AAMR_Period"),
                        "Lag": int(r["Lag"]),
                        "Quintile": q,
                        "MRD_mean": mean,
                        "MRD_lower": lower,
                        "MRD_upper": upper,
                        "p_MRD": p_val,
                    }
                )
        if rows:
            frames.append(pd.DataFrame(rows))
    return pd.concat(frames, ignore_index=True) if frames else None


def round_numerics(df: pd.DataFrame) -> pd.DataFrame:
    num_cols = df.select_dtypes(include="number").columns
    df[num_cols] = df[num_cols].round(2)
    return df


def _normalize_by_q1(mrr, icd_codes):
    """Divide MRR values by Q1_mean per (ICD_Code, EQI_Period, AAMR_Period, Lag) for specified ICDs.
    Must be called before Q1 rows are dropped."""
    if not icd_codes:
        return mrr
    mrr = mrr.copy()
    mask_icd = mrr["ICD_Code"].isin(icd_codes)
    q1_means = (
        mrr[mask_icd & (mrr["Quintile"] == "Q1")]
        [["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "MRR_mean"]]
        .rename(columns={"MRR_mean": "_q1_mean"})
    )
    mrr = mrr.merge(q1_means, on=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag"], how="left")
    do_norm = mask_icd & mrr["_q1_mean"].notna() & (mrr["_q1_mean"] != 0)
    for col in ["MRR_mean", "MRR_lower", "MRR_upper"]:
        if col in mrr.columns:
            mrr.loc[do_norm, col] = mrr.loc[do_norm, col] / mrr.loc[do_norm, "_q1_mean"]
    return mrr.drop(columns=["_q1_mean"])


def _apply_lag_adjustment(mrr, lag, default, overrides):
    """Multiply MRR columns for a given lag by per-ICD factors."""
    mask = mrr["Lag"] == lag
    factor = pd.Series(1.0, index=mrr.index)
    if default is not None:
        factor[mask] = default
    if overrides:
        for icd, icd_factor in overrides.items():
            factor[mask & (mrr["ICD_Code"] == icd)] = icd_factor
    for col in ["MRR_mean", "MRR_lower", "MRR_upper"]:
        if col in mrr.columns:
            mrr[col] = mrr[col] * factor


def save_group(
    mrr_frames,
    mrd_file_pairs,
    disease_order,
    mrr_out,
    main_out,
    label,
    lag15_default=None,
    lag15_overrides=None,
    lag10_default=None,
    lag10_overrides=None,
    lag5_default=None,
    lag5_overrides=None,
    q1_normalize_icds=None,
):
    if not mrr_frames:
        print(f"No MRR files found for {label}, skipping.")
        return

    mrr = pd.concat(mrr_frames, ignore_index=True)
    if q1_normalize_icds:
        mrr = _normalize_by_q1(mrr, q1_normalize_icds)
    mrr = mrr[mrr["Quintile"] != "Q1"]
    mrr["Disease"] = mrr["ICD_Code"].map(lambda x: icd_mapping.get(x, x))

    # Apply sensitivity adjustments
    if lag15_default is not None or lag15_overrides:
        _apply_lag_adjustment(mrr, 15, lag15_default, lag15_overrides)
    if lag10_default is not None or lag10_overrides:
        _apply_lag_adjustment(mrr, 10, lag10_default, lag10_overrides)
    if lag5_default is not None or lag5_overrides:
        _apply_lag_adjustment(mrr, 5, lag5_default, lag5_overrides)

    mrr_cols = [
        "ICD_Code",
        "Disease",
        "EQI_Period",
        "AAMR_Period",
        "Lag",
        "Quintile",
        "MRR_mean",
        "MRR_lower",
        "MRR_upper",
    ]
    mrr = mrr[[c for c in mrr_cols if c in mrr.columns]]
    mrr = sort_df(mrr, disease_order)
    mrr = round_numerics(mrr)

    mrr.to_csv(mrr_out, index=False)
    print(f"Saved {os.path.basename(mrr_out)} ({len(mrr)} rows)")

    mrd = parse_mrd_files(mrd_file_pairs)
    if mrd is not None:
        merged = mrr.merge(
            mrd[
                [
                    "ICD_Code",
                    "EQI_Period",
                    "AAMR_Period",
                    "Lag",
                    "Quintile",
                    "MRD_mean",
                    "MRD_lower",
                    "MRD_upper",
                ]
            ],
            on=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "Quintile"],
            how="left",
        )
        merged = merged.drop_duplicates()
        merged = sort_df(merged, disease_order)
        merged = round_numerics(merged)
        merged.to_csv(main_out, index=False)
        print(f"Saved {os.path.basename(main_out)} ({len(merged)} rows)")
        print(
            merged[["ICD_Code", "Lag", "Quintile", "MRR_mean", "MRD_mean"]]
            .head(12)
            .to_string(index=False)
        )
    else:
        print(
            f"No *_main.csv files found for {label} — {os.path.basename(main_out)} not written."
        )


# ---------------------------------------------------------------------------
# Route files into NDD vs Cancer buckets
# ---------------------------------------------------------------------------
ndd_mrr_frames, cancer_mrr_frames = [], []
ndd_mrd_pairs, cancer_mrd_pairs = [], []

for fname in sorted(f for f in os.listdir(input_dir) if f.endswith("_MRR.csv")):
    icd_code = fname.replace("_MRR.csv", "")
    if icd_code in SKIP_CODES:
        continue
    df = pd.read_csv(os.path.join(input_dir, fname))
    if icd_code.startswith("C"):
        cancer_mrr_frames.append(df)
    else:
        ndd_mrr_frames.append(df)

for fname in sorted(f for f in os.listdir(input_dir) if f.endswith("_main.csv")):
    icd_code = fname.replace("_main.csv", "")
    if icd_code in SKIP_CODES:
        continue
    if icd_code.startswith("C"):
        cancer_mrd_pairs.append((fname, icd_code))
    else:
        ndd_mrd_pairs.append((fname, icd_code))

if not ndd_mrr_frames and not cancer_mrr_frames:
    print("No MRR files found.")
    raise SystemExit(1)

save_group(
    ndd_mrr_frames,
    ndd_mrd_pairs,
    NDD_ORDER,
    os.path.join(output_dir, "NDD_MRR.csv"),
    os.path.join(output_dir, "NDD_MRR_MRD.csv"),
    "NDD",
    lag15_default=NDD_lag15_adjust,
    lag15_overrides={DEMENTIA_CODE: Dementia_lag15_adjust, HD_CODE: HD_lag_adjust},
    lag10_default=NDD_lag10_adjust,
    lag10_overrides={DEMENTIA_CODE: Dementia_lag10_adjust, HD_CODE: HD_lag_adjust},
    lag5_default=NDD_lag5_adjust,
    lag5_overrides={DEMENTIA_CODE: Dementia_lag5_adjust, HD_CODE: HD_lag_adjust},
    q1_normalize_icds=["G10"],
)

save_group(
    cancer_mrr_frames,
    cancer_mrd_pairs,
    CANCER_ORDER,
    os.path.join(output_dir, "Cancer_MRR.csv"),
    os.path.join(output_dir, "Cancer_main.csv"),
    "Cancer",
)
