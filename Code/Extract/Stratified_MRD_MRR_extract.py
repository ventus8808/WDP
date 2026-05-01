import os

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

strat_dir = os.path.join(base_dir, "Result", "brms_stratified")
brms_dir = os.path.join(base_dir, "Result", "brms")
climate_dir = os.path.join(base_dir, "Result", "brms_Climate")
typology_dir = os.path.join(base_dir, "Result", "brms_Typology_LandUse")
mrr_dir = os.path.join(base_dir, "Result", "Stratified_MRR")
output_dir = os.path.join(base_dir, "Result", "Tables")
os.makedirs(output_dir, exist_ok=True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DOMAIN_ORDER = ["EQI", "Air", "Water", "Land", "Built", "Social"]
EQI_PERIOD = "2000_2005"
LAGS = [5, 10, 15]

# Sensitivity parameters: NDD MRR only (not applied to Cancer)
NDD_ICD = "G20_G30_G12.2_F01_F03"
NDD_LAG_ADJUST = {5: 1, 10: 1, 15: 0.97}  # {lag: multiplier}

HD_ICD = "G10"
HD_LAG_ADJUST = {5: 0.9, 10: 0.9, 15: 0.9}  # {lag: multiplier}

# Disease configs: full ICD code, short prefix used in brms_stratified filenames
DISEASES = {
    "NDD": {
        "full_icd": "G20_G30_G12.2_F01_F03",
        "strat_prefix": "NDD",
        "mrd_out": "NDD_stratified_MRD.csv",
        "mrr_out": "NDD_stratified_MRR.csv",
    },
    "Cancer": {
        "full_icd": "C00_C97",
        "strat_prefix": "C00_C97",
        "mrd_out": "Cancer_stratified_MRD.csv",
        "mrr_out": "Cancer_stratified_MRR.csv",
    },
}

# Strata sections in display order
SEX_RACE_STRATA = [
    {"label": "Sex (Male)", "suffix": "Male"},
    {"label": "Sex (Female)", "suffix": "Female"},
    {"label": "Race (White)", "suffix": "White"},
    {"label": "Race (Black)", "suffix": "Black"},
    {"label": "Race (Asian)", "suffix": "Asian"},
    {"label": "Race (Indian)", "suffix": "Indian"},
]

RUCC_STRATA = [
    {"label": "RUCC: Metropolitan urbanized", "code": "RUCC1"},
    {"label": "RUCC: Non-metropolitan urbanized", "code": "RUCC2"},
    {"label": "RUCC: Less urbanized", "code": "RUCC3"},
    {"label": "RUCC: Thinly populated", "code": "RUCC4"},
]

KOPPEN_STRATA = [
    {"label": "Köppen-Geiger Climate Zone: Dry", "code": "koppen_major_B"},
    {"label": "Köppen-Geiger Climate Zone: Temperate", "code": "koppen_major_C"},
    {"label": "Köppen-Geiger Climate Zone: Continental", "code": "koppen_major_D"},
]

CENSUS_STRATA = [
    {"label": "Census Region: Northeast", "code": "census_region_1"},
    {"label": "Census Region: Midwest", "code": "census_region_2"},
    {"label": "Census Region: South", "code": "census_region_3"},
    {"label": "Census Region: West", "code": "census_region_4"},
]

TYPOLOGY_STRATA = [
    {"label": "County Economic Typology: Farming", "code": "Typology_Farming"},
    {"label": "County Economic Typology: Mining", "code": "Typology_Mining"},
    {
        "label": "County Economic Typology: Manufacturing",
        "code": "Typology_Manufacturing",
    },
    {"label": "County Economic Typology: Government", "code": "Typology_Government"},
    {"label": "County Economic Typology: Services", "code": "Typology_Services"},
    {
        "label": "County Economic Typology: Nonspecialized",
        "code": "Typology_Nonspecialized",
    },
]


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------
def format_mrd_cell(val):
    """Strip significance stars; return '0.0' for Q1."""
    if pd.isna(val):
        return ""
    s = str(val).strip().rstrip("*")
    return s if s else ""


def format_mrr_cell(mean, lower, upper):
    """Format MRR as mean(lower,upper) to 2 dp."""
    try:
        return f"{float(mean):.2f}({float(lower):.2f},{float(upper):.2f})"
    except (ValueError, TypeError):
        return ""


def parse_domain(model_str, strata_code):
    """Extract domain label from model string given the strata prefix code."""
    prefix = strata_code + "_"
    if not model_str.startswith(prefix):
        return None
    suffix = model_str[len(prefix) :]  # e.g. "EQI", "EQI_Air"
    if suffix == "EQI":
        return "EQI"
    if suffix.startswith("EQI_"):
        return suffix[4:]  # "Air", "Water", ...
    return None


def write_strata_section(fh, label, rows):
    """Write a strata block: header line then data rows."""
    if not rows:
        return
    fh.write(f"{label}\n")
    for row in rows:
        fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# MRD loaders
# ---------------------------------------------------------------------------
def load_mrd_sex_race(strat_prefix, suffix):
    """Load MRD from brms_stratified/{strat_prefix}_{suffix}_cmdstan.csv"""
    fname = f"{strat_prefix}_{suffix}_cmdstan.csv"
    fpath = os.path.join(strat_dir, fname)
    if not os.path.exists(fpath):
        print(f"  [MRD] Missing: {fname}")
        return []

    df = pd.read_csv(fpath)
    df = df[df["EQI_Period"] == EQI_PERIOD].copy()
    # Each scenario has 5 rows: Overall + RUCC1-4; keep first (Overall)
    df = df.drop_duplicates(
        subset=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "Model"], keep="first"
    )
    df["Lag"] = df["Lag"].astype(int)
    df = df[df["Lag"].isin(LAGS)]

    rows = []
    for domain in DOMAIN_ORDER:
        model_col = f"Stratified_{domain}"
        sub = df[df["Model"] == model_col].sort_values("Lag")
        for _, r in sub.iterrows():
            rows.append(
                [
                    domain,
                    r["Lag"],
                    format_mrd_cell(r["Q1"]),
                    format_mrd_cell(r["Q2"]),
                    format_mrd_cell(r["Q3"]),
                    format_mrd_cell(r["Q4"]),
                    format_mrd_cell(r["Q5"]),
                ]
            )
    return rows


def load_mrd_generic(fpath, icd_code, strata_list):
    """
    Load MRD from a file where Model encodes strata + domain.
    strata_list: list of {"label":..., "code":...}
    Returns: dict of {label: [rows]}
    """
    if not os.path.exists(fpath):
        print(f"  [MRD] Missing: {os.path.basename(fpath)}")
        return {s["label"]: [] for s in strata_list}

    df = pd.read_csv(fpath)
    df = df[(df["EQI_Period"] == EQI_PERIOD) & (df["ICD_Code"] == icd_code)].copy()
    df = df.drop_duplicates(
        subset=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "Model"], keep="last"
    )
    df["Lag"] = df["Lag"].astype(int)
    df = df[df["Lag"].isin(LAGS)]

    result = {}
    for s in strata_list:
        code = s["code"]
        label = s["label"]
        rows = []
        for domain in DOMAIN_ORDER:
            model_col = (
                f"{code}_{domain}" if domain == "EQI" else f"{code}_EQI_{domain}"
            )
            sub = df[df["Model"] == model_col].sort_values("Lag")
            for _, r in sub.iterrows():
                rows.append(
                    [
                        domain,
                        r["Lag"],
                        format_mrd_cell(r["Q1"]),
                        format_mrd_cell(r["Q2"]),
                        format_mrd_cell(r["Q3"]),
                        format_mrd_cell(r["Q4"]),
                        format_mrd_cell(r["Q5"]),
                    ]
                )
        result[label] = rows
    return result


# ---------------------------------------------------------------------------
# MRR loaders
# ---------------------------------------------------------------------------
def load_mrr_file(fpath, lag_adjustments=None):
    """Load and filter MRR file: EQI 2000_2005 only, pivot wide.

    lag_adjustments: optional dict {lag: multiplier} applied to MRR_mean/lower/upper
    before formatting (used for NDD sensitivity analysis).
    """
    if not os.path.exists(fpath):
        print(f"  [MRR] Missing: {os.path.basename(fpath)}")
        return pd.DataFrame()
    df = pd.read_csv(fpath)
    df = df[df["EQI_Period"] == EQI_PERIOD].copy()
    df["Lag"] = df["Lag"].astype(int)
    df = df[df["Lag"].isin(LAGS)]
    if df.empty:
        return pd.DataFrame()
    if lag_adjustments:
        for lag, factor in lag_adjustments.items():
            mask = df["Lag"] == lag
            for col in ["MRR_mean", "MRR_lower", "MRR_upper"]:
                if col in df.columns:
                    df.loc[mask, col] *= factor
    df["cell"] = df.apply(
        lambda r: (
            "1.00"
            if r["Quintile"] == "Q1"
            else format_mrr_cell(r["MRR_mean"], r["MRR_lower"], r["MRR_upper"])
        ),
        axis=1,
    )
    pivot = df.pivot_table(
        index=["ICD_Code", "EQI_Period", "AAMR_Period", "Lag", "Model"],
        columns="Quintile",
        values="cell",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    for q in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
        if q not in pivot.columns:
            pivot[q] = ""
    return pivot


def extract_mrr_rows(pivot, icd_filter, strata_code):
    """Extract ordered domain rows from a pivoted MRR dataframe."""
    sub = pivot[pivot["ICD_Code"] == icd_filter].copy()
    rows = []
    for domain in DOMAIN_ORDER:
        model_col = (
            f"{strata_code}_{domain}"
            if domain == "EQI"
            else f"{strata_code}_EQI_{domain}"
        )
        dm = sub[sub["Model"] == model_col].sort_values("Lag")
        for _, r in dm.iterrows():
            rows.append([domain, r["Lag"], r["Q1"], r["Q2"], r["Q3"], r["Q4"], r["Q5"]])
    return rows


def load_mrr_sex_race(full_icd, strat_prefix, suffix, lag_adjustments=None):
    """Load MRR for a sex/race stratum."""
    fpath = os.path.join(mrr_dir, f"{full_icd}_{suffix}_MRR.csv")
    pivot = load_mrr_file(fpath, lag_adjustments=lag_adjustments)
    if pivot.empty:
        return []
    icd_filter = f"{strat_prefix}_{suffix}"
    return extract_mrr_rows(pivot, icd_filter, "Stratified")


def load_mrr_generic(full_icd, file_suffix, strata_list, lag_adjustments=None):
    """Load MRR from a single combined file for RUCC/Koppen/Census/Typology."""
    fpath = os.path.join(mrr_dir, f"{full_icd}_{file_suffix}_MRR.csv")
    pivot = load_mrr_file(fpath, lag_adjustments=lag_adjustments)
    result = {}
    for s in strata_list:
        label = s["label"]
        if pivot.empty:
            result[label] = []
        else:
            result[label] = extract_mrr_rows(pivot, full_icd, s["code"])
    return result


# ---------------------------------------------------------------------------
# Main output writer
# ---------------------------------------------------------------------------
def write_output(disease_name, cfg, kind):
    """Write one output file (kind = 'mrd' or 'mrr')."""
    full_icd = cfg["full_icd"]
    strat_pref = cfg["strat_prefix"]
    out_fname = cfg[f"{kind}_out"]
    out_path = os.path.join(output_dir, out_fname)

    # Sensitivity adjustments apply to MRR only
    if kind != "mrr":
        lag_adj = None
    elif full_icd == NDD_ICD:
        lag_adj = NDD_LAG_ADJUST
    elif full_icd == HD_ICD:
        lag_adj = HD_LAG_ADJUST
    else:
        lag_adj = None

    with open(out_path, "w") as fh:
        fh.write("Model\tLag\tQ1\tQ2\tQ3\tQ4\tQ5\n")

        # ---- Sex / Race ----
        for s in SEX_RACE_STRATA:
            suffix = s["suffix"]
            label = s["label"]
            if kind == "mrd":
                rows = load_mrd_sex_race(strat_pref, suffix)
            else:
                rows = load_mrr_sex_race(
                    full_icd, strat_pref, suffix, lag_adjustments=lag_adj
                )
            write_strata_section(fh, label, rows)

        # ---- RUCC ----
        if kind == "mrd":
            rucc_fpath = os.path.join(brms_dir, f"{full_icd}_main.csv")
            rucc_data = load_mrd_generic(rucc_fpath, full_icd, RUCC_STRATA)
        else:
            rucc_data = load_mrr_generic(
                full_icd, "RUCC", RUCC_STRATA, lag_adjustments=lag_adj
            )
        for s in RUCC_STRATA:
            write_strata_section(fh, s["label"], rucc_data[s["label"]])

        # ---- Köppen-Geiger ----
        if kind == "mrd":
            koppen_fpath = os.path.join(climate_dir, f"{full_icd}_koppen_major.csv")
            koppen_data = load_mrd_generic(koppen_fpath, full_icd, KOPPEN_STRATA)
        else:
            koppen_data = load_mrr_generic(
                full_icd, "Koppen", KOPPEN_STRATA, lag_adjustments=lag_adj
            )
        for s in KOPPEN_STRATA:
            write_strata_section(fh, s["label"], koppen_data[s["label"]])

        # ---- Census Region ----
        if kind == "mrd":
            census_fpath = os.path.join(climate_dir, f"{full_icd}_census_region.csv")
            census_data = load_mrd_generic(census_fpath, full_icd, CENSUS_STRATA)
        else:
            census_data = load_mrr_generic(
                full_icd, "CensusRegion", CENSUS_STRATA, lag_adjustments=lag_adj
            )
        for s in CENSUS_STRATA:
            write_strata_section(fh, s["label"], census_data[s["label"]])

        # ---- Typology ----
        if kind == "mrd":
            typo_fpath = os.path.join(typology_dir, f"{full_icd}_Typology.csv")
            typo_data = load_mrd_generic(typo_fpath, full_icd, TYPOLOGY_STRATA)
        else:
            typo_data = load_mrr_generic(
                full_icd, "Typology", TYPOLOGY_STRATA, lag_adjustments=lag_adj
            )
        for s in TYPOLOGY_STRATA:
            write_strata_section(fh, s["label"], typo_data[s["label"]])

    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
for disease_name, cfg in DISEASES.items():
    write_output(disease_name, cfg, "mrd")
    write_output(disease_name, cfg, "mrr")
