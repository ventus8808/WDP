from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ─── Disease switches ──────────────────────────────────────────────────────────
CVD = 1  # I00_I99
CRD = 1  # J40_J47_J60_J70_J84_D86_C34
CLD = 1  # K70_K76_C22
CKD = 1  # N00_N29_C64_C65
Suicide = 1  # X60_X84_Y87.0
NDD = 1  # G20_G30_G12.2_F01_F03
Cancer = 1  # C00_C97

# ─── Projection switches ───────────────────────────────────────────────────────
Plane = 1
Mercator = 0

# ─── Font switches ─────────────────────────────────────────────────────────────
Font_Helvetica = 1
Font_Georgia = 0
# ──────────────────────────────────────────────────────────────────────────────

plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12
plt.rcParams["mathtext.fontset"] = "stix"

DISEASE_MAP = [
    ("CVD", "I00_I99"),
    ("CRD", "J40_J47_J60_J70_J84_D86_C34"),
    ("CLD", "K70_K76_C22"),
    ("CKD", "N00_N29_C64_C65"),
    ("Suicide", "X60_X84_Y87.0"),
    ("NDD", "G20_G30_G12.2_F01_F03"),
    ("Cancer", "C00_C97"),
]

CONTIGUOUS_STATES = [
    "01",
    "04",
    "05",
    "06",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
    "32",
    "33",
    "34",
    "35",
    "36",
    "37",
    "38",
    "39",
    "40",
    "41",
    "42",
    "44",
    "45",
    "46",
    "47",
    "48",
    "49",
    "50",
    "51",
    "53",
    "54",
    "55",
    "56",
]

OBLIQUE_CRS = (
    "+proj=omerc +lat_0=37 +lonc=-96 +alpha=1 +k=0.9996 "
    "+x_0=0 +y_0=0 +gamma=0 +ellps=WGS84 +units=m +no_defs"
)

AAMR_COLORS = {
    1: "#3a5dae",
    2: "#677eba",
    3: "#93a4cf",
    4: "#e8eaf0",
    5: "#e0bbc0",
    6: "#d88a91",
    7: "#be6e73",
    "No Data": "#cccccc",
}

PERCENTILE_LABELS = [
    "<5th",
    "5th-20th",
    "20th-40th",
    "40th-60th",
    "60th-80th",
    "80th-95th",
    ">95th",
]


def _load_counties(shapefile_path):
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    return counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()


def _bin_aamr(series):
    """Return (levels Series, bin_labels list) using up-to-7-level percentile bins."""
    percentiles = [0, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0]
    # Compute actual bin edges on non-null values; duplicates may reduce bin count
    _, bins = pd.qcut(series.dropna(), q=percentiles, retbins=True, duplicates="drop")
    bins[0] -= 1e-9  # ensure min value is included
    n = len(bins) - 1
    levels = pd.cut(series, bins=bins, labels=range(1, n + 1))
    bin_labels = []
    for i in range(n):
        plabel = PERCENTILE_LABELS[i] if i < len(PERCENTILE_LABELS) else f"Q{i + 1}"
        lo, hi = f"{bins[i]:.1f}", f"{bins[i + 1]:.1f}"
        if i == 0:
            bin_labels.append(f"{plabel} : < {hi} per 100k")
        elif i == n - 1:
            bin_labels.append(f"{plabel} : > {lo} per 100k")
        else:
            bin_labels.append(f"{plabel} : {lo} - {hi} per 100k")
    return levels, bin_labels


def _save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


def plot_aamr(counties_contiguous, shortname, icd, df_all, output_dir, fonts):
    df_icd = (
        df_all[df_all["Outcome"] == icd]
        .drop_duplicates(["COUNTY_FIPS", "Time_Period"])
        .copy()
    )
    if df_icd.empty:
        print(f"  No data for {shortname} ({icd})")
        return

    for time_period, df_period in df_icd.groupby("Time_Period"):
        df_period = df_period.copy()
        print(f"  {shortname} {time_period}...")

        if df_period["AAMR"].notna().sum() < 8:
            print(f"    Skipping: insufficient non-null AAMR values")
            continue

        levels, bin_labels = _bin_aamr(df_period["AAMR"])
        df_period["AAMR_Level"] = levels

        counties_merged = counties_contiguous.merge(
            df_period[["COUNTY_FIPS", "AAMR_Level"]], on="COUNTY_FIPS", how="left"
        )
        counties_merged["color"] = (
            counties_merged["AAMR_Level"]
            .astype(float)
            .map(AAMR_COLORS)
            .fillna(AAMR_COLORS["No Data"])
        )

        legend_elements = [
            mpatches.Patch(
                color=AAMR_COLORS.get(i + 1, AAMR_COLORS["No Data"]),
                label=bin_labels[i],
            )
            for i in range(len(bin_labels))
        ]
        stem = f"{shortname}_{time_period}"

        if Plane:
            state_boundaries = counties_merged.dissolve(by="STATEFP")
            for font in fonts:
                plt.rcParams["font.family"] = font
                fig, ax = plt.subplots(1, 1, figsize=(16, 10))
                counties_merged.plot(
                    color=counties_merged["color"],
                    linewidth=0,
                    edgecolor="none",
                    ax=ax,
                    zorder=1,
                )
                state_boundaries.boundary.plot(
                    ax=ax, color="#000000", linewidth=0.6, alpha=0.6, zorder=2
                )
                ax.set_axis_off()
                ax.legend(
                    handles=legend_elements,
                    bbox_to_anchor=(0.02, 0.02),
                    loc="lower left",
                    frameon=True,
                )
                _save(fig, output_dir / f"{stem}_Plane_{font}.png")

        if Mercator:
            counties_proj = counties_merged.to_crs(OBLIQUE_CRS)
            state_boundaries = counties_proj.dissolve(by="STATEFP")
            for font in fonts:
                plt.rcParams["font.family"] = font
                fig, ax = plt.subplots(1, 1, figsize=(16, 10))
                counties_proj.plot(
                    color=counties_proj["color"],
                    linewidth=0.1,
                    edgecolor="black",
                    ax=ax,
                )
                state_boundaries.boundary.plot(
                    ax=ax, color="black", linewidth=1.2, alpha=0.9
                )
                ax.set_axis_off()
                ax.legend(
                    handles=legend_elements,
                    bbox_to_anchor=(0.02, 0.02),
                    loc="lower left",
                    frameon=True,
                )
                _save(fig, output_dir / f"{stem}_Mercator_{font}.png")


def main():
    project_root = Path(__file__).resolve().parents[2]
    with open(project_root / "config.yaml") as f:
        config = yaml.safe_load(f)

    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    df_all = pd.read_csv(
        project_root / "Data/Processed/df.csv", dtype={"COUNTY_FIPS": str}
    )
    df_all["COUNTY_FIPS"] = df_all["COUNTY_FIPS"].str.zfill(5)

    counties_contiguous = _load_counties(config["data_sources"]["tiger"]["shapefile"])

    switches = {
        "CVD": CVD,
        "CRD": CRD,
        "CLD": CLD,
        "CKD": CKD,
        "Suicide": Suicide,
        "NDD": NDD,
        "Cancer": Cancer,
    }
    selected = [(name, icd) for name, icd in DISEASE_MAP if switches[name]]

    fonts = (["Helvetica"] if Font_Helvetica else []) + (
        ["Georgia"] if Font_Georgia else []
    )
    if not fonts:
        fonts = ["Helvetica"]

    for shortname, icd in selected:
        plot_aamr(counties_contiguous, shortname, icd, df_all, output_dir, fonts)


if __name__ == "__main__":
    main()
