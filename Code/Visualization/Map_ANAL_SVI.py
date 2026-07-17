from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ─── Run switches ─────────────────────────────────────────────────────────────
ANAL_CONTINUOUS = 1
SVI_CONTINUOUS = 1

Font_Helvetica = 1
Font_Georgia = 0

Plane = 1
Mercator = 0
# ──────────────────────────────────────────────────────────────────────────────

plt.rcParams["font.size"] = 12
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12

CONTIGUOUS_STATES = [
    "01", "04", "05", "06", "08", "09", "10", "11", "12", "13",
    "16", "17", "18", "19", "20", "21", "22", "23", "24", "25",
    "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
    "36", "37", "38", "39", "40", "41", "42", "44", "45", "46",
    "47", "48", "49", "50", "51", "53", "54", "55", "56",
]

OBLIQUE_CRS = (
    "+proj=omerc +lat_0=37 +lonc=-96 +alpha=1 +k=0.9996 "
    "+x_0=0 +y_0=0 +gamma=0 +ellps=WGS84 +units=m +no_defs"
)

ANAL_METRICS = {
    "popw_mean_rad": "Population-weighted mean radiance",
    "mean_rad": "Mean radiance",
    "sol": "Sum of lights",
    "lit_area_km2": "Lit area (km2)",
}

# Exposure windows used by the ANAL-AAMR lag design.
ANAL_PERIODS = [
    "2001-2005",
    "2006-2010",
    "2011-2015",
    "2016-2019",
    "2011-2014",
    "2006-2009",
]

# Custom low-to-high continuous palette:
# blue low values -> pale transition -> red high values -> purple highest values.
CONTINUOUS_COLORS = [
    "#006BAD",
    "#668FCA",
    "#5FA3CB",
    "#66A8CD",
    "#A5CDE2",
    "#D6DFEF",
    "#FCE8E6",
    "#FFC6BC",
    "#F8B9B8",
    "#C9CEFE",
]
CONTINUOUS_POSITIONS = [
    0.00,
    0.07,
    0.13,
    0.20,
    0.30,
    0.42,
    0.56,
    0.70,
    0.84,
    1.00,
]
CONTINUOUS_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "anal_blue_red_purple",
    list(zip(CONTINUOUS_POSITIONS, CONTINUOUS_COLORS)),
    N=256,
)
NO_DATA_COLOR = "#eeeeee"
COUNTY_EDGE_COLOR = "none"
COUNTY_LINEWIDTH = 0
STATE_EDGE_COLOR = "#404040"


def _load_counties(shapefile_path):
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    return counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()


def _save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


def _load_anal_svi(project_root):
    data_file = project_root / "Data/Processed/df_ANAL.csv"
    if not data_file.exists():
        raise FileNotFoundError(f"{data_file} not found")

    usecols = [
        "COUNTY_FIPS",
        "ANAL_Period",
        "popw_mean_rad",
        "mean_rad",
        "sol",
        "lit_area_km2",
        "SVI_Cont",
    ]
    df = pd.read_csv(data_file, usecols=usecols, dtype={"COUNTY_FIPS": str})
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)

    # df_ANAL has one row per county/outcome/lag pair. These fields are county
    # exposure attributes, so one row per county-period is sufficient for maps.
    return df.drop_duplicates(["COUNTY_FIPS", "ANAL_Period"]).reset_index(drop=True)


def _continuous_limits(series):
    vals = pd.to_numeric(series, errors="coerce").dropna()
    if vals.empty:
        return 0, 1
    lo = vals.quantile(0.01)
    hi = vals.quantile(0.99)
    if lo == hi:
        lo = vals.min()
        hi = vals.max()
    if lo == hi:
        hi = lo + 1
    return lo, hi


def _plot_continuous_map(
    counties,
    value_col,
    value_label,
    title,
    output_stem,
    fonts,
    output_dir,
):
    counties = counties.copy()
    counties[value_col] = pd.to_numeric(counties[value_col], errors="coerce")
    vmin, vmax = _continuous_limits(counties[value_col])

    cmap = CONTINUOUS_CMAP.copy()
    cmap.set_bad(NO_DATA_COLOR)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    if Plane:
        state_boundaries = counties.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties.plot(
                column=value_col,
                cmap=cmap,
                norm=norm,
                linewidth=COUNTY_LINEWIDTH,
                edgecolor=COUNTY_EDGE_COLOR,
                missing_kwds={"color": NO_DATA_COLOR},
                ax=ax,
            )
            state_boundaries.boundary.plot(
                ax=ax, color=STATE_EDGE_COLOR, linewidth=0.55, alpha=0.75
            )
            ax.set_axis_off()
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, shrink=0.65)
            cbar.set_label(value_label)
            plt.suptitle(title, y=0.82)
            _save(fig, output_dir / f"{output_stem}_Plane_{font}.png")

    if Mercator:
        counties_proj = counties.to_crs(OBLIQUE_CRS)
        state_boundaries = counties_proj.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_proj.plot(
                column=value_col,
                cmap=cmap,
                norm=norm,
                linewidth=COUNTY_LINEWIDTH,
                edgecolor=COUNTY_EDGE_COLOR,
                missing_kwds={"color": NO_DATA_COLOR},
                ax=ax,
            )
            state_boundaries.boundary.plot(
                ax=ax, color=STATE_EDGE_COLOR, linewidth=1.1, alpha=0.9
            )
            ax.set_axis_off()
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, shrink=0.65)
            cbar.set_label(value_label)
            plt.suptitle(title, y=0.82)
            _save(fig, output_dir / f"{output_stem}_Mercator_{font}.png")


def plot_anal_continuous(counties_contiguous, df, fonts, output_dir):
    for period in ANAL_PERIODS:
        period_df = df[df["ANAL_Period"] == period].copy()
        if period_df.empty:
            print(f"Warning: no ANAL rows for {period}")
            continue
        for metric, label in ANAL_METRICS.items():
            print(f"ANAL continuous map: {metric}, {period}...")
            plot_df = period_df[["COUNTY_FIPS", metric]]
            counties_merged = counties_contiguous.merge(plot_df, on="COUNTY_FIPS", how="left")
            _plot_continuous_map(
                counties_merged,
                metric,
                label,
                f"{label} ({period})",
                f"Map_ANAL_{metric}_{period}",
                fonts,
                output_dir,
            )


def plot_svi_continuous(counties_contiguous, df, fonts, output_dir):
    print("SVI continuous map...")
    svi = df.drop_duplicates("COUNTY_FIPS")[["COUNTY_FIPS", "SVI_Cont"]]
    counties_merged = counties_contiguous.merge(svi, on="COUNTY_FIPS", how="left")
    _plot_continuous_map(
        counties_merged,
        "SVI_Cont",
        "Continuous SVI",
        "Continuous Social Vulnerability Index",
        "Map_SVI_Continuous",
        fonts,
        output_dir,
    )


def main():
    project_root = Path(__file__).resolve().parents[2]
    with open(project_root / "config.yaml") as f:
        config = yaml.safe_load(f)

    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    counties_contiguous = _load_counties(config["data_sources"]["tiger"]["shapefile"])
    df = _load_anal_svi(project_root)

    fonts = (["Helvetica"] if Font_Helvetica else []) + (
        ["Georgia"] if Font_Georgia else []
    )
    if not fonts:
        fonts = ["Helvetica"]

    if ANAL_CONTINUOUS:
        plot_anal_continuous(counties_contiguous, df, fonts, output_dir)
    if SVI_CONTINUOUS:
        plot_svi_continuous(counties_contiguous, df, fonts, output_dir)


if __name__ == "__main__":
    main()
