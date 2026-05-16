from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ─── Run switches ─────────────────────────────────────────────────────────────
EQI_2000_2005 = 1  # EQI quintile map, period 2000-2005
EQI_2006_2010 = 0  # EQI quintile map, period 2006-2010
Delta = 1  # EQI change category map

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

EQI_COLORS = {
    1: "#2170b5",
    2: "#6baed6",
    3: "#c6dbef",
    4: "#fee0d2",
    5: "#fc9271",
    "No Data": "#cccccc",
}

CHANGE_COLORS = {
    "I": "#7494BC",  # Improved
    "S": "#EEEEEE",  # Stable
    "W": "#DE9F83",  # Worsened
}

PERIOD_LABELS = {"0005": "2000-2005", "0610": "2006-2010"}


def _load_counties(shapefile_path):
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    return counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()


def _save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


def plot_eqi_period(counties_contiguous, period, config, fonts, output_dir):
    print(f"EQI period {PERIOD_LABELS[period]}...")
    eqi_dir = (
        Path(__file__).resolve().parents[2]
        / config["data_sources"]["epa_eqi"]["processed"]
    )
    df_eqi = pd.read_csv(eqi_dir / f"EQI{period}.csv", usecols=["COUNTY_FIPS", "EQI"])
    df_eqi["COUNTY_FIPS"] = df_eqi["COUNTY_FIPS"].astype(str).str.zfill(5)

    counties_merged = counties_contiguous.merge(df_eqi, on="COUNTY_FIPS", how="left")
    counties_merged["color"] = (
        counties_merged["EQI"].map(EQI_COLORS).fillna(EQI_COLORS["No Data"])
    )

    legend_elements = [
        mpatches.Patch(color=EQI_COLORS[i], label=f"EQI {i}") for i in range(1, 6)
    ]
    title = f"EQI Distribution Map ({PERIOD_LABELS[period]})"

    if Plane:
        state_boundaries = counties_merged.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_merged.plot(
                color=counties_merged["color"], linewidth=0.1, edgecolor="black", ax=ax
            )
            state_boundaries.boundary.plot(
                ax=ax, color="black", linewidth=0.3, alpha=0.8
            )
            ax.set_axis_off()
            ax.legend(
                handles=legend_elements,
                bbox_to_anchor=(0.02, 0.02),
                loc="lower left",
                frameon=True,
            )
            plt.suptitle(title, y=0.82)
            _save(fig, output_dir / f"EQI_{period}_Plane_{font}.png")

    if Mercator:
        counties_proj = counties_merged.to_crs(OBLIQUE_CRS)
        state_boundaries = counties_proj.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_proj.plot(
                color=counties_proj["color"], linewidth=0.1, edgecolor="black", ax=ax
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
            plt.suptitle(title, y=0.82)
            _save(fig, output_dir / f"EQI_{period}_Mercator_{font}.png")


def plot_eqi_delta(counties_contiguous, project_root, fonts, output_dir):
    print("EQI Delta map...")
    data_file = project_root / "Data/Processed/df_Delta.csv"
    if not data_file.exists():
        print(f"Error: {data_file} not found")
        return

    df = pd.read_csv(data_file, dtype={"COUNTY_FIPS": str})
    df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
    if "EQI_Change_Category" not in df.columns:
        print("Error: EQI_Change_Category column not found")
        return

    df_uniq = df.drop_duplicates("COUNTY_FIPS")[["COUNTY_FIPS", "EQI_Change_Category"]]
    counties_merged = counties_contiguous.merge(df_uniq, on="COUNTY_FIPS", how="left")

    def get_color(cat):
        if pd.isna(cat):
            return CHANGE_COLORS["S"]
        return CHANGE_COLORS.get(str(cat).strip().upper()[0], CHANGE_COLORS["S"])

    counties_merged["color"] = counties_merged["EQI_Change_Category"].apply(get_color)

    legend_elements = [
        mpatches.Patch(color=CHANGE_COLORS["I"], label="Improved"),
        mpatches.Patch(color=CHANGE_COLORS["S"], label="Stable"),
        mpatches.Patch(color=CHANGE_COLORS["W"], label="Worsened"),
    ]

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
                title="EQI Change Status",
            )
            _save(fig, output_dir / f"Map_EQI_Delta_Plane_{font}.png")

    if Mercator:
        counties_proj = counties_merged.to_crs(OBLIQUE_CRS)
        state_boundaries = counties_proj.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_proj.plot(
                color=counties_proj["color"], linewidth=0.1, edgecolor="black", ax=ax
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
                title="EQI Change Status",
            )
            plt.suptitle("Distribution of EQI Change Categories", y=0.82)
            _save(fig, output_dir / f"Map_EQI_Delta_Mercator_{font}.png")


def main():
    project_root = Path(__file__).resolve().parents[2]
    with open(project_root / "config.yaml") as f:
        config = yaml.safe_load(f)

    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    counties_contiguous = _load_counties(config["data_sources"]["tiger"]["shapefile"])

    fonts = (["Helvetica"] if Font_Helvetica else []) + (
        ["Georgia"] if Font_Georgia else []
    )
    if not fonts:
        fonts = ["Helvetica"]

    if EQI_2000_2005:
        plot_eqi_period(counties_contiguous, "0005", config, fonts, output_dir)
    if EQI_2006_2010:
        plot_eqi_period(counties_contiguous, "0610", config, fonts, output_dir)
    if Delta:
        plot_eqi_delta(counties_contiguous, project_root, fonts, output_dir)


if __name__ == "__main__":
    main()
