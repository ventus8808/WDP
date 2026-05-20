from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ─── Type switches ─────────────────────────────────────────────────────────────
Cluster_EQI  = 1
Cluster_NLCD = 1

# ─── Projection switches ───────────────────────────────────────────────────────
Plane    = 1
Mercator = 0

# ─── Font switches ─────────────────────────────────────────────────────────────
Font_Helvetica = 1
Font_Georgia   = 0
# ──────────────────────────────────────────────────────────────────────────────

plt.rcParams["font.size"]       = 12
plt.rcParams["font.weight"]     = "normal"
plt.rcParams["axes.titlesize"]  = 16
plt.rcParams["axes.labelsize"]  = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["legend.fontsize"] = 12

CONTIGUOUS_STATES = [
    "01", "04", "05", "06", "08", "09", "10", "11", "12", "13", "16", "17",
    "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
    "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
    "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56",
]

OBLIQUE_CRS = (
    "+proj=omerc +lat_0=37 +lonc=-96 +alpha=1 +k=0.9996 "
    "+x_0=0 +y_0=0 +gamma=0 +ellps=WGS84 +units=m +no_defs"
)

# EQI clusters: k=3, sorted ascending by EQI quality (A=best environment)
EQI_CLUSTER_COLORS = {
    "A": "#2F7F4F",  # Low EQI burden
    "B": "#97c889",
    "C": "#E6EAB8",  # High EQI burden
    "No Data": "#cccccc",
}
EQI_CLUSTER_LABELS = {
    "A": "Cluster A (Low-burden)",
    "B": "Cluster B (Mixed-burden)",
    "C": "Cluster C (High-burden)",
}

# NLCD clusters: k=4, sorted ascending by intervention score (A=most natural)
NLCD_CLUSTER_COLORS = {
    "A": "#6AAA81",  # Natural
    "B": "#699DCB",  # Water-Sensitive
    "C": "#EFC085",  # Agricultural
    "D": "#E68785",  # Urban
    "No Data": "#cccccc",
}
NLCD_CLUSTER_LABELS = {
    "A": "Cluster A (Natural)",
    "B": "Cluster B (Water-sensitive)",
    "C": "Cluster C (Agricultural)",
    "D": "Cluster D (Urban)",
}


def _load_counties(shapefile_path):
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    return counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()


def _save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


def _plot(counties_contiguous, col, colors, labels, name, fonts, output_dir):
    counties_merged = counties_contiguous.copy()
    counties_merged["color"] = (
        counties_merged[col].map(colors).fillna(colors["No Data"])
    )
    state_boundaries = counties_merged.dissolve(by="STATEFP")

    legend_elements = [
        mpatches.Patch(color=colors[k], label=labels[k]) for k in labels
    ]

    if Plane:
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_merged.plot(
                color=counties_merged["color"],
                linewidth=0, edgecolor="none", ax=ax, zorder=1,
            )
            state_boundaries.boundary.plot(
                ax=ax, color="#000000", linewidth=0.3, alpha=0.6, zorder=2
            )
            ax.set_axis_off()
            ax.legend(handles=legend_elements, bbox_to_anchor=(0.02, 0.02),
                      loc="lower left", frameon=True)
            _save(fig, output_dir / f"Map_{name}_Plane_{font}.png")

    if Mercator:
        counties_proj = counties_merged.to_crs(OBLIQUE_CRS)
        state_boundaries_proj = counties_proj.dissolve(by="STATEFP")
        for font in fonts:
            plt.rcParams["font.family"] = font
            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            counties_proj.plot(
                color=counties_proj["color"],
                linewidth=0.1, edgecolor="black", ax=ax,
            )
            state_boundaries_proj.boundary.plot(
                ax=ax, color="black", linewidth=1.2, alpha=0.9
            )
            ax.set_axis_off()
            ax.legend(handles=legend_elements, bbox_to_anchor=(0.02, 0.02),
                      loc="lower left", frameon=True)
            _save(fig, output_dir / f"Map_{name}_Mercator_{font}.png")


def main():
    project_root = Path(__file__).resolve().parents[2]
    with open(project_root / "config.yaml") as f:
        config = yaml.safe_load(f)

    output_dir = project_root / "Result/Map"
    output_dir.mkdir(parents=True, exist_ok=True)

    counties_contiguous = _load_counties(config["data_sources"]["tiger"]["shapefile"])

    fonts = (["Helvetica"] if Font_Helvetica else []) + (["Georgia"] if Font_Georgia else [])
    if not fonts:
        fonts = ["Helvetica"]

    if Cluster_EQI:
        print("Cluster_EQI...")
        df = pd.read_csv(
            project_root / "Data/Processed/Cluster/Cluster_EQI.csv",
            dtype={"COUNTY_FIPS": str},
        )
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        counties_contiguous = counties_contiguous.merge(
            df[["COUNTY_FIPS", "Cluster_EQI"]], on="COUNTY_FIPS", how="left"
        )
        _plot(counties_contiguous, "Cluster_EQI",
              EQI_CLUSTER_COLORS, EQI_CLUSTER_LABELS, "Cluster_EQI", fonts, output_dir)

    if Cluster_NLCD:
        print("Cluster_NLCD...")
        df = pd.read_csv(
            project_root / "Data/Processed/Cluster/Cluster_NLCD.csv",
            dtype={"COUNTY_FIPS": str},
        )
        df["COUNTY_FIPS"] = df["COUNTY_FIPS"].str.zfill(5)
        # Re-load base counties to avoid leftover EQI column
        counties_nlcd = _load_counties(config["data_sources"]["tiger"]["shapefile"])
        counties_nlcd = counties_nlcd.merge(
            df[["COUNTY_FIPS", "Cluster_NLCD"]], on="COUNTY_FIPS", how="left"
        )
        _plot(counties_nlcd, "Cluster_NLCD",
              NLCD_CLUSTER_COLORS, NLCD_CLUSTER_LABELS, "Cluster_NLCD", fonts, output_dir)


if __name__ == "__main__":
    main()
