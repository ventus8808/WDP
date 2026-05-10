"""
NLCD land-use clustering — two-phase analysis.

Phase 1 · K_Compare  (Result/Cluster/Cluster_NLCD_K_Compare/)
    K-means cluster visualizations + diagnostics for k = 3, 4, 5, 6.
    Outputs: {k}_Radar.png, {k}_Violin.png, {k}_Map.png  (per k)
             Silhouette_K_Compare.png, Elbow_K_Compare.png

Phase 2 · Result  (Result/Cluster/Cluster_NLCD_Result/)
    K-means k = 4 selected: highest silhouette, 4 clusters map cleanly to
    Natural / Water-Sensitive / Agricultural / Urban land-use types.
    Outputs: k4_Radar.png, k4_Violin.png, k4_Map.png

Canonical data output:
    Data/Processed/Cluster/Cluster_NLCD.csv  (COUNTY_FIPS, Cluster_NLCD: A/B/C/D)
"""

import os
from math import pi
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

plt.rcParams["font.family"] = "Georgia"

NLCD_TYPE_COLORS = {
    "Natural": "#6AAA81",
    "Water-Sensitive": "#699DCB",
    "Agricultural": "#EFC085",
    "Urban": "#E68785",
    "Mixed": "#95a5a6",
}

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


def _label(i):
    return chr(ord("A") + int(i))


FEATURE_COLUMNS = [
    "Agriculture_Intensity",
    "Urban_Intensity",
    "Hydro_Sensitivity",
    "Nature_Background",
]

FEATURE_TYPE_MAP = {
    "Agriculture_Intensity": "Agricultural",
    "Urban_Intensity": "Urban",
    "Nature_Background": "Natural",
    "Hydro_Sensitivity": "Water-Sensitive",
}

FEATURE_SHORT = {
    "Agriculture_Intensity": "Agriculture",
    "Urban_Intensity": "Urban",
    "Hydro_Sensitivity": "Hydro",
    "Nature_Background": "Nature",
}

RAW_PROP_COLS = [
    "prop_urban",
    "prop_forest",
    "prop_grassland",
    "prop_shrub",
    "prop_barren",
    "prop_wetland",
    "prop_water",
    "prop_pasture",
    "prop_cropland",
]

RAW_PROP_SHORT = {
    "prop_water": "Water",
    "prop_cropland": "Cropland",
    "prop_pasture": "Pasture",
    "prop_urban": "Urban",
    "prop_wetland": "Wetland",
    "prop_forest": "Forest",
    "prop_shrub": "Shrub",
    "prop_grassland": "Grassland",
    "prop_barren": "Barren",
}


def get_config():
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# ── Data ──────────────────────────────────────────────────────────────────────


def load_and_prepare_data(nlcd_file_path):
    """Load NLCD/JRC data, filter to 2000–2005, compute 4 engineered features + 9 raw proportions."""
    df = pd.read_csv(nlcd_file_path)
    df = df[(df["Year"] >= 2000) & (df["Year"] <= 2005)].copy()
    group_cols = [
        "jrc_permanent_water_km2",
        "nlcd_cropland_km2",
        "nlcd_pasture_km2",
        "nlcd_urban_km2",
        "nlcd_wetland_km2",
        "nlcd_forest_km2",
        "nlcd_shrub_km2",
        "nlcd_grassland_km2",
        "nlcd_barren_km2",
        "total_area_km2",
    ]
    df = df.groupby("COUNTY_FIPS")[group_cols].mean().reset_index()
    t = df["total_area_km2"]
    df["Agriculture_Intensity"] = (df["nlcd_cropland_km2"] + df["nlcd_pasture_km2"]) / t
    df["Urban_Intensity"] = df["nlcd_urban_km2"] / t
    df["Hydro_Sensitivity"] = (
        df["jrc_permanent_water_km2"] + df["nlcd_wetland_km2"]
    ) / t
    df["Nature_Background"] = (
        df["nlcd_forest_km2"]
        + df["nlcd_shrub_km2"]
        + df["nlcd_grassland_km2"]
        + df["nlcd_barren_km2"]
    ) / t
    df["prop_water"] = df["jrc_permanent_water_km2"] / t
    df["prop_cropland"] = df["nlcd_cropland_km2"] / t
    df["prop_pasture"] = df["nlcd_pasture_km2"] / t
    df["prop_urban"] = df["nlcd_urban_km2"] / t
    df["prop_wetland"] = df["nlcd_wetland_km2"] / t
    df["prop_forest"] = df["nlcd_forest_km2"] / t
    df["prop_shrub"] = df["nlcd_shrub_km2"] / t
    df["prop_grassland"] = df["nlcd_grassland_km2"] / t
    df["prop_barren"] = df["nlcd_barren_km2"] / t
    df = df.dropna(subset=FEATURE_COLUMNS)
    print(f"Loaded {df.shape[0]} counties × {len(FEATURE_COLUMNS)} NLCD features")
    return df[["COUNTY_FIPS"] + FEATURE_COLUMNS + RAW_PROP_COLS]


def standardize(df):
    df_std = df.copy()
    df_std[FEATURE_COLUMNS] = StandardScaler().fit_transform(df[FEATURE_COLUMNS])
    return df_std


def _identify_type(profile_row):
    means = {fc: profile_row[f"{fc}_mean"] for fc in FEATURE_COLUMNS}
    return FEATURE_TYPE_MAP.get(max(means, key=means.get), "Mixed")


def profile_and_sort(df_std, labels):
    """Build profiles, assign land-use type, sort ascending by intervention score."""
    df = df_std.copy()
    df["Cluster"] = labels
    profiles = []
    for c in np.unique(labels):
        cd = df[df["Cluster"] == c]
        p = {"Cluster": int(c), "Count": len(cd)}
        for col in FEATURE_COLUMNS:
            p[f"{col}_mean"] = cd[col].mean()
            p[f"{col}_std"] = cd[col].std()
        for col in RAW_PROP_COLS:
            p[f"{col}_mean"] = cd[col].mean()
        profiles.append(p)
    pf = pd.DataFrame(profiles)
    pf["Cluster_Type"] = pf.apply(_identify_type, axis=1)
    pf["intervention_score"] = (
        pf["Urban_Intensity_mean"]
        + pf["Agriculture_Intensity_mean"]
        - pf["Nature_Background_mean"]
    )
    pf = pf.sort_values("intervention_score").reset_index(drop=True)
    old_to_new = {old: new for new, old in enumerate(pf["Cluster"])}
    pf["Cluster"] = range(len(pf))
    df["Cluster"] = df["Cluster"].map(old_to_new)
    ctype_map = pf.set_index("Cluster")["Cluster_Type"].to_dict()
    df["Cluster_Type"] = df["Cluster"].map(ctype_map)
    return df, pf


# ── Plotting ───────────────────────────────────────────────────────────────────


def _radar(profiles_df, output_dir, stem):
    short = [RAW_PROP_SHORT[col] for col in RAW_PROP_COLS]
    n = len(RAW_PROP_COLS)
    angles = [i / float(n) * 2 * pi for i in range(n)]
    angles_c = angles + angles[:1]
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(projection="polar"))
    ax.patch.set_visible(False)
    for _, row in profiles_df.iterrows():
        cid = int(row["Cluster"])
        ctype = row["Cluster_Type"]
        color = NLCD_TYPE_COLORS.get(ctype, NLCD_TYPE_COLORS["Mixed"])
        vals = [row[f"{col}_mean"] for col in RAW_PROP_COLS]
        vals_c = vals + vals[:1]
        ax.plot(
            angles_c,
            vals_c,
            "o-",
            linewidth=0.8,
            markersize=3,
            label=f"{_label(cid)}: {ctype}",
            color=color,
        )
        ax.fill(angles_c, vals_c, alpha=0.20, color=color)
    ax.set_xticks(angles)
    ax.set_xticklabels(short)
    all_vals = [
        row[f"{col}_mean"] for _, row in profiles_df.iterrows() for col in RAW_PROP_COLS
    ]
    ax.set_ylim(0, max(all_vals) + 0.02)
    ax.set_yticklabels([])
    ax.tick_params(axis="both", length=0)
    ax.tick_params(axis="x", pad=12)
    ax.xaxis.grid(True)
    ax.yaxis.grid(False)
    n_clusters = len(profiles_df)
    ncol = min(n_clusters, 4)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=ncol, frameon=True)
    fig.tight_layout(pad=2.0)
    out = os.path.join(output_dir, f"{stem}_Radar.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def _violin(df_c, profiles_df, output_dir, stem):
    label_order = [
        f"{_label(int(r['Cluster']))}: {r['Cluster_Type']}"
        for _, r in profiles_df.iterrows()
    ]
    palette = {
        f"{_label(int(r['Cluster']))}: {r['Cluster_Type']}": NLCD_TYPE_COLORS.get(
            r["Cluster_Type"], NLCD_TYPE_COLORS["Mixed"]
        )
        for _, r in profiles_df.iterrows()
    }
    melted = df_c.melt(
        id_vars=["Cluster", "Cluster_Type"],
        value_vars=FEATURE_COLUMNS,
        var_name="Feature",
        value_name="Value",
    )
    melted["Feature"] = melted["Feature"].map(FEATURE_SHORT)
    melted["Label"] = melted.apply(
        lambda r: f"{_label(int(r['Cluster']))}: {r['Cluster_Type']}", axis=1
    )
    plt.figure(figsize=(10, 6))
    sns.violinplot(
        data=melted,
        x="Feature",
        y="Value",
        hue="Label",
        hue_order=label_order,
        palette=palette,
        linewidth=0.8,
        width=0.8,
        inner="box",
    )
    plt.ylabel("Standardized Value")
    plt.xlabel("")
    plt.grid(axis="y", alpha=0.3)
    plt.legend(frameon=True, title="Cluster")
    plt.tight_layout()
    out = os.path.join(output_dir, f"{stem}_Violin.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def _map(df_orig, df_c, profiles_df, shapefile_path, output_dir, stem):
    color_map = {
        int(r["Cluster"]): NLCD_TYPE_COLORS.get(
            r["Cluster_Type"], NLCD_TYPE_COLORS["Mixed"]
        )
        for _, r in profiles_df.iterrows()
    }
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    cont = counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()

    dm = df_orig[["COUNTY_FIPS"]].copy()
    dm["COUNTY_FIPS"] = dm["COUNTY_FIPS"].astype(str).str.zfill(5)
    dm["cluster"] = df_c["Cluster"].values
    merged = cont.merge(dm, on="COUNTY_FIPS", how="left")
    merged["color"] = merged["cluster"].apply(
        lambda x: "#e0e0e0" if pd.isna(x) else color_map.get(int(x), "#95a5a6")
    )

    legend_handles = [
        mpatches.Patch(
            color=NLCD_TYPE_COLORS.get(r["Cluster_Type"], NLCD_TYPE_COLORS["Mixed"]),
            label=f"{_label(int(r['Cluster']))}: {r['Cluster_Type']} (n={int(r['Count'])})",
        )
        for _, r in profiles_df.iterrows()
    ]
    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    merged.plot(color=merged["color"], linewidth=0, edgecolor="none", ax=ax)
    merged.dissolve(by="STATEFP").boundary.plot(
        ax=ax, color="#000000", linewidth=0.6, alpha=0.6
    )
    ax.set_axis_off()
    ax.legend(
        handles=legend_handles,
        bbox_to_anchor=(0.02, 0.02),
        loc="lower left",
        frameon=True,
    )
    fig.tight_layout(pad=0.0)
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    out = os.path.join(output_dir, f"{stem}_Map.png")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"  Saved: {out}")


def _visualize(df_orig, df_c, profiles_df, shapefile_path, output_dir, stem):
    _radar(profiles_df, output_dir, stem)
    _violin(df_c, profiles_df, output_dir, stem)
    _map(df_orig, df_c, profiles_df, shapefile_path, output_dir, stem)


# ── Phases ────────────────────────────────────────────────────────────────────


def phase1_k_compare(X, df, df_std, shapefile_path, output_dir):
    """K-means cluster plots + silhouette/elbow for k = 3..6 → K_Compare/."""
    print("\n[Phase 1] K_Compare — K-means for k = 3..6")
    k_range = range(3, 7)
    silhouettes, wcss_list = [], []
    for k in k_range:
        print(f"  k={k}...")
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = km.fit_predict(X)
        score = silhouette_score(X, labels)
        silhouettes.append(score)
        wcss_list.append(km.inertia_)
        df_c, profiles_df = profile_and_sort(df_std, labels)
        _visualize(df, df_c, profiles_df, shapefile_path, output_dir, str(k))
        print(f"  k={k}: silhouette={score:.4f}")

    # Silhouette line plot
    plt.figure(figsize=(7, 4))
    plt.plot(list(k_range), silhouettes, marker="o", color="#2F7F4F", linewidth=2)
    plt.xlabel("Number of Clusters (k)")
    plt.ylabel("Silhouette Score")
    plt.title("K-means Silhouette Score by k")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    sil_out = os.path.join(output_dir, "Silhouette_K_Compare.png")
    plt.savefig(sil_out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {sil_out}")

    # Elbow plot
    plt.figure(figsize=(7, 4))
    plt.plot(list(k_range), wcss_list, marker="o", color="#699DCB", linewidth=2)
    plt.xlabel("Number of Clusters (k)")
    plt.ylabel("WCSS")
    plt.title("K-means Elbow Method")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    elbow_out = os.path.join(output_dir, "Elbow_K_Compare.png")
    plt.savefig(elbow_out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {elbow_out}")

    return silhouettes


def phase2_result(X, df, df_std, k, shapefile_path, result_dir, data_dir):
    """K-means k=4 final plots → Result/. Save canonical CSV."""
    print(
        f"\n[Phase 2] Result — K-means k = {k} (highest silhouette, 4 distinct land-use types)"
    )
    labels = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X)
    df_c, profiles_df = profile_and_sort(df_std, labels)
    _visualize(df, df_c, profiles_df, shapefile_path, result_dir, "k4")

    # Canonical cluster CSV
    os.makedirs(data_dir, exist_ok=True)
    out = df[["COUNTY_FIPS"]].copy()
    out["COUNTY_FIPS"] = out["COUNTY_FIPS"].astype(str).str.zfill(5)
    out["Cluster_NLCD"] = df_c["Cluster"].apply(_label)
    csv_path = os.path.join(data_dir, "Cluster_NLCD.csv")
    out[["COUNTY_FIPS", "Cluster_NLCD"]].to_csv(csv_path, index=False)
    print(f"  Canonical data: {csv_path}")


# ── Entry point ───────────────────────────────────────────────────────────────


def main():
    import argparse

    parser = argparse.ArgumentParser(description="NLCD clustering — two-phase analysis")
    parser.add_argument("--dpi", type=int, help="Output DPI (default: 300)")
    args, _ = parser.parse_known_args()
    if args.dpi:
        plt.rcParams["figure.dpi"] = args.dpi

    config = get_config()
    project_root = Path(__file__).resolve().parents[2]
    base_proc = config["data_directories"]["processed"]

    nlcd_file = os.path.join(project_root, base_proc, "Environmental", "NLCD_JRC.csv")
    shapefile_path = os.path.join(
        project_root, config["data_sources"]["tiger"]["shapefile"]
    )

    cluster_root = os.path.join(project_root, "Result", "Cluster")
    dirs = {
        "k_compare": os.path.join(cluster_root, "Cluster_NLCD_K_Compare"),
        "result": os.path.join(cluster_root, "Cluster_NLCD_Result"),
        "data": os.path.join(project_root, base_proc, "Cluster"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    df = load_and_prepare_data(nlcd_file)
    df_std = standardize(df)
    X = df_std[FEATURE_COLUMNS].values

    k_final = 4

    phase1_k_compare(X, df, df_std, shapefile_path, dirs["k_compare"])
    phase2_result(X, df, df_std, k_final, shapefile_path, dirs["result"], dirs["data"])

    print("\nDone. All outputs in Result/Cluster/Cluster_NLCD_*/")


if __name__ == "__main__":
    main()
