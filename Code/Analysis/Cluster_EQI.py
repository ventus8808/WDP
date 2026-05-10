"""
EQI clustering — three-phase analysis.

Phase 1 · K_Compare  (Result/Cluster/Cluster_EQI_K_Compare/)
    K-means clustering visualizations for k = 3, 4, 5, 6 (Georgia only).
    Outputs: KMeans_k{3-6}_G.png

Phase 2 · Method_Compare  (Result/Cluster/Cluster_EQI_Method_Compare/)
    Three methods at k = 3 (Georgia only) + silhouette scores k = 3..10.
    Outputs: {Method}_k3_G.png, Method_Compare.csv

Phase 3 · Result  (Result/Cluster/Cluster_EQI_Result/)
    K-means k = 3 selected: best overall silhouette, consistent interpretation.
    Outputs: KMeans_k3_{G,H}.png, Silhouette_{G,H}.png
             (scatter plot — all 3 methods × k = 3..10)

Canonical data output:
    Data/Processed/Cluster/Cluster_EQI.csv  (COUNTY_FIPS, Cluster_EQI: A/B/C)
"""

import os
import sys
from contextlib import contextmanager
from math import pi
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.image import imread as _imread
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from sklearn.cluster import AgglomerativeClustering, Birch, KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

sys.path.append(os.path.dirname(__file__))
from Cluster_Plot_Function import (
    FONT_FAMILIES,
    create_combined_with_labels,
    set_default_dpi,
    set_font_sizes,
)

CONTIGUOUS_STATES = [
    "01", "04", "05", "06", "08", "09", "10", "11", "12", "13", "16", "17",
    "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
    "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41",
    "42", "44", "45", "46", "47", "48", "49", "50", "51", "53", "54", "55", "56",
]

# Extended palette for up to 6 EQI clusters
EQI_COLORS = ["#2F7F4F", "#97c889", "#E6EAB8", "#699DCB", "#E68785", "#EFC085"]

# Colors for the 5-method silhouette scatter plot
METHOD_COLORS = {
    "KMeans":        "#2F7F4F",
    "Agglomerative": "#699DCB",
    "Birch":         "#EFC085",
}


METHODS = {
    "KMeans":        lambda X, k: KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(X),
    "Agglomerative": lambda X, k: AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(X),
    "Birch":         lambda X, k: Birch(n_clusters=k).fit_predict(X),
}


def get_config():
    config_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "config.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


@contextmanager
def _font(family):
    prev = plt.rcParams.get("font.family")
    plt.rcParams["font.family"] = family
    try:
        yield
    finally:
        if prev is not None:
            plt.rcParams["font.family"] = prev


def _fpath(base, family):
    if family == "Georgia":
        return base
    root, ext = os.path.splitext(base)
    return f"{root}_Helvetica{ext}"


def _suffix(family):
    return "G" if family == "Georgia" else "H"


def _rm(path):
    if os.path.exists(path):
        os.remove(path)


def _label(i):
    return chr(ord("A") + int(i))


def _colors(k):
    return EQI_COLORS[:k]


# ── Data ──────────────────────────────────────────────────────────────────────

def load_and_prepare_data(standard_file_path):
    df = pd.read_csv(standard_file_path)
    eqi_columns = ["EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]
    df_clean = df[["COUNTY_FIPS"] + eqi_columns].dropna().copy()
    print(f"Loaded {df_clean.shape[0]} counties × {len(eqi_columns)} EQI dimensions")
    return df_clean, eqi_columns


def standardize(df, eqi_columns):
    df_std = df.copy()
    df_std[eqi_columns] = StandardScaler().fit_transform(df[eqi_columns])
    return df_std


def profile_and_sort(df_std, eqi_columns, labels):
    """Build cluster profiles; sort ascending by EQI sum (better env → letter A)."""
    df = df_std.copy()
    df["Cluster"] = labels
    profiles = []
    for c in np.unique(labels):
        cd = df[df["Cluster"] == c]
        p = {"Cluster": int(c), "Count": len(cd)}
        for col in eqi_columns:
            p[f"{col}_mean"] = cd[col].mean()
            p[f"{col}_std"] = cd[col].std()
        profiles.append(p)
    pf = pd.DataFrame(profiles)
    pf["EQI_SUM"] = pf[[f"{col}_mean" for col in eqi_columns]].sum(axis=1)
    pf = pf.sort_values("EQI_SUM").reset_index(drop=True)
    old_to_new = {old: new for new, old in enumerate(pf["Cluster"])}
    pf["Cluster"] = range(len(pf))
    df["Cluster"] = df["Cluster"].map(old_to_new)
    return df, pf


# ── Plotting helpers ───────────────────────────────────────────────────────────

def _save_radar(profiles_df, eqi_columns, output_dir, stem, families=None):
    """Radar chart; families=None → Georgia only."""
    if families is None:
        families = ["Georgia"]
    k = len(profiles_df)
    cols = _colors(k)
    categories = [col.replace("EQI_", "") for col in eqi_columns]
    angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
    angles_c = angles + angles[:1]
    base_path = os.path.join(output_dir, f"{stem}_Radar.png")
    for family in families:
        with _font(family):
            fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(projection="polar"))
            ax.patch.set_visible(False)
            for _, row in profiles_df.iterrows():
                cid = int(row["Cluster"])
                vals = [row[f"{col}_mean"] for col in eqi_columns]
                vals_c = vals + vals[:1]
                ax.plot(angles_c, vals_c, "o-", linewidth=2,
                        label=f"Cluster {_label(cid)}", color=cols[cid])
                ax.fill(angles_c, vals_c, alpha=0.25, color=cols[cid])
            ax.set_xticks(angles)
            ax.set_xticklabels(categories)
            ax.set_ylim(-2, 1)
            ax.set_yticklabels([])
            ax.tick_params(axis="both", length=0)
            ax.tick_params(axis="x", pad=12)
            ax.grid(True)
            ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), frameon=True)
            fig.tight_layout(pad=1.0)
            plt.savefig(_fpath(base_path, family), dpi=300)
            plt.close()
    return base_path


def _save_box(df_clustered, eqi_columns, output_dir, stem, families=None):
    """Box plot; families=None → Georgia only."""
    if families is None:
        families = ["Georgia"]
    k = df_clustered["Cluster"].nunique()
    cols = _colors(k)
    melted = df_clustered.melt(
        id_vars=["Cluster"], value_vars=eqi_columns,
        var_name="EQI_Dimension", value_name="EQI_Value",
    )
    melted["EQI_Dimension"] = melted["EQI_Dimension"].str.replace("EQI_", "")
    melted["Cluster_Label"] = melted["Cluster"].apply(lambda x: f"Cluster {_label(x)}")
    hue_order = [f"Cluster {_label(i)}" for i in range(k)]
    palette = {f"Cluster {_label(i)}": cols[i] for i in range(k)}
    base_path = os.path.join(output_dir, f"{stem}_Box.png")
    for family in families:
        with _font(family):
            plt.figure(figsize=(9, 6))
            sns.boxplot(
                data=melted, x="EQI_Dimension", y="EQI_Value",
                hue="Cluster_Label", hue_order=hue_order, palette=palette,
                showfliers=False, linewidth=0.8, width=0.6,
            )
            plt.ylabel("EQI Value (Standardized)")
            plt.xlabel("")
            plt.grid(axis="x", alpha=0.3)
            plt.legend(frameon=True)
            plt.tight_layout()
            plt.savefig(_fpath(base_path, family), dpi=300)
            plt.close()
    return base_path


def _save_map(df_map, k_clusters, shapefile_path, output_dir, stem, families=None):
    """County choropleth; families=None → Georgia only."""
    if families is None:
        families = ["Georgia"]
    cols = _colors(k_clusters)
    color_map = {i: cols[i] for i in range(k_clusters)}
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]
    cont = counties[counties["STATEFP"].isin(CONTIGUOUS_STATES)].copy()
    merged = cont.merge(df_map[["COUNTY_FIPS", "cluster"]], on="COUNTY_FIPS", how="left")
    merged["color"] = merged["cluster"].apply(
        lambda x: "#d3d3d3" if pd.isna(x) else color_map.get(int(x), "#d3d3d3")
    )
    base_path = os.path.join(output_dir, f"{stem}_Map.png")
    for family in families:
        with _font(family):
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            merged.plot(color=merged["color"], linewidth=0, edgecolor="none", ax=ax)
            merged.dissolve(by="STATEFP").boundary.plot(
                ax=ax, color="#000000", linewidth=0.6, alpha=0.6
            )
            ax.set_axis_off()
            legend_handles = [
                mpatches.Patch(color=cols[i], label=f"Cluster {_label(i)}")
                for i in range(k_clusters)
            ]
            ax.legend(handles=legend_handles, bbox_to_anchor=(0.02, 0.02),
                      loc="lower left", frameon=True)
            fig.tight_layout(pad=0.0)
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            plt.savefig(_fpath(base_path, family), dpi=300)
            plt.close()
    return base_path


def _make_map_df(df, df_c):
    dm = df[["COUNTY_FIPS"]].copy()
    dm["COUNTY_FIPS"] = dm["COUNTY_FIPS"].astype(str).str.zfill(5)
    dm["cluster"] = df_c["Cluster"].values
    return dm


def _combine_georgia(box_path, radar_path, map_path, output_dir, final_stem):
    """Compose box+radar+map into a single labeled figure (Georgia, single-font)."""
    box_img   = _imread(box_path)
    radar_img = _imread(radar_path)
    map_img   = _imread(map_path)

    fig = plt.figure(figsize=(15, 15))
    top_h, bottom_h = 0.45, 0.55
    left_frac, right_frac = 0.6, 0.4

    axA = fig.add_axes([0,          bottom_h, left_frac,  top_h])
    axB = fig.add_axes([left_frac,  bottom_h, right_frac, top_h])
    axC = fig.add_axes([0,          0,        1.0,        bottom_h])

    for ax, img, lbl in [(axA, box_img, "(A)"), (axB, radar_img, "(B)"), (axC, map_img, "(C)")]:
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.02, 0.98, lbl, transform=ax.transAxes,
                fontsize=18, fontweight="bold", color="black", ha="left", va="top")

    out = os.path.join(output_dir, f"{final_stem}_G.png")
    plt.savefig(out, dpi=200)
    plt.close()
    print(f"  Saved: {out}")


def _combined_dual(box_path, radar_path, map_path, output_dir, method_name, k, final_stem):
    """Dual-font combined figure; cleans up intermediate subplot files."""
    create_combined_with_labels(
        box_path=box_path, radar_path=radar_path, map_path=map_path,
        output_dir=output_dir, method_name=method_name, k=k,
    )
    for family in FONT_FAMILIES:
        helv = "" if family == "Georgia" else "_Helvetica"
        src = os.path.join(output_dir, f"{method_name}_{k}_Combined_Labeled{helv}.png")
        dst = os.path.join(output_dir, f"{final_stem}_{_suffix(family)}.png")
        if os.path.exists(src):
            os.replace(src, dst)
            print(f"  Saved: {dst}")
    for p in [radar_path, box_path, map_path]:
        for family in FONT_FAMILIES:
            _rm(_fpath(p, family))


def _save_method_silhouette_scatter(method_silhouettes, k_range, output_dir):
    """Multi-line scatter plot: 5 methods × k=3..10."""
    for family in FONT_FAMILIES:
        with _font(family):
            fig, ax = plt.subplots(figsize=(9, 6))
            for method, scores in method_silhouettes.items():
                color = METHOD_COLORS.get(method, "#999999")
                ax.plot(list(k_range), scores, marker="o", label=method,
                        color=color, linewidth=2, markersize=6)
            ax.set_xlabel("Number of Clusters (k)")
            ax.set_ylabel("Silhouette Score")
            ax.set_title("Silhouette Score — All Methods, k = 3..10")
            ax.legend(frameon=True)
            ax.grid(True, alpha=0.3)
            ax.set_xticks(list(k_range))
            plt.tight_layout()
            out = os.path.join(output_dir, f"Silhouette_{_suffix(family)}.png")
            plt.savefig(out, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"  Saved: {out}")


def _run_and_visualize(method_name, cluster_fn, X, df, df_std, eqi_columns, k,
                       shapefile_path, output_dir, dual_font=False):
    """Cluster, visualize, clean intermediates. Returns (silhouette, df_clustered)."""
    labels = cluster_fn(X, k)
    score = silhouette_score(X, labels)
    df_c, profiles_df = profile_and_sort(df_std, eqi_columns, labels)
    stem = f"{method_name}_k{k}"
    families = FONT_FAMILIES if dual_font else ["Georgia"]

    radar_path = _save_radar(profiles_df, eqi_columns, output_dir, stem, families)
    box_path   = _save_box(df_c, eqi_columns, output_dir, stem, families)
    map_path   = _save_map(_make_map_df(df, df_c), k, shapefile_path, output_dir, stem, families)

    if dual_font:
        _combined_dual(box_path, radar_path, map_path, output_dir, method_name, k, stem)
    else:
        _combine_georgia(box_path, radar_path, map_path, output_dir, stem)
        _rm(box_path)
        _rm(radar_path)
        _rm(map_path)

    return score, df_c


# ── Phases ────────────────────────────────────────────────────────────────────

def phase1_k_compare(X, df, df_std, eqi_columns, shapefile_path, output_dir):
    """K-means cluster visualizations for k = 3..6 (Georgia only) → K_Compare/."""
    print("\n[Phase 1] K_Compare — K-means for k = 3..6 (Georgia)")
    for k in range(3, 7):
        print(f"  k={k}...")
        score, _ = _run_and_visualize(
            "KMeans", METHODS["KMeans"], X, df, df_std, eqi_columns,
            k, shapefile_path, output_dir, dual_font=False,
        )
        print(f"  k={k}: silhouette={score:.4f}")


def phase2_method_compare(X, df, df_std, eqi_columns, k_viz, k_sil_range, shapefile_path, output_dir):
    """5 methods: viz at k=k_viz (Georgia only) + silhouette k_sil_range. Returns silhouettes dict."""
    print(f"\n[Phase 2] Method_Compare — 5 methods: viz k={k_viz}, silhouette k={min(k_sil_range)}..{max(k_sil_range)}")
    method_silhouettes = {}
    rows = []
    k_list = list(k_sil_range)

    for method_name, cluster_fn in METHODS.items():
        print(f"  {method_name}...")
        silhouettes = []
        for kk in k_sil_range:
            labels = cluster_fn(X, kk)
            score = silhouette_score(X, labels)
            silhouettes.append(score)
            rows.append({"Method": method_name, "k": kk, "Silhouette_Score": score})
        method_silhouettes[method_name] = silhouettes
        # Visualization at the chosen k (single-font; re-runs clustering — fast for 3k counties)
        score_at_viz = silhouettes[k_list.index(k_viz)]
        _run_and_visualize(method_name, cluster_fn, X, df, df_std, eqi_columns,
                           k_viz, shapefile_path, output_dir, dual_font=False)
        print(f"  {method_name}: k={k_viz} silhouette={score_at_viz:.4f}")

    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "Method_Compare.csv"), index=False)
    return method_silhouettes


def phase3_result(X, df, df_std, eqi_columns, k, shapefile_path,
                  result_dir, data_dir, method_silhouettes, k_sil_range):
    """K-means k=3 final viz (dual-font) + method silhouette scatter → Result/."""
    print(f"\n[Phase 3] Result — K-means k={k} selected (best overall silhouette)")

    _, df_c = _run_and_visualize(
        "KMeans", METHODS["KMeans"], X, df, df_std, eqi_columns,
        k, shapefile_path, result_dir, dual_font=True,
    )
    _save_method_silhouette_scatter(method_silhouettes, k_sil_range, result_dir)

    os.makedirs(data_dir, exist_ok=True)
    out = df[["COUNTY_FIPS"]].copy()
    out["COUNTY_FIPS"] = out["COUNTY_FIPS"].astype(str).str.zfill(5)
    out["Cluster_EQI"] = df_c["Cluster"].apply(_label)
    csv_path = os.path.join(data_dir, "Cluster_EQI.csv")
    out[["COUNTY_FIPS", "Cluster_EQI"]].to_csv(csv_path, index=False)
    print(f"  Canonical data: {csv_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(description="EQI clustering — three-phase analysis")
    parser.add_argument("--dpi", type=int)
    parser.add_argument("--font-scale", type=float)
    args, _ = parser.parse_known_args()
    if args.dpi:
        set_default_dpi(args.dpi)
    if args.font_scale:
        set_font_sizes(scale=args.font_scale)

    config = get_config()
    project_root = Path(__file__).resolve().parents[2]
    base_proc = config["data_directories"]["processed"]

    standard_file  = os.path.join(project_root, base_proc, "EQI", "EQI0005_Standard.csv")
    shapefile_path = os.path.join(project_root, config["data_sources"]["tiger"]["shapefile"])

    cluster_root = os.path.join(project_root, "Result", "Cluster")
    dirs = {
        "k_compare": os.path.join(cluster_root, "Cluster_EQI_K_Compare"),
        "method":    os.path.join(cluster_root, "Cluster_EQI_Method_Compare"),
        "result":    os.path.join(cluster_root, "Cluster_EQI_Result"),
        "data":      os.path.join(project_root, base_proc, "Cluster"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)

    df, eqi_columns = load_and_prepare_data(standard_file)
    df_std = standardize(df, eqi_columns)
    X = df_std[eqi_columns].values

    k_final    = 3
    k_sil_range = range(3, 11)

    phase1_k_compare(X, df, df_std, eqi_columns, shapefile_path, dirs["k_compare"])
    method_silhouettes = phase2_method_compare(
        X, df, df_std, eqi_columns, k_final, k_sil_range, shapefile_path, dirs["method"]
    )
    phase3_result(
        X, df, df_std, eqi_columns, k_final, shapefile_path,
        dirs["result"], dirs["data"], method_silhouettes, k_sil_range,
    )

    print("\nDone. All outputs in Result/Cluster/Cluster_EQI_*/")


if __name__ == "__main__":
    main()
