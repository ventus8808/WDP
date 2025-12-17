import os
from math import pi

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib.image import imread

# Font size configuration with parameterization utilities
FONT_SIZES = {
    "base": 14,  # general text
    "title": 18,  # titles
    "axes": 14,  # axis labels
    "ticks": 14,  # tick labels
    "legend": 14,  # legend text
}


def _apply_font_sizes():
    plt.rcParams["font.family"] = "Georgia"
    plt.rcParams["font.weight"] = "normal"
    plt.rcParams["font.size"] = FONT_SIZES["base"]
    plt.rcParams["axes.titlesize"] = FONT_SIZES["title"]
    plt.rcParams["axes.labelsize"] = FONT_SIZES["axes"]
    plt.rcParams["xtick.labelsize"] = FONT_SIZES["ticks"]
    plt.rcParams["ytick.labelsize"] = FONT_SIZES["ticks"]
    plt.rcParams["legend.fontsize"] = FONT_SIZES["legend"]
    plt.rcParams["figure.titlesize"] = FONT_SIZES["title"]


def set_font_sizes(
    base=None, title=None, axes=None, ticks=None, legend=None, scale=None
):
    """
    Set global font sizes. Pass absolute sizes or a scale factor.
    Example:
      set_font_sizes(scale=1.2)  # enlarge all text by 20%
      set_font_sizes(base=16, title=20)  # set specific sizes
    """
    if scale is not None:
        factor = float(scale)
        for key in ("base", "title", "axes", "ticks", "legend"):
            FONT_SIZES[key] = max(1, int(round(FONT_SIZES[key] * factor)))
    if base is not None:
        FONT_SIZES["base"] = int(base)
    if title is not None:
        FONT_SIZES["title"] = int(title)
    if axes is not None:
        FONT_SIZES["axes"] = int(axes)
    if ticks is not None:
        FONT_SIZES["ticks"] = int(ticks)
    if legend is not None:
        FONT_SIZES["legend"] = int(legend)
    _apply_font_sizes()


def scale(factor: float):
    """Convenience: multiplicatively scale all font sizes."""
    set_font_sizes(scale=factor)


# Initialize rcParams with defaults
_apply_font_sizes()

# Global DPI and fixed-pixel targets (modifiable)
DEFAULT_DPI = 200  # change this to control global output DPI

# Fixed pixel targets
PIXELS_RADAR = (1000, 1000)
PIXELS_BOX = (2000, 1000)
PIXELS_MAP = (3000, 2000)
PIXELS_COMBINED = (3000, 3000)


def set_default_dpi(dpi: int):
    """Set global default DPI for all exported figures."""
    global DEFAULT_DPI
    DEFAULT_DPI = int(dpi)


def _figsize_from_pixels(w: int, h: int, dpi: int):
    """Convert pixel dimensions to matplotlib figsize in inches for a given DPI."""
    return (w / dpi, h / dpi)


# Centralized font families configuration
# Edit this list to add/remove font families; plotting functions will iterate it.
FONT_FAMILIES = ["Georgia", "Helvetica"]


def _font_variant_path(path: str, family: str):
    """Return output path variant for a given font family."""
    if family == "Georgia":
        return path
    root, ext = os.path.splitext(path)
    return f"{root}_Helvetica{ext}"


# Editable cluster palette (直接修改这些 hex 值来改变三类簇的颜色)
# Keys 0/1/2 对应 Cluster 0/1/2，'No Data' 用于缺失值
CLUSTER_PALETTE = {
    0: "#2F7F4F",  # Cluster 0 - 深绿（自然核心）
    1: "#97c889",  # Cluster 1 - 中绿（过渡）
    2: "#E6EAB8",  # Cluster 2 - 浅黄绿（人为）
    "No Data": "#d3d3d3",
}


def get_cluster_colors_list():
    """返回按索引排序的颜色列表（长度 3）。"""
    return [CLUSTER_PALETTE.get(i, "#000000") for i in range(3)]


# 供原有代码直接使用的列表；如果更改了 CLUSTER_PALETTE，
# 请调用 `refresh_cluster_colors()` 来同步此列表。
CLUSTER_COLORS = get_cluster_colors_list()


def set_cluster_palette_from_list(hex_list):
    """用长度为 3 的列表设置 Cluster 0/1/2 的颜色，然后刷新全局列表。

    例: set_cluster_palette_from_list(['#ff0000', '#00ff00', '#0000ff'])
    """
    if not isinstance(hex_list, (list, tuple)) or len(hex_list) != 3:
        raise ValueError("Provide a list/tuple of three hex color strings.")
    for i, h in enumerate(hex_list):
        CLUSTER_PALETTE[i] = h
    refresh_cluster_colors()


def set_cluster_palette_from_dict(palette_dict):
    """用 dict 更新 CLUSTER_PALETTE，支持 {0:hex,1:hex,2:hex,'No Data':hex}。"""
    CLUSTER_PALETTE.update(palette_dict)
    refresh_cluster_colors()


def refresh_cluster_colors():
    global CLUSTER_COLORS
    CLUSTER_COLORS = get_cluster_colors_list()


def get_color(cluster, color_map):
    """Get color for cluster"""
    if pd.isna(cluster):
        return color_map["No Data"]
    return color_map.get(int(cluster), color_map["No Data"])


def create_radar_chart(profiles_df, eqi_columns, output_dir, method_name, k):
    """Create radar chart for cluster profiles"""
    print("Creating radar chart...")

    # Prepare data
    categories = [col.replace("EQI_", "") for col in eqi_columns]
    # Base (Georgia) output path; non-Georgia variants will be suffixed automatically
    georgia_path = os.path.join(output_dir, f"{method_name}_{k}_Radar_Chart.png")

    # Render a radar chart for each configured font family
    for family in FONT_FAMILIES:
        prev_family = plt.rcParams.get("font.family", None)
        try:
            plt.rcParams["font.family"] = family

            fig, ax = plt.subplots(
                figsize=(5, 5),
                subplot_kw=dict(projection="polar"),
            )

            angles = [
                n / float(len(categories)) * 2 * pi for n in range(len(categories))
            ]
            angles += angles[:1]  # Close the loop

            # Remove circular background and keep a clean polar canvas
            ax.patch.set_visible(False)

            colors = CLUSTER_COLORS
            for _, row in profiles_df.iterrows():
                values = [row[f"{col}_mean"] for col in eqi_columns]
                values += values[:1]  # Close the loop
                cluster_id = int(row["Cluster"])

                ax.plot(
                    angles,
                    values,
                    "o-",
                    linewidth=2,
                    label=f"Cluster {cluster_id}",
                    color=colors[cluster_id],
                )
                ax.fill(angles, values, alpha=0.25, color=colors[cluster_id])

            # Place category labels on the outer rim
            ax.set_xticks(angles[:-1])
            ax.set_xticklabels(categories)
            # Keep radial gridlines but remove radial tick labels; hide tick marks
            ax.set_ylim(-2, 1)
            ax.grid(True)
            ax.set_yticklabels([])
            ax.tick_params(axis="both", which="both", length=0)

            # Emphasize spoke gridlines; hide concentric rings
            default_grid_color = plt.rcParams.get("grid.color", "0.8")
            for gl in ax.xaxis.get_gridlines():
                gl.set_visible(True)
                gl.set_linestyle("-")
                gl.set_alpha(0.9)
                gl.set_color(default_grid_color)
            for gl in ax.yaxis.get_gridlines():
                gl.set_visible(False)

            # Polar spine (outer circle) subtle black
            try:
                spine = ax.spines["polar"]
                spine.set_visible(True)
                spine.set_color("#000000")
                spine.set_linewidth(1.2)
                spine.set_alpha(0.9)
            except Exception:
                pass

            # Move category labels further out
            ax.tick_params(axis="x", pad=12)

            # Resolve per-font output path and save
            out_path = _font_variant_path(georgia_path, family)
            fig.tight_layout(pad=1.0)
            plt.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.12)
            plt.savefig(out_path, dpi=300)
            plt.close()
            print(f"Radar chart ({family}) saved to: {out_path}")
        finally:
            if prev_family is not None:
                plt.rcParams["font.family"] = prev_family

    # Return the Georgia (base) path for downstream compatibility
    return georgia_path


def create_box_plot(df, eqi_columns, output_dir, method_name, k):
    """Create box plot for EQI by cluster (exports one per FONT_FAMILIES)."""
    print("Creating box plot...")

    # Melt data for plotting
    melted_df = df.melt(
        id_vars=["Cluster"],
        value_vars=eqi_columns,
        var_name="EQI_Dimension",
        value_name="EQI_Value",
    )
    melted_df["EQI_Dimension"] = melted_df["EQI_Dimension"].str.replace("EQI_", "")
    melted_df["Cluster_Label"] = "Cluster " + melted_df["Cluster"].astype(str)

    # Base (Georgia) output path; non-Georgia variants will be suffixed automatically
    base_path = os.path.join(output_dir, f"{method_name}_{k}_Box_Plot.png")

    for family in FONT_FAMILIES:
        prev_family = plt.rcParams.get("font.family", None)
        try:
            plt.rcParams["font.family"] = family

            plt.figure(figsize=(9, 6))
            sns.boxplot(
                data=melted_df,
                x="EQI_Dimension",
                y="EQI_Value",
                hue="Cluster_Label",
                hue_order=["Cluster 0", "Cluster 1", "Cluster 2"],
                palette=CLUSTER_COLORS,
                showfliers=False,
                linewidth=0.8,
                width=0.6,
            )

            # Removed x-axis label for cleaner appearance
            plt.ylabel("EQI Value (Standardized)", fontsize=18)
            plt.xlabel(" ", fontsize=18)
            ax = plt.gca()
            leg = ax.get_legend()
            if leg is not None:
                leg.remove()
            ax.tick_params(axis="both", labelsize=18)
            ax.margins(x=0.08)
            # Grid settings aligned: hide horizontal, keep light vertical
            ax.yaxis.grid(False)
            ax.xaxis.grid(True, alpha=0.3)

            fig = plt.gcf()
            fig.tight_layout(pad=0.2)
            plt.subplots_adjust(left=0.10, right=0.90, top=0.92, bottom=0.12)

            out_path = _font_variant_path(base_path, family)
            plt.savefig(out_path, dpi=300)
            plt.close()
            print(f"Box plot ({family}) saved to: {out_path}")
        finally:
            if prev_family is not None:
                plt.rcParams["font.family"] = prev_family

    # Return the Georgia (base) path for downstream compatibility
    return base_path


def create_box_plot_helvetica(df, eqi_columns, output_dir, method_name, k):
    """Backward-compatible wrapper that returns the Helvetica-styled box plot path."""
    print("Creating box plot (Helvetica)...")
    base_path = create_box_plot(df, eqi_columns, output_dir, method_name, k)
    return _font_variant_path(base_path, "Helvetica")


def create_map_for_k(df_clusters, k, shapefile_path, output_dir, method_name):
    """Create and save choropleth map for given k and method for each font family; returns Georgia path."""
    print(f"Creating map for {method_name} k={k}...")

    # Load shapefile once
    counties = gpd.read_file(shapefile_path)
    counties["COUNTY_FIPS"] = counties["STATEFP"] + counties["COUNTYFP"]

    # Filter to contiguous US
    contiguous_states = [
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
    counties_contiguous = counties[counties["STATEFP"].isin(contiguous_states)].copy()

    # Merge with cluster data
    counties_merged = counties_contiguous.merge(
        df_clusters[["COUNTY_FIPS", f"cluster_{k}"]], on="COUNTY_FIPS", how="left"
    )

    # Assign discrete colors per cluster using the editable CLUSTER_PALETTE
    colors = get_cluster_colors_list()
    color_map = {i: colors[i] for i in range(k)}
    color_map["No Data"] = CLUSTER_PALETTE.get("No Data", colors[0])

    cluster_names = {i: f"Cluster {i}" for i in range(k)}
    cluster_names["No Data"] = "No Data"

    counties_merged["cluster_color"] = counties_merged[f"cluster_{k}"].apply(
        lambda x: get_color(x, color_map)
    )

    # Base (Georgia) output path; non-Georgia variants will be suffixed automatically
    base_path = os.path.join(output_dir, f"{method_name}_{k}_Cluster_Map.png")

    # Render one map per font family
    for family in FONT_FAMILIES:
        prev_family = plt.rcParams.get("font.family", None)
        try:
            plt.rcParams["font.family"] = family

            # Plot without county borders (edgecolor none), keep map clean
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            counties_merged.plot(
                color=counties_merged["cluster_color"],
                linewidth=0,  # remove county boundary lines
                edgecolor="none",
                ax=ax,
                zorder=1,
            )

            # Weak (muted) state boundaries
            state_boundaries = counties_merged.dissolve(by="STATEFP")
            state_boundaries.boundary.plot(
                ax=ax,
                color="#000000",  # muted gray
                linewidth=0.6,  # thin line
                alpha=0.6,
                zorder=2,
            )

            ax.set_axis_off()

            # Legend
            legend_elements = [
                mpatches.Patch(color=color_map[i], label=cluster_names[i])
                for i in range(k)
            ]
            ax.legend(
                handles=legend_elements,
                bbox_to_anchor=(0.02, 0.02),
                loc="lower left",
                frameon=True,
            )

            # Save per-font output
            out_path = _font_variant_path(base_path, family)
            fig.tight_layout(pad=0.0)
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            plt.savefig(out_path, dpi=300)
            print(f"Map ({family}) saved to: {out_path}")
            plt.close()
        finally:
            if prev_family is not None:
                plt.rcParams["font.family"] = prev_family

    # Return the Georgia (base) path for downstream compatibility
    return base_path


import os

import matplotlib.pyplot as plt
from matplotlib.image import imread


def create_combined_visualization(
    radar_path, box_path, map_path, output_dir, method_name, k
):
    """
    Create combined visualization:
    (A) Radar (top-left, 45% width)
    (B) Box (top-right, 55% width)
    (C) Map (bottom, 55% height)
    Layout proportions:
      - gap = 0.02 between A and B
      - top_h = 0.45 (top section)
      - bottom_h = 0.55 (bottom section)
      - left_w = 0.45 (A width)
      - right_w = 0.55 (B width)
    """
    print(f"Creating combined visualization for {method_name} k={k}...")

    # Load images
    radar_img = imread(radar_path)
    box_img = imread(box_path)
    map_img = imread(map_path)

    # Create figure
    fig = plt.figure(
        figsize=_figsize_from_pixels(
            PIXELS_COMBINED[0], PIXELS_COMBINED[1], DEFAULT_DPI
        )
    )

    # Layout parameters
    gap = 0.0  # no gap between A and B
    top_h = 1 / 3  # top section height (1000/3000)
    bottom_h = 2 / 3  # bottom section height (2000/3000)
    left_w = 2 / 3  # left width (2000/3000)
    right_w = 1 / 3  # right width (1000/3000)

    # Margins apply ONLY to top panels (A & B)
    margin_left_top = 0.0  # no inward margin for A
    margin_right_top = 0.0  # no inward margin for B

    # Compute top section usable width
    usable_w_top = 1.0 - margin_left_top - margin_right_top

    # Normalize A and B widths based on their proportions
    total_w = left_w + right_w
    left_frac = (left_w / total_w) * (usable_w_top - gap)
    right_frac = (right_w / total_w) * (usable_w_top - gap)

    # ---- Axes positions ----
    # [left, bottom, width, height]
    axA_pos = [margin_left_top, bottom_h + gap / 2, left_frac, top_h - gap / 2]  # (A)
    axB_pos = [
        margin_left_top + left_frac + gap,
        bottom_h + gap / 2,
        right_frac,
        top_h - gap / 2,
    ]  # (B)
    axC_pos = [0.0, 0.0, 1.0, bottom_h - gap / 2]

    # ---- (A) Box ----
    axA = fig.add_axes(axA_pos)
    axA.imshow(box_img)
    axA.axis("off")
    for spine in axA.spines.values():
        spine.set_visible(False)

    # ---- (B) Radar ----
    axB = fig.add_axes(axB_pos)
    axB.imshow(radar_img)
    axB.axis("off")
    for spine in axB.spines.values():
        spine.set_visible(False)

    # ---- (C) Map ----
    axC = fig.add_axes(axC_pos)
    axC.imshow(map_img)
    axC.axis("off")
    for spine in axC.spines.values():
        spine.set_visible(False)

    # Save figure
    combined_path = os.path.join(output_dir, f"{method_name}_{k}_Combined.png")
    plt.savefig(combined_path, dpi=DEFAULT_DPI)
    plt.close()
    print(f"✅ Combined visualization saved to: {combined_path}")


def create_combined_with_labels(
    box_path: str,
    radar_path: str,
    map_path: str,
    output_dir: str,
    method_name: str,
    k: int,
    dpi: int = DEFAULT_DPI,
):
    """
    Combine three images into a single figure and overlay corner labels for each font family:
      - (A) Box at top-left
      - (B) Radar at top-right
      - (C) Map at bottom (full width)

    This function iterates FONT_FAMILIES so that adding/removing a font only requires
    editing the list at the top of this file.
    """
    for family in FONT_FAMILIES:
        prev_family = plt.rcParams.get("font.family", None)
        try:
            plt.rcParams["font.family"] = family

            # Resolve per-font image paths
            box_img_path = _font_variant_path(box_path, family)
            radar_img_path = _font_variant_path(radar_path, family)
            map_img_path = _font_variant_path(map_path, family)

            # Load images
            box_img = imread(box_img_path)
            radar_img = imread(radar_img_path)
            map_img = imread(map_img_path)

            # Create figure canvas
            fig = plt.figure(
                figsize=_figsize_from_pixels(
                    PIXELS_COMBINED[0], PIXELS_COMBINED[1], dpi
                )
            )

            # Proportions (reduce map area, enlarge top row)
            gap = 0.0
            top_h = 0.45
            bottom_h = 0.55
            left_w = 3.0 / 5.0
            right_w = 2.0 / 5.0
            margin_left_top = 0.0
            margin_right_top = 0.0
            usable_w_top = 1.0 - margin_left_top - margin_right_top
            total_w = left_w + right_w
            left_frac = (left_w / total_w) * (usable_w_top - gap)
            right_frac = (right_w / total_w) * (usable_w_top - gap)

            # Axes positions
            axA_pos = [
                margin_left_top,
                bottom_h + gap / 2,
                left_frac,
                top_h - gap / 2,
            ]  # (A)
            axB_pos = [
                margin_left_top + left_frac + gap,
                bottom_h + gap / 2,
                right_frac,
                top_h - gap / 2,
            ]  # (B)
            axC_pos = [0.0, 0.0, 1.0, bottom_h - gap / 2]  # (C)

            # Draw images
            axA = fig.add_axes(axA_pos)
            axA.imshow(box_img)
            axA.axis("off")

            axB = fig.add_axes(axB_pos)
            axB.imshow(radar_img)
            axB.axis("off")

            axC = fig.add_axes(axC_pos)
            axC.imshow(map_img)
            axC.axis("off")

            # Overlay corner labels
            label_kwargs = dict(
                fontsize=FONT_SIZES["title"],
                fontweight="bold",
                color="black",
                ha="left",
                va="top",
            )
            axA.text(0.09, 0.98, "(A)", transform=axA.transAxes, **label_kwargs)
            axB.text(0.05, 0.98, "(B)", transform=axB.transAxes, **label_kwargs)
            axC.text(0.01, 0.99, "(C)", transform=axC.transAxes, **label_kwargs)

            # Save combined figure (suffix for non-default font)
            combined_labeled_path = os.path.join(
                output_dir,
                f"{method_name}_{k}_Combined_Labeled"
                f"{'' if family == 'Georgia' else '_Helvetica'}.png",
            )
            plt.savefig(combined_labeled_path, dpi=dpi)
            plt.close()
            print(
                f"✅ Combined labeled visualization ({family}) saved to: {combined_labeled_path}"
            )
        finally:
            if prev_family is not None:
                plt.rcParams["font.family"] = prev_family


def create_radar_chart_helvetica(profiles_df, eqi_columns, output_dir, method_name, k):
    """Create radar chart (Helvetica) for cluster profiles."""
    print("Creating radar chart (Helvetica)...")
    prev_family = plt.rcParams.get("font.family", None)
    try:
        plt.rcParams["font.family"] = "Helvetica"

        # Prepare data
        categories = [col.replace("EQI_", "") for col in eqi_columns]
        fig, ax = plt.subplots(
            figsize=(5, 5),
            subplot_kw=dict(projection="polar"),
        )

        angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
        angles += angles[:1]  # Close the loop

        colors = CLUSTER_COLORS
        for _, row in profiles_df.iterrows():
            values = [row[f"{col}_mean"] for col in eqi_columns]
            values += values[:1]  # Close the loop
            cluster_id = int(row["Cluster"])

            ax.plot(
                angles,
                values,
                "o-",
                linewidth=2,
                label=f"Cluster {cluster_id}",
                color=colors[cluster_id],
            )
            ax.fill(angles, values, alpha=0.25, color=colors[cluster_id])

        # Remove circular background and align styling with Georgia
        ax.patch.set_visible(False)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories)
        ax.set_ylim(-2, 1)

        # Grid and ticks to match Georgia style
        ax.grid(True)
        ax.set_yticklabels([])
        ax.tick_params(axis="both", which="both", length=0)

        # Emphasize spoke gridlines; hide concentric rings
        default_grid_color = plt.rcParams.get("grid.color", "0.8")
        for gl in ax.xaxis.get_gridlines():
            gl.set_visible(True)
            gl.set_linestyle("-")
            gl.set_alpha(0.9)
            gl.set_color(default_grid_color)
        for gl in ax.yaxis.get_gridlines():
            gl.set_visible(False)

        # Polar spine subtle
        try:
            spine = ax.spines["polar"]
            spine.set_visible(True)
            spine.set_color("#000000")
            spine.set_linewidth(1.2)
            spine.set_alpha(0.9)
        except Exception:
            pass

        # Move category labels further out
        ax.tick_params(axis="x", pad=12)
        # Restore theta (spoke) gridline color to default and make them more visible
        default_grid_color = plt.rcParams.get("grid.color", "0.8")
        x_gridlines = ax.xaxis.get_gridlines()
        for gl in x_gridlines:
            gl.set_visible(True)
            gl.set_linestyle("-")
            gl.set_alpha(0.9)
            gl.set_color(default_grid_color)

        # Hide concentric (radial) gridlines/rings
        y_gridlines = ax.yaxis.get_gridlines()
        for gl in y_gridlines:
            gl.set_visible(False)

        # Lighten outer circular axis (polar spine)
        try:
            spine = ax.spines["polar"]
            spine.set_visible(True)
            spine.set_color("#e8e8e8")
            spine.set_linewidth(1.2)
            spine.set_alpha(0.9)
        except Exception:
            pass

        # Move category labels further out
        ax.tick_params(axis="x", pad=18)

        radar_path = os.path.join(
            output_dir, f"{method_name}_{k}_Radar_Chart_Helvetica.png"
        )
        fig.tight_layout(pad=1.0)
        plt.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.12)
        plt.savefig(radar_path, dpi=300)
        plt.close()
        print(f"Radar chart (Helvetica) saved to: {radar_path}")
        return radar_path
    finally:
        if prev_family is not None:
            plt.rcParams["font.family"] = prev_family


def create_map_for_k_helvetica(df_clusters, k, shapefile_path, output_dir, method_name):
    """Backward-compatible wrapper that returns the Helvetica-styled map path."""
    base_path = create_map_for_k(
        df_clusters, k, shapefile_path, output_dir, method_name
    )
    return _font_variant_path(base_path, "Helvetica")


def create_combined_with_labels_helvetica(
    box_path: str,
    radar_path: str,
    map_path: str,
    output_dir: str,
    method_name: str,
    k: int,
    dpi: int = DEFAULT_DPI,
):
    """
    Combine three images into a single figure and overlay corner labels using Helvetica font:
      - (A) Box at top-left
      - (B) Radar at top-right
      - (C) Map at bottom (full width)
    """
    prev_family = plt.rcParams.get("font.family", None)
    try:
        plt.rcParams["font.family"] = "Helvetica"

        # Load images
        box_img = imread(box_path)
        radar_img = imread(radar_path)
        map_img = imread(map_path)

        # Create figure canvas
        fig = plt.figure(
            figsize=_figsize_from_pixels(PIXELS_COMBINED[0], PIXELS_COMBINED[1], dpi)
        )

        # Proportions (use the same as labeled version)
        gap = 0.0
        top_h = 0.45
        bottom_h = 0.55
        left_w = 3.0 / 5.0
        right_w = 2.0 / 5.0
        margin_left_top = 0.0
        margin_right_top = 0.0
        usable_w_top = 1.0 - margin_left_top - margin_right_top
        total_w = left_w + right_w
        left_frac = (left_w / total_w) * (usable_w_top - gap)
        right_frac = (right_w / total_w) * (usable_w_top - gap)

        # Axes positions
        axA_pos = [
            margin_left_top,
            bottom_h + gap / 2,
            left_frac,
            top_h - gap / 2,
        ]  # (A)
        axB_pos = [
            margin_left_top + left_frac + gap,
            bottom_h + gap / 2,
            right_frac,
            top_h - gap / 2,
        ]  # (B)
        axC_pos = [0.0, 0.0, 1.0, bottom_h - gap / 2]  # (C)

        # Draw images
        axA = fig.add_axes(axA_pos)
        axA.imshow(box_img)
        axA.axis("off")

        axB = fig.add_axes(axB_pos)
        axB.imshow(radar_img)
        axB.axis("off")

        axC = fig.add_axes(axC_pos)
        axC.imshow(map_img)
        axC.axis("off")

        # Overlay corner labels
        label_kwargs = dict(
            fontsize=FONT_SIZES["title"],
            fontweight="bold",
            color="black",
            ha="left",
            va="top",
        )
        axA.text(0.09, 0.98, "(A)", transform=axA.transAxes, **label_kwargs)
        axB.text(0.05, 0.98, "(B)", transform=axB.transAxes, **label_kwargs)
        axC.text(0.01, 0.99, "(C)", transform=axC.transAxes, **label_kwargs)

        # Save combined figure
        combined_labeled_path = os.path.join(
            output_dir, f"{method_name}_{k}_Combined_Labeled_Helvetica.png"
        )
        plt.savefig(combined_labeled_path, dpi=dpi)
        plt.close()
        print(
            f"✅ Combined labeled visualization (Helvetica) saved to: {combined_labeled_path}"
        )
    finally:
        if prev_family is not None:
            plt.rcParams["font.family"] = prev_family
