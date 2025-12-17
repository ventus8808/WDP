import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

# Configure matplotlib
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Garamond", "Times New Roman"],
        "font.size": 12,
    }
)

# Define paths
base_dir = Path(__file__).resolve().parents[2]
csv_path = base_dir / "Result" / "brms_heatmap" / "Top5Cancer.csv"
output_dir = base_dir / "Result" / "brms_heatmap"

# Load and filter data
df_all = pd.read_csv(csv_path)
df_all = df_all[
    (df_all["EQI_Period"] == "2000_2005")
    & (
        df_all["Model"].isin(
            ["EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]
        )
    )
]


# Parse Q values to extract estimate and CI
def parse_q(q_str):
    if pd.isna(q_str) or str(q_str) == "0.0":
        return 0.0, (0.0, 0.0)
    match = re.match(r"([-\d.]+)\(([^)]+)\)", str(q_str))
    if match:
        est = float(match.group(1))
        ci_low, ci_high = map(float, match.group(2).split(","))
        return est, (ci_low, ci_high)
    return float(q_str), (np.nan, np.nan)


# Parse all Q columns
for col in ["Q1", "Q2", "Q3", "Q4", "Q5"]:
    df_all[f"{col}_est"], df_all[f"{col}_ci"] = zip(*df_all[col].apply(parse_q))

# Create colormap based on all Q5 estimates across all lags
all_q5_values = df_all["Q5_est"].values
abs_max = max(abs(all_q5_values.min()), abs(all_q5_values.max()))
norm = plt.Normalize(vmin=-abs_max, vmax=abs_max)
cmap = ListedColormap(plt.cm.RdBu_r(np.linspace(0.1, 0.9, 256)))
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Define display configuration
outcomes = [
    "All-site\nCancer",
    "Lung\nCancer",
    "Colorectal\nCancer",
    "Breast\nCancer",
    "Pancreatic\nCancer",
    "Prostate\nCancer",
]
models = ["EQI", "Air", "Water", "Land", "Built", "Social"]
models_data = ["EQI", "EQI_Air", "EQI_Water", "EQI_Land", "EQI_Built", "EQI_Social"]
lags = [5, 10, 15]
qs = [2, 3, 4, 5]

# Create figure
fig, axes = plt.subplots(
    len(models), len(outcomes), figsize=(len(outcomes) * 2, len(models) * 2)
)

for i, (model_display, model_data) in enumerate(zip(models, models_data)):
    for j, outcome in enumerate(outcomes):
        ax = axes[i][j]

        # Build data matrix: rows=lags, cols=qs
        data = np.full((len(lags), len(qs)), np.nan)
        for r, lag in enumerate(lags):
            for c, q in enumerate(qs):
                row = df_all[
                    (df_all["Lag"] == lag)
                    & (df_all["Outcome"] == outcome.replace("\n", " "))
                    & (df_all["Model"] == model_data)
                ]
                if not row.empty:
                    data[r, c] = row[f"Q{q}_est"].values[0]

        # Plot heatmap
        ax.pcolormesh(data, cmap=cmap, norm=norm, edgecolors="white", linewidth=0.5)
        ax.set_aspect("equal")

        # Configure tick labels (only top-left subplot)
        if i == 0 and j == 0:
            ax.tick_params(
                axis="both", which="both", length=0, labelbottom=False, labelleft=False
            )
            for idx, q in enumerate(qs):
                ax.text(
                    0.45 + idx,
                    -0.15,
                    f"Q{q}",
                    ha="center",
                    va="top",
                    transform=ax.transData,
                    fontsize=12,
                )
            for idx, lag in enumerate(lags):
                ax.text(
                    4.05,
                    0.5 + idx,
                    f"{lag}y",
                    ha="left",
                    va="center",
                    transform=ax.transData,
                    fontsize=12,
                )
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        # Model labels on left
        if j == 0:
            ax.set_ylabel(model_display, fontsize=16, rotation=90, labelpad=10)

        # Cancer titles on top
        if i == 0:
            ax.set_title(outcome.replace("Cancer", "").strip(), fontsize=15)

        # Gray borders
        for spine in ax.spines.values():
            spine.set_edgecolor("gray")
            spine.set_linewidth(1)

# Add colorbar
cbar = fig.colorbar(
    sm, ax=axes, orientation="horizontal", fraction=0.08, pad=0.1, shrink=1
)
cbar.set_label("Point Estimate of Mortality Rate Difference", fontsize=15)

plt.subplots_adjust(
    left=0.15, right=0.9, top=0.9, bottom=0.2, wspace=0.25, hspace=-0.45
)

output_path = output_dir / "Top5_Cancer.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Color block map saved to {output_path}")
plt.close(fig)
