import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np

# Set font
plt.rcParams['font.family'] = 'Georgia'
plt.rcParams['font.size'] = 10

# Load the data
csv_path = '/Users/ventus/Repository/WDP/Result/brms_heatmap/brms_heatmap_2000_2005.csv'
df = pd.read_csv(csv_path)

# Filter for 5-year lag and specific models
df = df[(df['Lag'] == 5) & (df['Model'].isin(['RUCC1_EQI', 'RUCC2_EQI', 'RUCC3_EQI', 'RUCC4_EQI', 'EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']))]

# Function to parse Q values: extract estimate and CI
def parse_q(q_str):
    if pd.isna(q_str):
        return 0.0, (0.0, 0.0)
    q_str = str(q_str)
    if q_str == '0.0':
        return 0.0, (0.0, 0.0)
    match = re.match(r'([-\d.]+)\(([^)]+)\)', q_str)
    if match:
        est = float(match.group(1))
        ci_str = match.group(2)
        ci_low, ci_high = map(float, ci_str.split(','))
        return est, (ci_low, ci_high)
    else:
        return float(q_str), (np.nan, np.nan)

# Parse Q1-Q5
for col in ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']:
    df[f'{col}_est'], df[f'{col}_ci'] = zip(*df[col].apply(parse_q))

# Define outcomes and models in order
outcomes = ['All-site\nCancer', 'Lung\nCancer', 'Colorectal\nCancer', 'Breast\nCancer', 'Prostate\nCancer', 'Digestive\nSystem\nCancer']
models = ['RUCC1', 'RUCC2', 'RUCC3', 'RUCC4', 'EQI', 'Air', 'Water', 'Land', 'Built', 'Social']
models_data = ['RUCC1_EQI', 'RUCC2_EQI', 'RUCC3_EQI', 'RUCC4_EQI', 'EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']

# Define model groups with EQI first
models_all = ['EQI'] + [m for m in models if m != 'EQI']
models_rucc = ['EQI', 'RUCC1', 'RUCC2', 'RUCC3', 'RUCC4']
models_sub = ['EQI', 'Air', 'Water', 'Land', 'Built', 'Social']

model_groups = [
    ('rucc', models_rucc, ['EQI', 'RUCC1_EQI', 'RUCC2_EQI', 'RUCC3_EQI', 'RUCC4_EQI']),
    ('sub', models_sub, ['EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social'])
]

for group_name, models_group, models_data_group in model_groups:
    # Create figure - transposed, square subplots
    size_per_subplot = 3.5
    aspect_ratio = 0.8
    fig, axes = plt.subplots(len(models_group), len(outcomes), figsize=(len(outcomes)*size_per_subplot, len(models_group)*size_per_subplot * aspect_ratio), sharex=True)

    # Colormap for background: diverging for positive/negative effects
    q5_values = df['Q5_est'].values
    if len(q5_values) > 0:
        # Normalize to -1 to 1 for diverging colormap
        abs_max = max(abs(q5_values.min()), abs(q5_values.max()))
        norm = plt.Normalize(vmin=-abs_max, vmax=abs_max)
        cmap = plt.cm.RdBu_r  # Red for positive, blue for negative, white for neutral
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
    else:
        norm = None
        cmap = None
        sm = None

    for i, model_display in enumerate(models_group):
        # Calculate ylim for this model? Wait, now per model
        model_df = df[df['Model'] == models_data_group[i]]
        all_y = []
        for outcome in outcomes:
            row = model_df[model_df['Outcome'] == outcome.replace('\n', ' ')]
            if not row.empty:
                y_est = [0.0] + [row[f'Q{k}_est'].values[0] for k in [2,3,4,5]]
                y_ci_low = [0.0] + [row[f'Q{k}_ci'].values[0][0] for k in [2,3,4,5]]
                y_ci_high = [0.0] + [row[f'Q{k}_ci'].values[0][1] for k in [2,3,4,5]]
                all_y.extend(y_est + y_ci_low + y_ci_high)
        if all_y:
            max_abs = max(abs(min(all_y)), abs(max(all_y)))
            ylim_min = -max_abs * 0.5
            ylim_max = max_abs * 1.1
        else:
            ylim_min = -10
            ylim_max = 10

        for j, outcome in enumerate(outcomes):
            ax = axes[i][j]
            model_data = models_data_group[i]
            row = df[(df['Outcome'] == outcome.replace('\n', ' ')) & (df['Model'] == model_data)]
            if not row.empty:
                # Data points: Q1=0, Q2-Q5 estimates
                x = [1, 2, 3, 4, 5]
                y_est = [0.0] + [row[f'Q{k}_est'].values[0] for k in [2,3,4,5]]
                y_ci_low = [0.0] + [row[f'Q{k}_ci'].values[0][0] for k in [2,3,4,5]]
                y_ci_high = [0.0] + [row[f'Q{k}_ci'].values[0][1] for k in [2,3,4,5]]

                # Determine significance for Q5
                q5_ci_low = row['Q5_ci'].values[0][0]
                q5_ci_high = row['Q5_ci'].values[0][1]
                if q5_ci_low > 0:
                    line_color = '#333333'  # Dark gray for significant
                    fill_color = '#333333'
                    alpha = 0.3
                elif q5_ci_high < 0:
                    line_color = '#333333'  # Dark gray for significant
                    fill_color = '#333333'
                    alpha = 0.3
                else:
                    line_color = 'lightgray'  # Light gray for non-significant
                    fill_color = 'lightgray'
                    alpha = 0.1

                # Plot line
                ax.plot(x, y_est, color=line_color, linewidth=2)

                # Fill CI
                ax.fill_between(x, y_ci_low, y_ci_high, color=fill_color, alpha=alpha)

                # Set background color based on Q5_est
                if norm and cmap:
                    q5_val = row['Q5_est'].values[0]
                    color = cmap(norm(q5_val))
                    ax.set_facecolor(color)

                # Draw zero line
                ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)

                # Remove grid and ticks
                ax.grid(False)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])

                # Set ylim for this model
                ax.set_ylim(ylim_min, ylim_max)

            # Set row and column labels - transposed
            if j == 0:
                ax.set_ylabel(model_display, fontsize=16, rotation=90, ha='center', va='center', labelpad=10)
            if i == 0:
                ax.set_title(outcome.replace('\n', ' '), fontsize=14)  # Remove newline for title

    # Add shared Y label
    fig.text(0.1, 0.5, 'Mortality Rate of Difference (MRD) and 95% Confidence Intervals', va='center', rotation='vertical', fontsize=14)

    # Add colorbar
    if sm:
        cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.08)
        cbar.set_label('Effect Size of Q5 (The poorest environmental conditions)', fontsize=12)

    plt.subplots_adjust(left=0.15, right=0.85, top=0.9, bottom=0.1, wspace=0.2, hspace=0.2)
    plt.savefig(f'/Users/ventus/Repository/WDP/Result/brms_heatmap/Top5_Sparkline_Heatmap_2000_2005_Lag5_{group_name}_transposed.png', dpi=300, bbox_inches='tight')
    plt.close(fig)  # Close to free memory