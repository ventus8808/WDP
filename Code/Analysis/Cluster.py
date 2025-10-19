# Code/Analysis/Cluster.py

import pandas as pd
import numpy as np
import os
import yaml
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import statsmodels.api as sm
from math import pi
import seaborn as sns

def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config.yaml')
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_and_prepare_data(standard_file_path):
    """Load standardized EQI data and prepare for clustering"""
    print(f"Loading data from: {standard_file_path}")
    df = pd.read_csv(standard_file_path)
    
    # Select relevant columns
    eqi_columns = ['EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social']
    required_cols = ['COUNTY_FIPS', 'RUCC'] + eqi_columns
    
    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"Missing required columns. Expected: {required_cols}")
    
    df_clean = df[required_cols].dropna()
    # Filter valid RUCC values (1-4)
    df_clean = df_clean[df_clean['RUCC'].isin([1, 2, 3, 4])]
    print(f"Data loaded: {df_clean.shape[0]} counties, {len(eqi_columns)} EQI dimensions")
    
    return df_clean, eqi_columns

def compute_residuals_and_diagnostics(df, eqi_columns, output_dir):
    """Compute residuals after regressing out RUCC effects and save diagnostics"""
    print("Computing residuals and regression diagnostics...")
    
    residuals_df = df[['COUNTY_FIPS', 'RUCC']].copy()
    diagnostics = []
    
    for col in eqi_columns:
        # Prepare data for regression
        # Use statsmodels formula API for categorical variables
        import statsmodels.formula.api as smf
        formula = f'{col} ~ C(RUCC)'
        model = smf.ols(formula, data=df).fit()
        
        # Extract residuals
        residuals_df[f'{col}_residual'] = model.resid
        
        # Collect diagnostics
        diagnostics.append({
            'Domain': col,
            'R_squared': model.rsquared,
            'F_statistic': model.fvalue,
            'F_p_value': model.f_pvalue,
            'Intercept': model.params['Intercept'],
            'RUCC_coefficients': {k: v for k, v in model.params.items() if k.startswith('C(RUCC)')},
            'RUCC_p_values': {k: v for k, v in model.pvalues.items() if k.startswith('C(RUCC)')}
        })
    
    # Save diagnostics
    diagnostics_df = pd.DataFrame(diagnostics)
    diagnostics_path = os.path.join(output_dir, 'Regression_Diagnostics.csv')
    diagnostics_df.to_csv(diagnostics_path, index=False)
    print(f"Regression diagnostics saved to: {diagnostics_path}")
    
    # Standardize residuals
    residual_cols = [f'{col}_residual' for col in eqi_columns]
    scaler = StandardScaler()
    residuals_df[residual_cols] = scaler.fit_transform(residuals_df[residual_cols])
    
    return residuals_df, residual_cols

def evaluate_clusters(X, k_range, output_dir):
    """Evaluate different numbers of clusters using elbow method and silhouette score"""
    print("Evaluating optimal number of clusters...")
    
    wcss = []
    silhouette_scores = []
    
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(X)
        wcss.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(X, kmeans.labels_))
    
    # Elbow method plot
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.plot(k_range, wcss, 'bo-')
    plt.xlabel('Number of Clusters (K)')
    plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
    plt.title('Elbow Method')
    plt.grid(True)
    
    # Silhouette score plot
    plt.subplot(1, 2, 2)
    plt.plot(k_range, silhouette_scores, 'ro-')
    plt.xlabel('Number of Clusters (K)')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Analysis')
    plt.grid(True)
    
    plt.tight_layout()
    eval_plot_path = os.path.join(output_dir, 'Cluster_Evaluation.png')
    plt.savefig(eval_plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Cluster evaluation plots saved to: {eval_plot_path}")
    
    # Find optimal K
    optimal_k_elbow = k_range[np.argmin(np.diff(wcss))]  # Elbow point
    optimal_k_silhouette = k_range[np.argmax(silhouette_scores)]
    
    print(f"Optimal K by Elbow method: {optimal_k_elbow}")
    print(f"Optimal K by Silhouette score: {optimal_k_silhouette}")
    
    # Choose final K (prefer silhouette, but consider interpretability)
    final_k = optimal_k_silhouette if optimal_k_silhouette in range(3, 7) else optimal_k_elbow
    if final_k not in range(3, 7):
        final_k = 4  # Default fallback
    
    print(f"Selected final K: {final_k}")
    return final_k

def perform_clustering(X, n_clusters):
    """Perform K-Means clustering"""
    print(f"Performing K-Means clustering with {n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    return labels, kmeans

def profile_clusters(residuals_df, residual_cols, labels, output_dir):
    """Create cluster profiles"""
    print("Creating cluster profiles...")
    
    residuals_df['Cluster'] = labels
    profiles = []
    
    for cluster in np.unique(labels):
        cluster_data = residuals_df[residuals_df['Cluster'] == cluster]
        profile = {'Cluster': cluster, 'Count': len(cluster_data)}
        
        for col in residual_cols:
            profile[f'{col}_mean'] = cluster_data[col].mean()
            profile[f'{col}_std'] = cluster_data[col].std()
        
        profiles.append(profile)
    
    profiles_df = pd.DataFrame(profiles)
    profiles_path = os.path.join(output_dir, 'Cluster_Profiles.csv')
    profiles_df.to_csv(profiles_path, index=False)
    print(f"Cluster profiles saved to: {profiles_path}")
    
    return residuals_df, profiles_df

def create_radar_chart(profiles_df, residual_cols, output_dir):
    """Create radar chart for cluster profiles"""
    print("Creating radar chart...")
    
    # Prepare data
    categories = [col.replace('_residual', '').replace('EQI_', '') for col in residual_cols]
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection='polar'))
    
    angles = [n / float(len(categories)) * 2 * pi for n in range(len(categories))]
    angles += angles[:1]  # Close the loop
    
    for _, row in profiles_df.iterrows():
        values = [row[f'{col}_mean'] for col in residual_cols]
        values += values[:1]  # Close the loop
        
        ax.plot(angles, values, 'o-', linewidth=2, label=f'Cluster {int(row["Cluster"])}')
        ax.fill(angles, values, alpha=0.25)
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_ylim(-2, 2)  # Assuming standardized residuals
    ax.set_title('Cluster Profiles - Standardized Residuals', size=16, fontweight='bold')
    ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0))
    ax.grid(True)
    
    radar_path = os.path.join(output_dir, 'Cluster_Radar_Chart.png')
    plt.savefig(radar_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Radar chart saved to: {radar_path}")

def create_heatmap(profiles_df, residual_cols, output_dir):
    """Create heatmap for cluster profiles"""
    print("Creating heatmap...")
    
    # Prepare data for heatmap
    heatmap_data = profiles_df[[f'{col}_mean' for col in residual_cols]].T
    heatmap_data.columns = [f'Cluster {i}' for i in range(len(profiles_df))]
    heatmap_data.index = [col.replace('_residual', '').replace('EQI_', '') for col in residual_cols]
    
    plt.figure(figsize=(10, 6))
    plt.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto')
    plt.colorbar(label='Standardized Residual')
    plt.xticks(range(len(heatmap_data.columns)), heatmap_data.columns)
    plt.yticks(range(len(heatmap_data.index)), heatmap_data.index)
    plt.title('Cluster Profiles Heatmap', fontweight='bold')
    
    # Add value annotations
    for i in range(len(heatmap_data.index)):
        for j in range(len(heatmap_data.columns)):
            plt.text(j, i, f'{heatmap_data.iloc[i, j]:.2f}', 
                    ha='center', va='center', color='black', fontsize=8)
    
    heatmap_path = os.path.join(output_dir, 'Cluster_Heatmap.png')
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to: {heatmap_path}")

def create_swarm_plot(residuals_df, residual_cols, output_dir):
    """Create swarm plot for residuals by cluster"""
    print("Creating swarm plot...")
    
    # Melt data for plotting
    melted_df = residuals_df.melt(id_vars=['Cluster'], value_vars=residual_cols, 
                                  var_name='EQI_Dimension', value_name='Residual')
    melted_df['EQI_Dimension'] = melted_df['EQI_Dimension'].str.replace('_residual', '').str.replace('EQI_', '')
    
    plt.figure(figsize=(12, 8))
    sns.swarmplot(data=melted_df, x='EQI_Dimension', y='Residual', hue='Cluster', palette='Set1', dodge=True)
    plt.title('Swarm Plot of Residuals by Cluster and EQI Dimension')
    plt.xticks(rotation=45)
    plt.legend(title='Cluster')
    plt.grid(True, alpha=0.3)
    
    swarm_path = os.path.join(output_dir, 'Cluster_Swarm_Plot.png')
    plt.savefig(swarm_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Swarm plot saved to: {swarm_path}")

def create_violin_plot(residuals_df, residual_cols, output_dir):
    """Create violin plot for residuals by cluster"""
    print("Creating violin plot...")
    
    # Melt data for plotting
    melted_df = residuals_df.melt(id_vars=['Cluster'], value_vars=residual_cols, 
                                  var_name='EQI_Dimension', value_name='Residual')
    melted_df['EQI_Dimension'] = melted_df['EQI_Dimension'].str.replace('_residual', '').str.replace('EQI_', '')
    
    plt.figure(figsize=(12, 8))
    sns.violinplot(data=melted_df, x='EQI_Dimension', y='Residual', hue='Cluster', palette='Set1', split=True)
    plt.title('Violin Plot of Residuals by Cluster and EQI Dimension')
    plt.xticks(rotation=45)
    plt.legend(title='Cluster')
    plt.grid(True, alpha=0.3)
    
    violin_path = os.path.join(output_dir, 'Cluster_Violin_Plot.png')
    plt.savefig(violin_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Violin plot saved to: {violin_path}")

def create_box_plot(residuals_df, residual_cols, output_dir):
    """Create box plot for residuals by cluster"""
    print("Creating box plot...")
    
    # Melt data for plotting
    melted_df = residuals_df.melt(id_vars=['Cluster'], value_vars=residual_cols, 
                                  var_name='EQI_Dimension', value_name='Residual')
    melted_df['EQI_Dimension'] = melted_df['EQI_Dimension'].str.replace('_residual', '').str.replace('EQI_', '')
    
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=melted_df, x='EQI_Dimension', y='Residual', hue='Cluster', palette='Set1')
    plt.title('Box Plot of Residuals by Cluster and EQI Dimension')
    plt.xticks(rotation=45)
    plt.legend(title='Cluster')
    plt.grid(True, alpha=0.3)
    
    box_path = os.path.join(output_dir, 'Cluster_Box_Plot.png')
    plt.savefig(box_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Box plot saved to: {box_path}")

def create_raincloud_plot(residuals_df, residual_cols, output_dir):
    """Create raincloud plot (violin + swarm) for residuals by cluster"""
    print("Creating raincloud plot...")
    
    # Melt data for plotting
    melted_df = residuals_df.melt(id_vars=['Cluster'], value_vars=residual_cols, 
                                  var_name='EQI_Dimension', value_name='Residual')
    melted_df['EQI_Dimension'] = melted_df['EQI_Dimension'].str.replace('_residual', '').str.replace('EQI_', '')
    
    # Create subplots for each EQI dimension
    dimensions = melted_df['EQI_Dimension'].unique()
    n_dims = len(dimensions)
    fig, axes = plt.subplots(n_dims, 1, figsize=(10, 6 * n_dims), sharex=False)
    if n_dims == 1:
        axes = [axes]
    
    for i, dim in enumerate(dimensions):
        ax = axes[i]
        dim_data = melted_df[melted_df['EQI_Dimension'] == dim]
        
        # Violin plot
        sns.violinplot(data=dim_data, x='Cluster', y='Residual', ax=ax, palette='Set1', inner=None)
        # Swarm plot overlay
        sns.swarmplot(data=dim_data, x='Cluster', y='Residual', ax=ax, color='black', alpha=0.6, size=3)
        
        ax.set_title(f'Raincloud Plot: {dim}')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    raincloud_path = os.path.join(output_dir, 'Cluster_Raincloud_Plot.png')
    plt.savefig(raincloud_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Raincloud plot saved to: {raincloud_path}")

def save_final_results(original_df, residuals_df, output_dir):
    """Save final clustered dataset"""
    print("Saving final results...")
    
    # Merge back to original data
    final_df = original_df.merge(residuals_df[['COUNTY_FIPS', 'Cluster']], on='COUNTY_FIPS', how='left')
    
    # Ensure COUNTY_FIPS is properly formatted as 5-digit string
    final_df['COUNTY_FIPS'] = final_df['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    final_path = os.path.join(output_dir, 'EQI_Clusters.csv')
    final_df.to_csv(final_path, index=False)
    print(f"Final clustered data saved to: {final_path}")

def main():
    """Main clustering analysis pipeline"""
    print("Starting EQI Clustering Analysis...")
    
    # Configuration
    config = get_config()
    base_path = config['data_directories']['processed']
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    output_dir = os.path.join(project_root, 'Result', 'Cluster_Analysis')
    os.makedirs(output_dir, exist_ok=True)
    
    # File paths
    standard_file = os.path.join(base_path, 'EQI', 'EQI0005_Standard.csv')
    
    # Step 1: Load and prepare data
    df, eqi_columns = load_and_prepare_data(standard_file)
    
    # Step 2: Compute residuals and diagnostics
    residuals_df, residual_cols = compute_residuals_and_diagnostics(df, eqi_columns, output_dir)
    
    # Step 3: Evaluate optimal number of clusters
    k_range = range(2, 11)
    optimal_k = evaluate_clusters(residuals_df[residual_cols].values, k_range, output_dir)
    
    # Step 4: Perform clustering
    X = residuals_df[residual_cols].values
    labels, kmeans = perform_clustering(X, optimal_k)
    
    # Step 5: Profile clusters
    residuals_df, profiles_df = profile_clusters(residuals_df, residual_cols, labels, output_dir)
    
    # Step 6: Create visualizations
    create_radar_chart(profiles_df, residual_cols, output_dir)
    create_heatmap(profiles_df, residual_cols, output_dir)
    create_swarm_plot(residuals_df, residual_cols, output_dir)
    create_violin_plot(residuals_df, residual_cols, output_dir)
    create_box_plot(residuals_df, residual_cols, output_dir)
    create_raincloud_plot(residuals_df, residual_cols, output_dir)
    
    # Step 7: Save final results
    save_final_results(df, residuals_df, output_dir)
    
    print("Clustering analysis completed successfully!")
    print(f"All results saved to: {output_dir}")

if __name__ == '__main__':
    main()