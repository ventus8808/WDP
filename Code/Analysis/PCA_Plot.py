#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PCA Plotting Module

This script generates all visualizations for the PCA workflow.
It reads the pre-saved plotting data from the /Result/Figure_Original_Data directory
and creates high-quality scree plots and biplots.

This decouples the analysis from visualization, allowing for rapid iteration on plots.

Author: WDP Analysis Team (with Cursor AI)
Date: 2024-09-04
"""

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.font_manager as fm

# Clear font cache and force Georgia font
fm.fontManager.__init__()

# Force Georgia font settings with fallback
try:
    plt.rcParams.update({
        'font.family': 'Georgia',
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 16,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 14,
        'figure.titlesize': 20,
        'figure.dpi': 100,
        'savefig.dpi': 300
    })
    print("Georgia font configured successfully")
except Exception as e:
    print(f"Font configuration error: {e}")
    # Fallback to serif
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 16,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 14,
        'figure.titlesize': 20
    })

# Set seaborn to use the same font
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# --- Configuration and Setup ---
# 使用相对于项目根目录的路径
PROJECT_ROOT = Path(__file__).resolve().parents[2]
PLOT_DATA_DIR = PROJECT_ROOT / "Result/Figure_Original_Data"
FIGURES_DIR = PROJECT_ROOT / "Result/Figures/PCA_Analysis"

plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("colorblind")

# --- Plotting Functions ---

def create_scree_plot(data: pd.DataFrame, title: str, filename: Path):
    """Create and save a scree plot from eigenvalue data."""
    print(f"  Creating scree plot: {title}...")
    eigenvalues = data['Eigenvalue']
    n_components_kaiser = np.sum(eigenvalues > 1.0)

    plt.figure(figsize=(8, 5))
    sns.lineplot(x=data['Component'], y=eigenvalues, marker='o', color='navy', linestyle='-')
    plt.title(title, fontsize=16)
    plt.xlabel("Principal Component Number", fontsize=12)
    plt.ylabel("Eigenvalue", fontsize=12)
    plt.axhline(y=1, color='r', linestyle='--', label='Kaiser Criterion (Eigenvalue = 1)')
    
    if n_components_kaiser > 0:
        plt.axvline(x=n_components_kaiser, color='green', linestyle=':', 
                    label=f'{n_components_kaiser} Components Kept')
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

def create_loading_plot(loadings: pd.DataFrame, title: str, filename: Path):
    """Create a 1D bar plot for single-component PCA results."""
    print(f"  Creating 1D loading plot: {title}...")
    pc1_loadings = loadings[['PC1']].copy()
    # 修复类型检查问题：明确指定参数
    pc1_loadings = pc1_loadings.sort_values(by='PC1', ascending=False)  # type: ignore
    pc1_loadings.reset_index(inplace=True)
    pc1_loadings.rename(columns={'index': 'Variable'}, inplace=True)

    plt.figure(figsize=(10, 8))
    barplot = sns.barplot(x='PC1', y='Variable', data=pc1_loadings, hue='Variable', orient='h', legend=False)
    
    plt.xlabel('Loading on Principal Component 1', fontsize=12)
    plt.ylabel('Original Variable', fontsize=12)
    plt.title(f'{title} - PC1 Loadings', fontsize=16)
    plt.grid(True, axis='x', linestyle='--', alpha=0.6)
    plt.axvline(0, color='black', linewidth=1.5)

    # 不再显示数值标签，移除以下代码段（避免类型检查错误）
    # for i in barplot.patches:
    #     plt.text(i.get_width(), i.get_y() + i.get_height() / 2, f' {i.get_width():.3f}',
    #              va='center', ha='left' if i.get_width() >= 0 else 'right', fontsize=9)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# 删除3D绘图功能，已移除create_3d_plot函数

def create_all_2d_biplots(scores: pd.DataFrame, loadings: pd.DataFrame, title: str, base_filename: Path, analysis_name: str):
    """Create all possible 2D biplot combinations for 3+ component PCA."""
    n_components = len([col for col in loadings.columns if col.startswith('PC')])
    
    if n_components >= 3:
        print(f"  Creating all 2D biplot combinations for {n_components} components...")
        # Generate all pairwise combinations
        combinations = [(1, 2), (1, 3), (2, 3)]
        for pc1, pc2 in combinations:
            if pc1 <= n_components and pc2 <= n_components:
                # Use title-based filename: Environmental_Biplot_of_PC1_and_PC2.png
                filename = base_filename.parent / f"Environmental_Biplot_of_PC{pc1}_and_PC{pc2}.png"
                create_biplot_with_pcs(scores, loadings, title, filename, analysis_name, pc1, pc2)
    else:
        # For 2 or fewer components, create single biplot if create_biplot function exists
        # create_biplot(scores, loadings, title, base_filename, analysis_name)
        pass  # Skip if create_biplot is not defined

def create_biplot_with_pcs(scores: pd.DataFrame, loadings: pd.DataFrame, title: str, filename: Path, analysis_name: str, pc1_idx: int, pc2_idx: int):
    """Create a 2D PCA biplot for specified PC components."""
    print(f"    Creating 2D biplot: PC{pc1_idx} vs PC{pc2_idx}...")
    
    # Force font settings for this plot
    plt.rcParams.update({
        'font.family': 'Georgia',
        'font.size': 16,  # 增大基础字体
        'axes.titlesize': 20,  # 增大标题字体
        'axes.labelsize': 18,  # 增大轴标签字体
        'xtick.labelsize': 15,  # 增大刻度标签字体
        'ytick.labelsize': 15,
        'legend.fontsize': 16   # 增大图例字体
    })
    
    fig, ax = plt.subplots(figsize=(14, 14))
    
    # Column names in scores have analysis prefix
    pc1_col = f'{analysis_name}_PC{pc1_idx}'
    pc2_col = f'{analysis_name}_PC{pc2_idx}'
    pc1_loading = f'PC{pc1_idx}'
    pc2_loading = f'PC{pc2_idx}'
    
    # Plot sample scores
    sns.scatterplot(x=pc1_col, y=pc2_col, data=scores, alpha=0.2, color='dimgray', label='Samples (County-Year)', s=40, ax=ax)

    # Plot loading vectors and labels with improved positioning
    arrow_scale = np.abs(scores[[pc1_col, pc2_col]].values).max() / (np.abs(loadings[[pc1_loading, pc2_loading]].values).max() * 1.5)
    
    # Store label positions to avoid overlaps
    label_positions = []
    
    for i, var in enumerate(loadings.index):
        x, y = loadings[pc1_loading].iloc[i], loadings[pc2_loading].iloc[i]
        
        # Draw arrow
        ax.arrow(0, 0, x * arrow_scale, y * arrow_scale, 
                  color='darkred', head_width=0.08, head_length=0.1, linewidth=2, alpha=0.9)
        
        # Calculate label position with adaptive scaling
        label_x = x * arrow_scale * 1.15
        label_y = y * arrow_scale * 1.15
        
        # Adjust label position if it would be too close to plot edges
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        
        # Get plot dimensions for boundary checking
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]
        
        # Ensure label stays within plot boundaries (leave 10% margin)
        margin_x = x_range * 0.1
        margin_y = y_range * 0.1
        
        if label_x > xlim[1] - margin_x:
            label_x = xlim[1] - margin_x
        elif label_x < xlim[0] + margin_x:
            label_x = xlim[0] + margin_x
            
        if label_y > ylim[1] - margin_y:
            label_y = ylim[1] - margin_y
        elif label_y < ylim[0] + margin_y:
            label_y = ylim[0] + margin_y
        
        # Check for overlaps with existing labels and adjust if needed
        min_distance = 0.1 * min(x_range, y_range)  # Minimum distance between labels
        
        for existing_x, existing_y in label_positions:
            distance = ((label_x - existing_x)**2 + (label_y - existing_y)**2)**0.5
            if distance < min_distance:
                # Adjust position to avoid overlap
                angle = np.arctan2(label_y - existing_y, label_x - existing_x)
                label_x = existing_x + min_distance * np.cos(angle)
                label_y = existing_y + min_distance * np.sin(angle)
                
                # Re-check boundaries after adjustment
                label_x = max(xlim[0] + margin_x, min(xlim[1] - margin_x, label_x))
                label_y = max(ylim[0] + margin_y, min(ylim[1] - margin_y, label_y))
        
        label_positions.append((label_x, label_y))
        
        # Add text label with improved formatting - 将下划线替换为空格
        display_var = var.replace('_', ' ')
        ax.text(label_x, label_y, display_var, color='black', ha='center', va='center', 
                fontsize=16, fontweight='bold', fontfamily='Georgia',  # 调整环境图向量标签字号为16
                bbox=dict(facecolor='white', alpha=0.85, edgecolor='darkred', 
                         boxstyle='round,pad=0.15', linewidth=0.8))

    # Set labels and title with explicit font
    ax.set_xlabel(f'Loading on Principal Component {pc1_idx}', fontsize=16, fontweight='bold', fontfamily='Georgia')
    ax.set_ylabel(f'Loading on Principal Component {pc2_idx}', fontsize=16, fontweight='bold', fontfamily='Georgia')
    ax.set_title(f'Environmental Biplot of PC{pc1_idx} and PC{pc2_idx}', fontsize=18, fontweight='bold', fontfamily='Georgia')
    
    # Set solid line axes and remove grid
    ax.axhline(0, color='black', linestyle='-', linewidth=1.0)
    ax.axvline(0, color='black', linestyle='-', linewidth=1.0)
    ax.grid(False)  # Remove background grid
    
    # Expand axis limits to accommodate all labels
    current_xlim = ax.get_xlim()
    current_ylim = ax.get_ylim()
    
    # Add 15% padding to ensure labels fit within the plot
    x_range = current_xlim[1] - current_xlim[0]
    y_range = current_ylim[1] - current_ylim[0]
    
    ax.set_xlim(current_xlim[0] - 0.15 * x_range, current_xlim[1] + 0.15 * x_range)
    ax.set_ylim(current_ylim[0] - 0.15 * y_range, current_ylim[1] + 0.15 * y_range)
    
    # Set legend font
    legend = ax.legend(fontsize=14)
    for text in legend.get_texts():
        text.set_fontfamily('Georgia')
    
    # Set tick labels font explicitly
    for label in ax.get_xticklabels():
        label.set_fontfamily('Georgia')
    for label in ax.get_yticklabels():
        label.set_fontfamily('Georgia')
    
    ax.tick_params(axis='both', which='major', labelsize=13)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

def create_individual_loading_plots(loadings: pd.DataFrame, title: str, base_filename: Path):
    """Create individual loading plots for each principal component."""
    n_components = len([col for col in loadings.columns if col.startswith('PC')])
    
    for i in range(1, n_components + 1):
        pc_col = f'PC{i}'
        if pc_col in loadings.columns:
            print(f"  Creating loading plot for {pc_col}...")
            
            pc_data = loadings[[pc_col]].copy()
            # 按照loading的绝对值大小降序排列
            # 修复类型检查问题
            pc_data['abs_loading'] = pc_data[pc_col].abs()  # type: ignore
            pc_data = pc_data.sort_values(by='abs_loading', ascending=False)  # type: ignore
            pc_data.drop('abs_loading', axis=1, inplace=True)
            pc_data.reset_index(inplace=True)
            pc_data.rename(columns={'index': 'Variable'}, inplace=True)
            
            # Force font settings for this plot - 缩小B图字号
            plt.rcParams.update({
                'font.family': 'Georgia',
                'font.size': 11,  # 缩小基础字体
                'axes.titlesize': 16,  # 缩小标题字体
                'axes.labelsize': 13,  # 缩小轴标签字体
                'xtick.labelsize': 10,  # 缩小刻度标签字体
                'ytick.labelsize': 10
            })
            
            # 处理变量名：将下划线替换为空格
            pc_data['Variable'] = pc_data['Variable'].str.replace('_', ' ')
            
            fig, ax = plt.subplots(figsize=(10, 6))  # 减小高度，使图片不那么长
            
            # 使用灰色斜线填充的条形图
            barplot = sns.barplot(x=pc_col, y='Variable', data=pc_data, 
                                orient='h', legend=False, color='lightgray', ax=ax)
            
            # 添加斜线图案填充
            for patch in barplot.patches:
                patch.set_hatch('///')  # 斜线填充
                patch.set_edgecolor('gray')
                patch.set_linewidth(0.5)
            
            # Set font explicitly for all elements - 修改标题和轴标签，缩小B图字号
            ax.set_xlabel(f'Loading on Principal Component {i}', fontsize=13, fontweight='bold', fontfamily='Georgia')
            ax.set_ylabel('Original Variable', fontsize=13, fontweight='bold', fontfamily='Georgia')
            ax.set_title(f'Socioeconomic Biplot of PC1', fontsize=16, fontweight='bold', fontfamily='Georgia')
            ax.grid(True, axis='x', linestyle='--', alpha=0.6)
            ax.axvline(0, color='black', linewidth=2)
            
            # Set tick labels font explicitly
            for label in ax.get_xticklabels():
                label.set_fontfamily('Georgia')
            for label in ax.get_yticklabels():
                label.set_fontfamily('Georgia')
            
            ax.tick_params(axis='both', which='major', labelsize=10)  # 缩小B图刻度标签
            
            # Use title-based filename: Socioeconomic_Biplot_of_PC1.png
            filename = base_filename.parent / f"Socioeconomic_Biplot_of_PC1.png"
            plt.tight_layout()
            plt.savefig(filename, dpi=300)
            plt.close()

def create_combined_scree_plot(plot_data_dir: Path, filename: Path):
    """Create a single combined scree plot with both SVI and Climate analyses on the same axes."""
    print("  Creating combined scree plot...")
    
    # Force font settings for this plot
    plt.rcParams.update({
        'font.family': 'Georgia',
        'font.size': 14,
        'axes.titlesize': 20,
        'axes.labelsize': 16,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 14
    })
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # SVI Scree Plot - Purple color
    try:
        svi_data = pd.read_csv(plot_data_dir / 'SVI_scree_plot_data.csv')
        eigenvalues = svi_data['Eigenvalue']
        
        sns.lineplot(x=svi_data['Component'], y=eigenvalues, marker='o', color='purple', 
                    linestyle='-', linewidth=4, markersize=12, label='SVI Analysis', ax=ax)
        
    except FileNotFoundError:
        print("    Warning: SVI scree data not found")
    
    # Climate Scree Plot - Green color
    try:
        climate_data = pd.read_csv(plot_data_dir / 'Climate_scree_plot_data.csv')
        eigenvalues_climate = climate_data['Eigenvalue']
        
        sns.lineplot(x=climate_data['Component'], y=eigenvalues_climate, marker='s', color='green', 
                    linestyle='-', linewidth=4, markersize=12, label='Climate Analysis', ax=ax)
        
    except FileNotFoundError:
        print("    Warning: Climate scree data not found")
    
    # Add Kaiser criterion line - Darker red
    ax.axhline(y=1, color='darkred', linestyle='--', linewidth=3, alpha=0.8, label='Kaiser Criterion (Eigenvalue = 1)')
    
    # Set labels and title with explicit font
    ax.set_title('PCA Scree Plots Comparison', fontsize=20, fontweight='bold', pad=25, fontfamily='Georgia')
    ax.set_xlabel('Principal Component Number', fontsize=16, fontweight='bold', fontfamily='Georgia')
    ax.set_ylabel('Eigenvalue', fontsize=16, fontweight='bold', fontfamily='Georgia')
    
    # Set legend font
    legend = ax.legend(fontsize=14, frameon=True, fancybox=True, shadow=True)
    for text in legend.get_texts():
        text.set_fontfamily('Georgia')
        
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=1)
    
    # Improve aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    
    # Set tick labels font explicitly
    for label in ax.get_xticklabels():
        label.set_fontfamily('Georgia')
    for label in ax.get_yticklabels():
        label.set_fontfamily('Georgia')
    
    ax.tick_params(axis='both', which='major', labelsize=13)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_combined_pca_plot(figures_dir: Path, filename: Path):
    """Create a combined PCA plot with all subfigures arranged in specified layout."""
    print("  Creating combined PCA plot...")
    
    # Force font settings for this plot
    plt.rcParams.update({
        'font.family': 'Georgia',
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 11
    })
    
    # Define image files to load
    image_files = {
        'A': figures_dir / 'Scree_Plot.png',  # Scree图放左上角
        'B': figures_dir / 'Socioeconomic_Biplot_of_PC1.png',  # SVI图放右上角
        'C1': figures_dir / 'Environmental_Biplot_of_PC1_and_PC2.png',  # 环境图放下面
        'C2': figures_dir / 'Environmental_Biplot_of_PC1_and_PC3.png',
        'C3': figures_dir / 'Environmental_Biplot_of_PC2_and_PC3.png'
    }
    
    # Check if all required files exist
    missing_files = [label for label, path in image_files.items() if not path.exists()]
    if missing_files:
        print(f"    Warning: Missing image files for labels: {missing_files}")
        return
    
    # Load and display images
    import matplotlib.image as mpimg
    
    # Step 1: 创建下半部分 - C1, C2, C3拼接
    print("    Creating bottom row (C1, C2, C3)...")
    
    # 读取C图像
    img_c1 = mpimg.imread(str(image_files['C1']))
    img_c2 = mpimg.imread(str(image_files['C2']))
    img_c3 = mpimg.imread(str(image_files['C3']))
    
    # 计算拼接后的尺寸
    h1, w1 = img_c1.shape[:2]
    h2, w2 = img_c2.shape[:2]
    h3, w3 = img_c3.shape[:2]
    
    # 统一高度
    target_h = max(h1, h2, h3)
    
    # 调整图片尺寸（保持纵横比）
    from PIL import Image
    import numpy as np
    
    def resize_image_keep_ratio(img, target_height):
        h, w = img.shape[:2]
        ratio = w / h
        new_h = target_height
        new_w = int(target_height * ratio)
        img_pil = Image.fromarray((img * 255).astype(np.uint8))
        img_resized = img_pil.resize((new_w, new_h), Image.Resampling.LANCZOS)
        return np.array(img_resized) / 255.0
    
    img_c1_resized = resize_image_keep_ratio(img_c1, target_h)
    img_c2_resized = resize_image_keep_ratio(img_c2, target_h)
    img_c3_resized = resize_image_keep_ratio(img_c3, target_h)
    
    # 水平拼接C图（添加间距）
    gap_horizontal = 10  # C1、C2、C3之间的水平间距
    gap_ab = 30  # A和B之间的间距
    gap_vertical = 50  # 上下两排之间的间距
    
    # 创建带间距的底部拼接图
    def add_gap_between_images(images, gap):
        """在图片之间添加白色间距"""
        if len(images) == 1:
            return images[0]
        
        # 计算总宽度
        total_width = sum(img.shape[1] for img in images) + gap * (len(images) - 1)
        height = images[0].shape[0]  # 假设所有图片高度相同
        
        # 检查图片维度并创建相应的空白画布
        if len(images[0].shape) == 3:
            channels = images[0].shape[2]
            combined = np.ones((height, total_width, channels), dtype=images[0].dtype)
        else:
            combined = np.ones((height, total_width), dtype=images[0].dtype)
        
        # 逐个放置图片
        current_x = 0
        for img in images:
            img_width = img.shape[1]
            if len(img.shape) == 3 and len(combined.shape) == 3:
                combined[:, current_x:current_x + img_width, :] = img
            elif len(img.shape) == 2 and len(combined.shape) == 2:
                combined[:, current_x:current_x + img_width] = img
            else:
                # 处理维度不匹配的情况
                if len(img.shape) == 3 and len(combined.shape) == 2:
                    # 将RGB图像转换为灰度
                    img_gray = np.mean(img, axis=2)
                    combined[:, current_x:current_x + img_width] = img_gray
                elif len(img.shape) == 2 and len(combined.shape) == 3:
                    # 将灰度图像扩展为RGB
                    for c in range(combined.shape[2]):
                        combined[:, current_x:current_x + img_width, c] = img
            current_x += img_width + gap
        
        return combined
    
    bottom_row = add_gap_between_images([img_c1_resized, img_c2_resized, img_c3_resized], gap_horizontal)
    
    # Step 2: 创建上半部分 - A和B等高拼接（保持纵横比）
    print("    Creating top row (A and B)...")
    
    img_a = mpimg.imread(str(image_files['A']))
    img_b = mpimg.imread(str(image_files['B']))
    
    # 保持纵横比，使用较小的高度来避免压扁
    ha, wa = img_a.shape[:2]
    hb, wb = img_b.shape[:2]
    
    # 计算两个图片的纵横比
    ratio_a = wa / ha
    ratio_b = wb / hb
    
    # 使用较大的高度作为目标高度
    top_target_h = max(ha, hb)
    
    # 按照纵横比调整图片，保持比例不变（使用已定义的函数）
    
    img_a_resized = resize_image_keep_ratio(img_a, top_target_h)
    img_b_resized = resize_image_keep_ratio(img_b, top_target_h)
    
    # 让A和B图的总宽度与下面的C图宽度相等
    # 计算下面C图的总宽度
    bottom_total_w = bottom_row.shape[1]
    
    # 计算A和B当前的总宽度
    current_top_w = img_a_resized.shape[1] + img_b_resized.shape[1]
    
    # 如果当前宽度小于目标宽度，需要等比例放大A和B
    if current_top_w < bottom_total_w:
        scale_factor = bottom_total_w / current_top_w
        
        # 重新计算目标高度
        new_target_h = int(top_target_h * scale_factor)
        
        img_a_resized = resize_image_keep_ratio(img_a, new_target_h)
        img_b_resized = resize_image_keep_ratio(img_b, new_target_h)
    
    # 水平拼接A和B（添加间距）
    top_row = add_gap_between_images([img_a_resized, img_b_resized], gap_ab)
    
    # Step 3: 调整宽度使上下两行等宽（但保持纵横比）
    print("    Adjusting widths and combining...")
    
    top_w = top_row.shape[1]
    bottom_w = bottom_row.shape[1]
    target_w = max(top_w, bottom_w)
    
    # 不要强制调整宽度，而是通过添加空白边距来对齐
    def pad_to_width(img, target_width):
        h, w = img.shape[:2]
        if w >= target_width:
            return img
        
        # 计算需要添加的左右边距
        padding = target_width - w
        left_pad = padding // 2
        right_pad = padding - left_pad
        
        # 添加白色边距
        if len(img.shape) == 3:
            # 彩色图像
            padded = np.pad(img, ((0, 0), (left_pad, right_pad), (0, 0)), 
                          mode='constant', constant_values=1.0)  # 白色背景
        else:
            # 灰度图像
            padded = np.pad(img, ((0, 0), (left_pad, right_pad)), 
                          mode='constant', constant_values=1.0)  # 白色背景
        return padded
    
    # 使用边距对齐而不是拉伸
    if top_w != target_w:
        top_row = pad_to_width(top_row, target_w)
    if bottom_w != target_w:
        bottom_row = pad_to_width(bottom_row, target_w)
    
    # 垂直拼接上下两行（添加间距）
    def add_vertical_gap(top_img, bottom_img, gap):
        """在上下两张图之间添加白色间距"""
        h1, w1 = top_img.shape[:2]
        h2, w2 = bottom_img.shape[:2]
        target_width = max(w1, w2)
        total_height = h1 + h2 + gap
        
        # 检查图片维度
        if len(top_img.shape) == 3:
            channels = top_img.shape[2]
            combined = np.ones((total_height, target_width, channels), dtype=top_img.dtype)
        else:
            combined = np.ones((total_height, target_width), dtype=top_img.dtype)
        
        # 放置上部图片
        if len(combined.shape) == 3:
            combined[:h1, :w1, :] = top_img
            combined[h1 + gap:h1 + gap + h2, :w2, :] = bottom_img
        else:
            combined[:h1, :w1] = top_img
            combined[h1 + gap:h1 + gap + h2, :w2] = bottom_img
        
        return combined
    
    final_image = add_vertical_gap(top_row, bottom_row, gap_vertical)
    
    # Step 4: 创建matplotlib图形并添加标签和边框
    fig, ax = plt.subplots(figsize=(22, 14))  # 增大整体尺寸以适应间距
    ax.imshow(final_image)
    ax.axis('off')
    
    # 计算标签位置
    total_h, total_w = final_image.shape[:2]
    top_row_h = top_row.shape[0]
    
    # 计算每个图片的位置（考虑间距）
    # 上排A和B图的位置
    a_width = img_a_resized.shape[1]
    a_height = img_a_resized.shape[0]
    b_width = img_b_resized.shape[1]
    b_height = img_b_resized.shape[0]
    
    # 下排C1, C2, C3图的位置
    c1_width = img_c1_resized.shape[1]
    c1_height = img_c1_resized.shape[0]
    c2_width = img_c2_resized.shape[1]
    c2_height = img_c2_resized.shape[0]
    c3_width = img_c3_resized.shape[1]
    c3_height = img_c3_resized.shape[0]
    
    # 图片间距设置（用于拼接时的间距）
    gap = 10  # 图片之间的间距像素
    
    # 不再绘制边框，直接显示图片
    
    # A标签位置（左上角纯文字）
    ax.text(10, 30, 'A', 
           fontsize=22, fontweight='bold', fontfamily='Georgia', 
           color='black', ha='left', va='top')
    
    # B标签位置（B图左上角，考虑A和B之间的30像素间距）
    ax.text(a_width + gap_ab + 10, 30, 'B', 
           fontsize=22, fontweight='bold', fontfamily='Georgia', 
           color='black', ha='left', va='top')
    
    # C1, C2, C3标签位置（各自图片左上角，考虑垂直50像素间距）
    ax.text(10, top_row_h + gap_vertical + 30, 'C1', 
           fontsize=22, fontweight='bold', fontfamily='Georgia', 
           color='black', ha='left', va='top')
    
    ax.text(c1_width + gap_horizontal + 10, top_row_h + gap_vertical + 30, 'C2', 
           fontsize=22, fontweight='bold', fontfamily='Georgia', 
           color='black', ha='left', va='top')
    
    ax.text(c1_width + gap_horizontal + c2_width + gap_horizontal + 10, top_row_h + gap_vertical + 30, 'C3', 
           fontsize=22, fontweight='bold', fontfamily='Georgia', 
           color='black', ha='left', va='top')
    
    # 增加图片间距和外边距
    plt.tight_layout(pad=2.0)  # 增加外边距
    plt.subplots_adjust(wspace=0.15, hspace=0.25)  # 增加子图间距
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"    Combined PCA plot saved to: {filename}")

# --- Main Execution ---

def main():
    """Main function to find plot data and generate all visualizations."""
    print("\n" + "="*60)
    print("STARTING PCA PLOTTING PIPELINE")
    print("="*60)

    if not PLOT_DATA_DIR.exists():
        print(f"Error: Plot data directory not found at {PLOT_DATA_DIR}", file=sys.stderr)
        print("Please run the main PCA.py script first to generate the required data.", file=sys.stderr)
        return 1

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Figures will be saved to: {FIGURES_DIR}")

    # Find all unique analysis names (e.g., 'SVI', 'Climate')
    analysis_names = sorted(list(set([f.stem.split('_scree_plot_data')[0] for f in PLOT_DATA_DIR.glob('*_scree_plot_data.csv')])))

    if not analysis_names:
        print("No plot data found. Exiting.")
        return 0

    print(f"Found analyses to plot: {analysis_names}")

    for name in analysis_names:
        print(f"\n--- Generating plots for: {name.upper()} ---")
        try:
            # Load data
            loadings_data = pd.read_csv(PLOT_DATA_DIR / f'{name}_loadings_data.csv', index_col=0)
            scores_data = pd.read_csv(PLOT_DATA_DIR / f'{name}_scores_sample_data.csv')
            n_components = len([col for col in loadings_data.columns if col.startswith('PC')])
            
            if name.upper() == 'CLIMATE' and n_components >= 3:
                # For Climate with 3+ components: generate only 2D combinations (no 3D, no 1D)
                print(f"  Climate has {n_components} components - generating 2D biplots only...")
                
                # Create all 2D biplot combinations
                create_all_2d_biplots(scores_data, loadings_data, f'{name.upper()}', 
                                     FIGURES_DIR / f'{name}_biplot', name)
                
            elif name.upper() == 'SVI':
                # For SVI: generate individual loading plots (1D plots)
                create_individual_loading_plots(loadings_data, f'{name.upper()}', 
                                              FIGURES_DIR / f'{name}_loading')
                
            elif 'PC2' in loadings_data.columns:
                # For other analyses with 2+ components: standard biplot
                create_biplot_with_pcs(scores_data, loadings_data, f'{name.upper()}', 
                                     FIGURES_DIR / f'{name}_biplot.png', name, 1, 2)
            
            print(f"Successfully generated plots for {name.upper()}.")

        except FileNotFoundError as e:
            print(f"  Skipping {name}: Could not find all required data files. Missing {e.filename}", file=sys.stderr)
        except Exception as e:
            print(f"  An error occurred while plotting for {name}: {e}", file=sys.stderr)
    
    # Generate combined scree plot
    print("\n--- Creating Combined Scree Plot ---")
    create_combined_scree_plot(PLOT_DATA_DIR, FIGURES_DIR / 'Scree_Plot.png')
    
    # Generate combined PCA plot with all subfigures
    print("\n--- Creating Combined PCA Plot ---")
    create_combined_pca_plot(FIGURES_DIR, FIGURES_DIR / 'PCA_Plot.png')

    print("\n" + "="*60)
    print("PLOTTING PIPELINE COMPLETED")
    print("="*60)
    return 0

if __name__ == "__main__":
    sys.exit(main())


