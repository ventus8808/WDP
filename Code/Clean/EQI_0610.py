# /Code/Clean/EQI0610.py

import pandas as pd
import numpy as np
import os
import yaml

def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config.yaml')
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def clean_eqi_0610(source_path, output_path):
    """
    Cleans the 2006-2010 EPA EQI data.
    This involves selecting the pre-calculated quintile variables, renaming them
    to a standard format, and integrating the stratified indices.
    """
    print(f"Processing source file: {source_path}")

    # 1. 加载数据
    df = pd.read_csv(source_path)
    # 标准化列名，全部转为小写以便处理
    df.columns = df.columns.str.lower()

    # 2. 处理标识符
    df_processed = pd.DataFrame()
    df_processed['COUNTY_FIPS'] = df['stfips'].astype(str).str.replace(r'\\.0$', '', regex=True).str.zfill(5)
    
    # 修复RUCC类型转换问题，确保正确处理NaN值
    df['cat_rucc'] = pd.to_numeric(df['cat_rucc'], errors='coerce')
    # 用-1填充NaN值，避免后续条件判断失败
    df['cat_rucc'] = df['cat_rucc'].fillna(-1).astype(int)
    df_processed['RUCC'] = df['cat_rucc']

    # 3. 处理全国指数 (直接选择五分位数)
    national_cols_map = {
        'eqi_2jan2018_vc_5': 'EQI',
        'air_eqi_2jan2018_vc_5': 'EQI_air',
        'water_eqi_2jan2018_vc_5': 'EQI_water',
        'land_eqi_2jan2018_vc_5': 'EQI_land',
        'built_eqi_2jan2018_vc_5': 'EQI_built',
        'sociod_eqi_2jan2018_vc_5': 'EQI_Sociodemographic'  # 修改列名
    }
    for original_col, new_name in national_cols_map.items():
        if original_col in df.columns:
            df_processed[new_name] = df[original_col]

    # 4. 处理并整合城乡分层指数
    # 六个领域: EQI, air, water, land, built, sociod
    domains = ['EQI', 'air', 'water', 'land', 'built', 'sociod']
    
    for domain in domains:
        # 将sociod显示为Sociodemographic
        display_domain = 'Sociodemographic' if domain == 'sociod' else domain
        rucc_col_name = f'RUCC_EQI_{display_domain}' if domain != 'EQI' else 'RUCC_EQI'
        df_processed[rucc_col_name] = np.nan
        
        # 对于每个RUCC类别(1-4)，直接使用已有的五分位数值
        for i in range(1, 5):
            if domain == 'EQI':
                col_name = f'rucc{i}_eqi_2jan2018_vc_5'
            else:
                col_name = f'rucc{i}_{domain}_eqi_2jan2018_vc_5'
            
            # 检查列是否存在（考虑大小写）
            matching_cols = [col for col in df.columns if col.lower() == col_name.lower()]
            if matching_cols:
                col_name = matching_cols[0]
                mask = (df['cat_rucc'] == i) & (df[col_name].notna())
                if mask.sum() > 0:
                    df_processed.loc[mask, rucc_col_name] = df.loc[mask, col_name]

    # 5. 确保最终列顺序并保存
    final_columns = [
        'COUNTY_FIPS', 'RUCC', 'EQI', 'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic',
        'RUCC_EQI', 'RUCC_EQI_air', 'RUCC_EQI_water', 'RUCC_EQI_land', 'RUCC_EQI_built', 'RUCC_EQI_Sociodemographic'
    ]
    
    # 检查实际生成的列
    print("Generated columns:", list(df_processed.columns))
    
    # 确保所有需要的列都存在，不存在的列用NaN填充
    for col in final_columns:
        if col not in df_processed.columns:
            df_processed[col] = np.nan
    
    # 选择所有列
    df_final = df_processed[final_columns]
    
    # 将RUCC分层列转换为整数类型（允许NaN）
    rucc_columns = ['RUCC_EQI', 'RUCC_EQI_air', 'RUCC_EQI_water', 'RUCC_EQI_land', 'RUCC_EQI_built', 'RUCC_EQI_Sociodemographic']
    for col in rucc_columns:
        if col in df_final.columns:
            df_final[col] = df_final[col].astype('Int64')

    # 创建输出目录
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_final.to_csv(output_path, index=False)
    print(f"Successfully cleaned data and saved to: {output_path}")
    print(f"Final DataFrame shape: {df_final.shape}")
    print("Final columns:", df_final.columns.tolist())

if __name__ == '__main__':
    # 获取配置
    config = get_config()
    base_path = config['data_directories']['original']
    processed_path = config['data_directories']['processed']
    
    # 定义文件路径（使用相对路径）
    source_file = os.path.join(base_path, 'EPA EQI', '06_10_EQI.csv')
    output_file = os.path.join(processed_path, 'EQI', 'EQI0610.csv')
    
    clean_eqi_0610(source_file, output_file)