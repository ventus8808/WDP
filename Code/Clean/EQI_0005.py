# /Code/Clean/EQI0005.py

import pandas as pd
import numpy as np
import os
import yaml

def get_config():
    """Load configuration from config.yaml"""
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'config.yaml')
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def clean_eqi_0005(source_path, output_path):
    """
    Cleans the 2000-2005 EPA EQI data.
    This involves selecting key variables, calculating quintiles for the continuous
    index values, and integrating the stratified indices into a single set of columns.
    """
    print(f"Processing source file: {source_path}")

    # 1. 加载数据
    df = pd.read_csv(source_path)

    # 2. 处理标识符
    df_processed = pd.DataFrame()
    # 创建5位FIPS码，兼容'stfips'可能为数值或字符串的情况
    df_processed['COUNTY_FIPS'] = df['stfips'].astype(str).str.replace(r'\\.0$', '', regex=True).str.zfill(5)
    
    # 修复RUCC类型转换问题，确保正确处理NaN值
    df['cat_rucc'] = pd.to_numeric(df['cat_rucc'], errors='coerce')
    # 用-1填充NaN值，避免后续条件判断失败
    df['cat_rucc'] = df['cat_rucc'].fillna(-1).astype(int)
    df_processed['RUCC'] = df['cat_rucc']

    # 3. 处理全国指数 (计算五分位数)
    national_cols_map = {
        'EQI_22July2013': 'EQI',
        'air_EQI_22July2013': 'EQI_Air',
        'water_EQI_22July2013': 'EQI_Water',
        'land_EQI_22July2013': 'EQI_Land',
        'built_EQI_22July2013': 'EQI_Built',
        'sociod_EQI_22July2013': 'EQI_Social'  # 修改列名
    }
    for original_col, new_name in national_cols_map.items():
        # 使用qcut计算五分位数, labels=False得到0-4, +1得到1-5
        quintiles = pd.qcut(df[original_col], 5, labels=False, duplicates='drop')
        # 转换为Series以处理fillna等方法
        quintiles_series = pd.Series(quintiles)
        df_processed[new_name] = (quintiles_series.fillna(-1).astype(int) + 1).replace(0, np.nan)

    # 4. 处理并整合城乡分层指数
    # 六个领域: EQI, air, water, land, built, sociod
    domains = ['EQI', 'air', 'water', 'land', 'built', 'sociod']
    domain_display_map = {
        'EQI': 'EQI',
        'air': 'Air',
        'water': 'Water',
        'land': 'Land',
        'built': 'Built',
        'sociod': 'Social'
    }
    
    for domain in domains:
        display_domain = domain_display_map[domain]
        rucc_col_name = f'RUCC_EQI_{display_domain}' if domain != 'EQI' else 'RUCC_EQI'
        df_processed[rucc_col_name] = np.nan
        
        # 对于每个RUCC类别(1-4)，计算对应列的五分位数并填充
        for i in range(1, 5):
            if domain == 'EQI':
                col_name = f'RUCC{i}_EQI_22July2013'
            else:
                col_name = f'RUCC{i}_{domain}_EQI_22July2013'
            
            # 添加列名大小写不敏感匹配
            matching_cols = [col for col in df.columns if col.lower() == col_name.lower()]
            if matching_cols:
                col_name = matching_cols[0]
                mask = (df['cat_rucc'] == i) & (df[col_name].notna())
                if mask.sum() > 0:
                    domain_data = df.loc[mask, col_name]
                    # 直接计算五分位数，简化错误处理
                    quintiles = pd.qcut(domain_data, 5, labels=False, duplicates='drop')
                    quintiles_series = pd.Series(quintiles, index=domain_data.index)
                    df_processed.loc[mask, rucc_col_name] = quintiles_series.astype(int) + 1

    # 5. 确保最终列顺序并保存
    final_columns = [
        'COUNTY_FIPS', 'RUCC', 'EQI', 'EQI_Air', 'EQI_Water', 'EQI_Land', 'EQI_Built', 'EQI_Social',
        'RUCC_EQI', 'RUCC_EQI_Air', 'RUCC_EQI_Water', 'RUCC_EQI_Land', 'RUCC_EQI_Built', 'RUCC_EQI_Social'
    ]
    
    # 检查实际生成的列
    print("Generated columns:", list(df_processed.columns))
    
    # 确保所有需要的列都存在，不存在的列用NaN填充
    for col in final_columns:
        if col not in df_processed.columns:
            df_processed[col] = np.nan
    
    # 选择所有列
    df_final = df_processed[final_columns]
    
    # 修复：将RUCC分层列转换为整数类型（允许NaN）
    rucc_columns = ['RUCC_EQI', 'RUCC_EQI_Air', 'RUCC_EQI_Water', 'RUCC_EQI_Land', 'RUCC_EQI_Built', 'RUCC_EQI_Social']
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

    # 额外输出：计算标准化后的EQI
    df_standard = pd.DataFrame()
    df_standard['COUNTY_FIPS'] = df_processed['COUNTY_FIPS']
    df_standard['RUCC'] = df_processed['RUCC']

    # 全国指数标准化
    for original_col, new_name in national_cols_map.items():
        if original_col in df.columns:
            mean_val = df[original_col].mean()
            std_val = df[original_col].std()
            df_standard[new_name] = ((df[original_col] - mean_val) / std_val).round(4)

    # RUCC分层指数标准化
    for domain in domains:
        display_domain = domain_display_map[domain]
        rucc_col_name = f'RUCC_EQI_{display_domain}' if domain != 'EQI' else 'RUCC_EQI'
        df_standard[rucc_col_name] = np.nan
        
        for i in range(1, 5):
            if domain == 'EQI':
                col_name = f'RUCC{i}_EQI_22July2013'
            else:
                col_name = f'RUCC{i}_{domain}_EQI_22July2013'
            
            matching_cols = [col for col in df.columns if col.lower() == col_name.lower()]
            if matching_cols:
                col_name = matching_cols[0]
                mask = (df['cat_rucc'] == i) & (df[col_name].notna())
                if mask.sum() > 0:
                    domain_data = df.loc[mask, col_name]
                    mean_val = domain_data.mean()
                    std_val = domain_data.std()
                    df_standard.loc[mask, rucc_col_name] = ((domain_data - mean_val) / std_val).round(4)

    # 确保最终列顺序
    df_standard_final = df_standard[final_columns]
    
    # 保存标准化文件
    standard_output_path = output_path.replace('EQI0005.csv', 'EQI0005_Standard.csv')
    df_standard_final.to_csv(standard_output_path, index=False)
    print(f"Successfully saved standardized data to: {standard_output_path}")
    print(f"Standardized DataFrame shape: {df_standard_final.shape}")


if __name__ == '__main__':
    # 获取配置
    config = get_config()
    base_path = config['data_directories']['original']
    processed_path = config['data_directories']['processed']
    
    # 定义文件路径（使用相对路径）
    source_file = os.path.join(base_path, 'EPA EQI', '00_05_EQI.csv')
    output_file = os.path.join(processed_path, 'EQI', 'EQI0005.csv')
    
    clean_eqi_0005(source_file, output_file)