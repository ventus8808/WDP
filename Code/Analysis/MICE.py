#!/usr/bin/env python3
"""
使用MICE+PMM方法对EQI_AAMR_Point.csv中的AAMR和Smoking_Rate列进行插补
按疾病类型分层插补，每种疾病单独处理
只保留有完整EQI数据的记录
输出到EQI_AAMR_Point_MICE.csv
"""
import pandas as pd
import numpy as np
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import yaml
from pathlib import Path

# ---------------------- Config & Paths ----------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / 'config.yaml'

with CONFIG_PATH.open('r', encoding='utf-8') as f:
    CFG = yaml.safe_load(f)

# 输入和输出路径
INPUT_PATH = PROJECT_ROOT / "Data" / "Processed" / "df_EQI_AAMR" / "EQI_AAMR_Point.csv"
OUTPUT_PATH = PROJECT_ROOT / "Data" / "Processed" / "df_EQI_AAMR" / "EQI_AAMR_Point_MICE.csv"

# EQI相关列
EQI_COLS = ['RUCC','EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
            'RUCC_EQI','RUCC_EQI_Air','RUCC_EQI_Water','RUCC_EQI_Land','RUCC_EQI_Built','RUCC_EQI_Social']

def main():
    print("🧮 开始对EQI_AAMR_Point.csv进行分层MICE+PMM插补...")
    
    # 读取数据
    print(f"📂 读取数据: {INPUT_PATH}")
    df = pd.read_csv(INPUT_PATH)
    print(f"  原始数据形状: {df.shape}")
    
    # 保存原始的COUNTY_FIPS格式
    df['COUNTY_FIPS'] = df['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    # 删除缺少EQI数据的行
    eqi_cols = [col for col in EQI_COLS if col in df.columns]
    df_clean = df.dropna(subset=eqi_cols).copy()
    print(f"  删除缺失EQI数据的行后形状: {df_clean.shape}")
    
    # 获取需要插补的列
    columns_to_impute = ['AAMR', 'Smoking_Rate']
    
    # 检查这些列是否存在
    available_columns = [col for col in columns_to_impute if col in df_clean.columns]
    if not available_columns:
        print("  ⚠️ 没有找到需要插补的列")
        return
    
    print(f"  需要插补的列: {available_columns}")
    
    # 按疾病类型分层插补
    cancer_types = df_clean['Cancer_Type'].unique()
    print(f"  发现 {len(cancer_types)} 种癌症类型: {list(cancer_types)}")
    
    imputed_dfs = []
    
    for cancer_type in cancer_types:
        print(f"  🔄 正在处理癌症类型: {cancer_type}")
        cancer_df = df_clean[df_clean['Cancer_Type'] == cancer_type].copy()
        
        # 选择用于插补的特征列（数值型）
        feature_columns = []
        for col in cancer_df.columns:
            # 选择数值型列作为特征（除了需要插补的列）
            if col not in available_columns and np.issubdtype(cancer_df[col].dtype, np.number):
                feature_columns.append(col)
        
        # 创建完整的特征集（特征列 + 需要插补的列）
        all_features = feature_columns + available_columns
        
        # 提取用于插补的数据子集
        imputation_data = cancer_df[all_features].copy()
        
        # 使用MICE+PMM进行插补
        imputer = IterativeImputer(random_state=42, max_iter=10, sample_posterior=True)
        imputed_values = imputer.fit_transform(imputation_data)
        
        # 将插补后的数据放回原数据框
        imputation_data.iloc[:, :] = imputed_values
        
        # 更新原数据框中的插补列
        for col in available_columns:
            cancer_df[col] = imputation_data[col]
            
        imputed_dfs.append(cancer_df)
        print(f"    完成 {cancer_type} 的插补，处理了 {len(cancer_df)} 行数据")
    
    # 合并所有插补后的数据
    df_final = pd.concat(imputed_dfs, ignore_index=True)
    
    # 确保COUNTY_FIPS仍然是字符串格式，且有5位数字
    df_final['COUNTY_FIPS'] = df_final['COUNTY_FIPS'].astype(str).str.zfill(5)
    
    # 确保EQI相关列为整数类型
    eqi_columns = [col for col in EQI_COLS if col in df_final.columns]
    for col in eqi_columns:
        df_final[col] = pd.to_numeric(df_final[col], errors='coerce').astype('Int64')
    
    # 保持原数据的列顺序
    df_final = df_final[df.columns]
    
    # 保存结果
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(OUTPUT_PATH, index=False)
    print(f"  💾 保存插补后的数据到: {OUTPUT_PATH}")
    print(f"  最终数据形状: {df_final.shape}")

if __name__ == '__main__':
    main()