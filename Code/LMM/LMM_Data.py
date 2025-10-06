#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM数据处理模块 - 统一数据整合与预处理

主要功能:
1. 以COUNTY_FIPS为主键整合所有相关数据
2. 支持多个时间窗口的EQI和AAMR数据组合
3. 处理缺失值（直接删除策略）
4. 输出分析就绪的统一数据表

数据源:
- EQI数据: EQI0005.csv, EQI0610.csv
- AAMR数据: CDC_EQI_AAMR_2006_2010.csv, CDC_EQI_AAMR_2011_2015.csv, CDC_EQI_AAMR_2016_2020.csv
- 地理数据: Location.csv, Urbanization.csv

输出:
- 统一数据表: Data/df/EQI_LMM_Delete_df.csv
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import logging
from datetime import datetime

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ============ 路径配置 ============
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

if not CONFIG_PATH.exists():
    logger.error(f"配置文件不存在: {CONFIG_PATH}")
    sys.exit(1)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

class LMMDataIntegrator:
    """LMM数据整合器 - 基于新数据源的简化处理"""
    
    def __init__(self, use_mice: bool = False):
        """初始化数据整合器
        
        Args:
            use_mice: 是否使用MICE插补数据 (True=EQI_AAMR_Point_MICE.csv, False=EQI_AAMR_Point.csv)
        """
        self.project_root = PROJECT_ROOT
        self.use_mice = use_mice
        
        # 从config获取路径
        data_dirs = cfg.get("data_directories", {})
        
        # 设置输入和输出路径
        self.input_file = "EQI_AAMR_Point_MICE.csv" if use_mice else "EQI_AAMR_Point.csv"
        self.input_path = self.project_root / "Data" / "Processed" / "df_EQI_AAMR" / self.input_file
        
        # 输出路径 - 区分MICE和非MICE版本
        output_suffix = "_MI.csv" if use_mice else "_Delete_df.csv"
        self.output_path = self.project_root / data_dirs.get("df", "Data/df") / f"EQI_LMM{output_suffix}"
        
        # 确保输出目录存在
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 数据存储
        self.integrated_data = None  # 使用原来的变量名保持兼容性
        
        # 癌症类型映射
        self.cancer_types = {
            'C00_C97': 'All Cancers',
            'C15_C26': 'Digestive Organs',
            'C18_C21': 'Colorectal', 
            'C25': 'Pancreas',
            'C30_C39': 'Respiratory',
            'C34': 'Lung Cancer',
            'C50': 'Breast Cancer',
            'C51_C58': 'Female Genital',
            'C60_C63': 'Male Genital',
            'C61': 'Prostate',
            'C64_C68': 'Urinary Tract',
            'C76_C80': 'Ill-defined',
            'C81_C96': 'Lymphoid/Hematologic'
        }
        
    def load_all_data(self):
        """加载新的整合数据文件 - 保持原方法名"""
        logger.info("=== 加载整合数据文件 ===")
        logger.info(f"数据源: {self.input_path}")
        logger.info(f"使用MICE插补: {self.use_mice}")
        
        try:
            if not self.input_path.exists():
                logger.error(f"数据文件不存在: {self.input_path}")
                return False
            
            # 加载数据并直接赋值给integrated_data保持兼容性
            self.integrated_data = pd.read_csv(self.input_path)
            logger.info(f"数据加载成功: {self.integrated_data.shape}")
            
            # 数据预处理 - 确保数据类型正确
            # COUNTY_FIPS转换为字符串并补零
            self.integrated_data['COUNTY_FIPS'] = self.integrated_data['COUNTY_FIPS'].astype(str).str.zfill(5)
            
            # 添加必要的派生变量
            if 'Cancer_Description' not in self.integrated_data.columns:
                self.integrated_data['Cancer_Description'] = self.integrated_data['Cancer_Type'].map(self.cancer_types)
            
            if 'Analysis_Scenario' not in self.integrated_data.columns:
                # 构建分析场景名称
                eqi_period_map = {'2000-2005': '0005', '2006-2010': '0610'}
                eqi_period_str = self.integrated_data['EQI_Period'].map(eqi_period_map).fillna('0000')
                time_period_str = self.integrated_data['Time_Period'].str.replace('-', '_')
                self.integrated_data['Analysis_Scenario'] = 'EQI' + eqi_period_str + '_AAMR' + time_period_str
            
            if 'State_FIPS' not in self.integrated_data.columns:
                self.integrated_data['State_FIPS'] = self.integrated_data['COUNTY_FIPS'].str[:2]
            
            # 显示基本信息
            logger.info(f"列名: {list(self.integrated_data.columns)}")
            logger.info(f"时间段: {sorted(self.integrated_data['Time_Period'].unique())}")
            logger.info(f"癌症类型: {sorted(self.integrated_data['Cancer_Type'].unique())}")
            logger.info(f"滞后年数: {sorted(self.integrated_data['Lag_Years'].unique())}")
            
            return True
            
        except Exception as e:
            logger.error(f"数据加载失败: {e}")
            return False
    
    def integrate_all_scenarios(self):
        """整合所有场景数据 - 保持原方法名，数据已经是整合格式"""
        logger.info("=== 数据已是整合格式，跳过场景整合步骤 ===")
        # 数据已经在load_all_data中处理完毕
        return True
    
    def reshape_to_long_format(self):
        """转换为长格式 - 保持原方法名，数据已经是长格式"""
        logger.info("=== 数据已是长格式，跳过重塑步骤 ===")
        # 数据已经是长格式
        return True
    
    def remove_missing_data(self):
        """删除缺失数据 - 保持原方法名"""
        logger.info("=== 处理缺失数据 ===")
        
        if self.integrated_data is None:
            logger.error("整合数据不存在")
            return False
        
        try:
            original_count = len(self.integrated_data)
            
            # 显示缺失情况
            missing_summary = self.integrated_data.isnull().sum()
            missing_cols = missing_summary[missing_summary > 0]
            if len(missing_cols) > 0:
                logger.info("缺失数据情况:")
                for col, count in missing_cols.items():
                    logger.info(f"  {col}: {count:,} 个缺失值 ({count/original_count*100:.1f}%)")
            
            # 删除关键变量缺失的行
            key_vars = ['AAMR', 'EQI', 'COUNTY_FIPS', 'Smoking_Rate']
            available_vars = [var for var in key_vars if var in self.integrated_data.columns]
            
            logger.info(f"删除关键变量缺失的记录: {available_vars}")
            self.integrated_data = self.integrated_data.dropna(subset=available_vars)
            
            final_count = len(self.integrated_data)
            removed_count = original_count - final_count
            removal_rate = (removed_count / original_count) * 100 if original_count > 0 else 0
            
            logger.info(f"删除 {removed_count:,} 个缺失记录 ({removal_rate:.1f}%)")
            logger.info(f"保留 {final_count:,} 个完整记录")
            
            # 数据质量检查
            self._log_final_data_quality()
            
            return True
            
        except Exception as e:
            logger.error(f"缺失数据处理失败: {e}")
            return False
    

    
    def _log_final_data_quality(self):
        """记录最终数据质量统计"""
        logger.info("=== 最终数据质量报告 ===")
        
        if self.integrated_data is None:
            return
        
        # 基本统计
        logger.info(f"总记录数: {len(self.integrated_data):,}")
        logger.info(f"县数量: {self.integrated_data['COUNTY_FIPS'].nunique():,}")
        logger.info(f"州数量: {self.integrated_data['State_FIPS'].nunique()}")
        logger.info(f"癌症类型: {self.integrated_data['Cancer_Type'].nunique()}")
        logger.info(f"分析场景: {self.integrated_data['Analysis_Scenario'].nunique()}")
        
        # 按场景统计
        scenario_counts = self.integrated_data['Analysis_Scenario'].value_counts()
        logger.info("按场景记录数:")
        for scenario, count in scenario_counts.items():
            logger.info(f"  {scenario}: {count:,}")
        
        # 癌症类型统计
        cancer_counts = self.integrated_data['Cancer_Type'].value_counts()
        logger.info(f"癌症类型 (前5): {dict(cancer_counts.head())}")
        
        # AAMR质量统计
        aamr_stats = self.integrated_data['AAMR'].describe()
        logger.info(f"AAMR统计: 均值={aamr_stats['mean']:.2f}, 中位数={aamr_stats['50%']:.2f}, 范围=[{aamr_stats['min']:.2f}, {aamr_stats['max']:.2f}]")
        
        # 吸烟率质量统计
        if 'Smoking_Rate' in self.integrated_data.columns:
            smoking_stats = self.integrated_data['Smoking_Rate'].describe()
            logger.info(f"吸烟率统计: 均值={smoking_stats['mean']:.2f}%, 中位数={smoking_stats['50%']:.2f}%, 范围=[{smoking_stats['min']:.2f}%, {smoking_stats['max']:.2f}%]")
        
        # EQI分布统计
        if 'EQI' in self.integrated_data.columns:
            eqi_dist = self.integrated_data['EQI'].value_counts().sort_index()
            logger.info(f"EQI分布: {dict(eqi_dist)}")
        
        # RUCC分布统计
        if 'RUCC' in self.integrated_data.columns:
            rucc_dist = self.integrated_data['RUCC'].value_counts().sort_index()
            logger.info(f"RUCC分布: {dict(rucc_dist)}")
    
    def save_integrated_data(self):
        """保存整合后的数据 - 保持原方法名"""
        logger.info(f"保存整合数据到: {self.output_path}")
        
        if self.integrated_data is None:
            logger.error("没有数据可保存")
            return False
        
        try:
            # 保存数据
            self.integrated_data.to_csv(self.output_path, index=False)
            
            logger.info(f"数据保存成功: {self.output_path}")
            logger.info(f"文件大小: {self.output_path.stat().st_size / (1024*1024):.2f} MB")
            
            return True
            
        except Exception as e:
            logger.error(f"数据保存失败: {e}")
            return False
    
    def get_data_info(self):
        """获取数据信息"""
        if self.integrated_data is None:
            return None
            
        info = {
            'total_records': len(self.integrated_data),
            'counties': self.integrated_data['COUNTY_FIPS'].nunique(),
            'states': self.integrated_data['State_FIPS'].nunique(),
            'cancer_types': self.integrated_data['Cancer_Type'].nunique(),
            'scenarios': self.integrated_data['Analysis_Scenario'].nunique(),
            'columns': list(self.integrated_data.columns)
        }
        return info
    
    def get_integrated_data(self):
        """获取整合后的数据 - 保持原方法名"""
        return self.integrated_data
    
    def process_all(self):
        """执行完整的数据整合流程 - 保持原接口"""
        logger.info("=== 开始LMM数据整合流程 ===")
        
        steps = [
            ("加载所有数据", self.load_all_data),
            ("整合所有场景", self.integrate_all_scenarios), 
            ("转换为长格式", self.reshape_to_long_format),
            ("删除缺失数据", self.remove_missing_data),
            ("保存整合数据", self.save_integrated_data)
        ]
        
        for step_name, step_func in steps:
            logger.info(f"执行步骤: {step_name}")
            if not step_func():
                logger.error(f"步骤失败: {step_name}")
                return False
        
        logger.info("=== 数据整合流程完成! ===")
        logger.info(f"输出文件: {self.output_path}")
        
        return True


def main():
    """主函数 - 数据整合演示"""
    import argparse
    
    parser = argparse.ArgumentParser(description="EQI LMM 数据整合")
    parser.add_argument("--mice", action="store_true", help="使用MICE插补数据")
    args = parser.parse_args()
    
    print("=== EQI LMM 数据整合 ===")
    print(f"使用MICE插补: {args.mice}")
    
    integrator = LMMDataIntegrator(use_mice=args.mice)
    
    if integrator.process_all():
        print(f"\n数据整合成功完成!")
        print(f"输出文件: {integrator.output_path}")
        
        # 显示数据预览
        data = integrator.get_integrated_data()
        if data is not None:
            print(f"\n数据预览:")
            print(f"形状: {data.shape}")
            print(f"列名: {list(data.columns)}")
            print(f"\n前5行:")
            print(data.head())
    else:
        print("数据整合失败!")


if __name__ == "__main__":
    main()