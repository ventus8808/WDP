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
- 统一数据表: /Users/ventus/Repository/WDP/Data/df/EQI_LMM_Delete_df.csv
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
    """LMM数据整合器 - 统一数据整合和预处理"""
    
    def __init__(self):
        """初始化数据整合器"""
        self.project_root = PROJECT_ROOT
        self.output_path = self.project_root / "Data" / "df" / "EQI_LMM_Delete_df.csv"
        
        # 确保输出目录存在
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 数据存储
        self.eqi_data = {}  # 存储不同时期的EQI数据
        self.aamr_data = {}  # 存储不同时期的AAMR数据
        self.smoking_data = None  # 吸烟率协变量数据
        self.location_data = None
        self.urbanization_data = None
        self.integrated_data = None
        
        # 分析场景定义
        self.analysis_scenarios = [
            {
                "name": "EQI0005_AAMR2006_2010",
                "eqi_period": "0005",
                "aamr_period": "2006_2010",
                "lag_years": 5,
                "description": "EQI 2000-2005 → AAMR 2006-2010 (5年滞后)"
            },
            {
                "name": "EQI0005_AAMR2011_2015", 
                "eqi_period": "0005",
                "aamr_period": "2011_2015",
                "lag_years": 10,
                "description": "EQI 2000-2005 → AAMR 2011-2015 (10年滞后)"
            },
            {
                "name": "EQI0610_AAMR2011_2015",
                "eqi_period": "0610", 
                "aamr_period": "2011_2015",
                "lag_years": 5,
                "description": "EQI 2006-2010 → AAMR 2011-2015 (5年滞后)"
            },
            {
                "name": "EQI0610_AAMR2016_2020",
                "eqi_period": "0610",
                "aamr_period": "2016_2020", 
                "lag_years": 10,
                "description": "EQI 2006-2010 → AAMR 2016-2020 (10年滞后)"
            }
        ]
        
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
        
    def _setup_data_paths(self):
        """设置数据路径"""
        base_processed = self.project_root / "Data" / "Processed"
        
        paths = {
            'eqi_0005': base_processed / "EQI" / "EQI0005.csv",
            'eqi_0610': base_processed / "EQI" / "EQI0610.csv",
            'aamr_2006_2010': base_processed / "CDC" / "CDC_EQI_AAMR_2006_2010.csv",
            'aamr_2011_2015': base_processed / "CDC" / "CDC_EQI_AAMR_2011_2015.csv", 
            'aamr_2016_2020': base_processed / "CDC" / "CDC_EQI_AAMR_2016_2020.csv",
            'smoking': base_processed / "Smoking" / "County_Smoking_EQI.csv",
            'location': base_processed / "CDC" / "Location.csv",
            'urbanization': base_processed / "CDC" / "Urbanization.csv"
        }
        
        return paths
    
    def load_all_data(self):
        """加载所有必需的数据文件"""
        logger.info("=== 开始加载所有数据文件 ===")
        
        data_paths = self._setup_data_paths()
        
        try:
            # 加载EQI数据
            logger.info("加载EQI数据...")
            for period in ["0005", "0610"]:
                path_key = f"eqi_{period}"
                if path_key in data_paths:
                    logger.info(f"  加载 EQI{period}: {data_paths[path_key]}")
                    self.eqi_data[period] = pd.read_csv(data_paths[path_key])
                    logger.info(f"    形状: {self.eqi_data[period].shape}")
            
            # 加载AAMR数据
            logger.info("加载AAMR数据...")
            for period in ["2006_2010", "2011_2015", "2016_2020"]:
                path_key = f"aamr_{period}"
                if path_key in data_paths:
                    logger.info(f"  加载 AAMR {period}: {data_paths[path_key]}")
                    self.aamr_data[period] = pd.read_csv(data_paths[path_key])
                    logger.info(f"    形状: {self.aamr_data[period].shape}")
            
            # 加载吸烟率数据
            logger.info("加载吸烟率协变量数据...")
            self.smoking_data = pd.read_csv(data_paths['smoking'])
            logger.info(f"  吸烟率数据: {self.smoking_data.shape}")
            
            # 加载地理数据
            logger.info("加载地理和城市化数据...")
            self.location_data = pd.read_csv(data_paths['location'])
            logger.info(f"  Location数据: {self.location_data.shape}")
            
            # 加载城市化数据（选择2018年作为代表年份）
            urbanization_full = pd.read_csv(data_paths['urbanization'])
            self.urbanization_data = urbanization_full[
                urbanization_full['Year'] == 2018
            ].copy()
            logger.info(f"  Urbanization数据 (2018年): {self.urbanization_data.shape}")
            
            logger.info("所有数据文件加载完成!")
            return True
            
        except Exception as e:
            logger.error(f"数据加载失败: {e}")
            return False
    
    def integrate_scenario_data(self, scenario):
        """整合单个分析场景的数据"""
        logger.info(f"整合场景数据: {scenario['description']}")
        
        try:
            # 获取EQI和AAMR数据
            eqi_period = scenario['eqi_period']
            aamr_period = scenario['aamr_period']
            
            if eqi_period not in self.eqi_data:
                logger.error(f"EQI数据不存在: {eqi_period}")
                return None
                
            if aamr_period not in self.aamr_data:
                logger.error(f"AAMR数据不存在: {aamr_period}")
                return None
            
            eqi_df = self.eqi_data[eqi_period].copy()
            aamr_df = self.aamr_data[aamr_period].copy()
            
            # 从AAMR开始合并
            merged_df = aamr_df.copy()
            logger.info(f"  起始AAMR数据: {len(merged_df)} 行")
            
            # 合并EQI数据
            merged_df = pd.merge(
                merged_df,
                eqi_df,
                on='COUNTY_FIPS',
                how='left',
                suffixes=('', '_eqi')
            )
            logger.info(f"  合并EQI后: {len(merged_df)} 行")
            
            # 合并吸烟率数据（根据EQI时期选择对应的吸烟率）
            if self.smoking_data is not None:
                smoking_col = f"{eqi_period}_SR"  # 0005_SR 或 0610_SR
                if smoking_col in self.smoking_data.columns:
                    smoking_subset = self.smoking_data[['COUNTY_FIPS', smoking_col]].copy()
                    smoking_subset = smoking_subset.rename(columns={smoking_col: 'Smoking_Rate'})
                    merged_df = pd.merge(
                        merged_df,
                        smoking_subset,
                        on='COUNTY_FIPS',
                        how='left'
                    )
                    logger.info(f"  合并吸烟率后: {len(merged_df)} 行")
                else:
                    logger.warning(f"吸烟率列 {smoking_col} 不存在")
            
            # 跳过Location和Urbanization数据合并，因为不需要这些列
            
            # 添加场景标识
            merged_df['Analysis_Scenario'] = scenario['name']
            merged_df['Lag_Years'] = scenario['lag_years']
            merged_df['EQI_Period'] = eqi_period
            merged_df['AAMR_Period'] = aamr_period
            
            return merged_df
            
        except Exception as e:
            logger.error(f"场景数据整合失败: {e}")
            return None
    
    def integrate_all_scenarios(self):
        """整合所有分析场景的数据"""
        logger.info("=== 开始整合所有分析场景 ===")
        
        all_scenario_data = []
        
        for scenario in self.analysis_scenarios:
            scenario_df = self.integrate_scenario_data(scenario)
            
            if scenario_df is not None:
                all_scenario_data.append(scenario_df)
                logger.info(f"场景 {scenario['name']}: {len(scenario_df)} 行")
            else:
                logger.warning(f"场景 {scenario['name']} 整合失败")
        
        if not all_scenario_data:
            logger.error("没有成功整合任何场景数据")
            return False
        
        # 合并所有场景数据
        self.integrated_data = pd.concat(all_scenario_data, ignore_index=True)
        logger.info(f"总整合数据: {len(self.integrated_data)} 行")
        
        return True
    
    def reshape_to_long_format(self):
        """将AAMR数据从宽格式转换为长格式"""
        logger.info("转换数据为长格式...")
        
        if self.integrated_data is None:
            logger.error("整合数据不存在")
            return False
        
        try:
            # 识别AAMR列
            aamr_columns = [col for col in self.integrated_data.columns if col.startswith('AAMR_')]
            logger.info(f"发现 {len(aamr_columns)} 个AAMR列")
            
            # 识别非AAMR列（用作id_vars）
            id_vars = [col for col in self.integrated_data.columns if not col.startswith('AAMR_')]
            
            # 执行melt操作
            long_data = pd.melt(
                self.integrated_data,
                id_vars=id_vars,
                value_vars=aamr_columns,
                var_name='Cancer_Type_Raw',
                value_name='AAMR'
            )
            
            # 清理Cancer_Type列名
            long_data['Cancer_Type'] = long_data['Cancer_Type_Raw'].str.replace('AAMR_', '')
            long_data = long_data.drop(columns=['Cancer_Type_Raw'])
            
            # 转换AAMR为数值类型
            long_data['AAMR'] = pd.to_numeric(long_data['AAMR'], errors='coerce')
            
            # 数据类型转换和标准化
            # COUNTY_FIPS转换为字符串
            long_data['COUNTY_FIPS'] = long_data['COUNTY_FIPS'].astype(str).str.zfill(5)
            
            # EQI相关变量转换为整型
            eqi_columns = ['RUCC', 'EQI', 'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 
                          'EQI_Sociodemographic', 'RUCC_EQI', 'RUCC_EQI_air', 'RUCC_EQI_water', 
                          'RUCC_EQI_land', 'RUCC_EQI_built', 'RUCC_EQI_Sociodemographic']
            
            for col in eqi_columns:
                if col in long_data.columns:
                    long_data[col] = pd.to_numeric(long_data[col], errors='coerce').astype('Int64')
            
            # 转换吸烟率为数值类型
            if 'Smoking_Rate' in long_data.columns:
                long_data['Smoking_Rate'] = pd.to_numeric(long_data['Smoking_Rate'], errors='coerce')
            
            # 删除不需要的列
            columns_to_drop = ['Census_Region', 'Census_Division', 'Urbanization_Code', 
                             'Urbanization_Type', 'RUCC_Category']
            long_data = long_data.drop(columns=[col for col in columns_to_drop if col in long_data.columns])
            
            # 添加癌症类型描述
            long_data['Cancer_Description'] = long_data['Cancer_Type'].map(self.cancer_types)
            
            # 添加地理标识
            long_data['State_FIPS'] = long_data['COUNTY_FIPS'].str[:2]
            
            self.integrated_data = long_data
            
            logger.info(f"长格式数据: {len(self.integrated_data)} 行")
            logger.info(f"缺失AAMR值: {self.integrated_data['AAMR'].isna().sum()} 个")
            logger.info(f"癌症类型数: {self.integrated_data['Cancer_Type'].nunique()} 种")
            logger.info(f"分析场景数: {self.integrated_data['Analysis_Scenario'].nunique()} 个")
            
            return True
            
        except Exception as e:
            logger.error(f"数据格式转换失败: {e}")
            return False
    
    def remove_missing_data(self):
        """删除缺失数据（直接删除策略）"""
        logger.info("删除缺失数据...")
        
        if self.integrated_data is None:
            logger.error("整合数据不存在")
            return False
        
        try:
            original_count = len(self.integrated_data)
            
            # 删除关键变量缺失的行
            key_vars = ['AAMR', 'EQI', 'COUNTY_FIPS', 'Smoking_Rate']
            available_vars = [var for var in key_vars if var in self.integrated_data.columns]
            
            self.integrated_data = self.integrated_data.dropna(subset=available_vars)
            
            final_count = len(self.integrated_data)
            removed_count = original_count - final_count
            removal_rate = (removed_count / original_count) * 100
            
            logger.info(f"删除 {removed_count:,} 个缺失记录 ({removal_rate:.1f}%)")
            logger.info(f"保留 {final_count:,} 个完整记录")
            
            # 数据质量检查
            self._log_final_data_quality()
            
            return True
            
        except Exception as e:
            logger.error(f"删除缺失数据失败: {e}")
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
        """保存整合后的数据"""
        logger.info(f"保存整合数据到: {self.output_path}")
        
        if self.integrated_data is None:
            logger.error("没有数据可保存")
            return False
        
        try:
            # 保存数据
            self.integrated_data.to_csv(self.output_path, index=False)
            
            # 生成数据字典
            self._generate_data_dictionary()
            
            logger.info(f"数据保存成功: {self.output_path}")
            logger.info(f"文件大小: {self.output_path.stat().st_size / (1024*1024):.2f} MB")
            
            return True
            
        except Exception as e:
            logger.error(f"数据保存失败: {e}")
            return False
    
    def _generate_data_dictionary(self):
        """生成数据字典"""
        dict_path = self.output_path.parent / "EQI_LMM_Data_Dictionary.txt"
        
        try:
            with open(dict_path, 'w', encoding='utf-8') as f:
                f.write("EQI LMM 统一数据表 - 数据字典\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"数据文件: {self.output_path.name}\n\n")
                
                f.write("数据概述:\n")
                f.write(f"- 总记录数: {len(self.integrated_data):,}\n")
                f.write(f"- 列数: {len(self.integrated_data.columns)}\n")
                f.write(f"- 县数量: {self.integrated_data['COUNTY_FIPS'].nunique():,}\n")
                f.write(f"- 癌症类型: {self.integrated_data['Cancer_Type'].nunique()}\n")
                f.write(f"- 分析场景: {self.integrated_data['Analysis_Scenario'].nunique()}\n\n")
                
                f.write("分析场景:\n")
                for i, scenario in enumerate(self.analysis_scenarios, 1):
                    f.write(f"{i}. {scenario['name']}: {scenario['description']}\n")
                f.write("\n")
                
                f.write("列变量说明:\n")
                column_descriptions = {
                    'COUNTY_FIPS': '县FIPS代码 (主键, 字符串)',
                    'State': '州缩写',
                    'RUCC': '农村-城市连续码 (1-4整型)',
                    'EQI': '环境质量指数总分 (1-5五分位数整型)',
                    'EQI_air': '空气质量指数 (1-5整型)',
                    'EQI_water': '水质量指数 (1-5整型)', 
                    'EQI_land': '土地质量指数 (1-5整型)',
                    'EQI_built': '建成环境质量指数 (1-5整型)',
                    'EQI_Sociodemographic': '社会人口学质量指数 (1-5整型)',
                    'RUCC_EQI': 'RUCC分层EQI总分 (1-5整型)',
                    'RUCC_EQI_air': 'RUCC分层空气质量指数 (1-5整型)',
                    'RUCC_EQI_water': 'RUCC分层水质量指数 (1-5整型)',
                    'RUCC_EQI_land': 'RUCC分层土地质量指数 (1-5整型)',
                    'RUCC_EQI_built': 'RUCC分层建成环境质量指数 (1-5整型)',
                    'RUCC_EQI_Sociodemographic': 'RUCC分层社会人口学质量指数 (1-5整型)',
                    'Analysis_Scenario': '分析场景标识',
                    'Lag_Years': '滞后年数 (5或10年)',
                    'EQI_Period': 'EQI时间段 (0005或0610)',
                    'AAMR_Period': 'AAMR时间段',
                    'Cancer_Type': '癌症类型ICD代码',
                    'Cancer_Description': '癌症类型描述',
                    'AAMR': '年龄调整癌症死亡率',
                    'Smoking_Rate': '县级吸烟率 (%) - 与EQI时期匹配',
                    'State_FIPS': '2位州FIPS代码'
                }
                
                for col in self.integrated_data.columns:
                    desc = column_descriptions.get(col, '待补充描述')
                    f.write(f"- {col}: {desc}\n")
                
                f.write("\n癌症类型映射:\n")
                for code, desc in self.cancer_types.items():
                    f.write(f"- {code}: {desc}\n")
            
            logger.info(f"数据字典已保存: {dict_path}")
            
        except Exception as e:
            logger.warning(f"数据字典生成失败: {e}")
    
    def get_integrated_data(self):
        """获取整合后的数据"""
        return self.integrated_data
    
    def process_all(self):
        """执行完整的数据整合流程"""
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
    print("=== EQI LMM 数据整合 ===")
    
    integrator = LMMDataIntegrator()
    
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