#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM数据处理模块 - 原始数据加载处理读取整合

主要功能:
1. 多源数据加载和整合 (AAMR, EQI, Location, Urbanization) 
2. 支持普通分析和滞后效应分析的数据准备
3. 宽格式到长格式数据转换
4. 智能缺失值插补
5. 数据质量检查和验证

输入数据源:
- CDC AAMR数据 (各时间周期)
- EPA EQI数据 (2000-2005, 2006-2010)
- CDC地理和城市化数据
- 控制变量数据

输出:
- 标准化的长格式分析数据集
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import warnings
import logging

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

class LMMDataProcessor:
    """LMM数据处理器 - 统一数据加载和处理接口"""
    
    def __init__(self, eqi_period="0610", aamr_period="2016_2020"):
        """
        初始化数据处理器 - 横断面分析
        
        参数:
            eqi_period: EQI时间跨度 ("0005"=2000-2005合并, "0610"=2006-2010合并)
            aamr_period: AAMR时间跨度 ("2006_2010", "2011_2015", "2016_2020"合并)
            
        注: 所有数据都是时间跨度合并的横断面数据，不是从多年份中选择
        """
        self.eqi_period = eqi_period
        self.aamr_period = aamr_period
        self.project_root = PROJECT_ROOT
        
        # 滞后场景定义
        self.lag_scenarios = {
            ("0005", "2006_2010"): "5年滞后",
            ("0005", "2011_2015"): "10年滞后", 
            ("0610", "2011_2015"): "5年滞后",
            ("0610", "2016_2020"): "10年滞后"
        }
        
        # 数据存储
        self.aamr_data = None
        self.eqi_data = None
        self.location_data = None
        self.urbanization_data = None
        self.final_data = None
        
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
        
        # EQI变量定义（支持滞后分析）
        self.eqi_variables = {
            "EQI": "总环境质量指数",
            "EQI_air": "空气质量",
            "EQI_water": "水质量",
            "EQI_land": "土地质量", 
            "EQI_built": "建成环境质量",
            "EQI_Sociodemographic": "社会人口学质量"
        }
        
    def _setup_data_paths(self):
        """设置数据路径"""
        base_processed = self.project_root / "Data" / "Processed"
        
        paths = {
            'aamr': base_processed / "CDC" / f"CDC_EQI_AAMR_{self.aamr_period}.csv",
            'eqi': base_processed / "EQI" / f"EQI{self.eqi_period}.csv", 
            'location': base_processed / "CDC" / "Location.csv",
            'urbanization': base_processed / "CDC" / "Urbanization.csv"
        }
        
        return paths
    
    def get_scenario_description(self):
        """获取当前分析场景的描述"""
        return f"横断面分析 - EQI {self.eqi_period} 时间跨度, AAMR {self.aamr_period} 时间跨度"
    
    def load_data(self):
        """加载所有必需的数据文件"""
        logger.info("开始加载数据文件...")
        logger.info(f"分析场景: {self.get_scenario_description()}")
        
        data_paths = self._setup_data_paths()
        
        try:
            # 加载AAMR数据
            logger.info(f"加载AAMR数据: {data_paths['aamr']}")
            self.aamr_data = pd.read_csv(data_paths['aamr'])
            logger.info(f"AAMR数据形状: {self.aamr_data.shape}")
            
            # 加载EQI数据
            logger.info(f"加载EQI数据: {data_paths['eqi']}")
            self.eqi_data = pd.read_csv(data_paths['eqi'])
            # 标准化FIPS列名
            if 'COUNTY_FIPS' in self.eqi_data.columns:
                self.eqi_data['FIPS'] = self.eqi_data['COUNTY_FIPS']
            logger.info(f"EQI数据形状: {self.eqi_data.shape}")
            
            # 加载Location数据
            logger.info(f"加载Location数据: {data_paths['location']}")
            self.location_data = pd.read_csv(data_paths['location'])
            if 'COUNTY_FIPS' in self.location_data.columns:
                self.location_data['FIPS'] = self.location_data['COUNTY_FIPS']
            logger.info(f"Location数据形状: {self.location_data.shape}")
            
            # 加载Urbanization数据
            logger.info(f"加载Urbanization数据: {data_paths['urbanization']}")
            urbanization_full = pd.read_csv(data_paths['urbanization'])
            
            # 根据AAMR时间周期选择对应的城市化数据年份
            target_years = self._get_urbanization_years()
            best_year = self._select_best_year(urbanization_full, target_years)
            
            self.urbanization_data = urbanization_full[
                urbanization_full['Year'] == best_year
            ].copy()
            
            if 'COUNTY_FIPS' in self.urbanization_data.columns:
                self.urbanization_data['FIPS'] = self.urbanization_data['COUNTY_FIPS']
            
            logger.info(f"选择年份{best_year}的Urbanization数据: {self.urbanization_data.shape}")
            logger.info("所有数据文件加载完成!")
            return True
            
        except Exception as e:
            logger.error(f"数据加载失败: {e}")
            return False
    
    def _get_urbanization_years(self):
        """根据AAMR时间周期获取对应的城市化数据年份"""
        year_mapping = {
            "2016_2020": [2018, 2017, 2019, 2020, 2016],
            "2011_2015": [2013, 2014, 2012, 2015, 2011], 
            "2006_2010": [2008, 2009, 2007, 2010, 2006]
        }
        return year_mapping.get(self.aamr_period, [2018, 2020, 2019, 2017])
    
    def _select_best_year(self, urbanization_full, target_years):
        """选择覆盖县数最多的年份"""
        best_year = None
        best_coverage = 0
        
        for year in target_years:
            year_data = urbanization_full[urbanization_full['Year'] == year]
            if not year_data.empty:
                coverage = len(year_data)
                logger.info(f"  年份 {year}: {coverage} 个县")
                if coverage > best_coverage:
                    best_coverage = coverage
                    best_year = year
        
        if best_year is None:
            available_years = sorted(urbanization_full['Year'].unique(), reverse=True)
            best_year = available_years[0]
            logger.warning(f"使用最新可用年份: {best_year}")
            
        return best_year
    
    def reshape_aamr_data(self):
        """将AAMR数据从宽格式转换为长格式"""
        logger.info("开始AAMR数据格式转换 (宽 -> 长)...")
        
        if self.aamr_data is None:
            logger.error("AAMR数据未加载")
            return False
        
        try:
            # 识别AAMR列
            aamr_columns = [col for col in self.aamr_data.columns if col.startswith('AAMR_')]
            logger.info(f"发现 {len(aamr_columns)} 个AAMR列: {aamr_columns[:5]}...")
            
            # 执行melt操作
            aamr_long = pd.melt(
                self.aamr_data,
                id_vars=['FIPS'],
                value_vars=aamr_columns,
                var_name='Cancer_Type_Raw',
                value_name='AAMR'
            )
            
            # 清理Cancer_Type列名
            aamr_long['Cancer_Type'] = aamr_long['Cancer_Type_Raw'].str.replace('AAMR_', '')
            aamr_long = aamr_long.drop(columns=['Cancer_Type_Raw'])
            
            # 转换AAMR为数值类型
            aamr_long['AAMR'] = pd.to_numeric(aamr_long['AAMR'], errors='coerce')
            
            # 添加地理标识
            aamr_long['FIPS_str'] = aamr_long['FIPS'].astype(str).str.zfill(5)
            aamr_long['State_FIPS'] = aamr_long['FIPS_str'].str[:2]
            
            # 存储长格式数据
            self.aamr_long = aamr_long
            
            logger.info(f"长格式AAMR数据形状: {self.aamr_long.shape}")
            logger.info(f"缺失AAMR值数量: {self.aamr_long['AAMR'].isna().sum()}")
            logger.info(f"癌症类型数量: {self.aamr_long['Cancer_Type'].nunique()}")
            
            return True
            
        except Exception as e:
            logger.error(f"AAMR数据转换失败: {e}")
            return False
    
    def merge_all_data(self):
        """合并所有数据源"""
        logger.info("开始合并所有数据源...")
        
        if any(data is None for data in [self.aamr_long, self.eqi_data, 
                                        self.location_data, self.urbanization_data]):
            logger.error("存在未加载的数据")
            return False
        
        try:
            # 从长格式AAMR开始
            merged_data = self.aamr_long.copy()
            initial_rows = len(merged_data)
            logger.info(f"起始数据: {initial_rows} 行")
            
            # 准备EQI数据列（横断面分析包含总EQI和五大领域）
            eqi_cols = ['FIPS', 'EQI', 'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic', 'RUCC']
            # 只保留数据中实际存在的列
            eqi_cols = [col for col in eqi_cols if col in self.eqi_data.columns]
            
            # 合并EQI数据
            merged_data = pd.merge(
                merged_data, 
                self.eqi_data[eqi_cols], 
                on='FIPS', 
                how='left'
            )
            logger.info(f"合并EQI后: {len(merged_data)} 行")
            
            # 合并Location数据
            merged_data = pd.merge(
                merged_data,
                self.location_data[['FIPS', 'Census_Region', 'Census_Division']],
                on='FIPS',
                how='left'
            )
            logger.info(f"合并Location后: {len(merged_data)} 行")
            
            # 合并Urbanization数据
            merged_data = pd.merge(
                merged_data,
                self.urbanization_data[['FIPS', 'Urbanization_Code', 'Urbanization_Type']],
                on='FIPS',
                how='left'
            )
            logger.info(f"合并Urbanization后: {len(merged_data)} 行")
            
            # 添加RUCC分类（支持滞后分析）
            if 'RUCC' in merged_data.columns:
                merged_data['RUCC_Category'] = merged_data['RUCC'].apply(
                    lambda x: f'RUCC{int(x)}' if pd.notna(x) and x in [1, 2, 3, 4] else np.nan
                )
            
            # 数据质量检查
            self._log_data_quality(merged_data)
            
            self.final_data = merged_data
            return True
            
        except Exception as e:
            logger.error(f"数据合并失败: {e}")
            return False
    
    def _log_data_quality(self, data):
        """记录数据质量统计"""
        logger.info("数据质量检查:")
        logger.info(f"  总记录数: {len(data):,}")
        logger.info(f"  县数量: {data['FIPS'].nunique()}")
        logger.info(f"  州数量: {data['State_FIPS'].nunique()}")
        logger.info(f"  癌症类型数: {data['Cancer_Type'].nunique()}")
        
        # 缺失值统计
        key_vars = ['AAMR', 'EQI', 'Census_Region', 'Urbanization_Type']
        for var in key_vars:
            if var in data.columns:
                missing_pct = data[var].isna().mean() * 100
                logger.info(f"  {var}缺失: {missing_pct:.1f}%")
    
    def prepare_categorical_variables(self):
        """准备分类变量和参照组"""
        logger.info("准备分类变量...")
        
        if self.final_data is None:
            logger.error("合并数据不存在")
            return False
        
        try:
            # EQI分类变量 (1-5, 以1为参照)
            self.final_data['EQI_factor'] = pd.Categorical(
                self.final_data['EQI'], 
                categories=[1, 2, 3, 4, 5], 
                ordered=True
            )
            
            # 为横断面分析准备所有EQI变量的分类版本
            eqi_vars_to_categorize = ['EQI', 'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic']
            for eqi_var in eqi_vars_to_categorize:
                if eqi_var in self.final_data.columns:
                    # 检查是否有足够的分位数
                    unique_vals = self.final_data[eqi_var].dropna().unique()
                    if len(unique_vals) >= 5:
                        self.final_data[f'{eqi_var}_factor'] = pd.Categorical(
                            self.final_data[eqi_var],
                            categories=[1, 2, 3, 4, 5],
                            ordered=True
                        )
            
            # Census Region分类变量
            if 'Census_Region' in self.final_data.columns:
                region_categories = self.final_data['Census_Region'].dropna().unique()
                self.final_data['Census_Region_factor'] = pd.Categorical(
                    self.final_data['Census_Region'],
                    categories=sorted(region_categories)
                )
            
            # Urbanization Type分类变量
            if 'Urbanization_Type' in self.final_data.columns:
                urban_categories = self.final_data['Urbanization_Type'].dropna().unique()
                self.final_data['Urbanization_Type_factor'] = pd.Categorical(
                    self.final_data['Urbanization_Type'],
                    categories=sorted(urban_categories)
                )
            
            logger.info("分类变量设置完成")
            self._log_categorical_summary()
            
            return True
            
        except Exception as e:
            logger.error(f"分类变量准备失败: {e}")
            return False
    
    def _log_categorical_summary(self):
        """记录分类变量摘要"""
        if 'EQI' in self.final_data.columns:
            eqi_dist = self.final_data['EQI'].value_counts().sort_index()
            logger.info(f"  EQI分布: {dict(eqi_dist)}")
        
        if 'Census_Region' in self.final_data.columns:
            region_count = self.final_data['Census_Region'].nunique()
            logger.info(f"  Census Region数量: {region_count}")
        
        if 'Urbanization_Type' in self.final_data.columns:
            urban_count = self.final_data['Urbanization_Type'].nunique()
            logger.info(f"  Urbanization Type数量: {urban_count}")
    
    def remove_missing_data(self):
        """删除缺失值，不进行插值"""
        logger.info("删除缺失数据...")
        
        if self.final_data is None:
            logger.error("最终数据不存在")
            return False
        
        try:
            original_count = len(self.final_data)
            
            # 删除AAMR缺失的行
            self.final_data = self.final_data.dropna(subset=['AAMR']).copy()
            
            final_count = len(self.final_data)
            removed_count = original_count - final_count
            
            logger.info(f"删除 {removed_count} 个缺失AAMR的记录")
            logger.info(f"保留 {final_count} 个完整记录")
            
            return True
            
        except Exception as e:
            logger.error(f"删除缺失数据失败: {e}")
            return False
    
    def validate_final_data(self):
        """验证最终数据质量"""
        logger.info("验证最终数据质量...")
        
        if self.final_data is None:
            logger.error("最终数据不存在")
            return False
        
        try:
            # 基本统计
            logger.info("=== 数据质量报告 ===")
            logger.info(f"总记录数: {len(self.final_data):,}")
            logger.info(f"癌症类型数: {self.final_data['Cancer_Type'].nunique()}")
            logger.info(f"县数量: {self.final_data['FIPS'].nunique()}")
            logger.info(f"州数量: {self.final_data['State_FIPS'].nunique()}")
            
            # AAMR质量检查
            if 'AAMR' in self.final_data.columns:
                valid_aamr = self.final_data['AAMR'].notna().sum()
                aamr_range = [self.final_data['AAMR'].min(), self.final_data['AAMR'].max()]
                aamr_mean = self.final_data['AAMR'].mean()
                logger.info(f"AAMR质量: 有效数{valid_aamr:,}, 范围{aamr_range}, 均值{aamr_mean:.2f}")
            
            # 异常值检查
            warnings_count = 0
            if 'AAMR' in self.final_data.columns:
                negative_aamr = (self.final_data['AAMR'] < 0).sum()
                if negative_aamr > 0:
                    logger.warning(f"发现 {negative_aamr} 个负AAMR值")
                    warnings_count += 1
            
            if 'EQI' in self.final_data.columns:
                invalid_eqi = (~self.final_data['EQI'].isin([1, 2, 3, 4, 5])).sum()
                if invalid_eqi > 0:
                    logger.warning(f"发现 {invalid_eqi} 个无效EQI值")
                    warnings_count += 1
            
            if warnings_count == 0:
                logger.info("数据质量验证通过!")
            else:
                logger.warning(f"数据质量验证完成，发现{warnings_count}个警告")
            
            return True
            
        except Exception as e:
            logger.error(f"数据验证失败: {e}")
            return False
    
    def get_analysis_ready_data(self, cancer_type=None, stratum=None, stratum_variable=None):
        """
        获取分析就绪的数据
        
        参数:
            cancer_type: 特定癌症类型，None表示所有类型
            stratum: 分层变量值，None表示不分层
            stratum_variable: 分层变量名（如'Urbanization_Type', 'RUCC_Category'）
            
        返回:
            分析用数据框
        """
        if self.final_data is None:
            logger.error("最终数据不存在")
            return None
        
        data = self.final_data.copy()
        
        # 按癌症类型筛选
        if cancer_type:
            data = data[data['Cancer_Type'] == cancer_type]
        
        # 按分层变量筛选
        if stratum and stratum_variable:
            if stratum_variable in data.columns:
                data = data[data[stratum_variable] == stratum]
            else:
                logger.warning(f"分层变量 {stratum_variable} 不存在")
        elif stratum:
            # 向后兼容：自动检测分层类型
            for col in ['Urbanization_Type', 'RUCC_Category', 'Census_Region']:
                if col in data.columns and stratum in data[col].values:
                    data = data[data[col] == stratum]
                    break
            else:
                logger.warning(f"未找到分层值: {stratum}")
        
        # 移除关键变量缺失的记录
        required_cols = ['AAMR', 'EQI', 'State_FIPS']
        available_required = [col for col in required_cols if col in data.columns]
        data = data.dropna(subset=available_required)
        
        # 重置索引
        data = data.reset_index(drop=True)
        
        if len(data) > 0:
            logger.debug(f"分析数据准备: {len(data)}行, 癌症={cancer_type}, 分层={stratum}")
        else:
            logger.warning(f"分析数据为空: 癌症={cancer_type}, 分层={stratum}")
        
        return data
    
    def get_available_cancer_types(self):
        """获取可用的癌症类型列表"""
        if self.final_data is None:
            return []
        return sorted(self.final_data['Cancer_Type'].unique())
    
    def get_available_strata(self, stratum_variable='Urbanization_Type'):
        """获取可用的分层值列表"""
        if self.final_data is None or stratum_variable not in self.final_data.columns:
            return []
        return sorted(self.final_data[stratum_variable].dropna().unique())
    
    def save_processed_data(self, output_path=None):
        """保存处理后的数据"""
        if output_path is None:
            output_dir = self.project_root / "Data" / "Processed" / "EQI"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # 生成文件名（横断面分析）
            filename = f"LMM_CrossSectional_Data_{self.eqi_period}_{self.aamr_period}.csv"
            
            output_path = output_dir / filename
        
        if self.final_data is not None:
            self.final_data.to_csv(output_path, index=False)
            logger.info(f"处理后数据已保存到: {output_path}")
            return str(output_path)
        else:
            logger.error("没有数据可保存")
            return None
    
    def process_all(self):
        """执行完整的数据处理流程"""
        logger.info(f"=== 开始LMM数据处理流程 ===")
        logger.info(f"EQI时间跨度: {self.eqi_period}")
        logger.info(f"AAMR时间跨度: {self.aamr_period}")
        logger.info(f"场景: 横断面分析")
        
        steps = [
            ("加载数据", self.load_data),
            ("转换AAMR格式", self.reshape_aamr_data),
            ("合并数据源", self.merge_all_data),
            ("准备分类变量", self.prepare_categorical_variables),
            ("删除缺失数据", self.remove_missing_data),
            ("验证数据质量", self.validate_final_data)
        ]
        
        for step_name, step_func in steps:
            logger.info(f"执行步骤: {step_name}")
            if not step_func():
                logger.error(f"步骤失败: {step_name}")
                return False
        
        logger.info("=== 数据处理流程完成! ===")
        return True


def main():
    """主函数 - 演示用法"""
    # 横断面分析示例
    print("=== 横断面分析示例 ===")
    processor = LMMDataProcessor(eqi_period="0610", aamr_period="2016_2020")
    
    if processor.process_all():
        output_file = processor.save_processed_data()
        print(f"横断面分析数据已处理完成: {output_file}")
    



if __name__ == "__main__":
    main()