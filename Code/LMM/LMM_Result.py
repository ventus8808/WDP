#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM结果输出模块 - 标准化表格生成

主要功能:
1. 从LMM模型结果生成标准化表格
2. 按照指定格式输出系数、置信区间和显著性标记
3. 支持多癌症类型和多模型的结果整理
4. 生成4个结果表（2时间跨度×2滞后期）

输出格式:
ICD_Code | Model        | Q1   | Q2              | Q3              | Q4              | Q5
C00_C97  | EQI          | 0.00 | -2.80(95%CI)*   | -1.28(95%CI)    | -14.70(95%CI)*  | ...

显著性标记:
- ***: p < 0.001
- **: p < 0.01  
- *: p < 0.05
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Optional
from datetime import datetime
from statsmodels.stats.multitest import multipletests

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LMMResultFormatter:
    """LMM结果格式化器"""
    
    def __init__(self):
        """初始化结果格式化器"""
        
        # 显著性标记配置
        self.significance_markers = {
            0.001: "***",
            0.01: "**",
            0.05: "*"
        }
        
        # 多重校正标记配置
        self.fdr_significance_markers = {
            0.001: "***†",
            0.01: "**†", 
            0.05: "*†"
        }
        
        # 模型顺序配置
        self.model_order = [
            'EQI',
            'EQI_air', 
            'EQI_water',
            'EQI_land',
            'EQI_built',
            'EQI_Sociodemographic',
            'RUCC1_RUCC_EQI',
            'RUCC1_RUCC_EQI_air',
            'RUCC1_RUCC_EQI_water', 
            'RUCC1_RUCC_EQI_land',
            'RUCC1_RUCC_EQI_built',
            'RUCC1_RUCC_EQI_Sociodemographic',
            'RUCC2_RUCC_EQI',
            'RUCC2_RUCC_EQI_air',
            'RUCC2_RUCC_EQI_water',
            'RUCC2_RUCC_EQI_land', 
            'RUCC2_RUCC_EQI_built',
            'RUCC2_RUCC_EQI_Sociodemographic',
            'RUCC3_RUCC_EQI',
            'RUCC3_RUCC_EQI_air',
            'RUCC3_RUCC_EQI_water',
            'RUCC3_RUCC_EQI_land',
            'RUCC3_RUCC_EQI_built', 
            'RUCC3_RUCC_EQI_Sociodemographic',
            'RUCC4_RUCC_EQI',
            'RUCC4_RUCC_EQI_air',
            'RUCC4_RUCC_EQI_water',
            'RUCC4_RUCC_EQI_land',
            'RUCC4_RUCC_EQI_built',
            'RUCC4_RUCC_EQI_Sociodemographic'
        ]
        
        # 简化模型名称映射
        self.model_name_mapping = {
            'RUCC1_RUCC_EQI': 'RUCC1_EQI',
            'RUCC1_RUCC_EQI_air': 'RUCC1_EQI_air',
            'RUCC1_RUCC_EQI_water': 'RUCC1_EQI_water',
            'RUCC1_RUCC_EQI_land': 'RUCC1_EQI_land',
            'RUCC1_RUCC_EQI_built': 'RUCC1_EQI_built',
            'RUCC1_RUCC_EQI_Sociodemographic': 'RUCC1_EQI_Sociodemographic',
            'RUCC2_RUCC_EQI': 'RUCC2_EQI',
            'RUCC2_RUCC_EQI_air': 'RUCC2_EQI_air',
            'RUCC2_RUCC_EQI_water': 'RUCC2_EQI_water',
            'RUCC2_RUCC_EQI_land': 'RUCC2_EQI_land',
            'RUCC2_RUCC_EQI_built': 'RUCC2_EQI_built',
            'RUCC2_RUCC_EQI_Sociodemographic': 'RUCC2_EQI_Sociodemographic',
            'RUCC3_RUCC_EQI': 'RUCC3_EQI',
            'RUCC3_RUCC_EQI_air': 'RUCC3_EQI_air',
            'RUCC3_RUCC_EQI_water': 'RUCC3_EQI_water',
            'RUCC3_RUCC_EQI_land': 'RUCC3_EQI_land',
            'RUCC3_RUCC_EQI_built': 'RUCC3_EQI_built',
            'RUCC3_RUCC_EQI_Sociodemographic': 'RUCC3_EQI_Sociodemographic',
            'RUCC4_RUCC_EQI': 'RUCC4_EQI',
            'RUCC4_RUCC_EQI_air': 'RUCC4_EQI_air',
            'RUCC4_RUCC_EQI_water': 'RUCC4_EQI_water',
            'RUCC4_RUCC_EQI_land': 'RUCC4_EQI_land',
            'RUCC4_RUCC_EQI_built': 'RUCC4_EQI_built',
            'RUCC4_RUCC_EQI_Sociodemographic': 'RUCC4_EQI_Sociodemographic'
        }
        
        # 时间场景映射
        self.scenario_mapping = {
            'EQI0005_AAMR2006_2010': {
                'eqi_period': '2000_2005',
                'aamr_period': '2006_2010',
                'lag': 5
            },
            'EQI0005_AAMR2011_2015': {
                'eqi_period': '2000_2005', 
                'aamr_period': '2011_2015',
                'lag': 10
            },
            'EQI0610_AAMR2011_2015': {
                'eqi_period': '2006_2010',
                'aamr_period': '2011_2015',
                'lag': 5
            },
            'EQI0610_AAMR2016_2020': {
                'eqi_period': '2006_2010',
                'aamr_period': '2016_2020',
                'lag': 10
            }
        }
    
    def apply_multiple_correction(self, model_results: Dict, method: str = 'fdr_bh') -> Dict:
        """
        对所有p值应用多重校正
        
        参数:
            model_results: 原始模型结果
            method: 多重校正方法 ('fdr_bh', 'bonferroni', 'holm', 'sidak')
            
        返回:
            包含校正后p值的模型结果
        """
        logger.info(f"开始应用多重校正: {method}")
        
        # 收集所有p值
        all_pvalues = []
        pvalue_locations = []  # 记录p值位置信息
        
        corrected_results = model_results.copy()
        
        # 遍历所有结果收集p值
        for scenario_name, scenario_data in model_results.get('scenario_results', {}).items():
            for cancer_type, cancer_data in scenario_data.get('cancer_results', {}).items():
                for model_name, model_result in cancer_data.items():
                    if model_result and 'coefficients' in model_result:
                        for quintile, coef_info in model_result['coefficients'].items():
                            if 'p_value' in coef_info and not pd.isna(coef_info['p_value']):
                                # 排除参照组(Q1通常为0)
                                if not (coef_info['coefficient'] == 0.0 and 
                                       coef_info['lower_ci'] == 0.0 and 
                                       coef_info['upper_ci'] == 0.0):
                                    all_pvalues.append(coef_info['p_value'])
                                    pvalue_locations.append({
                                        'scenario': scenario_name,
                                        'cancer': cancer_type, 
                                        'model': model_name,
                                        'quintile': quintile
                                    })
        
        if not all_pvalues:
            logger.warning("未找到有效的p值，跳过多重校正")
            return corrected_results
        
        logger.info(f"收集到 {len(all_pvalues)} 个p值进行校正")
        
        # 应用多重校正
        try:
            reject, pvals_corrected, alpha_sidak, alpha_bonf = multipletests(
                all_pvalues, 
                method=method,
                alpha=0.05
            )
            
            logger.info(f"多重校正完成: {method}")
            logger.info(f"校正前显著结果: {sum(np.array(all_pvalues) < 0.05)}")
            logger.info(f"校正后显著结果: {sum(reject)}")
            
        except Exception as e:
            logger.error(f"多重校正失败: {e}")
            return corrected_results
        
        # 将校正后的p值放回结果中
        for i, location in enumerate(pvalue_locations):
            scenario_name = location['scenario']
            cancer_type = location['cancer']
            model_name = location['model']
            quintile = location['quintile']
            
            # 添加校正后的p值
            coef_path = (corrected_results['scenario_results'][scenario_name]
                        ['cancer_results'][cancer_type][model_name]['coefficients'][quintile])
            
            coef_path['p_value_corrected'] = pvals_corrected[i]
            coef_path['significant_corrected'] = reject[i]
        
        return corrected_results
        
    def format_coefficient(self, coef_info: Dict, use_corrected: bool = False) -> str:
        """
        格式化单个系数及其置信区间
        
        参数:
            coef_info: 包含coefficient, lower_ci, upper_ci, p_value的字典
            use_corrected: 是否使用校正后的p值
            
        返回:
            格式化字符串
        """
        if pd.isna(coef_info['coefficient']):
            return ""
        
        coef = coef_info['coefficient']
        lower_ci = coef_info['lower_ci']
        upper_ci = coef_info['upper_ci']
        
        # 选择使用原始p值还是校正后p值
        if use_corrected and 'p_value_corrected' in coef_info:
            p_value = coef_info['p_value_corrected']
            significance_markers = self.fdr_significance_markers
        else:
            p_value = coef_info['p_value']
            significance_markers = self.significance_markers
        
        # Q1参照组特殊处理
        if coef == 0.0 and lower_ci == 0.0 and upper_ci == 0.0:
            return "0.00"
        
        # 格式化系数
        coef_str = f"{coef:.2f}"
        
        # 格式化置信区间
        ci_str = f"({lower_ci:.2f}, {upper_ci:.2f})"
        
        # 添加显著性标记
        sig_marker = ""
        if not pd.isna(p_value):
            for threshold, marker in significance_markers.items():
                if p_value < threshold:
                    sig_marker = marker
                    break
        
        # 组合结果
        if coef == 0.0:
            return "0.00"
        else:
            return f"{coef_str}{ci_str}{sig_marker}"
    
    def extract_model_results(self, model_results: Dict, cancer_type: str, scenario: str) -> Dict:
        """
        从模型结果中提取特定癌症类型和场景的结果
        
        参数:
            model_results: 完整模型结果
            cancer_type: 癌症类型
            scenario: 分析场景
            
        返回:
            模型结果字典
        """
        extracted_results = {}
        
        # 检查路径是否存在
        if 'scenario_results' not in model_results:
            return extracted_results
            
        if scenario not in model_results['scenario_results']:
            return extracted_results
            
        scenario_data = model_results['scenario_results'][scenario]
        if 'cancer_results' not in scenario_data:
            return extracted_results
            
        if cancer_type not in scenario_data['cancer_results']:
            return extracted_results
        
        cancer_data = scenario_data['cancer_results'][cancer_type]
        
        # 提取结果并重新组织
        for model_key, model_result in cancer_data.items():
            if model_result and 'coefficients' in model_result:
                extracted_results[model_key] = model_result['coefficients']
        
        return extracted_results
    
    def create_result_table(self, model_results: Dict, cancer_types: List[str], scenarios: List[str], use_corrected: bool = False) -> pd.DataFrame:
        """
        创建结果表格
        
        参数:
            model_results: 完整模型结果
            cancer_types: 癌症类型列表
            scenarios: 场景列表
            use_corrected: 是否使用校正后的p值
            
        返回:
            结果表格DataFrame
        """
        all_rows = []
        
        for scenario in scenarios:
            # 获取场景信息
            scenario_info = self.scenario_mapping.get(scenario, {})
            eqi_period = scenario_info.get('eqi_period', 'Unknown')
            aamr_period = scenario_info.get('aamr_period', 'Unknown') 
            lag = scenario_info.get('lag', 0)
            
            for cancer_type in cancer_types:
                # 提取该癌症类型在该场景下的结果
                scenario_results = self.extract_model_results(model_results, cancer_type, scenario)
                
                if not scenario_results:
                    logger.warning(f"没有找到结果: {scenario} - {cancer_type}")
                    continue
                
                # 按模型顺序处理
                for model_name in self.model_order:
                    if model_name in scenario_results:
                        coeffs = scenario_results[model_name]
                        
                        # 简化模型名称
                        display_model_name = self.model_name_mapping.get(model_name, model_name)
                        
                        # 创建行数据
                        row = {
                            'ICD_Code': cancer_type,
                            'EQI_Period': eqi_period,
                            'AAMR_Period': aamr_period,
                            'Lag': lag,
                            'Model': display_model_name,
                            'Q1': self.format_coefficient(coeffs.get('Q1', {}), use_corrected),
                            'Q2': self.format_coefficient(coeffs.get('Q2', {}), use_corrected),
                            'Q3': self.format_coefficient(coeffs.get('Q3', {}), use_corrected),
                            'Q4': self.format_coefficient(coeffs.get('Q4', {}), use_corrected),
                            'Q5': self.format_coefficient(coeffs.get('Q5', {}), use_corrected)
                        }
                        
                        all_rows.append(row)
        
        # 创建DataFrame
        if all_rows:
            df = pd.DataFrame(all_rows)
            # 重新排列列顺序，ICD_Code放第一列
            df = df[['ICD_Code', 'EQI_Period', 'AAMR_Period', 'Lag', 'Model', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']]
        else:
            # 创建空表格
            df = pd.DataFrame(columns=['ICD_Code', 'EQI_Period', 'AAMR_Period', 'Lag', 'Model', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5'])
        
        return df
    
    # 删除lag_specific_tables方法，因为不再需要
    
    def save_results(self, model_results: Dict, cancer_types: List[str], output_dir: str, apply_correction: bool = True) -> Dict[str, str]:
        """
        保存结果表格到文件 - 按癌症类型分别输出
        
        参数:
            model_results: 完整模型结果
            cancer_types: 癌症类型列表
            output_dir: 输出目录
            apply_correction: 是否应用多重校正
            
        返回:
            保存的文件路径字典
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        saved_files = {}
        all_scenarios = list(self.scenario_mapping.keys())
        
        # 如果需要，应用多重校正
        corrected_results = model_results
        if apply_correction:
            corrected_results = self.apply_multiple_correction(model_results, method='fdr_bh')
        
        # 为每个癌症类型创建结果表
        for cancer_type in cancer_types:
            # 保存原始结果（无校正）
            cancer_table = self.create_result_table(model_results, [cancer_type], all_scenarios, use_corrected=False)
            
            if not cancer_table.empty:
                cancer_path = output_path / f"LMM_{cancer_type}_Results.csv"
                cancer_table.to_csv(cancer_path, index=False)
                saved_files[f'cancer_{cancer_type}'] = str(cancer_path)
                logger.info(f"{cancer_type} 结果表已保存: {cancer_path}")
            else:
                logger.warning(f"{cancer_type} 无有效结果数据")
            
            # 如果应用了多重校正，保存校正后结果
            if apply_correction:
                cancer_table_corrected = self.create_result_table(corrected_results, [cancer_type], all_scenarios, use_corrected=True)
                
                if not cancer_table_corrected.empty:
                    cancer_path_corrected = output_path / f"LMM_{cancer_type}_Results_FDR.csv"
                    cancer_table_corrected.to_csv(cancer_path_corrected, index=False)
                    saved_files[f'cancer_{cancer_type}_fdr'] = str(cancer_path_corrected)
                    logger.info(f"{cancer_type} FDR校正结果已保存: {cancer_path_corrected}")
        
        return saved_files
    
    # 删除_generate_summary_file方法，因为不再需要txt输出
    
    def create_display_table(self, model_results: Dict, cancer_types: List[str], scenario: str) -> pd.DataFrame:
        """
        创建用于显示的简化表格
        
        参数:
            model_results: 模型结果
            cancer_types: 癌症类型列表
            scenario: 特定场景
            
        返回:
            显示表格
        """
        table = self.create_result_table(model_results, cancer_types, [scenario])
        
        if not table.empty:
            # 保留新的列结构
            return table
        else:
            return pd.DataFrame(columns=['EQI_Period', 'AAMR_Period', 'Lag', 'ICD_Code', 'Model', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5'])


def main():
    """主函数 - 演示用法"""
    print("=== LMM结果格式化测试 ===")
    
    # 创建模拟结果数据进行测试
    mock_results = {
        'analysis_info': {
            'total_scenarios': 4,
            'total_cancer_types': 1,
            'data_source': '/test/data.csv',
            'analysis_timestamp': '2025-10-01 17:42:32'
        },
        'scenario_results': {
            'EQI0005_AAMR2006_2010': {
                'cancer_results': {
                    'C00_C97': {
                        'EQI': {
                            'coefficients': {
                                'Q1': {'coefficient': 0.0, 'lower_ci': 0.0, 'upper_ci': 0.0, 'p_value': np.nan},
                                'Q2': {'coefficient': -2.80, 'lower_ci': -5.20, 'upper_ci': -0.40, 'p_value': 0.023},
                                'Q3': {'coefficient': -1.28, 'lower_ci': -3.80, 'upper_ci': 1.24, 'p_value': 0.320},
                                'Q4': {'coefficient': -14.70, 'lower_ci': -18.20, 'upper_ci': -11.20, 'p_value': 0.0001},
                                'Q5': {'coefficient': -8.50, 'lower_ci': -12.10, 'upper_ci': -4.90, 'p_value': 0.008}
                            }
                        }
                    }
                }
            }
        }
    }
    
    # 创建格式化器
    formatter = LMMResultFormatter()
    
    # 测试表格创建
    cancer_types = ['C00_C97']
    table = formatter.create_display_table(mock_results, cancer_types, 'EQI0005_AAMR2006_2010')
    
    print("测试结果表格:")
    print(table.to_string(index=False))
    
    print(f"\n格式化器创建成功!")
    print(f"支持的模型数: {len(formatter.model_order)}")
    print(f"支持的场景数: {len(formatter.scenario_mapping)}")
    print(f"简化的模型名称: {len(formatter.model_name_mapping)}")


if __name__ == "__main__":
    main()