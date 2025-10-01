#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM结果提取和输出模块 (简化版)

主要功能:
1. 从模型结果提取系数和置信区间
2. 生成标准化结果表格 (CSV格式)
3. 支持主效应模型和领域探索模型的结果处理

输出格式:
- 详细结果表: 长格式CSV，包含所有统计信息
- 汇总结果表: 宽格式CSV，便于比较
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
import logging
from datetime import datetime
import warnings

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LMMResultProcessor:
    """LMM结果处理器 - 统一的结果提取和输出接口"""
    
    def __init__(self, project_root: Path = None, analysis_type: str = "standard"):
        """
        初始化结果处理器
        
        参数:
            project_root: 项目根目录路径
            analysis_type: 分析类型
        """
        if project_root is None:
            project_root = Path(__file__).resolve().parents[2]
        
        self.project_root = project_root
        self.analysis_type = analysis_type
        
        # 设置输出目录
        self.output_dir = project_root / "Result" / "EQI_GLMM"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储所有结果
        self.all_results = {}
        
        # 显著性标记配置
        self.significance_markers = {
            0.001: "***",
            0.01: "**", 
            0.05: "*"
        }
    
    def save_results(self, results: Dict[str, Any], analysis_type: str = "Overall", 
                    timestamp: str = None) -> None:
        """
        保存分析结果到CSV文件
        
        参数:
            results: 模型结果字典
            analysis_type: 分析类型描述
            timestamp: 时间戳
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        try:
            # 转换结果为DataFrame
            detailed_df = self._create_detailed_results_table(results)
            
            if len(detailed_df) > 0:
                # 保存详细结果
                detailed_file = self.output_dir / f"Detailed_Results_{analysis_type}_{timestamp}.csv"
                detailed_df.to_csv(detailed_file, index=False)
                logger.info(f"保存详细结果: {detailed_file}")
                
                # 保存汇总结果 
                summary_df = self._create_summary_results_table(detailed_df)
                if len(summary_df) > 0:
                    summary_file = self.output_dir / f"Summary_Results_{analysis_type}_{timestamp}.csv" 
                    summary_df.to_csv(summary_file, index=False)
                    logger.info(f"保存汇总结果: {summary_file}")
            
        except Exception as e:
            logger.error(f"保存结果失败: {e}")
    
    def _create_detailed_results_table(self, results: Dict[str, Any]) -> pd.DataFrame:
        """
        创建详细结果表
        
        参数:
            results: 模型结果字典
            
        返回:
            详细结果DataFrame
        """
        try:
            all_rows = []
            
            for model_key, model_result in results.items():
                # 解析模型类型和癌症类型
                if 'main_effect_' in model_key:
                    model_type = 'Main Effect'
                    cancer_type = model_key.replace('main_effect_', '')
                elif 'domain_exploration_' in model_key:
                    model_type = 'Domain Exploration'
                    cancer_type = model_key.replace('domain_exploration_', '')
                else:
                    continue
                
                # 提取EQI系数
                if 'eqi_coefficients' in model_result:
                    for param_name, coef_info in model_result['eqi_coefficients'].items():
                        
                        # 解析EQI变量和分位数
                        eqi_var, quartile = self._parse_eqi_parameter(param_name)
                        
                        row = {
                            'Model_Type': model_type,
                            'Cancer_Type': cancer_type,
                            'EQI_Variable': eqi_var,
                            'Quartile': quartile,
                            'Parameter': param_name,
                            'Coefficient': coef_info.get('coefficient', np.nan),
                            'PValue': coef_info.get('pvalue', np.nan),
                            'CI_Lower': coef_info.get('ci_lower', np.nan),
                            'CI_Upper': coef_info.get('ci_upper', np.nan),
                            'Significant': coef_info.get('significant', False)
                        }
                        
                        # 添加模型信息
                        if 'model_info' in model_result:
                            row.update({
                                'N_Obs': model_result['model_info'].get('n_obs', np.nan),
                                'AIC': model_result['model_info'].get('aic', np.nan),
                                'BIC': model_result['model_info'].get('bic', np.nan)
                            })
                        
                        all_rows.append(row)
            
            return pd.DataFrame(all_rows)
            
        except Exception as e:
            logger.error(f"创建详细结果表失败: {e}")
            return pd.DataFrame()
    
    def _create_summary_results_table(self, detailed_df: pd.DataFrame) -> pd.DataFrame:
        """
        创建汇总结果表 (宽格式)
        
        参数:
            detailed_df: 详细结果DataFrame
            
        返回:
            汇总结果DataFrame
        """
        try:
            if len(detailed_df) == 0:
                return pd.DataFrame()
            
            # 创建汇总表 - 以癌症类型为行，EQI变量和分位数为列
            summary_rows = []
            
            for (model_type, cancer_type), group in detailed_df.groupby(['Model_Type', 'Cancer_Type']):
                
                row = {
                    'Model_Type': model_type,
                    'Cancer_Type': cancer_type
                }
                
                # 为每个EQI变量和分位数创建列
                for _, result in group.iterrows():
                    eqi_var = result['EQI_Variable']
                    quartile = result['Quartile']
                    
                    col_prefix = f"{eqi_var}_{quartile}"
                    
                    row[f"{col_prefix}_Coef"] = result['Coefficient']
                    row[f"{col_prefix}_PVal"] = result['PValue']
                    row[f"{col_prefix}_Sig"] = "***" if result['PValue'] < 0.001 else "**" if result['PValue'] < 0.01 else "*" if result['PValue'] < 0.05 else ""
                
                # 添加模型拟合信息
                if 'N_Obs' in group.columns:
                    row['N_Obs'] = group['N_Obs'].iloc[0]
                if 'AIC' in group.columns:
                    row['AIC'] = group['AIC'].iloc[0]
                
                summary_rows.append(row)
            
            return pd.DataFrame(summary_rows)
            
        except Exception as e:
            logger.error(f"创建汇总结果表失败: {e}")
            return pd.DataFrame()
    
    def _parse_eqi_parameter(self, param_name: str) -> Tuple[str, str]:
        """
        解析EQI参数名称，提取变量和分位数
        
        参数:
            param_name: 参数名称，如 "C(EQI, Treatment(1))[T.2]"
            
        返回:
            (eqi_variable, quartile) 元组
        """
        try:
            # 示例: "C(EQI_air, Treatment(1))[T.2]"
            # 注意：必须先匹配长名称，再匹配短名称，避免EQI被提前匹配
            if 'EQI_Sociodemographic' in param_name:
                eqi_var = 'EQI_Sociodemographic'
            elif 'EQI_air' in param_name:
                eqi_var = 'EQI_air'
            elif 'EQI_water' in param_name:
                eqi_var = 'EQI_water'
            elif 'EQI_land' in param_name:
                eqi_var = 'EQI_land'
            elif 'EQI_built' in param_name:
                eqi_var = 'EQI_built'
            elif 'EQI' in param_name:
                eqi_var = 'EQI'
            else:
                eqi_var = 'Unknown'
            
            # 提取分位数 - 修改匹配规则处理 "[T.2]" 和 "[T.2.0]" 格式
            if '[T.2' in param_name:  # 匹配 [T.2] 或 [T.2.0]
                quartile = 'Q2'
            elif '[T.3' in param_name:  # 匹配 [T.3] 或 [T.3.0]
                quartile = 'Q3' 
            elif '[T.4' in param_name:  # 匹配 [T.4] 或 [T.4.0]
                quartile = 'Q4'
            elif '[T.5' in param_name:  # 匹配 [T.5] 或 [T.5.0]
                quartile = 'Q5'
            else:
                quartile = 'Unknown'
            
            return eqi_var, quartile
            
        except Exception as e:
            logger.error(f"解析参数名称失败: {param_name}, {e}")
            return 'Unknown', 'Unknown'
    
    def generate_analysis_report(self, timestamp: str = None) -> str:
        """
        生成简单的分析报告
        
        参数:
            timestamp: 时间戳
            
        返回:
            报告文件路径
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        try:
            report_lines = [
                "LMM分析报告",
                "=" * 50,
                f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                f"分析类型: {self.analysis_type}",
                "",
                "分析内容:",
                "1. 主效应模型 (Main Effect Model)",
                "   - 总体EQI vs 癌症死亡率",
                "   - 模型: AAMR ~ C(EQI) + 控制变量 + (1|State)",
                "",
                "2. 领域探索模型 (Domain Exploration Model)", 
                "   - 五大环境领域 vs 癌症死亡率",
                "   - 模型: AAMR ~ C(EQI_air) + C(EQI_water) + ... + 控制变量 + (1|State)",
                "",
                "结果文件:",
                "- Detailed_Results_*.csv: 详细结果表",
                "- Summary_Results_*.csv: 汇总结果表",
                "",
                "注: 详见CSV文件获取具体统计结果"
            ]
            
            report_content = "\n".join(report_lines)
            
            # 保存报告
            report_file = self.output_dir / f"Analysis_Report_{timestamp}.txt"
            with report_file.open("w", encoding="utf-8") as f:
                f.write(report_content)
            
            logger.info(f"分析报告已生成: {report_file}")
            return str(report_file)
            
        except Exception as e:
            logger.error(f"生成报告失败: {e}")
            return ""
    
    def create_compact_table(self, results: Dict[str, Any], eqi_period: str, aamr_period: str, 
                            timestamp: str = None) -> pd.DataFrame:
        """
        创建简洁的结果表格
        
        参数:
            results: 模型结果字典
            eqi_period: EQI时间段 ("0005" 或 "0610")
            aamr_period: AAMR时间段 ("2006_2010" 或 "2016_2020")
            timestamp: 时间戳
            
        返回:
            简洁的结果DataFrame
        """
        try:
            if timestamp is None:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            compact_rows = []
            
            for model_key, model_result in results.items():
                # 解析模型信息
                if 'main_effect_' in model_key:
                    analysis_model = 'EQI_Main_Effect'
                    cancer_type = model_key.replace('main_effect_', '')
                elif 'domain_exploration_' in model_key:
                    analysis_model = 'EQI_Domain_Exploration'
                    cancer_type = model_key.replace('domain_exploration_', '')
                else:
                    continue
                
                if 'eqi_coefficients' in model_result:
                    # 创建行数据，先初始化基本信息
                    row = {
                        'Cancer_Type': cancer_type,
                        'Model': analysis_model,
                        'Q1': '0.00 (ref)'  # Q1作为参考组
                    }
                    
                    # 处理Q2-Q5的系数
                    for param_name, coef_info in model_result['eqi_coefficients'].items():
                        # 只处理总体EQI参数（主效应模型）或所有EQI参数（领域探索模型）
                        eqi_var, quartile = self._parse_eqi_parameter(param_name)
                        
                        if quartile in ['Q2', 'Q3', 'Q4', 'Q5']:
                            # 格式化系数和置信区间
                            coef = coef_info['coefficient']
                            ci_lower = coef_info['ci_lower']
                            ci_upper = coef_info['ci_upper']
                            pvalue = coef_info['pvalue']
                            
                            # 确定显著性标记
                            sig_marker = ""
                            for threshold, marker in sorted(self.significance_markers.items()):
                                if pvalue < threshold:
                                    sig_marker = marker
                                    break
                            
                            # 格式化为 "系数(95%CI下限, 95%CI上限)显著性标记"
                            formatted_value = f"{coef:.2f}({ci_lower:.2f}, {ci_upper:.2f}){sig_marker}"
                            
                            # 对于领域探索模型，在分位数前加上EQI变量名
                            if analysis_model == 'EQI_Domain_Exploration' and eqi_var != 'EQI':
                                col_name = f"{eqi_var}_{quartile}"
                            else:
                                col_name = quartile
                            
                            row[col_name] = formatted_value
                    
                    # 确保基本的Q2-Q5列都存在（对于主效应模型）
                    if analysis_model == 'EQI_Main_Effect':
                        for q in ['Q2', 'Q3', 'Q4', 'Q5']:
                            if q not in row:
                                row[q] = "N/A"
                    
                    compact_rows.append(row)
            
            if not compact_rows:
                logger.warning("没有找到适合的模型结果用于生成简洁表格")
                return pd.DataFrame()
            
            # 创建DataFrame并排序
            compact_df = pd.DataFrame(compact_rows)
            
            # 根据模型类型确定列顺序
            base_cols = ['Cancer_Type', 'Model', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5']
            all_cols = list(compact_df.columns)
            
            # 保持基本列在前，其他列在后
            ordered_cols = [col for col in base_cols if col in all_cols]
            ordered_cols += [col for col in all_cols if col not in base_cols]
            
            compact_df = compact_df[ordered_cols]
            compact_df = compact_df.sort_values(['Cancer_Type', 'Model'])
            
            # 保存表格
            filename = f"Compact_EQI{eqi_period}_AAMR{aamr_period}_Results_{timestamp}.csv"
            filepath = self.output_dir / filename
            compact_df.to_csv(filepath, index=False)
            logger.info(f"简洁表格已保存: {filepath}")
            
            return compact_df
            
        except Exception as e:
            logger.error(f"创建简洁表格失败: {e}")
            return pd.DataFrame()
    
    def save_all_compact_tables(self, all_results: Dict[str, Dict[str, Any]], timestamp: str = None) -> List[str]:
        """
        保存所有4个简洁表格（4种年份组合）
        
        参数:
            all_results: 包含所有年份组合分析结果的字典
            timestamp: 时间戳
            
        返回:
            保存的文件路径列表
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        saved_files = []
        
        # 生成4个表格 - 4种年份组合
        period_combinations = [
            ("0005", "2006_2010", "EQI2000-2005 vs AAMR2006-2010"),
            ("0005", "2016_2020", "EQI2000-2005 vs AAMR2016-2020"),
            ("0610", "2006_2010", "EQI2006-2010 vs AAMR2006-2010"), 
            ("0610", "2016_2020", "EQI2006-2010 vs AAMR2016-2020")
        ]
        
        for eqi_period, aamr_period, description in period_combinations:
            # 构造结果键
            result_key = f"EQI{eqi_period}_AAMR{aamr_period}"
            
            # 获取对应的结果数据
            if result_key in all_results:
                results = all_results[result_key]
                
                compact_df = self.create_compact_table(results, eqi_period, aamr_period, timestamp)
                
                if not compact_df.empty:
                    filename = f"Compact_EQI{eqi_period}_AAMR{aamr_period}_Results_{timestamp}.csv"
                    filepath = self.output_dir / filename
                    saved_files.append(str(filepath))
                    logger.info(f"✅ {description}: {filename}")
                else:
                    logger.warning(f"❌ {description}: 无数据")
            else:
                logger.warning(f"❌ {description}: 结果数据不存在 (键: {result_key})")
        
        return saved_files