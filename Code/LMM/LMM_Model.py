#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM模型构建模块 - 多层线性回归分析

主要功能:
1. 固定斜率随机截距多层线性回归模型
2. 支持总体EQI分析和RUCC分层EQI分析
3. 包含吸烟率作为协变量控制
4. 州级随机效应处理

模型设计:
- Level 1 (County): 固定效应 - EQI五分位数 + 吸烟率
- Level 2 (State): 随机效应 - 随机截距 (1|State)

输出:
- 回归系数及95%置信区间
- 显著性检验结果
- 模型拟合统计
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional
import warnings

# 统计建模包
try:
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    from statsmodels.stats.outliers_influence import variance_inflation_factor
except ImportError:
    print("请安装statsmodels: pip install statsmodels")
    sys.exit(1)

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 抑制统计建模的警告信息
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

class LMMAnalyzer:
    """LMM多层回归分析器"""
    
    def __init__(self, data_path: str):
        """
        初始化分析器
        
        参数:
            data_path: 数据文件路径
        """
        self.data_path = Path(data_path)
        self.data = None
        
        # 模型变量配置
        self.eqi_variables = [
            'EQI',                    # 总体EQI
            'EQI_air',               # 空气质量
            'EQI_water',             # 水质量
            'EQI_land',              # 土地质量
            'EQI_built',             # 建成环境质量
            'EQI_Sociodemographic'   # 社会人口学质量
        ]
        
        self.rucc_eqi_variables = [
            'RUCC_EQI',                    # RUCC分层总体EQI
            'RUCC_EQI_air',               # RUCC分层空气质量
            'RUCC_EQI_water',             # RUCC分层水质量
            'RUCC_EQI_land',              # RUCC分层土地质量
            'RUCC_EQI_built',             # RUCC分层建成环境质量
            'RUCC_EQI_Sociodemographic'   # RUCC分层社会人口学质量
        ]
        
        # 多EQI域联合模型变量
        self.multi_eqi_variables = [
            'EQI_air',               # 空气质量
            'EQI_water',             # 水质量
            'EQI_land',              # 土地质量
            'EQI_built',             # 建成环境质量
            'EQI_Sociodemographic'   # 社会人口学质量
        ]
        
        # 模型结果存储
        self.model_results = {}
        self.analysis_summary = {}
        
    def load_data(self) -> bool:
        """加载分析数据"""
        logger.info(f"加载分析数据: {self.data_path}")
        
        try:
            self.data = pd.read_csv(self.data_path)
            logger.info(f"数据加载成功: {self.data.shape}")
            
            # 数据类型转换
            self.data['COUNTY_FIPS'] = self.data['COUNTY_FIPS'].astype(str)
            if 'State_FIPS' in self.data.columns:
                self.data['State_FIPS'] = self.data['State_FIPS'].astype(str)
            elif 'State' in self.data.columns:
                # 如果没有State_FIPS但有State列，使用State作为分组变量
                self.data['State_FIPS'] = self.data['State'].astype(str)
            
            # 确保EQI变量为分类变量
            for var in self.eqi_variables + self.rucc_eqi_variables + self.multi_eqi_variables:
                if var in self.data.columns:
                    self.data[f'{var}_factor'] = pd.Categorical(
                        self.data[var], 
                        categories=[1, 2, 3, 4, 5], 
                        ordered=True
                    )
            
            logger.info("数据预处理完成")
            return True
            
        except Exception as e:
            logger.error(f"数据加载失败: {e}")
            return False
    
    def get_analysis_data(self, scenario: str, cancer_type: str, rucc_filter: Optional[int] = None) -> pd.DataFrame:
        """
        获取特定分析的数据子集
        
        参数:
            scenario: 分析场景 (如 'EQI0005_AAMR2006_2010')
            cancer_type: 癌症类型 (如 'C00_C97')
            rucc_filter: RUCC过滤条件 (1-4, None表示不过滤)
            
        返回:
            分析用数据框
        """
        # 解析scenario格式 (如 'EQI0005_AAMR2006_2010')
        parts = scenario.split('_')
        eqi_period = f"{parts[0][3:7]}_{parts[0][7:11]}"  # EQI0005 -> 2000_2005
        aamr_period = f"{parts[1][4:8]}_{parts[1][8:12]}"  # AAMR2006 -> 2006_2010
        
        # 筛选条件
        mask = (self.data['EQI_Period'] == eqi_period) & (self.data['Time_Period'] == aamr_period) & (self.data['Cancer_Type'] == cancer_type)
        
        if rucc_filter is not None:
            mask = mask & (self.data['RUCC'] == rucc_filter)
        
        analysis_data = self.data[mask].copy()
        
        # 检查数据量
        if len(analysis_data) < 50:
            logger.warning(f"数据量过小: {len(analysis_data)} 记录")
        
        # 删除关键变量缺失的记录
        before_count = len(analysis_data)
        required_cols = ['AAMR', 'Smoking_Rate', 'State_FIPS']
        # 检查State_FIPS是否存在，如果不存在就用State
        if 'State_FIPS' not in analysis_data.columns and 'State' in analysis_data.columns:
            required_cols = ['AAMR', 'Smoking_Rate', 'State']
        analysis_data = analysis_data.dropna(subset=required_cols)
        after_count = len(analysis_data)
        
        if before_count > after_count:
            logger.info(f"删除 {before_count - after_count} 个缺失记录")
        
        return analysis_data
    
    def fit_mixed_model(self, data: pd.DataFrame, eqi_var: str, model_name: str) -> Dict:
        """
        拟合混合效应模型（单一EQI变量模型）
        
        参数:
            data: 分析数据
            eqi_var: EQI变量名
            model_name: 模型名称
            
        返回:
            模型结果字典
        """
        try:
            # 检查EQI变量有效性
            eqi_factor_var = f'{eqi_var}_factor'
            if eqi_factor_var not in data.columns:
                logger.error(f"EQI因子变量不存在: {eqi_factor_var}")
                return None
            
            # 检查EQI分布
            eqi_dist = data[eqi_var].value_counts().sort_index()
            if len(eqi_dist) < 3:
                logger.warning(f"EQI变量 {eqi_var} 分类不足: {dict(eqi_dist)}")
                return None
            
            # 构建公式 - 固定斜率随机截距模型
            # Level 1: EQI五分位数 + 吸烟率 (固定效应)
            # Level 2: 州随机截距 (随机效应)
            formula = f"AAMR ~ C({eqi_factor_var}, Treatment(reference=1)) + Smoking_Rate"
            
            # 拟合混合效应模型
            logger.debug(f"拟合模型: {formula}")
            model = smf.mixedlm(formula, data, groups=data['State_FIPS'])
            result = model.fit(method='lbfgs', maxiter=1000)
            
            # 提取系数和置信区间
            coefficients = self._extract_coefficients(result, eqi_var)
            
            # 模型诊断信息
            model_info = {
                'model_name': model_name,
                'eqi_variable': eqi_var,
                'sample_size': len(data),
                'counties': data['COUNTY_FIPS'].nunique(),
                'states': data['State_FIPS'].nunique(),
                'aic': result.aic,
                'bic': result.bic,
                'log_likelihood': result.llf,
                'random_effect_var': result.scale if hasattr(result, 'scale') else None,
                'converged': result.converged if hasattr(result, 'converged') else True
            }
            
            return {
                'coefficients': coefficients,
                'model_info': model_info,
                'raw_result': result
            }
            
        except Exception as e:
            logger.error(f"模型拟合失败 ({model_name}): {e}")
            return None
    
    def fit_multi_eqi_model(self, data: pd.DataFrame, model_name: str) -> Dict:
        """
        拟合多EQI域联合模型
        
        参数:
            data: 分析数据
            model_name: 模型名称
            
        返回:
            模型结果字典
        """
        try:
            # 检查所有需要的变量是否存在
            required_vars = [f'{var}_factor' for var in self.multi_eqi_variables]
            missing_vars = [var for var in required_vars if var not in data.columns]
            if missing_vars:
                logger.error(f"缺少必要的变量: {missing_vars}")
                return None
            
            # 检查变量分布
            for var in self.multi_eqi_variables:
                var_dist = data[var].value_counts().sort_index()
                if len(var_dist) < 3:
                    logger.warning(f"EQI变量 {var} 分类不足: {dict(var_dist)}")
            
            # 构建多EQI域联合模型公式
            # 包含所有5个EQI细分域的五分位数 + 吸烟率 + 州随机截距
            eqi_factors = [f'C({var}_factor, Treatment(reference=1))' for var in self.multi_eqi_variables]
            formula = f"AAMR ~ {' + '.join(eqi_factors)} + Smoking_Rate"
            
            # 拟合混合效应模型
            logger.debug(f"拟合多EQI模型: {formula}")
            model = smf.mixedlm(formula, data, groups=data['State_FIPS'])
            result = model.fit(method='lbfgs', maxiter=1000)
            
            # 提取系数和置信区间
            coefficients = self._extract_multi_eqi_coefficients(result)
            
            # 模型诊断信息
            model_info = {
                'model_name': model_name,
                'eqi_variables': self.multi_eqi_variables,
                'sample_size': len(data),
                'counties': data['COUNTY_FIPS'].nunique(),
                'states': data['State_FIPS'].nunique(),
                'aic': result.aic,
                'bic': result.bic,
                'log_likelihood': result.llf,
                'random_effect_var': result.scale if hasattr(result, 'scale') else None,
                'converged': result.converged if hasattr(result, 'converged') else True
            }
            
            return {
                'coefficients': coefficients,
                'model_info': model_info,
                'raw_result': result
            }
            
        except Exception as e:
            logger.error(f"多EQI模型拟合失败 ({model_name}): {e}")
            return None
    
    def _extract_coefficients(self, result, eqi_var: str) -> Dict:
        """提取模型系数和置信区间"""
        coefficients = {}
        
        # 获取系数表
        coef_table = result.summary().tables[1]
        params = result.params
        conf_int = result.conf_int()
        pvalues = result.pvalues
        
        # 提取EQI系数 (Q2-Q5, Q1为参照组)
        for q in range(2, 6):
            param_name = f'C({eqi_var}_factor, Treatment(reference=1))[T.{q}]'
            
            if param_name in params.index:
                coef = params[param_name]
                lower_ci = conf_int.loc[param_name, 0]
                upper_ci = conf_int.loc[param_name, 1] 
                p_value = pvalues[param_name]
                
                coefficients[f'Q{q}'] = {
                    'coefficient': coef,
                    'lower_ci': lower_ci,
                    'upper_ci': upper_ci,
                    'p_value': p_value,
                    'significant': p_value < 0.05
                }
            else:
                # 如果该分位数不存在，设为缺失
                coefficients[f'Q{q}'] = {
                    'coefficient': np.nan,
                    'lower_ci': np.nan,
                    'upper_ci': np.nan,
                    'p_value': np.nan,
                    'significant': False
                }
        
        # Q1始终为参照组 (coefficient = 0)
        coefficients['Q1'] = {
            'coefficient': 0.0,
            'lower_ci': 0.0,
            'upper_ci': 0.0,
            'p_value': np.nan,
            'significant': False
        }
        
        # 吸烟率系数
        smoking_param = 'Smoking_Rate'
        if smoking_param in params.index:
            coef = params[smoking_param]
            lower_ci = conf_int.loc[smoking_param, 0]
            upper_ci = conf_int.loc[smoking_param, 1]
            p_value = pvalues[smoking_param]
            
            coefficients['Smoking_Rate'] = {
                'coefficient': coef,
                'lower_ci': lower_ci,
                'upper_ci': upper_ci,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        
        return coefficients
    
    def _extract_multi_eqi_coefficients(self, result) -> Dict:
        """提取多EQI域联合模型系数和置信区间"""
        coefficients = {}
        
        # 获取系数表
        params = result.params
        conf_int = result.conf_int()
        pvalues = result.pvalues
        
        # 为每个EQI域提取系数 (Q2-Q5, Q1为参照组)
        for eqi_var in self.multi_eqi_variables:
            for q in range(1, 6):  # 包括Q1作为参照组
                if q == 1:
                    # Q1始终为参照组 (coefficient = 0)
                    coefficients[f'{eqi_var}_Q{q}'] = {
                        'coefficient': 0.0,
                        'lower_ci': 0.0,
                        'upper_ci': 0.0,
                        'p_value': np.nan,
                        'significant': False
                    }
                else:
                    param_name = f'C({eqi_var}_factor, Treatment(reference=1))[T.{q}]'
                    
                    if param_name in params.index:
                        coef = params[param_name]
                        lower_ci = conf_int.loc[param_name, 0]
                        upper_ci = conf_int.loc[param_name, 1] 
                        p_value = pvalues[param_name]
                        
                        coefficients[f'{eqi_var}_Q{q}'] = {
                            'coefficient': coef,
                            'lower_ci': lower_ci,
                            'upper_ci': upper_ci,
                            'p_value': p_value,
                            'significant': p_value < 0.05
                        }
                    else:
                        # 如果该分位数不存在，设为缺失
                        coefficients[f'{eqi_var}_Q{q}'] = {
                            'coefficient': np.nan,
                            'lower_ci': np.nan,
                            'upper_ci': np.nan,
                            'p_value': np.nan,
                            'significant': False
                        }
        
        # 吸烟率系数
        smoking_param = 'Smoking_Rate'
        if smoking_param in params.index:
            coef = params[smoking_param]
            lower_ci = conf_int.loc[smoking_param, 0]
            upper_ci = conf_int.loc[smoking_param, 1]
            p_value = pvalues[smoking_param]
            
            coefficients['Smoking_Rate'] = {
                'coefficient': coef,
                'lower_ci': lower_ci,
                'upper_ci': upper_ci,
                'p_value': p_value,
                'significant': p_value < 0.05
            }
        
        return coefficients
    
    def run_scenario_analysis(self, scenario: str, cancer_types: List[str]) -> Dict:
        """
        运行单个场景的完整分析
        
        参数:
            scenario: 分析场景
            cancer_types: 癌症类型列表
            
        返回:
            场景分析结果
        """
        logger.info(f"=== 开始场景分析: {scenario} ===")
        
        scenario_results = {
            'scenario': scenario,
            'cancer_results': {}
        }
        
        for cancer_type in cancer_types:
            logger.info(f"分析癌症类型: {cancer_type}")
            
            # 获取分析数据
            analysis_data = self.get_analysis_data(scenario, cancer_type)
            
            if len(analysis_data) < 50:
                logger.warning(f"跳过 {cancer_type}: 数据量不足 ({len(analysis_data)} 记录)")
                continue
            
            cancer_results = {}
            
            # 1. 总体EQI分析
            logger.info("  进行总体EQI分析...")
            for eqi_var in self.eqi_variables:
                if eqi_var in analysis_data.columns:
                    model_name = f"{scenario}_{cancer_type}_{eqi_var}"
                    result = self.fit_mixed_model(analysis_data, eqi_var, model_name)
                    if result:
                        cancer_results[eqi_var] = result
            
            # 2. 多EQI域联合模型分析
            logger.info("  进行多EQI域联合模型分析...")
            # 检查是否所有需要的变量都存在
            required_multi_vars = [var for var in self.multi_eqi_variables if var in analysis_data.columns]
            if len(required_multi_vars) == len(self.multi_eqi_variables):
                model_name = f"{scenario}_{cancer_type}_Multi_EQI"
                result = self.fit_multi_eqi_model(analysis_data, model_name)
                if result:
                    cancer_results['Multi_EQI'] = result
            else:
                missing_vars = set(self.multi_eqi_variables) - set(required_multi_vars)
                logger.warning(f"跳过多EQI模型: 缺少变量 {missing_vars}")
            
            # 3. RUCC分层EQI分析 (只分析RUCC1-RUCC4)
            logger.info("  进行RUCC分层分析...")
            for rucc in [1, 2, 3, 4]:
                rucc_data = self.get_analysis_data(scenario, cancer_type, rucc_filter=rucc)
                
                if len(rucc_data) < 30:  # 分层数据的最小样本量要求更低
                    logger.debug(f"跳过 RUCC{rucc}: 数据量不足 ({len(rucc_data)} 记录)")
                    continue
                
                for rucc_eqi_var in self.rucc_eqi_variables:
                    if rucc_eqi_var in rucc_data.columns:
                        model_name = f"{scenario}_{cancer_type}_RUCC{rucc}_{rucc_eqi_var}"
                        result = self.fit_mixed_model(rucc_data, rucc_eqi_var, model_name)
                        if result:
                            key = f"RUCC{rucc}_{rucc_eqi_var}"
                            cancer_results[key] = result
            
            if cancer_results:
                scenario_results['cancer_results'][cancer_type] = cancer_results
                logger.info(f"  {cancer_type} 分析完成: {len(cancer_results)} 个模型")
            else:
                logger.warning(f"  {cancer_type} 无有效模型")
        
        logger.info(f"场景 {scenario} 分析完成")
        return scenario_results
    
    def run_full_analysis(self, scenarios: List[str], cancer_types: List[str]) -> Dict:
        """
        运行完整的LMM分析
        
        参数:
            scenarios: 分析场景列表
            cancer_types: 癌症类型列表
            
        返回:
            完整分析结果
        """
        logger.info("=== 开始完整LMM分析 ===")
        
        all_results = {
            'analysis_info': {
                'total_scenarios': len(scenarios),
                'total_cancer_types': len(cancer_types),
                'data_source': str(self.data_path),
                'analysis_timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
            },
            'scenario_results': {}
        }
        
        for scenario in scenarios:
            scenario_result = self.run_scenario_analysis(scenario, cancer_types)
            all_results['scenario_results'][scenario] = scenario_result
        
        # 保存结果到实例变量
        self.model_results = all_results
        
        # 生成分析摘要
        self._generate_analysis_summary()
        
        logger.info("=== LMM分析完成 ===")
        return all_results
    
    def _generate_analysis_summary(self):
        """生成分析摘要统计"""
        summary = {
            'total_models': 0,
            'successful_models': 0,
            'failed_models': 0,
            'scenarios_analyzed': 0,
            'cancer_types_analyzed': set(),
            'model_convergence': {},
            'sample_sizes': []
        }
        
        if 'scenario_results' in self.model_results:
            summary['scenarios_analyzed'] = len(self.model_results['scenario_results'])
            
            for scenario, scenario_data in self.model_results['scenario_results'].items():
                if 'cancer_results' in scenario_data:
                    for cancer, cancer_models in scenario_data['cancer_results'].items():
                        summary['cancer_types_analyzed'].add(cancer)
                        
                        for model_name, model_result in cancer_models.items():
                            summary['total_models'] += 1
                            if model_result and 'model_info' in model_result:
                                summary['successful_models'] += 1
                                
                                # 收集样本量信息
                                if 'sample_size' in model_result['model_info']:
                                    summary['sample_sizes'].append(model_result['model_info']['sample_size'])
                                
                                # 收集收敛信息
                                if 'converged' in model_result['model_info']:
                                    converged = model_result['model_info']['converged']
                                    if converged not in summary['model_convergence']:
                                        summary['model_convergence'][converged] = 0
                                    summary['model_convergence'][converged] += 1
                            else:
                                summary['failed_models'] += 1
        
        # 转换集合为列表
        summary['cancer_types_analyzed'] = list(summary['cancer_types_analyzed'])
        
        # 计算样本量统计
        if summary['sample_sizes']:
            summary['sample_size_stats'] = {
                'min': min(summary['sample_sizes']),
                'max': max(summary['sample_sizes']),
                'mean': np.mean(summary['sample_sizes']),
                'median': np.median(summary['sample_sizes'])
            }
        
        self.analysis_summary = summary
        
        # 记录摘要信息
        logger.info("=== 分析摘要 ===")
        logger.info(f"总模型数: {summary['total_models']}")
        logger.info(f"成功模型: {summary['successful_models']}")
        logger.info(f"失败模型: {summary['failed_models']}")
        logger.info(f"分析场景: {summary['scenarios_analyzed']}")
        logger.info(f"癌症类型: {len(summary['cancer_types_analyzed'])}")
        
        if 'sample_size_stats' in summary:
            stats = summary['sample_size_stats']
            logger.info(f"样本量: 最小={stats['min']}, 最大={stats['max']}, 均值={stats['mean']:.0f}")
    
    def get_model_results(self) -> Dict:
        """获取模型结果"""
        return self.model_results
    
    def get_analysis_summary(self) -> Dict:
        """获取分析摘要"""
        return self.analysis_summary


def main():
    """主函数 - 演示用法"""
    print("=== LMM多层回归分析 ===")
    
    # 数据路径
    data_path = "Data/df/EQI_LMM_Delete_df.csv"
    
    # 创建分析器
    analyzer = LMMAnalyzer(data_path)
    
    # 加载数据
    if not analyzer.load_data():
        print("数据加载失败!")
        return
    
    # 定义分析参数
    scenarios = [
        'EQI0005_AAMR2006_2010',
        'EQI0005_AAMR2011_2015', 
        'EQI0610_AAMR2011_2015',
        'EQI0610_AAMR2016_2020'
    ]
    
    # 选择主要癌症类型进行测试
    cancer_types = ['C00_C97']  # 先只测试总癌症
    
    # 运行分析
    results = analyzer.run_full_analysis(scenarios, cancer_types)
    
    # 显示结果
    print(f"\n分析完成!")
    print(f"分析摘要: {analyzer.get_analysis_summary()}")


if __name__ == "__main__":
    main()