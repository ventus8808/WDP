#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM模型拟合与计算模块 (简化版)

主要功能:
1. 主效应模型 (The Main Effect Model) - 总体EQI vs 癌症死亡率
2. 领域探索模型 (The Domain Exploration Model) - 五大环境领域 vs 癌症死亡率

模型公式:
- 主效应模型: AAMR ~ C(EQI) + C(Census_Region) + C(Urbanization_Type) + (1|State)
- 领域探索模型: AAMR ~ C(EQI_air) + C(EQI_water) + C(EQI_land) + C(EQI_built) + C(EQI_sociod) 
                      + C(Census_Region) + C(Urbanization_Type) + (1|State)

输出: 系数估计和统计检验结果
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import warnings
import logging
from typing import Dict, List, Optional, Any

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 忽略statsmodels警告
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', message='.*convergence.*')

class LMMAnalyzer:
    """线性混合效应模型分析器 - 实现主效应模型和领域探索模型"""
    
    def __init__(self):
        """初始化LMM分析器"""
        self.models = {}  # 存储拟合的模型
        self.results = {}  # 存储模型结果
        
        # 控制变量
        self.control_vars = ['Census_Region', 'Urbanization_Type']
        
        # 五大环境领域
        self.environmental_domains = ['air', 'water', 'land', 'built', 'Sociodemographic']
        
    def fit_main_effect_model(self, data: pd.DataFrame, cancer_type: str) -> Optional[Any]:
        """
        拟合主效应模型 (The Main Effect Model)
        
        研究问题: "总体的、累积的环境质量与特定癌症的死亡率是否存在关联？"
        
        模型方程: AAMR ~ C(EQI) + C(Census_Region) + C(Urbanization_Type) + (1 | State)
        
        参数:
            data: 分析数据
            cancer_type: 癌症类型
            
        返回:
            拟合的模型对象
        """
        try:
            # 检查必需的变量 - 数据是长格式，癌症类型在Cancer_Type列中
            required_vars = ['AAMR', 'Cancer_Type', 'EQI'] + self.control_vars + ['State_FIPS']
            missing_vars = [var for var in required_vars if var not in data.columns]
            if missing_vars:
                logger.error(f"缺少必需变量: {missing_vars}")
                return None
            
            # 筛选特定癌症类型的数据
            cancer_data = data[data['Cancer_Type'] == cancer_type].copy()
            if len(cancer_data) == 0:
                logger.error(f"没有找到癌症类型 {cancer_type} 的数据")
                return None
            
            # 准备数据 - 删除缺失值
            model_data = cancer_data[required_vars].dropna().copy()
            
            if len(model_data) < 50:  # 最小样本量检查
                logger.warning(f"样本量过小: {len(model_data)}")
                return None
            
            # 构建模型公式 - 使用AAMR作为因变量，EQI=1作为参考组
            formula = f"AAMR ~ C(EQI, Treatment(1)) + C(Census_Region) + C(Urbanization_Type)"
            
            # 拟合混合效应模型
            print(f"    拟合主效应模型: {cancer_type}")
            model = smf.mixedlm(formula, model_data, groups=model_data['State_FIPS']).fit()
            
            model_key = f"main_effect_{cancer_type}"
            self.models[model_key] = model
            self.results[model_key] = self._extract_coefficients(model)
            
            return model
            
        except Exception as e:
            logger.error(f"主效应模型拟合失败 ({cancer_type}): {e}")
            return None
    
    def fit_domain_exploration_model(self, data: pd.DataFrame, cancer_type: str) -> Optional[Any]:
        """
        拟合领域探索模型 (The Domain Exploration Model)
        
        研究问题: "如果总体环境质量有关联，那么是空气、水、土地、建筑还是社会人口学
                   这五个领域中的哪一个或哪几个在起主导作用？"
        
        模型方程: AAMR ~ C(EQI_air) + C(EQI_water) + C(EQI_land) + C(EQI_built) + C(EQI_sociod) 
                      + C(Census_Region) + C(Urbanization_Type) + (1 | State)
        
        参数:
            data: 分析数据
            cancer_type: 癌症类型
            
        返回:
            拟合的模型对象
        """
        try:
            # 检查必需的变量 - 数据是长格式
            eqi_domains = [f'EQI_{domain}' for domain in self.environmental_domains]
            required_vars = ['AAMR', 'Cancer_Type'] + eqi_domains + self.control_vars + ['State_FIPS']
            
            missing_vars = [var for var in required_vars if var not in data.columns]
            if missing_vars:
                logger.error(f"缺少必需变量: {missing_vars}")
                return None
            
            # 筛选特定癌症类型的数据
            cancer_data = data[data['Cancer_Type'] == cancer_type].copy()
            if len(cancer_data) == 0:
                logger.error(f"没有找到癌症类型 {cancer_type} 的数据")
                return None
            
            # 准备数据 - 删除缺失值
            model_data = cancer_data[required_vars].dropna().copy()
            
            if len(model_data) < 50:  # 最小样本量检查
                logger.warning(f"样本量过小: {len(model_data)}")
                return None
            
            # 构建模型公式 - 所有五个领域同时进入模型，都以1作为参考组
            domain_terms = []
            for domain in self.environmental_domains:
                domain_terms.append(f"C(EQI_{domain}, Treatment(1))")
            
            formula = ("AAMR ~ " + 
                      " + ".join(domain_terms) + 
                      " + C(Census_Region) + C(Urbanization_Type)")
            
            # 拟合混合效应模型
            print(f"    拟合领域探索模型: {cancer_type}")
            model = smf.mixedlm(formula, model_data, groups=model_data['State_FIPS']).fit()
            
            model_key = f"domain_exploration_{cancer_type}"
            self.models[model_key] = model
            self.results[model_key] = self._extract_coefficients(model)
            
            return model
            
        except Exception as e:
            logger.error(f"领域探索模型拟合失败 ({cancer_type}): {e}")
            return None
    
    def fit_multiple_models(self, data: pd.DataFrame, cancer_types: List[str]) -> Dict[str, Any]:
        """
        为多种癌症类型拟合两种模型
        
        参数:
            data: 分析数据
            cancer_types: 癌症类型列表
            
        返回:
            所有模型结果的字典
        """
        all_results = {}
        
        for cancer_type in cancer_types:
            print(f"  正在分析癌症类型: {cancer_type}")
            
            # 拟合主效应模型
            main_model = self.fit_main_effect_model(data, cancer_type)
            if main_model:
                all_results[f"main_effect_{cancer_type}"] = self.results[f"main_effect_{cancer_type}"]
            
            # 拟合领域探索模型
            domain_model = self.fit_domain_exploration_model(data, cancer_type)  
            if domain_model:
                all_results[f"domain_exploration_{cancer_type}"] = self.results[f"domain_exploration_{cancer_type}"]
        
        return all_results
    
    def _extract_coefficients(self, model) -> Dict[str, Any]:
        """
        从模型中提取EQI系数和统计信息
        
        参数:
            model: 拟合的模型对象
            
        返回:
            系数字典
        """
        try:
            # 获取模型摘要
            params = model.params
            conf_int = model.conf_int()
            pvalues = model.pvalues
            
            # 提取EQI相关系数
            eqi_results = {}
            
            for param_name in params.index:
                # 寻找EQI相关的参数 - 修改匹配规则
                # 匹配形如 "C(EQI, Treatment(1))[T.2]" 或 "C(EQI_air, Treatment(1))[T.2]" 的参数
                if 'EQI' in param_name and ('[T.' in param_name):
                    
                    eqi_results[param_name] = {
                        'coefficient': params[param_name],
                        'pvalue': pvalues[param_name],
                        'ci_lower': conf_int.iloc[params.index.get_loc(param_name), 0],
                        'ci_upper': conf_int.iloc[params.index.get_loc(param_name), 1],
                        'significant': pvalues[param_name] < 0.05
                    }
            
            # 添加模型拟合信息
            model_info = {
                'n_obs': int(model.nobs),
                'aic': float(model.aic),
                'bic': float(model.bic),
                'log_likelihood': float(model.llf)
            }
            
            return {
                'eqi_coefficients': eqi_results,
                'model_info': model_info
            }
            
        except Exception as e:
            logger.error(f"系数提取失败: {e}")
            return {}