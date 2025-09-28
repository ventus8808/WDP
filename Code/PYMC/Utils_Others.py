#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PyMC辅助函数模块
其他可能用到的实用函数
Author: WDP Analysis Team
Date: 2025-09-26
"""

import numpy as np
import pandas as pd
from pathlib import Path
import yaml
from typing import Dict, List, Optional, Union, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')


def load_config(config_path: Optional[Union[str, Path]] = None) -> Dict:
    """
    加载项目配置文件
    
    Parameters
    ----------
    config_path : str or Path, optional
        配置文件路径，默认使用项目根目录的config.yaml
        
    Returns
    -------
    Dict
        配置字典
    """
    if config_path is None:
        project_root = Path(__file__).resolve().parents[2]
        config_path = project_root / "config.yaml"
    
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"配置文件不存在: {config_path}")
    
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def get_pymc_config(config_path: Optional[Union[str, Path]] = None) -> Dict:
    """
    获取PyMC分析配置
    
    Parameters
    ----------
    config_path : str or Path, optional
        配置文件路径
        
    Returns
    -------
    Dict
        PyMC配置字典
    """
    config = load_config(config_path)
    return config.get('pymc_analysis', {})


def get_sampling_config(mode: str = 'test', 
                       config_path: Optional[Union[str, Path]] = None) -> Dict:
    """
    获取MCMC采样配置
    
    Parameters
    ----------
    mode : str
        采样模式 ('test' 或 'production')
    config_path : str or Path, optional
        配置文件路径
        
    Returns
    -------
    Dict
        采样配置字典
    """
    pymc_config = get_pymc_config(config_path)
    sampling_config = pymc_config.get('sampling', {})
    
    if mode not in sampling_config:
        print(f"⚠️  采样模式 '{mode}' 不存在，使用默认配置")
        return {
            'draws': 1000,
            'tune': 500,
            'chains': 2,
            'cores': 2,
            'target_accept': 0.8
        }
    
    return sampling_config[mode]


def validate_disease_code(disease_code: str, 
                         config_path: Optional[Union[str, Path]] = None) -> bool:
    """
    验证疾病编码是否有效
    
    Parameters
    ----------
    disease_code : str
        疾病编码
    config_path : str or Path, optional
        配置文件路径
        
    Returns
    -------
    bool
        是否有效
    """
    pymc_config = get_pymc_config(config_path)
    analysis_config = pymc_config.get('analysis', {})
    valid_diseases = analysis_config.get('disease_codes', ['C81-C96', 'C50', 'C34'])
    
    return disease_code in valid_diseases


def validate_model_type(model_type: str, 
                       config_path: Optional[Union[str, Path]] = None) -> bool:
    """
    验证模型类型是否有效
    
    Parameters
    ----------
    model_type : str
        模型类型
    config_path : str or Path, optional
        配置文件路径
        
    Returns
    -------
    bool
        是否有效
    """
    pymc_config = get_pymc_config(config_path)
    models = pymc_config.get('models', {})
    
    return model_type in models


def parse_compound_list(compound_input: str) -> List[str]:
    """
    解析化合物输入参数
    
    Parameters
    ----------
    compound_input : str
        化合物输入字符串，支持多种格式
        
    Returns
    -------
    List[str]
        化合物列表
    """
    # 移除空格并分割
    compounds = [c.strip() for c in compound_input.split(',')]
    
    # 处理特殊关键词
    processed_compounds = []
    for compound in compounds:
        if compound.upper() == 'ALL':
            # 返回所有可用化合物（需要从数据中获取）
            processed_compounds.append('ALL')
        elif compound.upper() == 'TEST':
            # 测试化合物
            processed_compounds.extend(['24D', 'Atrazine', 'Glyphosate'])
        else:
            processed_compounds.append(compound)
    
    return list(set(processed_compounds))  # 去重


def parse_model_list(model_input: str) -> List[str]:
    """
    解析模型类型输入参数
    
    Parameters
    ----------
    model_input : str
        模型输入字符串
        
    Returns
    -------
    List[str]
        模型类型列表
    """
    models = [m.strip().upper() for m in model_input.split(',')]
    
    # 验证模型类型
    valid_models = ['M0', 'M1', 'M2', 'M3']
    validated_models = []
    
    for model in models:
        if model in valid_models:
            validated_models.append(model)
        else:
            print(f"⚠️  忽略无效的模型类型: {model}")
    
    return validated_models if validated_models else ['M0']


def parse_lag_years(lag_input: Union[str, int, List]) -> List[int]:
    """
    解析滞后年份输入参数
    
    Parameters
    ----------
    lag_input : str, int, or list
        滞后年份输入
        
    Returns
    -------
    List[int]
        滞后年份列表
    """
    if isinstance(lag_input, int):
        return [lag_input]
    elif isinstance(lag_input, list):
        return [int(x) for x in lag_input]
    elif isinstance(lag_input, str):
        return [int(x.strip()) for x in lag_input.split(',')]
    else:
        return [5]  # 默认值


def create_output_filename(disease: str, compound: str, 
                         timestamp: Optional[str] = None) -> str:
    """
    创建输出文件名
    
    Parameters
    ----------
    disease : str
        疾病编码
    compound : str
        化合物名称
    timestamp : str, optional
        时间戳
        
    Returns
    -------
    str
        文件名
    """
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    return f"{disease}_{compound}_Results.csv"


def setup_logging(level: str = 'INFO') -> None:
    """
    设置日志系统
    
    Parameters
    ----------
    level : str
        日志级别
    """
    import logging
    
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def plot_trace(trace: az.InferenceData, var_names: Optional[List[str]] = None,
               figsize: Tuple[int, int] = (12, 8)) -> None:
    """
    绘制MCMC轨迹图
    
    Parameters
    ----------
    trace : az.InferenceData
        MCMC采样结果
    var_names : List[str], optional
        要绘制的变量名
    figsize : Tuple[int, int]
        图像大小
    """
    if var_names is None:
        var_names = ['beta_exposure', 'intercept']
    
    az.plot_trace(trace, var_names=var_names, figsize=figsize)
    plt.tight_layout()
    plt.show()


def plot_posterior(trace: az.InferenceData, var_names: Optional[List[str]] = None,
                  figsize: Tuple[int, int] = (12, 6)) -> None:
    """
    绘制后验分布图
    
    Parameters
    ----------
    trace : az.InferenceData
        MCMC采样结果
    var_names : List[str], optional
        要绘制的变量名
    figsize : Tuple[int, int]
        图像大小
    """
    if var_names is None:
        var_names = ['beta_exposure', 'intercept']
    
    az.plot_posterior(trace, var_names=var_names, figsize=figsize)
    plt.tight_layout()
    plt.show()


def calculate_relative_risk(beta_samples: np.ndarray, 
                          exposure_change: float) -> Dict:
    """
    计算相对风险及其置信区间
    
    Parameters
    ----------
    beta_samples : np.ndarray
        回归系数的后验样本
    exposure_change : float
        暴露变化量
        
    Returns
    -------
    Dict
        相对风险统计量
    """
    rr_samples = np.exp(beta_samples * exposure_change)
    
    return {
        'rr_mean': np.mean(rr_samples),
        'rr_median': np.median(rr_samples),
        'rr_ci': np.percentile(rr_samples, [2.5, 97.5]),
        'rr_prob_gt_1': np.mean(rr_samples > 1.0)
    }


def compare_models(traces: Dict[str, az.InferenceData]) -> pd.DataFrame:
    """
    比较多个模型的拟合度
    
    Parameters
    ----------
    traces : Dict[str, az.InferenceData]
        模型名称到采样结果的映射
        
    Returns
    -------
    pd.DataFrame
        模型比较表
    """
    comparison_data = []
    
    for model_name, trace in traces.items():
        try:
            waic = az.waic(trace)
            loo = az.loo(trace)
            
            comparison_data.append({
                'Model': model_name,
                'WAIC': float(waic.waic),
                'WAIC_SE': float(waic.se) if hasattr(waic, 'se') else np.nan,
                'LOO': float(loo.loo),
                'LOO_SE': float(loo.se) if hasattr(loo, 'se') else np.nan
            })
        except Exception as e:
            print(f"⚠️  模型 {model_name} 比较失败: {e}")
            comparison_data.append({
                'Model': model_name,
                'WAIC': np.nan,
                'WAIC_SE': np.nan,
                'LOO': np.nan,
                'LOO_SE': np.nan
            })
    
    df = pd.DataFrame(comparison_data)
    
    # 计算Delta WAIC和Delta LOO
    if len(df) > 0 and not df['WAIC'].isna().all():
        best_waic = df['WAIC'].min()
        df['Delta_WAIC'] = df['WAIC'] - best_waic
        
        best_loo = df['LOO'].min()
        df['Delta_LOO'] = df['LOO'] - best_loo
    
    return df.sort_values('WAIC', na_last=True)


def check_data_availability(disease_code: str, compound: str, 
                           config_path: Optional[Union[str, Path]] = None) -> Dict:
    """
    检查指定分析所需数据的可用性
    
    Parameters
    ----------
    disease_code : str
        疾病编码
    compound : str
        化合物名称
    config_path : str or Path, optional
        配置文件路径
        
    Returns
    -------
    Dict
        数据可用性报告
    """
    from Utils_Data import WDPDataLoader
    
    try:
        loader = WDPDataLoader(config_path)
        
        report = {
            'disease_code': disease_code,
            'compound': compound,
            'data_available': True,
            'issues': []
        }
        
        # 检查疾病数据
        try:
            mortality_df = loader.load_mortality_data(disease_code)
            report['mortality_records'] = len(mortality_df)
            report['mortality_counties'] = mortality_df['COUNTY_FIPS'].nunique()
        except Exception as e:
            report['data_available'] = False
            report['issues'].append(f"疾病数据加载失败: {e}")
        
        # 检查协变量数据
        try:
            covariate_df = loader.load_covariate_data()
            report['covariate_records'] = len(covariate_df)
        except Exception as e:
            report['data_available'] = False
            report['issues'].append(f"协变量数据加载失败: {e}")
        
        # 检查农药数据
        try:
            pesticide_df = loader.load_pesticide_data("Weight")
            
            # 使用与Utils_Data相同的匹配逻辑
            matching_cols = []
            
            # 1. 直接匹配
            if compound in pesticide_df.columns:
                matching_cols.append(compound)
            
            # 2. 化学编号匹配 (如 24D -> chem24 或 cat24)
            if compound.upper().endswith('D'):
                chem_num = compound[:-1]
                matching_cols.extend([col for col in pesticide_df.columns 
                                    if f"chem{chem_num}_" in col.lower() or f"cat{chem_num}_" in col.lower()])
            
            # 3. 模糊匹配
            if not matching_cols:
                matching_cols.extend([col for col in pesticide_df.columns 
                                    if compound.lower() in col.lower()])
            
            if matching_cols:
                report['pesticide_available'] = True
                # 优先选择avg估计值
                avg_cols = [col for col in matching_cols if 'avg' in col]
                report['pesticide_column'] = avg_cols[0] if avg_cols else matching_cols[0]
            else:
                report['pesticide_available'] = False
                report['issues'].append(f"未找到化合物 '{compound}' 的数据")
        except Exception as e:
            report['data_available'] = False
            report['issues'].append(f"农药数据加载失败: {e}")
        
        # 检查空间数据
        try:
            adj_matrix, county_fips = loader.load_spatial_adjacency()
            report['spatial_counties'] = len(county_fips)
        except Exception as e:
            report['data_available'] = False
            report['issues'].append(f"空间邻接数据加载失败: {e}")
        
        return report
        
    except Exception as e:
        return {
            'disease_code': disease_code,
            'compound': compound,
            'data_available': False,
            'issues': [f"数据检查失败: {e}"]
        }


def print_analysis_summary(model_data: Dict, trace: az.InferenceData) -> None:
    """
    打印分析摘要信息
    
    Parameters
    ----------
    model_data : Dict
        模型数据字典
    trace : az.InferenceData
        MCMC采样结果
    """
    print(f"\n" + "="*60)
    print(f"WDP PyMC分析摘要")
    print(f"="*60)
    print(f"疾病编码: {model_data['disease_code']}")
    print(f"化合物: {model_data['compound']}")
    print(f"模型类型: {model_data['model_type']}")
    print(f"滞后年份: {model_data['lag_years']}")
    print(f"测量类型: {model_data['measure_type']}")
    print(f"-"*60)
    print(f"数据统计:")
    print(f"  县数: {model_data['n_counties']}")
    print(f"  年数: {model_data['n_years']}")
    print(f"  观测点: {model_data['n_obs_points']}")
    print(f"  审查点: {model_data['n_cens_points']}")
    print(f"  协变量数: {model_data['n_covariates']}")
    
    if model_data['n_covariates'] > 0:
        print(f"  协变量: {', '.join(model_data['covariate_names'])}")
    
    print(f"-"*60)
    print(f"采样统计:")
    
    # 计算采样统计
    posterior = trace.posterior
    n_chains = posterior.dims['chain']
    n_draws = posterior.dims['draw']
    
    print(f"  链数: {n_chains}")
    print(f"  每链样本数: {n_draws}")
    print(f"  总样本数: {n_chains * n_draws}")
    
    # 收敛统计
    try:
        rhat = az.rhat(trace)
        ess = az.ess(trace)
        
        max_rhat = float(rhat.max()) if hasattr(rhat, 'max') else np.nan
        min_ess = float(ess.min()) if hasattr(ess, 'min') else np.nan
        
        print(f"  最大R̂: {max_rhat:.4f}")
        print(f"  最小ESS: {min_ess:.0f}")
        
        if max_rhat > 1.01:
            print(f"  ⚠️  警告: 部分参数可能未收敛")
        else:
            print(f"  ✅ 收敛良好")
            
    except Exception as e:
        print(f"  ⚠️  收敛统计计算失败: {e}")
    
    print(f"="*60)


if __name__ == "__main__":
    # 测试辅助函数
    print("测试WDP PyMC辅助函数...")
    
    # 测试配置加载
    try:
        config = get_pymc_config()
        print(f"✅ 配置加载成功")
        
        # 测试采样配置
        test_config = get_sampling_config('test')
        print(f"✅ 测试采样配置: {test_config}")
        
        # 测试验证函数
        print(f"疾病编码C81-C96有效: {validate_disease_code('C81-C96')}")
        print(f"模型类型M0有效: {validate_model_type('M0')}")
        
        # 测试解析函数
        compounds = parse_compound_list("24D, Atrazine, TEST")
        print(f"解析化合物: {compounds}")
        
        models = parse_model_list("M0,M1,M2,M3")
        print(f"解析模型: {models}")
        
        print("✅ 辅助函数测试成功！")
        
    except Exception as e:
        print(f"❌ 辅助函数测试失败: {e}")
        import traceback
        traceback.print_exc()