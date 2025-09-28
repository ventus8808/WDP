#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PyMC数据处理工具模块
数据检查，数据加载，准备好所有传入模型的，需要计算的数据
Author: WDP Analysis Team
Date: 2025-09-26
"""

import pandas as pd
import numpy as np
from pathlib import Path
import yaml
import geopandas as gpd
from scipy import sparse
import networkx as nx
from typing import Tuple, Dict, Optional, Union, List
from Utils_Others import get_pymc_config
import warnings
warnings.filterwarnings('ignore')


class WDPDataLoader:
    """WDP项目数据加载和预处理类"""
    
    def __init__(self, config_path: Optional[Union[str, Path]] = None):
        """
        初始化数据加载器
        
        Parameters
        ----------
        config_path : str or Path, optional
            配置文件路径，默认使用项目根目录的config.yaml
        """
        # 统一使用中心配置加载（单一事实来源）
        if config_path is None:
            project_root = Path(__file__).resolve().parents[2]
            config_path = project_root / "config.yaml"
        self.config_path = Path(config_path)
        self.project_root = self.config_path.parent
        # 仅保留pymc_analysis部分作为本模块配置
        self.pymc_config = get_pymc_config(self.config_path)
    
    def get_data_path(self, key: str, **kwargs) -> Path:
        """
        获取数据文件路径
        
        Parameters
        ----------
        key : str
            配置文件中的数据路径键
        **kwargs : dict
            格式化参数
            
        Returns
        -------
        Path
            数据文件绝对路径
        """
        data_files = self.pymc_config.get('data_files', {})

        if key not in data_files:
            raise KeyError(f"数据路径键 '{key}' 未在配置中找到")

        path_template = data_files[key]

        # 格式化路径模板
        if kwargs:
            path_template = path_template.format(**kwargs)

        return self.project_root / path_template
    
    def load_mortality_data(self, disease_code: str) -> pd.DataFrame:
        """
        加载死亡率数据（支持审查数据）
        
        Parameters
        ----------
        disease_code : str
            疾病编码，如'C81-C96'
            
        Returns
        -------
        pd.DataFrame
            死亡率数据框，包含必需的列
        """
        # 获取CDC数据路径
        cdc_file = self.get_data_path('cdc_data_template', disease_code=disease_code)
        
        if not cdc_file.exists():
            raise FileNotFoundError(f"疾病数据文件不存在: {cdc_file}")
        
        print(f"加载疾病数据: {cdc_file}")
        df = pd.read_csv(cdc_file)
        
        # 验证必需的列
        required_cols = ['COUNTY_FIPS', 'Year', 'Deaths_Type', 'Population']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"缺少必需的列: {missing_cols}")
        
        # 处理审查数据
        df = self._process_censored_data(df)
        
        # 数据类型转换（安全处理NaN值）
        df['COUNTY_FIPS'] = pd.to_numeric(df['COUNTY_FIPS'], errors='coerce').astype('int64')
        df['Year'] = pd.to_numeric(df['Year'], errors='coerce').astype('int32')
        df['Population'] = pd.to_numeric(df['Population'], errors='coerce').fillna(0).astype('int32')
        
        print(f"加载了 {len(df)} 条记录，覆盖 {df['COUNTY_FIPS'].nunique()} 个县，{df['Year'].nunique()} 个年份")
        
        return df
    
    def _process_censored_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        处理审查数据，生成观测和审查变量
        
        Parameters
        ----------
        df : pd.DataFrame
            原始数据框
            
        Returns
        -------
        pd.DataFrame
            处理后的数据框
        """
        # 复制数据框避免修改原始数据
        df = df.copy()
        
        # 为观测数据和审查数据创建列
        if 'Deaths_Observed' not in df.columns:
            df['Deaths_Observed'] = np.nan
        if 'Deaths_Censored_Lower' not in df.columns:
            df['Deaths_Censored_Lower'] = np.nan
        if 'Deaths_Censored_Upper' not in df.columns:
            df['Deaths_Censored_Upper'] = np.nan
        
        # 处理观测数据
        observed_mask = df['Deaths_Type'] == 'observed'
        if observed_mask.any():
            # 检查是否有Deaths列，如果没有则使用Deaths_Observed列
            if 'Deaths' in df.columns:
                df.loc[observed_mask, 'Deaths_Observed'] = df.loc[observed_mask, 'Deaths'].fillna(0)
            elif 'Deaths_Observed' in df.columns:
                # 如果Deaths_Observed列已存在且有数据，保持不变
                pass
            else:
                print("⚠️  警告: 既没有Deaths列也没有Deaths_Observed列数据")
        
        # 处理审查数据
        censored_mask = df['Deaths_Type'] == 'censored'
        if censored_mask.any():
            df.loc[censored_mask, 'Deaths_Censored_Lower'] = 1
            df.loc[censored_mask, 'Deaths_Censored_Upper'] = 9
        
        return df
    
    def load_covariate_data(self) -> pd.DataFrame:
        """
        加载协变量数据（PCA主成分）
        
        Returns
        -------
        pd.DataFrame
            协变量数据框，包含PCA分数
        """
        pca_file = self.get_data_path('pca_covariates')
        
        if not pca_file.exists():
            raise FileNotFoundError(f"协变量文件不存在: {pca_file}")
        
        print(f"加载协变量数据: {pca_file}")
        df = pd.read_csv(pca_file)
        
        # 验证必需的列
        required_cols = ['COUNTY_FIPS', 'Year']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"协变量数据缺少必需的列: {missing_cols}")
        
        # 检查PCA分数列
        pca_cols = [col for col in df.columns if 'PC' in col]
        if not pca_cols:
            print("⚠️  警告: 未找到PCA分数列")
        else:
            print(f"找到PCA分数列: {pca_cols}")
        
        return df
    
    def load_spatial_adjacency(self) -> Tuple[sparse.csr_matrix, np.ndarray]:
        """
        加载县级空间邻接矩阵
        
        Returns
        -------
        Tuple[sparse.csr_matrix, np.ndarray]
            邻接矩阵(稀疏格式)和县FIPS代码数组
        """
        adj_file = self.get_data_path('adjacency_data')
        
        if not adj_file.exists():
            raise FileNotFoundError(f"邻接数据文件不存在: {adj_file}")
        
        print(f"加载邻接数据: {adj_file}")
        df = pd.read_csv(adj_file)
        
        # 验证列名
        if 'county_from' in df.columns and 'county_to' in df.columns:
            county_from_col = 'county_from'
            county_to_col = 'county_to'
        elif 'COUNTY_FIPS' in df.columns and 'NEIGHBOR_FIPS' in df.columns:
            county_from_col = 'COUNTY_FIPS'
            county_to_col = 'NEIGHBOR_FIPS'
        else:
            raise ValueError(f"邻接数据文件缺少必需的列结构")
        
        # 获取所有唯一的县FIPS
        all_counties = set(df[county_from_col].unique()) | set(df[county_to_col].unique())
        county_fips = np.array(sorted(all_counties))
        n_counties = len(county_fips)
        
        # 创建县FIPS到索引的映射
        county_to_idx = {fips: i for i, fips in enumerate(county_fips)}
        
        # 构建邻接矩阵
        row_indices = []
        col_indices = []
        
        for _, row in df.iterrows():
            from_idx = county_to_idx[row[county_from_col]]
            to_idx = county_to_idx[row[county_to_col]]
            
            row_indices.extend([from_idx, to_idx])
            col_indices.extend([to_idx, from_idx])
        
        # 创建稀疏邻接矩阵
        data = np.ones(len(row_indices))
        adj_matrix = sparse.csr_matrix((data, (row_indices, col_indices)), 
                                       shape=(n_counties, n_counties))
        
        print(f"构建邻接矩阵: {n_counties}×{n_counties}，{len(df)} 条邻接关系")
        
        return adj_matrix, county_fips
    
    def load_pesticide_data(self, measure_type: str = "Weight") -> pd.DataFrame:
        """
        加载农药使用数据
        
        Parameters
        ----------
        measure_type : str
            测量类型，"Weight"或"Density"
            
        Returns
        -------
        pd.DataFrame
            农药使用数据框
        """
        if measure_type == "Weight":
            pesticide_file = self.get_data_path('pesticide_data')
        elif measure_type == "Density":
            pesticide_file = self.get_data_path('pesticide_density_data')
        else:
            raise ValueError(f"不支持的测量类型: {measure_type}")
        
        if not pesticide_file.exists():
            raise FileNotFoundError(f"农药数据文件不存在: {pesticide_file}")
        
        print(f"加载农药数据: {pesticide_file}")
        df = pd.read_csv(pesticide_file)
        
        return df
    
    def calculate_lagged_exposure(self, pesticide_df: pd.DataFrame, 
                                 compound: str, lag_years: int = 5) -> Tuple[pd.DataFrame, str, str]:
        """
        计算滞后暴露
        
        Parameters
        ----------
        pesticide_df : pd.DataFrame
            农药使用数据
        compound : str
            化合物名称或列名模式
        lag_years : int
            滞后年份
            
        Returns
        -------
        pd.DataFrame
            包含滞后暴露的数据框
        """
        # 查找匹配的化合物列
        if compound in pesticide_df.columns:
            compound_col = compound
        else:
            # 尝试多种匹配模式
            potential_cols = []
            
            # 1. 直接匹配
            potential_cols.extend([col for col in pesticide_df.columns 
                                 if compound.lower() == col.lower()])
            
            # 2. 化学编号匹配 (如 24D -> chem24 或 cat24)
            if compound.upper().endswith('D'):
                chem_num = compound[:-1]
                potential_cols.extend([col for col in pesticide_df.columns 
                                     if f"chem{chem_num}_" in col.lower() or f"cat{chem_num}_" in col.lower()])
            
            # 3. 模糊匹配
            if not potential_cols:
                potential_cols.extend([col for col in pesticide_df.columns 
                                     if compound.lower() in col.lower()])
            
            if not potential_cols:
                raise ValueError(f"未找到匹配的化合物列: {compound}")
            
            # 优先选择avg估计值
            avg_cols = [col for col in potential_cols if 'avg' in col]
            compound_col = avg_cols[0] if avg_cols else potential_cols[0]
            print(f"使用化合物列: {compound_col}")
        # 估计类型推断
        col_lower = compound_col.lower()
        if 'avg' in col_lower:
            estimate_type = 'avg'
        elif 'mean' in col_lower:
            estimate_type = 'mean'
        elif 'median' in col_lower:
            estimate_type = 'median'
        else:
            estimate_type = 'value'
        
        # 标准化关键列名（Year 与 COUNTY_FIPS）
        year_candidates = ['Year', 'YEAR', 'year', 'CalendarYear', 'calendar_year', 'yr']
        county_candidates = ['COUNTY_FIPS', 'county_fips', 'FIPS', 'FIPS_COUNTY', 'CountyFIPS', 'county']

        year_col = next((c for c in year_candidates if c in pesticide_df.columns), None)
        county_col = next((c for c in county_candidates if c in pesticide_df.columns), None)

        if year_col is None:
            raise KeyError(f"未找到年份列。可接受的列名: {year_candidates}；实际列: {list(pesticide_df.columns)}")
        if county_col is None:
            raise KeyError(f"未找到县FIPS列。可接受的列名: {county_candidates}；实际列: {list(pesticide_df.columns)}")

        print(f"识别到年份列: {year_col}，县列: {county_col}")

        # 创建滞后暴露数据并统一列名
        exposure_col_name = f'{compound}_lag{lag_years}'
        lagged_df = pesticide_df[[county_col, year_col, compound_col]].copy()
        lagged_df = lagged_df.rename(columns={county_col: 'COUNTY_FIPS', year_col: 'Year', compound_col: exposure_col_name})
        
        # 年份为数值，排序
        lagged_df['Year'] = pd.to_numeric(lagged_df['Year'], errors='coerce')
        lagged_df = lagged_df.sort_values(by=['COUNTY_FIPS', 'Year'])

        # 将年份整体平移 lag_years，用于与死亡年份对齐
        lagged_df['Year'] = lagged_df['Year'] + lag_years
        
        # 计算按县的滚动平均（仅基于历史数据，不包含当前年）
        rolling_avg = (
            lagged_df
            .groupby('COUNTY_FIPS')[exposure_col_name]
            .rolling(window=lag_years, min_periods=1)
            .mean()
            .shift(1)  # 确保只使用历史数据
            .reset_index(level=0, drop=True)
        )
        lagged_df[exposure_col_name] = rolling_avg
        
        print(f"计算滞后暴露: {compound} -> lag{lag_years}年")
        
        return lagged_df, compound_col, estimate_type
    
    def get_model_covariates(self, model_type: str) -> List[str]:
        """
        获取指定模型的协变量列表
        
        Parameters
        ----------
        model_type : str
            模型类型 (M0, M1, M2, M3)
            
        Returns
        -------
        List[str]
            协变量列名列表
        """
        models = self.pymc_config.get('models', {})
        
        if model_type not in models:
            raise ValueError(f"未知的模型类型: {model_type}")
        
        model_config = models[model_type]
        covariates = model_config.get('covariates', [])
        
        # 根据PCA诊断结果映射协变量名称
        covariate_mapping = {
            'SVI_std': 'SVI_PC1',
            'Climate1_std': 'ENV_PC1', 
            'Climate2_std': 'ENV_PC2'
        }
        
        actual_covariates = []
        for covar in covariates:
            if covar in covariate_mapping:
                actual_covariates.append(covariate_mapping[covar])
            else:
                actual_covariates.append(covar)
        
        return actual_covariates
    
    def prepare_model_data(self, disease_code: str, compound: str, 
                          model_type: str = 'M0', lag_years: int = 5,
                          measure_type: str = 'Weight') -> Dict:
        """
        准备完整的模型数据
        
        Parameters
        ----------
        disease_code : str
            疾病编码
        compound : str
            化合物名称
        model_type : str
            模型类型
        lag_years : int
            滞后年份
        measure_type : str
            农药测量类型
            
        Returns
        -------
        Dict
            包含所有模型数据的字典
        """
        print(f"\n=== 准备模型数据 ===")
        print(f"疾病: {disease_code}")
        print(f"化合物: {compound}")
        print(f"模型: {model_type}")
        print(f"滞后: {lag_years}年")
        
        # 1. 加载死亡率数据
        mortality_df = self.load_mortality_data(disease_code)
        
        # 2. 加载协变量数据
        covariate_df = self.load_covariate_data()
        
        # 3. 加载空间邻接矩阵
        adj_matrix, county_fips = self.load_spatial_adjacency()
        
        # 4. 加载农药数据并计算滞后暴露
        pesticide_df = self.load_pesticide_data(measure_type)
        lagged_exposure_df, selected_exposure_column, estimate_type = self.calculate_lagged_exposure(
            pesticide_df, compound, lag_years
        )
        
        # 5. 合并数据
        print("合并数据...")
        
        # 合并死亡率和协变量数据（使用inner以避免隐式缺失带来的偏倚）
        initial_rows = len(mortality_df)
        merged_df = pd.merge(mortality_df, covariate_df,
                           on=['COUNTY_FIPS', 'Year'], how='inner')
        print(f"合并协变量后: {initial_rows} -> {len(merged_df)} 行 (排除 {initial_rows - len(merged_df)} 行)")

        # 合并滞后暴露数据（inner）
        prev_rows = len(merged_df)
        merged_df = pd.merge(merged_df, lagged_exposure_df,
                           on=['COUNTY_FIPS', 'Year'], how='inner')
        print(f"合并暴露后: {prev_rows} -> {len(merged_df)} 行 (排除 {prev_rows - len(merged_df)} 行)")
        
        # 6. 过滤有效数据
        # 确保在空间邻接矩阵中的县
        valid_counties = set(county_fips)
        merged_df = merged_df[merged_df['COUNTY_FIPS'].isin(valid_counties)]
        
        # 移除缺失暴露数据的记录
        exposure_col = f'{compound}_lag{lag_years}'
        merged_df = merged_df.dropna(subset=[exposure_col])
        
        if len(merged_df) == 0:
            raise ValueError(f"合并后无有效数据记录")
        
        # 7. 创建索引映射
        county_to_idx = {fips: i for i, fips in enumerate(county_fips)}
        merged_df['county_idx'] = merged_df['COUNTY_FIPS'].map(county_to_idx)
        
        years = sorted(merged_df['Year'].unique())
        year_to_idx = {year: i for i, year in enumerate(years)}
        merged_df['time_idx'] = merged_df['Year'].map(year_to_idx)
        
        # 8. 分离观测数据和审查数据
        merged_df['is_censored'] = merged_df['Deaths_Type'] == 'censored'
        observed_mask = ~merged_df['is_censored']
        censored_mask = merged_df['is_censored']
        
        obs_data = merged_df[observed_mask].copy()
        cens_data = merged_df[censored_mask].copy()
        
        # 9. 准备协变量矩阵
        model_covariates = self.get_model_covariates(model_type)
        
        if model_covariates:
            missing_covars = [col for col in model_covariates if col not in merged_df.columns]
            if missing_covars:
                print(f"⚠️  警告: 缺少协变量列 {missing_covars}，将被忽略")
                model_covariates = [col for col in model_covariates if col in merged_df.columns]
            
            if model_covariates:
                X = merged_df[model_covariates].values
            else:
                X = np.empty((len(merged_df), 0))
        else:
            X = np.empty((len(merged_df), 0))
        
        # 10. 对数变换和标准化暴露
        exposure_values = merged_df[exposure_col].values
        
        # 更合理的零值处理
        non_zero_mask = exposure_values > 0
        if np.any(non_zero_mask):
            min_non_zero = np.min(exposure_values[non_zero_mask])
            offset = min_non_zero / 10  # 使用最小非零值的1/10作为偏移
        else:
            offset = 1e-6  # 安全回退
        log_exposure = np.log(exposure_values + offset)
        
        log_exposure_mean = float(np.mean(log_exposure))
        log_exposure_std = float(np.std(log_exposure)) if np.std(log_exposure) > 0 else 1.0
        log_exposure_stdized = (log_exposure - log_exposure_mean) / log_exposure_std
        
        # 读取农药类别映射（可选）
        category = 'Unknown'
        try:
            mapping_path = self.get_data_path('pesticide_mapping')
            if mapping_path.exists():
                mapping_df = pd.read_csv(mapping_path)

                # 支持两种列名风格：旧(Compound/Category) 与 新(compound_name/category1_name)
                col_compound = None
                col_category = None
                lower_cols = {c.lower(): c for c in mapping_df.columns}

                if 'compound' in lower_cols and 'category' in lower_cols:
                    col_compound = lower_cols['compound']
                    col_category = lower_cols['category']
                elif 'compound_name' in lower_cols and 'category1_name' in lower_cols:
                    col_compound = lower_cols['compound_name']
                    col_category = lower_cols['category1_name']

                if col_compound and col_category:
                    def normalize_name(s: str) -> str:
                        s = str(s).strip().lower()
                        # 常见别名归一化
                        aliases = {
                            '24d': '2,4-d',
                            '2,4d': '2,4-d',
                            '2,4-d': '2,4-d',
                            'glyphosate': 'glyphosate',
                            'atrazine': 'atrazine',
                        }
                        return aliases.get(s, s)

                    target = normalize_name(compound)
                    mapping_df['_norm'] = mapping_df[col_compound].astype(str).map(normalize_name)
                    matched = mapping_df[mapping_df['_norm'] == target]
                    if len(matched) > 0:
                        category = str(matched.iloc[0][col_category])
        except Exception as e:
            print(f"⚠️  农药类别映射读取失败: {e}")

        # 构建最终数据字典
        model_data = {
            # 基本信息
            'disease_code': disease_code,
            'compound': compound,
            'model_type': model_type,
            'lag_years': lag_years,
            'measure_type': measure_type,
            'estimate_type': estimate_type,
            'category': category,
            
            # 原始数据
            'full_data': merged_df,
            'observed_data': obs_data,
            'censored_data': cens_data,
            
            # 空间结构
            'adj_matrix': adj_matrix,
            'county_fips': county_fips,
            'n_counties': len(county_fips),
            
            # 时间结构
            'years': years,
            'n_years': len(years),
            
            # 暴露数据
            'exposure_raw': exposure_values,
            'exposure_log': log_exposure,
            'exposure_log_mean': log_exposure_mean,
            'exposure_log_std': log_exposure_std,
            'exposure_log_stdized': log_exposure_stdized,
            'exposure_col': exposure_col,
            'selected_exposure_column': selected_exposure_column,
            
            # 观测变量
            'y_obs': obs_data['Deaths_Observed'].values if len(obs_data) > 0 else np.array([]),
            'n_obs': obs_data['Population'].values if len(obs_data) > 0 else np.array([]),
            'exposure_obs': log_exposure[observed_mask] if observed_mask.any() else np.array([]),
            'county_obs': obs_data['county_idx'].values if len(obs_data) > 0 else np.array([]),
            'time_obs': obs_data['time_idx'].values if len(obs_data) > 0 else np.array([]),
            
            # 审查变量
            'y_cens_lower': cens_data['Deaths_Censored_Lower'].values if len(cens_data) > 0 else np.array([]),
            'y_cens_upper': cens_data['Deaths_Censored_Upper'].values if len(cens_data) > 0 else np.array([]),
            'n_cens': cens_data['Population'].values if len(cens_data) > 0 else np.array([]),
            'exposure_cens': log_exposure[censored_mask] if censored_mask.any() else np.array([]),
            'county_cens': cens_data['county_idx'].values if len(cens_data) > 0 else np.array([]),
            'time_cens': cens_data['time_idx'].values if len(cens_data) > 0 else np.array([]),
            
            # 协变量
            'X': X,
            'covariate_names': model_covariates,
            'n_covariates': X.shape[1],
            
            # 数据计数
            'n_obs_points': len(obs_data),
            'n_cens_points': len(cens_data),
            'n_total_points': len(merged_df)
        }
        
        # 打印数据摘要
        print(f"\n=== 数据摘要 ===")
        print(f"县数: {model_data['n_counties']}")
        print(f"年数: {model_data['n_years']}")
        print(f"协变量数: {model_data['n_covariates']}")
        print(f"观测点数: {model_data['n_obs_points']}")
        print(f"审查点数: {model_data['n_cens_points']}")
        print(f"协变量: {model_data['covariate_names']}")
        
        return model_data


if __name__ == "__main__":
    # 测试数据加载
    loader = WDPDataLoader()
    
    try:
        # 测试准备模型数据
        model_data = loader.prepare_model_data(
            disease_code="C81-C96",
            compound="24D",  # 假设的化合物名称
            model_type="M1",
            lag_years=5,
            measure_type="Weight"
        )
        
        print("\n✅ 数据加载测试成功！")
        
    except Exception as e:
        print(f"❌ 数据加载测试失败: {e}")
        import traceback
        traceback.print_exc()