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
        
        # --- START: 修改建议 ---
        # 为关键列指定更节省内存的类型
        dtype_spec = {
            'COUNTY_FIPS': str, # 始终以字符串形式读取FIPS以保留前导零
            'Year': 'int16',
            'Population': 'float64', # 先读为浮点数以处理可能的空值
            'Deaths_Type': 'category', # 如果类型不多，category很高效
        }
        df = pd.read_csv(cdc_file, dtype=dtype_spec)
        # --- END: 修改建议 ---
        
        # 验证必需的列
        required_cols = ['COUNTY_FIPS', 'Year', 'Deaths_Type', 'Population']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            raise ValueError(f"缺少必需的列: {missing_cols}")
        
        # 处理审查数据
        df = self._process_censored_data(df)
        
        # 将FIPS转换为可空整数
        df['COUNTY_FIPS'] = pd.to_numeric(df['COUNTY_FIPS'], errors='coerce').astype('Int64')
        df['Year'] = pd.to_numeric(df['Year'], errors='coerce').astype('int16')
        # 将Population转换为可空整数
        df['Population'] = df['Population'].astype('Int32')
        
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

    # --- Mapping helpers ---
    def _load_mapping_df(self) -> Optional[pd.DataFrame]:
        try:
            mapping_path = self.get_data_path('pesticide_mapping')
            if mapping_path.exists():
                return pd.read_csv(mapping_path)
        except Exception as e:
            print(f"⚠️  农药mapping读取失败: {e}")
        return None

    @staticmethod
    def _norm_name(s: str) -> str:
        # 小写 + 去除空格、逗号、连字符和下划线
        import re
        return re.sub(r"[\s,_\-]", "", str(s).strip().lower())

    def _resolve_from_mapping(self, compound_input: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """根据mapping将输入名字解析为 (compound_id, compound_display, category1_name)"""
        mapping_df = self._load_mapping_df()
        if mapping_df is None:
            return None, None, None

        lower_cols = {c.lower(): c for c in mapping_df.columns}
        col_cid = lower_cols.get('compound_id')
        col_cname = lower_cols.get('compound_name')
        col_cat = lower_cols.get('category1_name') or lower_cols.get('category')
        if not (col_cid and col_cname):
            return None, None, None

        # 规范化匹配
        norm_target = self._norm_name(compound_input)
        mapping_df['_norm'] = mapping_df[col_cname].astype(str).map(self._norm_name)
        matched = mapping_df[mapping_df['_norm'] == norm_target]
        if len(matched) == 0:
            # 回退：尝试直接按ID匹配
            if compound_input.isdigit():
                matched = mapping_df[mapping_df[col_cid].astype(str) == compound_input]
        if len(matched) == 0:
            return None, None, None

        cid = str(matched.iloc[0][col_cid])
        cname = str(matched.iloc[0][col_cname])
        cat = str(matched.iloc[0][col_cat]) if col_cat in mapping_df.columns else 'Unknown'
        return cid, cname, cat
    
    def calculate_lagged_exposure(self, pesticide_df: pd.DataFrame,
                                 compound: str, lag_years: int,
                                 estimate_type: str) -> Tuple[pd.DataFrame, str, str]:
        """
        计算滞后暴露 (v2.2，支持 estimate_type 选择 min/avg/max)。
        
        接受的 `compound` 格式: "2", "cat21", "Atrazine"
        接受的 `estimate_type` 格式: "min", "avg", "max"
        """
        compound_input = str(compound).strip().lower()
        exposure_type = ''
        exposure_id = ''

        # 1. 解析输入，确定是 chem 还是 cat
        if compound_input.startswith('cat'):
            exposure_type = 'cat'
            exposure_id = compound_input[3:]
            if not exposure_id.isdigit():
                raise ValueError(f"类别格式错误: '{compound}'. 应为 'cat' + 数字 (例如 'cat21').")
        elif compound_input.isdigit():
            exposure_type = 'chem'
            exposure_id = compound_input
        else:
            cid, _, _ = self._resolve_from_mapping(compound_input)
            if cid is None:
                raise ValueError(f"无法从 mapping.csv 中解析化合物名称: '{compound}'")
            exposure_type = 'chem'
            exposure_id = cid
        
        # 2. 根据 estimate_type 构造后缀候选项
        selected_estimate = estimate_type.lower()
        if selected_estimate in ['avg', 'mean', 'median']:
            suffixes = ['_avg', '_mean', '_median']
            actual_estimate_type = 'avg'
        elif selected_estimate == 'min':
            suffixes = ['_min']
            actual_estimate_type = 'min'
        elif selected_estimate == 'max':
            suffixes = ['_max']
            actual_estimate_type = 'max'
        else:
            raise ValueError(f"不支持的 estimate_type: '{estimate_type}'. 请使用 'min', 'avg', 或 'max'.")

        candidates = [f"{exposure_type}{exposure_id}{suffix}" for suffix in suffixes]
        
        compound_col = None
        for c in candidates:
            if c in pesticide_df.columns:
                compound_col = c
                break
        
        if compound_col is None:
            raise ValueError(f"在 PNSP.csv 中未找到与 '{compound}' (estimate='{estimate_type}') 匹配的数据列。尝试查找: {candidates}")

        print(f"输入 '{compound}' (estimate='{estimate_type}') 被解析为列: '{compound_col}'")

        # 4. 标准化关键列名
        year_candidates = ['Year', 'YEAR', 'year', 'CalendarYear', 'calendar_year', 'yr']
        county_candidates = ['COUNTY_FIPS', 'county_fips', 'FIPS', 'FIPS_COUNTY', 'CountyFIPS', 'county']

        year_col = next((c for c in year_candidates if c in pesticide_df.columns), None)
        county_col = next((c for c in county_candidates if c in pesticide_df.columns), None)

        if year_col is None:
            raise KeyError(f"未找到年份列。可接受的列名: {year_candidates}；实际列: {list(pesticide_df.columns)}")
        if county_col is None:
            raise KeyError(f"未找到县FIPS列。可接受的列名: {county_candidates}；实际列: {list(pesticide_df.columns)}")

        # --- 修正后的滞后计算逻辑 ---
        print(f"识别到年份列: {year_col}，县列: {county_col}")
        print(f"正在计算 {lag_years} 年滞后暴露，请稍候...")

        # 确保数据类型正确并按县、年排序
        df = pesticide_df[[county_col, year_col, compound_col]].copy()
        df[year_col] = pd.to_numeric(df[year_col], errors='coerce')
        df = df.sort_values(by=[county_col, year_col])

        # 创建一个包含所有县和完整年份范围的 MultiIndex
        all_counties = df[county_col].unique()
        year_range = range(df[year_col].min(), df[year_col].max() + 1 + lag_years)
        multi_index = pd.MultiIndex.from_product([all_counties, year_range], names=[county_col, year_col])

        # 扩展 DataFrame 以填补缺失的年份，用0填充暴露值
        df_full = df.set_index([county_col, year_col]).reindex(multi_index).fillna(0).reset_index()

        # 使用 rolling 计算过去 lag_years 年的平均值
        # closed='left' 确保窗口是 [t-lag_years, t-1]
        df_full['rolling_avg'] = df_full.groupby(county_col)[compound_col].transform(
            lambda x: x.rolling(window=lag_years, closed='left').mean()
        )

        # 生成最终的滞后暴露列名
        exposure_col_name = f'{compound}_{actual_estimate_type}_lag{lag_years}'
        df_full = df_full.rename(columns={'rolling_avg': exposure_col_name})

        # 滞后暴露的年份应该与死亡年份对齐，所以将年份加 lag_years 是不正确的。
        # 正确的做法是直接在原始年份上计算历史暴露。
        # 我们需要的是 `mortality_year` 对应的 `exposure_year` 的历史平均。
        # 因此，我们直接返回带有 `exposure_col_name` 的 df_full，在合并时，
        # `mortality_df` 中的 `Year` 会自动匹配 `df_full` 中计算好的 `Year` 的滞后值。

        # 统一列名
        lagged_df = df_full[[county_col, year_col, exposure_col_name]].copy()
        lagged_df = lagged_df.rename(columns={county_col: 'COUNTY_FIPS', year_col: 'Year'})

        # 移除全为NA的行（这些是窗口期不足导致的）
        lagged_df = lagged_df.dropna(subset=[exposure_col_name])
        
        print(f"计算滞后暴露: {compound} -> 列 {compound_col} -> lag{lag_years}年")
        
        return lagged_df, compound_col, actual_estimate_type
    
    def get_model_covariates(self, model_type: str) -> List[str]:
        """
        获取指定模型的协变量"显示名称"列表（保持与config一致，包含交互表达式）。
        """
        models = self.pymc_config.get('models', {})
        if model_type not in models:
            raise ValueError(f"未知的模型类型: {model_type}")
        model_config = models[model_type]
        covariates = model_config.get('covariates', [])
        return list(covariates)
    
    def prepare_model_data(self, disease_code: str, compound: str, 
                          model_type: str = 'M0', lag_years: int = 5,
                          measure_type: str = 'Weight',
                          estimate_type: str = 'avg') -> Dict: # <--- 修改
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
        lagged_exposure_df, selected_exposure_column, estimate_type_found = self.calculate_lagged_exposure(
            pesticide_df, compound, lag_years, estimate_type # <--- 修改
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
        # 注意：calculate_lagged_exposure 生成的列名包含估计类型 (min/avg/max)
        exposure_col = f"{compound}_{estimate_type_found}_lag{lag_years}"
        merged_df = merged_df.dropna(subset=[exposure_col])

        # 移除人口为缺失或非正数的记录（偏移项需要正人口）
        merged_df = merged_df.copy()
        merged_df['Population'] = pd.to_numeric(merged_df['Population'], errors='coerce')
        before_pop = len(merged_df)
        merged_df = merged_df[merged_df['Population'] > 0]
        after_pop = len(merged_df)
        if after_pop < before_pop:
            print(f"过滤人口<=0或缺失: {before_pop} -> {after_pop} 行 (排除 {before_pop - after_pop} 行)")
        
        # 根据模型需要的协变量，移除NaN行
        model_covariates_display_early = self.get_model_covariates(model_type)
        covariate_mapping_early = {
            'SVI_std': 'SVI_PC1',
            'Climate1_std': 'ENV_PC1',
            'Climate2_std': 'ENV_PC2',
            'Climate3_std': 'ENV_PC3'
        }
        required_cov_cols = []
        for covar_disp in model_covariates_display_early:
            name = str(covar_disp).strip()
            if ' * ' in name and name.endswith('exposure'):
                base_name = name.split('*')[0].strip()
                actual_col = covariate_mapping_early.get(base_name, base_name)
                if actual_col in merged_df.columns:
                    required_cov_cols.append(actual_col)
            else:
                actual_col = covariate_mapping_early.get(name, name)
                if actual_col in merged_df.columns:
                    required_cov_cols.append(actual_col)
        if required_cov_cols:
            before_cov = len(merged_df)
            merged_df = merged_df.dropna(subset=required_cov_cols)
            after_cov = len(merged_df)
            if after_cov < before_cov:
                print(f"过滤协变量缺失: {before_cov} -> {after_cov} 行 (排除 {before_cov - after_cov} 行)")

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

        # 确保计数为整数类型（Poisson支持集为非负整数）
        if 'Deaths_Observed' in obs_data.columns:
            obs_data['Deaths_Observed'] = pd.to_numeric(obs_data['Deaths_Observed'], errors='coerce').fillna(0).round().astype('int64')
        if 'Deaths_Censored_Lower' in cens_data.columns:
            cens_data['Deaths_Censored_Lower'] = pd.to_numeric(cens_data['Deaths_Censored_Lower'], errors='coerce').fillna(1).round().astype('int64')
        if 'Deaths_Censored_Upper' in cens_data.columns:
            cens_data['Deaths_Censored_Upper'] = pd.to_numeric(cens_data['Deaths_Censored_Upper'], errors='coerce').fillna(9).round().astype('int64')

        # 9. 对数变换和标准化暴露（先于构建设计矩阵，以便交互项使用）
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

        # 10. 准备协变量矩阵（支持交互项）
        model_covariates_display = self.get_model_covariates(model_type)
        covariate_mapping = {
            'SVI_std': 'SVI_PC1',
            'Climate1_std': 'ENV_PC1', 
            'Climate2_std': 'ENV_PC2',
            'Climate3_std': 'ENV_PC3'
        }
        X_cols = []
        kept_display_names = []

        for covar_disp in model_covariates_display:
            name = str(covar_disp).strip()
            if ' * ' in name and name.endswith('exposure'):
                # 交互项: <base> * exposure
                base_name = name.split('*')[0].strip()
                actual_col = covariate_mapping.get(base_name, base_name)
                if actual_col not in merged_df.columns:
                    print(f"⚠️  交互项基础列缺失，跳过: {base_name} -> {actual_col}")
                    continue
                base_vec = merged_df[actual_col].values.astype(float)
                inter_vec = base_vec * log_exposure_stdized  # 使用标准化后的暴露
                X_cols.append(inter_vec.reshape(-1, 1))
                kept_display_names.append(name)
            else:
                # 主效应
                actual_col = covariate_mapping.get(name, name)
                if actual_col not in merged_df.columns:
                    print(f"⚠️  主效应列缺失，跳过: {name} -> {actual_col}")
                    continue
                X_cols.append(merged_df[actual_col].values.astype(float).reshape(-1, 1))
                kept_display_names.append(name)

        if len(X_cols) > 0:
            X = np.hstack(X_cols)
        else:
            X = np.empty((len(merged_df), 0))
        
        # 读取农药类别映射 & 显示名（优先使用mapping）
        category = 'Unknown'
        compound_display = compound
        cid, cname, cat = self._resolve_from_mapping(compound)
        if cname:
            compound_display = cname
        if cat:
            category = cat

        # 构建最终数据字典
        model_data = {
            # 基本信息
            'disease_code': disease_code,
            'compound': compound,
            'compound_display': compound_display,
            'model_type': model_type,
            'lag_years': lag_years,
            'measure_type': measure_type,
            'estimate_type': estimate_type_found, # 使用函数返回的实际找到的类型
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
            'covariate_names': kept_display_names,
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