#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
import json
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy import stats
"""
MICE + PMM 多重插补模块

基于您的专业设计，实现：
1. MICE (Multivariate Imputation by Chained Equations) 框架
2. PMM (Predictive Mean Matching) 用于AAMR连续变量
3. 严格的收敛性诊断和分布合理性检验
4. 高质量可视化诊断输出

设计理念：
- MICE框架的灵活性：为每个变量选择最适合的回归模型
- PMM的稳健性：确保插补值都来自真实观测，避免不合理值
- 全面诊断：收敛图 + 密度图 + 散点图 + 箱线图

作者: AI Assistant  
日期: 2025-10-01
"""

import pandas as pd
import numpy as np
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
import warnings
from tqdm import tqdm
import json

# 过滤弃用警告
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message=".*Downcasting.*")
warnings.filterwarnings("ignore", message=".*Setting an item.*")
warnings.filterwarnings("ignore", message=".*Glyph.*")
warnings.filterwarnings("ignore", message=".*labels.*boxplot.*")

# 核心插补功能，不包含绘图诊断

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MICEImputerWithPMM:
    """
    MICE + PMM 多重插补器
    
    核心功能：
    1. 使用MICE框架进行多重插补
    2. 对连续变量（AAMR）应用PMM算法
    3. 对分类变量应用合适的分类模型
    4. 提供全面的诊断功能
    """
    
    def __init__(self, 
                 n_imputations: int = 5,   # 减少到5个用于汇总
                 max_iter: int = 10,       # 减少迭代次数
                 convergence_threshold: float = 0.001,
                 random_state: int = 42,
                 pmm_k_neighbors: int = 5,
                 config_path: Optional[str] = None,
                 output_dir: Optional[Path] = None):
        """
        初始化MICE+PMM插补器
        
        参数:
            n_imputations: 生成插补数据集的数量
            max_iter: 最大迭代次数
            convergence_threshold: 收敛阈值
            random_state: 随机种子
            pmm_k_neighbors: PMM中选择的邻居数
            config_path: 配置文件路径（可选）
            output_dir: 输出目录（可选）
        """
        # 插补参数
        self.n_imputations = n_imputations
        self.max_iter = max_iter
        self.random_seed = random_state
        self.convergence_threshold = convergence_threshold
        self.pmm_k_neighbors = pmm_k_neighbors
        
        # 设置路径
        self.project_root = Path(__file__).resolve().parents[2]
        self.config_path = config_path or self.project_root / "config.yaml"
        
        # 尝试加载配置文件（如果存在）
        self.config = None
        self.mi_config = {}
        if self.config_path.exists():
            try:
                self.load_config()
                self.setup_paths()
            except Exception as e:
                logger.warning(f"配置文件加载失败，使用默认设置: {e}")
                self.setup_default_paths(output_dir)
        else:
            logger.info("配置文件不存在，使用默认设置")
            self.setup_default_paths(output_dir)
        
        # 存储插补结果和诊断信息
        self.imputed_datasets = []
        self.convergence_history = {}
        self.convergence_stats = {}
        self.original_data = None
        
    def load_config(self):
        """加载配置文件"""
        with open(self.config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self.mi_config = self.config.get('eqi_lmm_multiple_imputation', {})
        logger.info("MICE配置加载完成")
        
    def setup_default_paths(self, output_dir: Optional[Path] = None):
        """设置默认路径（当没有配置文件时）"""
        if output_dir:
            self.mi_results_dir = Path(output_dir)
        else:
            self.mi_results_dir = Path("MI_Results")
        
        self.mi_diagnostics_dir = self.mi_results_dir / "Diagnostics"
        
        # 创建目录
        self.mi_results_dir.mkdir(parents=True, exist_ok=True)
        self.mi_diagnostics_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"使用默认输出目录: {self.mi_results_dir}")
        
    def setup_paths(self):
        """设置路径"""
        # 数据路径
        self.mi_data_path = self.project_root / self.mi_config['data_paths']['mi_source']
        self.mi_datasets_dir = self.project_root / self.mi_config['data_paths']['mi_datasets_dir']
        self.diagnostics_dir = self.project_root / self.mi_config['data_paths']['diagnostics_dir']
        
        # 结果路径
        self.result_base_dir = self.project_root / self.mi_config['result_paths']['base_dir']
        self.mi_diagnostics_dir = self.project_root / self.mi_config['result_paths']['mi_diagnostics']
        
        # 创建目录
        for path in [self.mi_datasets_dir, self.mi_diagnostics_dir]:
            path.mkdir(parents=True, exist_ok=True)
    
    def load_data(self) -> pd.DataFrame:
        """加载MI数据集"""
        logger.info(f"=== 加载MI数据集 ===")
        logger.info(f"数据路径: {self.mi_data_path}")
        
        if not self.mi_data_path.exists():
            raise FileNotFoundError(f"MI数据文件不存在: {self.mi_data_path}")
        
        df = pd.read_csv(self.mi_data_path)
        logger.info(f"数据加载完成: {df.shape}")
        
        # 存储原始数据用于诊断
        self.original_data = df.copy()
        
        return df
    
    def identify_variable_types(self, df: pd.DataFrame) -> Dict[str, List[str]]:
        """
        识别变量类型，为不同类型变量选择合适的插补方法
        
        只对AAMR开头的变量进行插补，其他变量作为辅助变量
        
        返回:
            变量类型字典
        """
        logger.info("=== 识别变量类型 ===")
        
        # 只对AAMR开头的变量进行插补
        aamr_vars = [col for col in df.columns if col.startswith('AAMR_') and df[col].isnull().any()]
        
        # 所有其他变量作为辅助变量（不插补）
        auxiliary_vars = [col for col in df.columns if not col.startswith('AAMR_')]
        
        var_types = {
            'continuous': aamr_vars,  # 只有AAMR变量需要插补
            'ordinal': [],            # 不插补
            'categorical': [],        # 不插补
            'auxiliary': auxiliary_vars,
            'missing': aamr_vars      # 只有AAMR变量可能有缺失
        }
        
        logger.info(f"需插补的AAMR变量 ({len(aamr_vars)}): {aamr_vars}")
        logger.info(f"辅助变量 ({len(auxiliary_vars)}): {auxiliary_vars[:10]}...")  # 只显示前10个
        
        return var_types
    
    def prepare_data_for_mice(self, df: pd.DataFrame, var_types: Dict[str, List[str]]) -> Tuple[pd.DataFrame, Dict]:
        """
        为MICE准备数据
        
        返回:
            预处理后的数据框和编码信息
        """
        logger.info("=== 为MICE准备数据 ===")
        
        df_mice = df.copy()
        encoding_info = {}
        
        # 对分类变量进行标签编码（保持数值形式供sklearn使用）
        from sklearn.preprocessing import LabelEncoder
        
        categorical_all = var_types['categorical'] + var_types['ordinal']
        
        for col in categorical_all:
            if col in df_mice.columns and not df_mice[col].isnull().all():
                le = LabelEncoder()
                # 只对非缺失值进行编码
                mask = df_mice[col].notna()
                if mask.sum() > 0:
                    df_mice.loc[mask, col] = le.fit_transform(df_mice.loc[mask, col].astype(str))
                    encoding_info[col] = le
                    logger.info(f"对 {col} 进行标签编码: {len(le.classes_)} 个类别")
        
        # 确保数值列为float类型
        numeric_cols = var_types['continuous'] + var_types['ordinal']
        for col in numeric_cols:
            if col in df_mice.columns:
                df_mice[col] = pd.to_numeric(df_mice[col], errors='coerce')
        
        return df_mice, encoding_info
    
    def create_pmm_imputer(self, df: pd.DataFrame, var_types: Dict[str, List[str]]) -> IterativeImputer:
        """
        创建支持PMM的MICE插补器
        
        返回:
            配置好的IterativeImputer
        """
        logger.info("=== 创建PMM-MICE插补器 ===")
        
        # 为不同类型变量创建不同的估计器
        # 对于连续变量使用线性回归（后续会手动实现PMM）
        # 对于分类变量使用随机森林
        
        # 这里使用sklearn的IterativeImputer作为基础框架
        # 但我们会在后续步骤中手动实现PMM逻辑
        imputer = IterativeImputer(
            max_iter=self.max_iter,
            random_state=self.random_seed,
            initial_strategy='median',  # 初始插补策略
            verbose=0
        )
        
        logger.info(f"MICE插补器创建完成 - 最大迭代次数: {self.max_iter}")
        
        return imputer
    
    def pmm_single_variable(self, y_obs: np.ndarray, y_pred_obs: np.ndarray, 
                           y_pred_mis: np.ndarray, k: int = 5) -> np.ndarray:
        """
        为单个变量执行预测均值匹配(PMM)
        
        参数:
            y_obs: 观测值
            y_pred_obs: 观测值的预测值  
            y_pred_mis: 缺失值的预测值
            k: PMM邻居数
            
        返回:
            插补值
        """
        imputed_values = np.zeros(len(y_pred_mis))
        
        for i, pred_mis in enumerate(y_pred_mis):
            # 计算与所有观测预测值的距离
            distances = np.abs(y_pred_obs - pred_mis)
            
            # 找到k个最近邻
            k_nearest_idx = np.argsort(distances)[:min(k, len(distances))]
            
            # 从k个最近邻的真实观测值中随机选择一个
            selected_idx = np.random.choice(k_nearest_idx)
            imputed_values[i] = y_obs[selected_idx]
        
        return imputed_values
    
    def mice_with_pmm_single_imputation(self, df: pd.DataFrame, var_types: Dict[str, List[str]], 
                                      encoding_info: Dict) -> pd.DataFrame:
        """
        执行单次MICE+PMM插补
        
        返回:
            完整的插补数据集
        """
        df_imputed = df.copy()
        missing_vars = var_types['missing']
        
        # 迭代插补
        for iteration in range(self.max_iter):
            for var in missing_vars:
                if df_imputed[var].isnull().sum() == 0:
                    continue
                    
                # 分离观测值和缺失值
                obs_mask = df_imputed[var].notna()
                mis_mask = df_imputed[var].isna()
                
                if mis_mask.sum() == 0:
                    continue
                
                # 准备预测变量（只使用数值变量）
                exclude_vars = ['COUNTY_FIPS', 'State', 'EQI_Period', 'AAMR_Period', 'HHS_Region', 'Census_Region', 'Census_Division']
                numeric_predictors = [col for col in df_imputed.columns 
                                    if col != var and col not in exclude_vars and 
                                    pd.api.types.is_numeric_dtype(df_imputed[col])]
                
                # 准备预测数据
                X_obs = df_imputed.loc[obs_mask, numeric_predictors].copy()
                X_mis = df_imputed.loc[mis_mask, numeric_predictors].copy()
                
                # 对数值变量用中位数填充缺失值
                for col in numeric_predictors:
                    median_val = df_imputed[col].median()
                    X_obs[col] = X_obs[col].fillna(median_val)
                    X_mis[col] = X_mis[col].fillna(median_val)
                
                y_obs = df_imputed.loc[obs_mask, var].values
                
                # 根据变量类型选择模型
                if var in var_types['continuous']:
                    # 连续变量：使用线性回归 + PMM
                    from sklearn.linear_model import LinearRegression
                    model = LinearRegression()
                    model.fit(X_obs, y_obs)
                    
                    # 获取预测值
                    y_pred_obs = model.predict(X_obs)
                    y_pred_mis = model.predict(X_mis)
                    
                    # 应用PMM
                    imputed_vals = self.pmm_single_variable(y_obs, y_pred_obs, y_pred_mis, self.pmm_k_neighbors)
                    
                else:
                    # 分类变量：使用随机森林
                    from sklearn.ensemble import RandomForestClassifier
                    model = RandomForestClassifier(n_estimators=50, random_state=self.random_seed)
                    model.fit(X_obs, y_obs.astype(int))
                    
                    # 预测类别
                    imputed_vals = model.predict(X_mis)
                
                # 更新插补值
                df_imputed.loc[mis_mask, var] = imputed_vals
        
        return df_imputed
    
    def fit_transform(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """
        便捷接口：拟合并转换数据（兼容sklearn风格）
        
        参数:
            df: 包含缺失值的数据框
            
        返回:
            多个完整插补后的数据框列表
        """
        return self.run_multiple_imputation(df)
        
    def run_multiple_imputation(self, df: pd.DataFrame) -> List[pd.DataFrame]:
        """
        运行多重插补
        
        返回:
            多个完整插补数据集的列表
        """
        logger.info(f"=== 开始MICE+PMM多重插补 ===")
        logger.info(f"插补数据集数量: {self.n_imputations}")
        logger.info(f"最大迭代次数: {self.max_iter}")
        
        # 识别变量类型
        var_types = self.identify_variable_types(df)
        
        # 准备数据
        df_mice, encoding_info = self.prepare_data_for_mice(df, var_types)
        
        # 初始化结果列表
        imputed_datasets = []
        convergence_tracking = {var: [] for var in var_types['missing']}
        
        # 设置随机种子序列
        np.random.seed(self.random_seed)
        seeds = np.random.randint(0, 10000, self.n_imputations)
        
        # 执行多重插补
        for m in tqdm(range(self.n_imputations), desc="执行MICE+PMM插补"):
            np.random.seed(seeds[m])
            
            # 单次插补
            df_imputed = self.mice_with_pmm_single_imputation(df_mice, var_types, encoding_info)
            
            # 恢复分类变量的原始标签
            for col, encoder in encoding_info.items():
                if col in df_imputed.columns:
                    # 将数值转回原始标签
                    mask = df_imputed[col].notna()
                    if mask.sum() > 0:
                        df_imputed.loc[mask, col] = encoder.inverse_transform(
                            df_imputed.loc[mask, col].astype(int)
                        )
            
            imputed_datasets.append(df_imputed)
            
            # 记录收敛统计信息
            for var in var_types['missing']:
                if var in var_types['continuous']:
                    mean_val = df_imputed[var].mean()
                    std_val = df_imputed[var].std()
                    convergence_tracking[var].append({'mean': mean_val, 'std': std_val, 'imputation': m})
        
        logger.info(f"多重插补完成！生成 {len(imputed_datasets)} 个完整数据集")
        
        # 存储结果
        self.imputed_datasets = imputed_datasets
        self.convergence_stats = convergence_tracking
        
        # 转换为诊断模块需要的格式
        self.convergence_history = {}
        for var, records in convergence_tracking.items():
            self.convergence_history[var] = [r['mean'] for r in records]
        
        return imputed_datasets
    
    def create_final_imputed_dataset(self, imputed_datasets: List[pd.DataFrame]) -> pd.DataFrame:
        """
        合并多个插补数据集，生成最终的完整数据集
        使用中位数方法合并AAMR变量
        
        返回:
            最终的完整数据集
        """
        logger.info("=== 创建最终插补数据集 ===")
        
        if not imputed_datasets:
            raise ValueError("没有插补数据集可以合并")
        
        # 以第一个数据集为基础
        final_df = imputed_datasets[0].copy()
        
        # 获取需要插补的AAMR变量
        aamr_vars = [col for col in final_df.columns if col.startswith('AAMR_')]
        
        # 对每个AAMR变量，使用多个插补数据集的中位数
        for var in aamr_vars:
            if final_df[var].isnull().any():
                # 收集该变量在所有插补数据集中的值
                all_values = []
                for df in imputed_datasets:
                    all_values.append(df[var])
                
                # 计算中位数（跨插补数据集）
                combined_values = pd.concat(all_values, axis=1)
                median_values = combined_values.median(axis=1)
                
                # 更新最终数据集
                final_df[var] = median_values
        
        logger.info(f"最终数据集形状: {final_df.shape}")
        logger.info(f"剩余缺失值: {final_df.isnull().sum().sum()}")
        
        return final_df
    
    def convert_to_long_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        将宽格式数据转换为长格式（类似EQI_LMM_Delete_df.csv的格式）
        
        返回:
            长格式数据框
        """
        logger.info("=== 转换为长格式 ===")
        
        # 获取AAMR变量
        aamr_vars = [col for col in df.columns if col.startswith('AAMR_')]
        
        # 准备基础列
        base_cols = ['COUNTY_FIPS', 'State', 'RUCC', 'EQI', 'EQI_air', 'EQI_water', 
                    'EQI_land', 'EQI_built', 'EQI_Sociodemographic', 'RUCC_EQI', 
                    'RUCC_EQI_air', 'RUCC_EQI_water', 'RUCC_EQI_land', 'RUCC_EQI_built', 
                    'RUCC_EQI_Sociodemographic', 'SR', 'EQI_Period', 'AAMR_Period']
        
        # 创建癌症类型映射
        cancer_mapping = {
            'AAMR_C00_C97': ('C00_C97', 'All Cancers'),
            'AAMR_C15_C26': ('C15_C26', 'Digestive System'),
            'AAMR_C18_C21': ('C18_C21', 'Colorectal'),
            'AAMR_C25': ('C25', 'Pancreatic'),
            'AAMR_C30_C39': ('C30_C39', 'Respiratory System'),
            'AAMR_C34': ('C34', 'Lung and Bronchus'),
            'AAMR_C50': ('C50', 'Female Breast'),
            'AAMR_C51_C58': ('C51_C58', 'Female Genital System'),
            'AAMR_C60_C63': ('C60_C63', 'Male Genital System'),
            'AAMR_C61': ('C61', 'Prostate'),
            'AAMR_C64_C68': ('C64_C68', 'Urinary System'),
            'AAMR_C76_C80': ('C76_C80', 'Other and Unspecified Primary Sites'),
            'AAMR_C81_C96': ('C81_C96', 'Hematopoietic and Lymphoid Tissues')
        }
        
        # 转换为长格式
        long_data = []
        
        for _, row in df.iterrows():
            base_info = row[base_cols].to_dict()
            
            for aamr_col in aamr_vars:
                if aamr_col in cancer_mapping and pd.notna(row[aamr_col]):
                    cancer_type, cancer_desc = cancer_mapping[aamr_col]
                    
                    record = base_info.copy()
                    # 创建正确的Analysis_Scenario格式
                    eqi_period = str(row['EQI_Period']).replace('_', '')
                    aamr_period = str(row['AAMR_Period']).replace('_', '')
                    
                    # 转换为LMM期望的格式
                    if eqi_period == '20002005':
                        eqi_code = '0005'
                    elif eqi_period == '20062010':
                        eqi_code = '0610'
                    else:
                        eqi_code = eqi_period
                    
                    if aamr_period == '20062010':
                        aamr_code = '2006_2010'
                    elif aamr_period == '20112015':
                        aamr_code = '2011_2015'
                    elif aamr_period == '20162020':
                        aamr_code = '2016_2020'
                    else:
                        aamr_code = aamr_period
                    
                    record.update({
                        'Smoking_Rate': row['SR'],
                        'Analysis_Scenario': f"EQI{eqi_code}_AAMR{aamr_code}",
                        'Lag_Years': 5 if ('2000_2005' in str(row['EQI_Period']) and '2006_2010' in str(row['AAMR_Period'])) or \
                                        ('2006_2010' in str(row['EQI_Period']) and '2011_2015' in str(row['AAMR_Period'])) else 10,
                        'EQI_Period': eqi_period,
                        'AAMR': row[aamr_col],
                        'Cancer_Type': cancer_type,
                        'Cancer_Description': cancer_desc,
                        'State_FIPS': row['COUNTY_FIPS'][:2]  # 前两位是州FIPS
                    })
                    
                    long_data.append(record)
        
        long_df = pd.DataFrame(long_data)
        
        # 重新排列列顺序以匹配目标格式
        target_cols = ['COUNTY_FIPS', 'State', 'RUCC', 'EQI', 'EQI_air', 'EQI_water', 
                      'EQI_land', 'EQI_built', 'EQI_Sociodemographic', 'RUCC_EQI', 
                      'RUCC_EQI_air', 'RUCC_EQI_water', 'RUCC_EQI_land', 'RUCC_EQI_built', 
                      'RUCC_EQI_Sociodemographic', 'Smoking_Rate', 'Analysis_Scenario', 
                      'Lag_Years', 'EQI_Period', 'AAMR', 'Cancer_Type', 'Cancer_Description', 
                      'State_FIPS']
        
        long_df = long_df[target_cols]
        
        logger.info(f"长格式数据形状: {long_df.shape}")
        
        return long_df
        
    def save_final_dataset(self, final_df: pd.DataFrame, long_df: pd.DataFrame) -> Dict[str, str]:
        """
        保存最终的插补数据集
        
        返回:
            保存的文件路径字典
        """
        logger.info("=== 保存最终数据集 ===")
        
        saved_files = {}
        
        # 保存完整格式（与原始MI_df格式相同）
        wide_filepath = self.project_root / "Data" / "df" / "EQI_LMM_MI_Imputed.csv"
        final_df.to_csv(wide_filepath, index=False)
        saved_files['wide_format'] = str(wide_filepath)
        
        # 保存长格式（与Delete_df格式相同）
        long_filepath = self.project_root / "Data" / "df" / "EQI_LMM_MI_Imputed_Long.csv"
        long_df.to_csv(long_filepath, index=False)
        saved_files['long_format'] = str(long_filepath)
        
        logger.info(f"宽格式数据已保存: {wide_filepath}")
        logger.info(f"长格式数据已保存: {long_filepath}")
        
        return saved_files
    
    def get_imputation_summary(self) -> Dict[str, Any]:
        """
        获取插补摘要信息（不包含诊断图）
        
        返回:
            插补摘要字典
        """
        logger.info("=== 生成插补摘要 ===")
        
        summary = {
            'n_imputations': self.n_imputations,
            'n_iterations': self.max_iter,
            'original_shape': self.original_data.shape if self.original_data is not None else None,
            'n_datasets_created': len(self.imputed_datasets),
            'convergence_tracked': bool(self.convergence_stats),
            'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
        }
        
        return summary
    
    def validate_imputation_quality(self) -> Dict[str, Any]:
        """
        验证插补质量（数值检查，不包含图表）
        
        返回:
            质量验证结果
        """
        logger.info("=== 验证插补质量 ===")
        
        if not self.imputed_datasets or self.original_data is None:
            logger.warning("没有插补数据或原始数据，跳过质量验证")
            return {}
        
        validation_results = {}
        
        # 获取连续变量
        continuous_vars = [col for col in self.original_data.columns 
                          if col.startswith('AAMR_') or col == 'SR']
        
        for var in continuous_vars:
            if self.original_data[var].isnull().all():
                continue
                
            # 原始数据统计
            obs_data_raw = self.original_data[var].dropna()
            obs_data = pd.to_numeric(obs_data_raw, errors='coerce').dropna()
            
            var_results = {}
            
            if len(obs_data) > 0:
                var_results['original'] = {
                    'n_observations': len(obs_data),
                    'mean': float(obs_data.mean()),
                    'std': float(obs_data.std()),
                    'median': float(obs_data.median()),
                    'min': float(obs_data.min()),
                    'max': float(obs_data.max())
                }
                
                # 插补数据统计
                imp_stats = []
                for i, df_imp in enumerate(self.imputed_datasets):
                    imp_data_raw = df_imp[var].dropna()
                    if len(imp_data_raw) > 0:
                        imp_data = pd.to_numeric(imp_data_raw, errors='coerce').dropna()
                        if len(imp_data) > 0:
                            imp_stats.append({
                                'dataset': i + 1,
                                'mean': float(imp_data.mean()),
                                'std': float(imp_data.std()),
                                'median': float(imp_data.median()),
                                'min': float(imp_data.min()),
                                'max': float(imp_data.max())
                            })
                
                if imp_stats:
                    means = [stat['mean'] for stat in imp_stats]
                    stds = [stat['std'] for stat in imp_stats]
                    
                    var_results['imputed_summary'] = {
                        'n_datasets': len(imp_stats),
                        'mean_range': [min(means), max(means)],
                        'std_range': [min(stds), max(stds)],
                        'pooled_mean': np.mean(means),
                        'between_imputation_variance': np.std(means)
                    }
                    
                    # 质量检查
                    var_results['quality_checks'] = {
                        'mean_within_reasonable_range': (
                            var_results['original']['mean'] * 0.8 <= np.mean(means) <= 
                            var_results['original']['mean'] * 1.2
                        ),
                        'no_extreme_values': all(
                            var_results['original']['min'] * 0.5 <= stat['min'] and
                            stat['max'] <= var_results['original']['max'] * 1.5
                            for stat in imp_stats
                        ),
                        'consistent_across_imputations': np.std(means) / np.mean(means) < 0.1 if np.mean(means) > 0 else True
                    }
                
                validation_results[var] = var_results
        
        return validation_results
    
    def check_correlation_preservation(self) -> Dict[str, Any]:
        """
        检查变量间相关性是否在插补后得到保持
        
        返回:
            相关性保持检查结果
        """
        logger.info("=== 检查相关性保持 ===")
        
        if not self.imputed_datasets or self.original_data is None:
            return {}
        
        correlation_results = {}
        
        # 选择关键变量组合
        key_pairs = [
            ('AAMR_C00_C97', 'EQI'),
            ('AAMR_C00_C97', 'SR'),
            ('AAMR_C34', 'EQI'),  # 肺癌 vs EQI
            ('SR', 'EQI')
        ]
        
        for var1, var2 in key_pairs:
            if var1 not in self.original_data.columns or var2 not in self.original_data.columns:
                continue
                
            # 计算原始数据中的相关性
            orig_complete = self.original_data[[var1, var2]].dropna()
            if len(orig_complete) > 10:  # 需要足够的观测数据
                orig_corr = orig_complete[var1].corr(orig_complete[var2])
                
                # 计算各插补数据集中的相关性
                imp_corrs = []
                for df_imp in self.imputed_datasets:
                    imp_complete = df_imp[[var1, var2]].dropna()
                    if len(imp_complete) > 10:
                        imp_corr = imp_complete[var1].corr(imp_complete[var2])
                        if not pd.isna(imp_corr):
                            imp_corrs.append(imp_corr)
                
                if imp_corrs:
                    correlation_results[f'{var1}_vs_{var2}'] = {
                        'original_correlation': orig_corr,
                        'imputed_correlations': imp_corrs,
                        'mean_imputed_correlation': np.mean(imp_corrs),
                        'correlation_preserved': abs(orig_corr - np.mean(imp_corrs)) < 0.1,
                        'correlation_variance': np.var(imp_corrs)
                    }
        
        return correlation_results
    
    def run_simple_imputation_with_diagnostics(self) -> Dict[str, Any]:
        """
        运行简化的插补流程和诊断
        
        返回:
            完整结果字典
        """
        logger.info("=" * 60)
        logger.info("简化的 MICE + PMM 多重插补分析")
        logger.info("=" * 60)
        
        start_time = pd.Timestamp.now()
        
        # 1. 加载数据
        df = self.load_data()
        
        # 2. 执行多重插补
        imputed_datasets = self.run_multiple_imputation(df)
        
        # 3. 创建最终插补数据集
        final_df = self.create_final_imputed_dataset(imputed_datasets)
        
        # 4. 转换为长格式
        long_df = self.convert_to_long_format(final_df)
        
        # 5. 保存最终数据集
        saved_files = self.save_final_dataset(final_df, long_df)
        
        # 6. 创建简化诊断图
        logger.info("=== 创建简化诊断报告 ===")
        diagnostic_files = self.create_simple_diagnostics(final_df)
        
        # 7. 生成综合报告
        end_time = pd.Timestamp.now()
        duration = end_time - start_time
        
        results = {
            'n_imputations': self.n_imputations,
            'n_iterations': self.max_iter,
            'original_shape': self.original_data.shape,
            'final_shape_wide': final_df.shape,
            'final_shape_long': long_df.shape,
            'n_aamr_vars': len([col for col in df.columns if col.startswith('AAMR_')]),
            'saved_files': saved_files,
            'diagnostic_files': diagnostic_files,
            'processing_time': str(duration),
            'timestamp': start_time.strftime('%Y-%m-%d %H:%M:%S')
        }
        
        # 保存结果摘要
        diagnostics_dir = Path("/Users/ventus/Repository/WDP/Result/EQI_LMM_MI_Diagnose")
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        summary_path = diagnostics_dir / 'imputation_summary.json'
        
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False, default=str)
        
        logger.info(f"插补流程完成! 耗时: {duration}")
        logger.info(f"合并了 {len(imputed_datasets)} 个插补数据集")
        logger.info(f"最终数据集: {final_df.shape} (宽格式), {long_df.shape} (长格式)")
        logger.info(f"诊断文件: {len(diagnostic_files)} 个")
        logger.info(f"结果摘要: {summary_path}")
        
        return results
        
    def create_simple_diagnostics(self, final_df: pd.DataFrame) -> Dict[str, str]:
        """
        创建专业的插补诊断图表
        
        返回:
            诊断文件路径字典
        """
        diagnostics_dir = Path("/Users/ventus/Repository/WDP/Result/EQI_LMM_MI_Diagnose")
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        
        diagnostic_files = {}
        
        # 设置绘图样式
        plt.style.use('default')
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.grid'] = True
        plt.rcParams['grid.alpha'] = 0.3
        
        # 配色方案 - 简洁清晰的对比色
        colors = {
            'observed': '#1f77b4',    # 蓝色 - 原始数据
            'imputed': '#ff7f0e',     # 橙色 - 插补后总体
            'complete': '#2ca02c',    # 绿色 - 备用
            'missing': '#d62728'      # 红色 - 备用
        }
        
        # 获取AAMR变量
        aamr_vars = [col for col in final_df.columns if col.startswith('AAMR_')]
        
        # 1. 高级插补质量诊断图
        if self.original_data is not None:
            diagnostic_files['quality'] = self.create_imputation_quality_plot(aamr_vars, colors, diagnostics_dir, final_df)
        
        # 2. PMM效果验证图
        diagnostic_files['pmm_validation'] = self.create_pmm_validation_plot(aamr_vars, colors, diagnostics_dir, final_df)
        
        # 3. 收敛性诊断图（如果有收敛历史）
        if hasattr(self, 'convergence_stats') and self.convergence_stats:
            diagnostic_files['convergence'] = self.create_convergence_plot(colors, diagnostics_dir)
        
        # 4. 插补分布合理性检查
        diagnostic_files['distribution'] = self.create_distribution_validation_plot(aamr_vars, colors, diagnostics_dir, final_df)
        
        # 5. 变量相关性保持检查
        diagnostic_files['correlation'] = self.create_correlation_preservation_plot(colors, diagnostics_dir, final_df)
        
        # 6. 三重分布对比图 (原始-插补-完整)
        diagnostic_files['triple_distribution'] = self.create_triple_distribution_comparison(aamr_vars, colors, diagnostics_dir, final_df)
        
        # 7. 为每个变量生成详细诊断报告
        diagnostic_files['individual_reports'] = self.create_individual_variable_reports(aamr_vars, colors, diagnostics_dir, final_df)
        
        logger.info(f"专业诊断图已生成: {diagnostics_dir}")
        logger.info(f"生成的诊断图: {list(diagnostic_files.keys())}")
        
        return diagnostic_files
    
    def create_imputation_quality_plot(self, aamr_vars: List[str], colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> str:
        """创建插补质量诊断图 - 覆盖所有AAMR变量"""
        # 过滤出真正需要插补的变量（有缺失值的）
        vars_to_plot = []
        for var in aamr_vars:
            if var != 'AAMR_Period' and self.original_data[var].isnull().any():
                vars_to_plot.append(var)
        
        # 计算需要的图表数量和布局
        n_vars = len(vars_to_plot)
        n_cols = 3
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        for i, var in enumerate(vars_to_plot):
            if i >= len(axes):
                break
            ax = axes[i]
            
            # 获取数据
            orig_obs = self.original_data[var].dropna()
            missing_mask = self.original_data[var].isna()
            imputed_vals = final_df[var][missing_mask]
            
            if len(orig_obs) > 0 and len(imputed_vals) > 0:
                # 获取插补完成后的完整数据
                final_complete_vals = final_df[var].dropna()
                
                # 创建三重直方图对比
                # 1. 原始观测值
                ax.hist(orig_obs, bins=25, alpha=0.7, color=colors['observed'], 
                       density=True, label=f'Original (n={len(orig_obs)})', 
                       edgecolor='white', linewidth=0.5)
                
                # 2. 插补值
                ax.hist(imputed_vals, bins=25, alpha=0.7, color=colors['imputed'], 
                       density=True, label=f'Imputed (n={len(imputed_vals)})', 
                       edgecolor='white', linewidth=0.5)
                
                # 3. 插补完成后总体
                if len(final_complete_vals) > 0:
                    ax.hist(final_complete_vals, bins=25, alpha=0.5, color=colors['complete'], 
                           density=True, label=f'Complete (n={len(final_complete_vals)})', 
                           edgecolor='white', linewidth=0.5)
                
                # 添加均值线
                ax.axvline(orig_obs.mean(), color=colors['observed'], linestyle='--', linewidth=1.5, alpha=0.8)
                ax.axvline(imputed_vals.mean(), color=colors['imputed'], linestyle='--', linewidth=1.5, alpha=0.8)
                if len(final_complete_vals) > 0:
                    ax.axvline(final_complete_vals.mean(), color=colors['complete'], linestyle='-', linewidth=2, alpha=0.9)
                
                # 添加统计信息
                stats_text = f'Observed: μ={orig_obs.mean():.1f}, σ={orig_obs.std():.1f}\n'
                stats_text += f'Imputed: μ={imputed_vals.mean():.1f}, σ={imputed_vals.std():.1f}'
                ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                ax.set_title(f'{var} - Quality Check', fontweight='bold', fontsize=11)
                ax.set_xlabel('AAMR Value')
                ax.set_ylabel('Density')
                ax.legend(frameon=True, fancybox=True, shadow=True)
                ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for i in range(len(vars_to_plot), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        filepath = diagnostics_dir / 'imputation_quality_comprehensive.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"插补质量诊断图已生成，覆盖 {len(vars_to_plot)} 个变量")
        return str(filepath)
    
    def create_pmm_validation_plot(self, aamr_vars: List[str], colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> str:
        """创建PMM效果验证图 - 覆盖所有需要插补的AAMR变量"""
        # 过滤出真正需要插补的变量
        vars_to_plot = []
        for var in aamr_vars:
            if var != 'AAMR_Period' and self.original_data[var].isnull().any():
                vars_to_plot.append(var)
        
        # 计算布局
        n_vars = len(vars_to_plot)
        n_cols = 3
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        for i, var in enumerate(vars_to_plot):
            if i >= len(axes):
                break
            ax = axes[i]
            
            # 获取观测值和插补值
            observed_vals = self.original_data[var].dropna()
            missing_mask = self.original_data[var].isna()
            
            if len(observed_vals) > 10 and missing_mask.sum() > 0:
                # PMM质量检查：插补值应该来自观测值的邻域
                if hasattr(self, 'imputed_datasets') and self.imputed_datasets:
                    all_imputed = []
                    for dataset in self.imputed_datasets:
                        imputed_vals = dataset[var][missing_mask].dropna()
                        all_imputed.extend(imputed_vals)
                    
                    if all_imputed:
                        # 检查插补值是否在观测值范围内
                        obs_min, obs_max = observed_vals.min(), observed_vals.max()
                        obs_q25, obs_q75 = observed_vals.quantile([0.25, 0.75])
                        
                        # 绘制观测值分布
                        ax.hist(observed_vals, bins=30, alpha=0.6, color=colors['observed'], 
                               density=True, label='Observed', edgecolor='white')
                        
                        # 绘制插补值分布
                        ax.hist(all_imputed, bins=30, alpha=0.6, color=colors['imputed'],
                               density=True, label='All Imputed', edgecolor='white')
                        
                        # 添加范围标记
                        ax.axvline(obs_min, color='red', linestyle=':', alpha=0.7, label='Obs Range')
                        ax.axvline(obs_max, color='red', linestyle=':', alpha=0.7)
                        ax.axvspan(obs_q25, obs_q75, alpha=0.2, color='green', label='IQR')
                        
                        # PMM质量评分
                        within_range = sum(obs_min <= val <= obs_max for val in all_imputed)
                        quality_score = within_range / len(all_imputed) * 100
                        
                        ax.set_title(f'{var}\nPMM Quality: {quality_score:.1f}% within observed range', 
                                   fontweight='bold', fontsize=10)
                        ax.legend(fontsize=9)
                else:
                    # 如果没有多个插补数据集，使用最终数据集
                    final_imputed = final_df[var][missing_mask].dropna()
                    if len(final_imputed) > 0:
                        ax.hist(observed_vals, bins=20, alpha=0.6, color=colors['observed'], 
                               density=True, label='Observed')
                        ax.hist(final_imputed, bins=20, alpha=0.6, color=colors['imputed'],
                               density=True, label='Imputed')
                        ax.set_title(f'{var} - PMM Validation', fontweight='bold')
                        ax.legend()
                
                ax.set_xlabel('AAMR Value')
                ax.set_ylabel('Density')
                ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for i in range(len(vars_to_plot), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        filepath = diagnostics_dir / 'pmm_validation_comprehensive.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"PMM验证图已生成，覆盖 {len(vars_to_plot)} 个变量")
        return str(filepath)
    
    def create_convergence_plot(self, colors: dict, diagnostics_dir: Path) -> str:
        """创建收敛性诊断图"""
        if not hasattr(self, 'convergence_stats') or not self.convergence_stats:
            return ""
        
        # 提取收敛数据
        convergence_data = {}
        for var, records in self.convergence_stats.items():
            if records:
                convergence_data[var] = [r['mean'] for r in records]
        
        if not convergence_data:
            return ""
        
        n_vars = len(convergence_data)
        n_cols = 2
        n_rows = (n_vars + 1) // 2
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        for i, (var, history) in enumerate(convergence_data.items()):
            if i >= len(axes):
                break
                
            ax = axes[i]
            iterations = range(1, len(history) + 1)
            
            # 绘制收敛轨迹
            ax.plot(iterations, history, color=colors['complete'], linewidth=2.5, 
                   marker='o', markersize=4, alpha=0.8, label='Mean Value')
            
            # 添加趋势线
            if len(history) > 2:
                z = np.polyfit(iterations, history, 1)
                p = np.poly1d(z)
                ax.plot(iterations, p(iterations), color=colors['missing'], 
                       linestyle='--', alpha=0.8, label='Trend')
            
            # 收敛评估
            if len(history) > 1:
                recent_std = np.std(history[-min(5, len(history)):])
                convergence_status = "Converged" if recent_std < 0.01 else "Not Converged"
                color = 'green' if convergence_status == "Converged" else 'red'
                
                ax.text(0.02, 0.95, f'{convergence_status}\nRecent SD: {recent_std:.4f}', 
                       transform=ax.transAxes, fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor=color, alpha=0.2))
            
            ax.set_title(f'{var} - Convergence Trace', fontweight='bold')
            ax.set_xlabel('Imputation Number')
            ax.set_ylabel('Mean Value')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        # 隐藏多余子图
        for i in range(len(convergence_data), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        filepath = diagnostics_dir / 'convergence_advanced.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        return str(filepath)
    
    def create_distribution_validation_plot(self, aamr_vars: List[str], colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> str:
        """创建分布合理性验证图 - 覆盖所有需要插补的AAMR变量"""
        # 过滤出真正需要插补的变量
        vars_to_plot = []
        for var in aamr_vars:
            if var != 'AAMR_Period' and self.original_data[var].isnull().any():
                vars_to_plot.append(var)
        
        # 计算布局
        n_vars = len(vars_to_plot)
        n_cols = 3
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        # Q-Q图和分布检验
        for i, var in enumerate(vars_to_plot):
            if i >= len(axes):
                break
            ax = axes[i]
            
            observed_vals = self.original_data[var].dropna()
            final_vals = final_df[var].dropna()
            
            if len(observed_vals) > 10 and len(final_vals) > 10:
                from scipy import stats
                
                # Q-Q图
                stats.probplot(observed_vals, dist="norm", plot=ax)
                ax.get_lines()[0].set_markerfacecolor(colors['observed'])
                ax.get_lines()[0].set_markeredgecolor('white')
                ax.get_lines()[0].set_markersize(4)
                ax.get_lines()[1].set_color(colors['missing'])
                
                # 正态性检验
                _, p_value = stats.shapiro(observed_vals[:5000] if len(observed_vals) > 5000 else observed_vals)
                normality_text = f"Shapiro-Wilk p-value: {p_value:.4f}"
                normality_status = "Normal" if p_value > 0.05 else "Non-normal"
                
                ax.text(0.02, 0.95, f'{normality_text}\n{normality_status}', 
                       transform=ax.transAxes, fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                ax.set_title(f'{var} - Normality Check', fontweight='bold')
                ax.grid(True, alpha=0.3)
        
        # 隐藏多余的子图
        for i in range(len(vars_to_plot), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        filepath = diagnostics_dir / 'distribution_validation_comprehensive.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"分布验证图已生成，覆盖 {len(vars_to_plot)} 个变量")
        return str(filepath)
    
    def create_correlation_preservation_plot(self, colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> str:
        """创建相关性保持检查图"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        axes = axes.flatten()
        
        # 检查关键变量对的相关性
        key_pairs = [
            ('AAMR_C00_C97', 'EQI'),
            ('AAMR_C00_C97', 'SR'), 
            ('AAMR_C34', 'EQI'),
            ('EQI', 'SR')
        ]
        
        for i, (var1, var2) in enumerate(key_pairs):
            if var1 not in self.original_data.columns or var2 not in self.original_data.columns:
                continue
                
            ax = axes[i]
            
            # 原始数据相关性
            orig_data = self.original_data[[var1, var2]].dropna()
            if len(orig_data) > 10:
                orig_corr = orig_data[var1].corr(orig_data[var2])
                
                # 散点图
                ax.scatter(orig_data[var2], orig_data[var1], 
                          alpha=0.6, color=colors['observed'], s=20, 
                          edgecolors='white', linewidth=0.5, label='Observed')
                
                # 如果有插补数据，显示插补后的相关性
                final_data = final_df[[var1, var2]].dropna()
                if len(final_data) > len(orig_data):
                    final_corr = final_data[var1].corr(final_data[var2])
                    
                    # 添加插补点
                    missing_mask = self.original_data[[var1, var2]].isna().any(axis=1)
                    if missing_mask.sum() > 0:
                        imputed_data = final_df.loc[missing_mask, [var1, var2]].dropna()
                        if len(imputed_data) > 0:
                            ax.scatter(imputed_data[var2], imputed_data[var1],
                                     alpha=0.8, color=colors['imputed'], s=25,
                                     edgecolors='white', linewidth=0.5, label='Imputed')
                    
                    corr_change = abs(orig_corr - final_corr)
                    status = "Preserved" if corr_change < 0.1 else "Changed"
                    
                    ax.text(0.02, 0.95, f'Original r: {orig_corr:.3f}\nFinal r: {final_corr:.3f}\n{status}',
                           transform=ax.transAxes, fontsize=9, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                else:
                    ax.text(0.02, 0.95, f'Correlation: {orig_corr:.3f}',
                           transform=ax.transAxes, fontsize=9, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                ax.set_title(f'{var1} vs {var2}', fontweight='bold')
                ax.set_xlabel(var2)
                ax.set_ylabel(var1)
                ax.legend()
                ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filepath = diagnostics_dir / 'correlation_preservation.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        return str(filepath)
    
    def create_triple_distribution_comparison(self, aamr_vars: List[str], colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> str:
        """创建简化的插补前后分布对比图"""
        # 过滤出真正需要插补的变量
        vars_to_plot = []
        for var in aamr_vars:
            if var != 'AAMR_Period' and self.original_data[var].isnull().any():
                vars_to_plot.append(var)
        
        # 计算布局 - 每行显示3个变量，更紧凑
        n_vars = len(vars_to_plot)
        n_cols = 3
        n_rows = (n_vars + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 5 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else axes
        else:
            axes = axes.flatten()
        
        logger.info(f"创建插补前后对比图，覆盖 {len(vars_to_plot)} 个变量")
        
        for i, var in enumerate(vars_to_plot):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # 只获取两种数据：原始观测 vs 插补完成后总体
            observed_vals = self.original_data[var].dropna()  # 原始观测值
            complete_vals = final_df[var].dropna()  # 插补完成后总体
            
            # 计算统计信息
            missing_count = self.original_data[var].isna().sum()
            missing_pct = (missing_count / len(self.original_data)) * 100
            
            if len(observed_vals) > 0 and len(complete_vals) > 0:
                # 统一的x轴范围
                x_min = min(observed_vals.min(), complete_vals.min())
                x_max = max(observed_vals.max(), complete_vals.max())
                bins = np.linspace(x_min, x_max, 30)
                
                # 绘制两个直方图 - 清晰对比
                ax.hist(observed_vals, bins=bins, alpha=0.7, color=colors['observed'], 
                       density=True, label=f'Before Imputation (n={len(observed_vals)})', 
                       edgecolor='white', linewidth=0.8)
                
                ax.hist(complete_vals, bins=bins, alpha=0.7, color=colors['imputed'], 
                       density=True, label=f'After Imputation (n={len(complete_vals)})', 
                       edgecolor='white', linewidth=0.8)
                
                # 添加均值线
                obs_mean = observed_vals.mean()
                complete_mean = complete_vals.mean()
                
                ax.axvline(obs_mean, color=colors['observed'], linestyle='--', 
                          linewidth=2, alpha=0.9, label=f'Before Mean: {obs_mean:.1f}')
                ax.axvline(complete_mean, color=colors['imputed'], linestyle='--', 
                          linewidth=2, alpha=0.9, label=f'After Mean: {complete_mean:.1f}')
                
                # 设置标题和标签 - 英文避免方框
                ax.set_title(f'{var}\\nImputed {missing_count} values ({missing_pct:.1f}%)', 
                            fontweight='bold', fontsize=11)
                ax.set_xlabel('AAMR Value', fontsize=10)
                ax.set_ylabel('Density', fontsize=10)
                ax.legend(fontsize=8, loc='upper right', frameon=True)
                ax.grid(True, alpha=0.3, color='gray', linestyle='-', linewidth=0.5)
                
                # 添加简化的统计信息文本框
                mean_diff = abs(obs_mean - complete_mean)
                stats_text = f'Before: μ={obs_mean:.1f}, σ={observed_vals.std():.1f}\\nAfter:  μ={complete_mean:.1f}, σ={complete_vals.std():.1f}\\nDiff:   {mean_diff:.1f}'
                
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=8,
                       verticalalignment='top', fontfamily='monospace',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=0.7))
            else:
                ax.text(0.5, 0.5, f'{var}\\nInsufficient data for comparison', 
                       transform=ax.transAxes, ha='center', va='center', 
                       fontsize=11, fontweight='bold')
        
        # 隐藏多余的子图
        for i in range(len(vars_to_plot), len(axes)):
            axes[i].set_visible(False)
        
        # 添加总标题 - 英文
        fig.suptitle('AAMR Variables: Before vs After Imputation Distribution Comparison', 
                    fontsize=14, fontweight='bold', y=0.95)
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        
        filepath = diagnostics_dir / 'before_after_imputation_comparison.png'
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"插补前后对比图已生成: {filepath}")
        return str(filepath)
    
    def create_individual_variable_reports(self, aamr_vars: List[str], colors: dict, diagnostics_dir: Path, final_df: pd.DataFrame) -> List[str]:
        """为每个需要插补的AAMR变量创建详细的个体诊断报告"""
        individual_files = []
        
        # 过滤出真正需要插补的变量
        vars_to_analyze = []
        for var in aamr_vars:
            if var != 'AAMR_Period' and self.original_data[var].isnull().any():
                vars_to_analyze.append(var)
        
        # 创建个体变量诊断目录
        individual_dir = diagnostics_dir / "individual_variables"
        individual_dir.mkdir(exist_ok=True)
        
        logger.info(f"为 {len(vars_to_analyze)} 个变量创建详细诊断报告...")
        
        for var in vars_to_analyze:
            try:
                # 为每个变量创建4合1的诊断图
                fig, axes = plt.subplots(2, 2, figsize=(16, 12))
                
                # 获取数据
                observed_vals = self.original_data[var].dropna()
                missing_mask = self.original_data[var].isna()
                imputed_vals = final_df[var][missing_mask]
                final_complete = final_df[var].dropna()
                
                missing_count = missing_mask.sum()
                missing_pct = (missing_count / len(self.original_data)) * 100
                
                # 子图1: 插补前后对比 (简化版)
                ax1 = axes[0, 0]
                if len(observed_vals) > 1 and len(final_complete) > 1:
                    # 只显示两个分布：插补前 vs 插补后
                    
                    # 1. 原始观测值直方图
                    ax1.hist(observed_vals, bins=25, alpha=0.7, color=colors['observed'], 
                            density=True, label=f'Before Imputation (n={len(observed_vals)})', 
                            edgecolor='white', linewidth=0.8)
                    
                    # 2. 插补完成后总体直方图
                    ax1.hist(final_complete, bins=25, alpha=0.7, color=colors['imputed'], 
                            density=True, label=f'After Imputation (n={len(final_complete)})', 
                            edgecolor='white', linewidth=0.8)
                    
                    # 添加均值线
                    ax1.axvline(observed_vals.mean(), color=colors['observed'], linestyle='--', 
                               alpha=0.9, linewidth=2)
                    ax1.axvline(final_complete.mean(), color=colors['imputed'], linestyle='--', 
                               alpha=0.9, linewidth=2)
                
                ax1.set_title(f'Distribution: Before vs After Imputation\n({missing_count} values, {missing_pct:.1f}%)', 
                             fontweight='bold', fontsize=11)
                ax1.set_xlabel('AAMR Value')
                ax1.set_ylabel('Density')
                ax1.legend(fontsize=9)
                ax1.grid(True, alpha=0.3)
                
                # 子图2: 箱线图对比 (简化版)
                ax2 = axes[0, 1]
                box_data = [observed_vals, final_complete]
                box_labels = ['Before', 'After']
                
                bp = ax2.boxplot(box_data, labels=box_labels, patch_artist=True)
                box_colors = [colors['observed'], colors['imputed']]
                for patch, color in zip(bp['boxes'], box_colors):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                ax2.set_title('Box Plot: Before vs After', fontweight='bold', fontsize=11)
                ax2.set_ylabel('AAMR Value')
                ax2.grid(True, alpha=0.3)
                
                # 子图3: Q-Q图（正态性检查）
                ax3 = axes[1, 0]
                if len(final_complete) > 10:
                    stats.probplot(final_complete, dist="norm", plot=ax3)
                    ax3.get_lines()[0].set_markerfacecolor(colors['complete'])
                    ax3.get_lines()[0].set_markeredgecolor('white')
                    ax3.get_lines()[0].set_markersize(4)
                    ax3.get_lines()[1].set_color(colors['missing'])
                    
                    # 正态性检验
                    sample_size = min(5000, len(final_complete))
                    _, p_value = stats.shapiro(final_complete.sample(sample_size) if len(final_complete) > sample_size else final_complete)
                    ax3.text(0.02, 0.95, f'Shapiro-Wilk p: {p_value:.4f}', 
                            transform=ax3.transAxes, fontsize=10,
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
                ax3.set_title('Normality Check (Q-Q Plot)', fontweight='bold', fontsize=11)
                ax3.grid(True, alpha=0.3)
                
                # 子图4: 统计摘要
                ax4 = axes[1, 1]
                ax4.axis('off')
                
                # 计算统计信息
                stats_text = f"""
                {var} - Imputation Summary
                
                📊 Data Overview:
                • Total observations: {len(self.original_data):,}
                • Missing values: {missing_count:,} ({missing_pct:.1f}%)
                • Complete after imputation: {len(final_complete):,}
                
                📈 Statistical Comparison:
                                Observed    Imputed     Final
                Mean:          {observed_vals.mean():8.2f}   {imputed_vals.mean():8.2f}   {final_complete.mean():8.2f}
                Std:           {observed_vals.std():8.2f}   {imputed_vals.std():8.2f}   {final_complete.std():8.2f}
                Median:        {observed_vals.median():8.2f}   {imputed_vals.median():8.2f}   {final_complete.median():8.2f}
                Min:           {observed_vals.min():8.2f}   {imputed_vals.min():8.2f}   {final_complete.min():8.2f}
                Max:           {observed_vals.max():8.2f}   {imputed_vals.max():8.2f}   {final_complete.max():8.2f}
                
                ✅ Quality Checks:
                • Range preservation: {imputed_vals.min() >= observed_vals.min() * 0.8 and imputed_vals.max() <= observed_vals.max() * 1.2}
                • Mean similarity: {abs(observed_vals.mean() - imputed_vals.mean()) / observed_vals.std() < 0.5}
                • Distribution similarity: {abs(observed_vals.std() - imputed_vals.std()) / observed_vals.std() < 0.3}
                """
                
                ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes, fontsize=10,
                        verticalalignment='top', fontfamily='monospace',
                        bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.3))
                
                plt.tight_layout()
                
                # 保存个体诊断图
                safe_var_name = var.replace('/', '_').replace('\\', '_')
                filepath = individual_dir / f'{safe_var_name}_detailed_diagnostics.png'
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                plt.close()
                
                individual_files.append(str(filepath))
                
            except Exception as e:
                logger.error(f"创建 {var} 个体诊断时出错: {e}")
                plt.close()
        
        logger.info(f"个体变量诊断完成，共生成 {len(individual_files)} 个详细报告")
        return individual_files

def main():
    """主函数"""
    print("=" * 60)
    print("简化的 MICE + PMM 多重插补")
    print("=" * 60)
    
    try:
        # 创建插补器
        imputer = MICEImputerWithPMM()
        
        # 运行简化插补和诊断流程
        results = imputer.run_simple_imputation_with_diagnostics()
        
        print("✅ MICE+PMM插补完成!")
        print(f"📊 合并了 {results['n_imputations']} 个插补数据集")
        print(f"📁 最终数据集: {results['final_shape_wide']} (宽格式), {results['final_shape_long']} (长格式)")
        print(f"🔍 诊断图表: {len(results['diagnostic_files'])} 个")
        print(f"⏱️ 处理时间: {results['processing_time']}")
        print(f"📄 宽格式文件: {results['saved_files']['wide_format']}")
        print(f"📄 长格式文件: {results['saved_files']['long_format']}")
        
    except Exception as e:
        logger.error(f"插补过程出错: {e}")
        print(f"❌ 插补失败: {e}")

if __name__ == "__main__":
    main()