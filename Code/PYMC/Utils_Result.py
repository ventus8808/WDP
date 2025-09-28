#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PyMC结果提取模块 (最终版)
提取计算出的各种结果输出到表，按化合物-疾病格式输出
Author: WDP Analysis Team
Date: 2025-09-28
"""

import numpy as np
import pandas as pd
import arviz as az
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional
import warnings
warnings.filterwarnings('ignore')


class BYM2ResultExtractor:
    """BYM2模型结果提取和格式化器"""

    def __init__(self, output_dir: Optional[Path] = None):
        """初始化结果提取器"""
        if output_dir is None:
            project_root = Path(__file__).resolve().parents[2]
            output_dir = project_root / "Result" / "PyMC_Results"

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def extract_exposure_effects(self, trace: az.InferenceData, model_data: Dict) -> Dict:
        """从后验中计算需要的效应量指标"""
        beta_samples = trace.posterior['beta_exposure'].values.flatten()
        exposure_raw = model_data['exposure_raw']
        exposure_log = model_data['exposure_log']

        # 1) RR per SD (log-scale)
        exposure_std = np.std(exposure_log)
        rr_per_sd_samples = np.exp(beta_samples * exposure_std)

        # 2) RR per IQR (raw scale, log transformed diff)
        p25, p75 = np.percentile(exposure_raw, [25, 75])
        log_ratio_iqr = np.log(p75 + 1e-6) - np.log(p25 + 1e-6)
        rr_per_iqr_samples = np.exp(beta_samples * log_ratio_iqr)

        # 3) Quartile contrasts Q2/Q3/Q4 vs Q1
        p50 = np.percentile(exposure_raw, 50)
        q1 = exposure_raw[exposure_raw <= p25]
        q2 = exposure_raw[(exposure_raw > p25) & (exposure_raw <= p50)]
        q3 = exposure_raw[(exposure_raw > p50) & (exposure_raw <= p75)]
        q4 = exposure_raw[exposure_raw > p75]

        # Fallback avoid empty slices
        def safe_mean(x):
            return float(np.mean(x)) if x.size > 0 else float(np.percentile(exposure_raw, 12.5))

        mean_q1 = safe_mean(q1)
        mean_q2 = safe_mean(q2)
        mean_q3 = safe_mean(q3)
        mean_q4 = safe_mean(q4)

        log_ratio_q2_q1 = np.log(mean_q2 + 1e-6) - np.log(mean_q1 + 1e-6)
        log_ratio_q3_q1 = np.log(mean_q3 + 1e-6) - np.log(mean_q1 + 1e-6)
        log_ratio_q4_q1 = np.log(mean_q4 + 1e-6) - np.log(mean_q1 + 1e-6)

        rr_q2_vs_q1_samples = np.exp(beta_samples * log_ratio_q2_q1)
        rr_q3_vs_q1_samples = np.exp(beta_samples * log_ratio_q3_q1)
        rr_q4_vs_q1_samples = np.exp(beta_samples * log_ratio_q4_q1)

        # 4) Bayesian two-sided p-value for beta
        p_value = 2 * min(np.mean(beta_samples > 0), np.mean(beta_samples < 0))
        significance = "***" if p_value < 0.001 else "**" if p_value < 0.01 else "*" if p_value < 0.05 else ""

        return {
            "rr_per_sd": rr_per_sd_samples,
            "rr_per_iqr": rr_per_iqr_samples,
            "rr_q2_vs_q1": rr_q2_vs_q1_samples,
            "rr_q3_vs_q1": rr_q3_vs_q1_samples,
            "rr_q4_vs_q1": rr_q4_vs_q1_samples,
            "p_value": p_value,
            "significance": significance,
        }

    def create_result_row(self, trace: az.InferenceData, model: object, model_data: Dict) -> Dict:
        """创建结果表的一行数据, 包含格式化和诊断。 (v2.0, 支持协变量和交互项)"""
        effects = self.extract_exposure_effects(trace, model_data)

        def format_rr(samples: np.ndarray) -> str:
            mean = float(np.mean(samples))
            ci_l, ci_u = np.percentile(samples, [2.5, 97.5])
            return f"{mean:.3f} ({ci_l:.3f}, {ci_u:.3f})"

        def format_coef(summary_row: pd.Series) -> str:
            # 兼容不同版本的HDI列名
            mean = float(summary_row.get('mean', np.nan))
            hdi_l = float(summary_row.get('hdi_3%', summary_row.get('hdi_2.5%', np.nan)))
            hdi_u = float(summary_row.get('hdi_97%', summary_row.get('hdi_97.5%', np.nan)))
            return f"{mean:.3f} ({hdi_l:.3f}, {hdi_u:.3f})"

        # 诊断 (主暴露变量)
        try:
            summary_exposure = az.summary(trace, var_names=['beta_exposure'])
            r_hat = float(summary_exposure.loc['beta_exposure', 'r_hat'])
            ess_bulk = float(summary_exposure.loc['beta_exposure', 'ess_bulk'])
        except Exception:
            r_hat, ess_bulk = np.nan, np.nan

        # 拟合指标
        try:
            waic = az.waic(trace, scale="deviance")
            waic_val = float(waic.waic)
        except Exception:
            waic_val = np.nan

        # 基础结果行
        result_row = {
            'Timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'Disease': model_data['disease_code'],
            'Exposure': model_data.get('compound_display', model_data['compound']),
            'Category': model_data.get('category', 'Unknown'),
            'Measure': model_data['measure_type'],
            'Estimate': model_data.get('estimate_type', 'avg'),
            'Lag': model_data['lag_years'],
            'Model': model_data['model_type'],

            'RR_per_SD': format_rr(effects['rr_per_sd']),
            'RR_per_IQR': format_rr(effects['rr_per_iqr']),
            'RR_Q1_vs_Q1': "1.0 (Reference)",
            'RR_Q2_vs_Q1': format_rr(effects['rr_q2_vs_q1']),
            'RR_Q3_vs_Q1': format_rr(effects['rr_q3_vs_q1']),
            'RR_Q4_vs_Q1': format_rr(effects['rr_q4_vs_q1']),

            'P_Value': f"{effects['p_value']:.3f}{effects['significance']}",
            'R_hat': f"{r_hat:.4f}" if not np.isnan(r_hat) else 'NA',
            'ESS_bulk': int(ess_bulk) if not np.isnan(ess_bulk) else 'NA',
            'WAIC': f"{waic_val:.2f}" if not np.isnan(waic_val) else 'NA',

            'N_Counties': model_data['n_counties'],
            'N_Records': model_data['n_total_points'],
            'Status_Message': 'SUCCESS'
        }

        # --- 新增核心逻辑：提取协变量和交互项的效应 ---
        covariate_names = model_data.get('covariate_names', [])
        if covariate_names:
            try:
                summary_cov = az.summary(trace, var_names=['beta_covariates'])
            except Exception:
                summary_cov = None
            for i, name in enumerate(covariate_names):
                # 交互项名称规范化: "A * B" -> "A_B"
                norm_name = name.replace(' * ', '_').replace('*', '_')
                if summary_cov is not None and f'beta_covariates[{i}]' in summary_cov.index:
                    summary_row = summary_cov.loc[f'beta_covariates[{i}]']
                else:
                    summary_row = pd.Series({'mean': np.nan, 'hdi_3%': np.nan, 'hdi_97%': np.nan})

                if "*" in name:
                    # 对于交互项，直接报告原始系数
                    col_name = f"Coef_{norm_name}"
                    result_row[col_name] = format_coef(summary_row)
                else:
                    # 对于主效应协变量，报告 RR per SD = exp(beta)
                    col_name = f"RR_per_SD_{name}"
                    try:
                        beta_samples = trace.posterior['beta_covariates'].sel(beta_covariates_dim_0=i).values.flatten()
                    except Exception:
                        # 兼容不同维度名称
                        beta_arr = trace.posterior['beta_covariates'].values
                        beta_samples = beta_arr[..., i].reshape(-1)
                    rr_samples = np.exp(beta_samples)
                    result_row[col_name] = format_rr(rr_samples)

        return result_row

    def save_results(self, result_rows: list, disease_code: str, compound: str) -> Path:
        """将结果保存到CSV文件，使用动态扩展的列顺序。"""
        filename = f"{disease_code}_{compound}_Results.csv"
        output_file = self.output_dir / filename
        results_df = pd.DataFrame(result_rows)

        # --- 修改列顺序逻辑 ---
        # 核心固定列
        core_columns = [
            'Timestamp', 'Disease', 'Exposure', 'Lag', 'Model',
            'RR_per_SD', 'RR_per_IQR',
            'RR_Q1_vs_Q1', 'RR_Q2_vs_Q1', 'RR_Q3_vs_Q1', 'RR_Q4_vs_Q1',
            'P_Value'
        ]
        
        # 动态生成的效应列 (协变量和交互项)
        effect_columns = sorted([col for col in results_df.columns if col.startswith(('RR_per_SD_', 'Coef_'))])

        # 诊断和元数据列
        diag_columns = [
            'R_hat', 'ESS_bulk', 'WAIC', 'N_Counties', 'N_Records', 
            'Status_Message', 'Category', 'Measure', 'Estimate'
        ]

        # 组合成最终顺序
        final_column_order = core_columns + effect_columns + diag_columns
        
        # 确保所有预期的列都存在，不存在的用NA填充
        for col in final_column_order:
            if col not in results_df.columns:
                results_df[col] = 'NA'
        
        # 过滤掉不在最终顺序里的列，并按最终顺序排序
        results_df = results_df[[col for col in final_column_order if col in results_df.columns]]

        if output_file.exists():
            try:
                existing_df = pd.read_csv(output_file)
                
                # --- START: 去重逻辑 (保持不变) ---
                unique_keys = ['Disease', 'Exposure', 'Lag', 'Model', 'Measure', 'Estimate']
                # 确保新旧DF都有这些键
                results_df_keys = results_df[unique_keys].apply(tuple, axis=1)
                existing_keys = set(existing_df[unique_keys].apply(tuple, axis=1))
                
                new_rows_mask = ~results_df_keys.isin(existing_keys)
                results_df_to_append = results_df[new_rows_mask]
                
                if results_df_to_append.empty:
                    print(f"结果已存在，无需追加: {output_file}")
                    return output_file
                # --- END: 去重逻辑 ---

                combined_df = pd.concat([existing_df, results_df_to_append], ignore_index=True)
            except Exception as e:
                print(f"⚠️  追加结果失败，将覆盖文件。错误: {e}")
                combined_df = results_df
            combined_df.to_csv(output_file, index=False)
            print(f"结果已追加到: {output_file}")
        else:
            results_df.to_csv(output_file, index=False)
            print(f"结果已保存到: {output_file}")

        return output_file

    def process_single_analysis(self, trace: az.InferenceData, model: object, model_data: Dict) -> Path:
        """处理单个分析的完整结果提取和保存。"""
        print("\n=== 结果提取 ===")
        result_row = self.create_result_row(trace, model, model_data)

        print("\n=== 分析结果摘要 ===")
        print(f"  RR (per SD):     {result_row['RR_per_SD']}")
        print(f"  RR (Q4 vs Q1):   {result_row['RR_Q4_vs_Q1']}")
        print(f"  诊断 (R-hat):    {result_row['R_hat']}")
        print(f"  诊断 (ESS):      {result_row['ESS_bulk']}")

        output_file = self.save_results([result_row], model_data['disease_code'], model_data['compound'])
        return output_file

    def create_summary_table(self, results_dir: Optional[Path] = None) -> pd.DataFrame:
        """创建所有结果的汇总表"""
        if results_dir is None:
            results_dir = self.output_dir

        result_files = list(results_dir.glob("*_Results.csv"))
        if not result_files:
            print("⚠️  未找到结果文件")
            return pd.DataFrame()

        all_results = []
        for file in result_files:
            try:
                df = pd.read_csv(file)
                all_results.append(df)
            except Exception as e:
                print(f"⚠️  读取文件失败 {file}: {e}")

        if not all_results:
            return pd.DataFrame()

        summary_df = pd.concat(all_results, ignore_index=True)
        summary_file = results_dir / "All_Results_Summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"汇总表保存到: {summary_file}")

        return summary_df


def test_result_extraction():
    """测试结果提取功能"""
    from Utils_Data import WDPDataLoader
    from Utils_Model import BYM2ModelFitter
    
    print("测试结果提取...")
    
    try:
        # 准备测试数据
        loader = WDPDataLoader()
        model_data = loader.prepare_model_data(
            disease_code="C81-C96",
            compound="24D",
            model_type="M0",
            lag_years=5
        )
        
        # 拟合模型
        fitter = BYM2ModelFitter(sampling_config={
            'draws': 100,
            'tune': 50, 
            'chains': 1,
            'cores': 1
        })
        
        model, trace = fitter.run_analysis(model_data)
        
        # 提取结果
        extractor = BYM2ResultExtractor()
        output_file = extractor.process_single_analysis(trace, model, model_data)
        print(f"✅ 结果提取测试成功！输出文件: {output_file}")

    except Exception as e:
        print(f"❌ 结果提取测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # 运行测试
    test_result_extraction()