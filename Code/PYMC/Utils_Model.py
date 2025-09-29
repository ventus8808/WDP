#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PyMC模型拟合模块
BYM2空间时间模型实现，支持M0-M3多种配置
Author: WDP Analysis Team
Date: 2025-09-26
"""
import pymc as pm
import numpy as np
import pytensor.tensor as pt
import arviz as az
from scipy import sparse
from typing import Dict, Optional, Tuple
import os
import warnings
warnings.filterwarnings('ignore')


class BYM2ModelFitter:
    """BYM2贝叶斯空间时间模型拟合器"""
    
    def __init__(self, sampling_config: Optional[Dict] = None):
        """
        初始化模型拟合器
        
        Parameters
        ----------
        sampling_config : dict, optional
            采样配置参数
        """
        # 默认采样配置
        self.sampling_config = sampling_config or {
            'draws': 1000,
            'tune': 500,
            'chains': 2,
            'cores': 2,
            'target_accept': 0.8
        }
        
    def build_bym2_model(self, model_data: Dict) -> pm.Model:
        """
        构建BYM2空间时间模型
        
        Parameters
        ----------
        model_data : dict
            从Utils_Data准备的模型数据字典
            
        Returns
        -------
        pm.Model
            PyMC模型对象
        """
        # 提取数据
        n_counties = model_data['n_counties']
        n_years = model_data['n_years']
        n_covariates = model_data['n_covariates']
        adj_matrix = model_data['adj_matrix']
        
        # 观测数据
        y_obs = model_data['y_obs']
        n_obs = model_data['n_obs']
        exposure_obs = model_data['exposure_obs']
        county_obs = model_data['county_obs']
        time_obs = model_data['time_obs']
        
        # 审查数据  
        y_cens_lower = model_data['y_cens_lower']
        y_cens_upper = model_data['y_cens_upper']
        n_cens = model_data['n_cens']
        exposure_cens = model_data['exposure_cens']
        county_cens = model_data['county_cens']
        time_cens = model_data['time_cens']
        
        # 协变量数据（与 full_data 行顺序对齐的设计矩阵）
        X = model_data['X']
        
        print(f"构建BYM2模型 - {model_data['model_type']}")
        print(f"县数: {n_counties}, 年数: {n_years}, 协变量数: {n_covariates}")
        print(f"观测点: {len(y_obs)}, 审查点: {len(y_cens_lower)}")
        
        # Validate input data for NaNs or invalid values
        assert not np.any(np.isnan(model_data['y_obs'])), "y_obs contains NaN values"
        assert not np.any(np.isnan(model_data['y_cens_lower'])), "y_cens_lower contains NaN values"
        assert not np.any(np.isnan(model_data['y_cens_upper'])), "y_cens_upper contains NaN values"
        assert np.all(model_data['y_obs'] >= 0), "y_obs contains negative values"
        assert np.all(model_data['y_cens_lower'] >= 0), "y_cens_lower contains negative values"
        assert np.all(model_data['y_cens_upper'] >= 0), "y_cens_upper contains negative values"
        
        # Debugging data preparation
        print("Debugging data preparation...")
        print(f"y_obs: {model_data['y_obs'][:10]}")
        print(f"y_cens_lower: {model_data['y_cens_lower'][:10]}")
        print(f"y_cens_upper: {model_data['y_cens_upper'][:10]}")
        print(f"n_obs: {model_data['n_obs']}")
        print(f"n_cens: {model_data['n_cens']}")
        
        # 使用偏移项方案后，不需要对 n_obs / n_cens 做归一化；保留原始规模
        # （之前的归一化会与使用原始人口规模的似然定义不一致，导致数值问题）
        
        with pm.Model() as model:
            
            # === 回归系数 ===
            
            # 截距
            intercept = pm.Normal('intercept', mu=0, sigma=10)
            
            # 暴露效应
            beta_exposure = pm.Normal('beta_exposure', mu=0, sigma=1)
            
            # 协变量效应
            if n_covariates > 0:
                beta_covariates = pm.Normal('beta_covariates', 
                                          mu=0, sigma=1, 
                                          shape=n_covariates)
            
            # === 空间随机效应 (BYM2) ===
            # BYM2: 结构化(CAR) + 非结构化(IID) 的混合，并通过sigma_spatial缩放
            sigma_spatial = pm.HalfNormal('sigma_spatial', sigma=1.0)
            rho = pm.Beta('rho', alpha=1.0, beta=1.0)  # 0: IID, 1: CAR

            # 结构化空间效应 (标准的ICAR)
            adj_no_diag = adj_matrix.copy().tocsr()
            adj_no_diag.setdiag(0)
            adj_no_diag.eliminate_zeros()
            # 使用PyMC的标准ICAR实现，默认alpha=1.0通常更稳定
            u_spatial_raw = pm.CAR('u_spatial_raw', mu=pt.zeros(n_counties), W=adj_no_diag, tau=1.0, alpha=0.999, shape=n_counties)  # alpha 需在(-1,1)
            u_spatial = u_spatial_raw - pt.mean(u_spatial_raw) # 中心化约束仍然是好的实践

            # 非结构化空间效应（标准正态），并中心化
            v_spatial_raw = pm.Normal('v_spatial_raw', mu=0, sigma=1.0, shape=n_counties)
            v_spatial = v_spatial_raw - pt.mean(v_spatial_raw)

            # BYM2 组合场
            phi_spatial = pm.Deterministic(
                'phi_spatial',
                sigma_spatial * (pt.sqrt(rho) * u_spatial + pt.sqrt(1.0 - rho) * v_spatial)
            )
            
            # === 时间随机效应 (RW1) ===
            
            sigma_temporal = pm.Exponential('sigma_temporal', lam=1.0)
            
            # 随机游走 (RW1)
            if n_years > 1:
                # 增量
                delta_temporal = pm.Normal('delta_temporal', mu=0, sigma=sigma_temporal, shape=n_years - 1)
                # 初始点
                phi_temporal_init = pm.Normal('phi_temporal_init', mu=0, sigma=sigma_temporal)
                # 累积求和得到随机游走
                phi_temporal = pm.Deterministic('phi_temporal',
                    pt.concatenate([[phi_temporal_init], phi_temporal_init + pt.cumsum(delta_temporal)])
                )
            else:
                # 如果只有一年
                delta_temporal = pt.as_tensor_variable(np.array([0.0]))
                phi_temporal = pm.Normal('phi_temporal', mu=0, sigma=sigma_temporal, shape=1)
            
            # === 线性预测子与似然函数 ===

            # 准备索引和输入向量
            idx_all = model_data['full_data']['county_idx'].values
            idx_year_all = model_data['full_data']['time_idx'].values
            exposure = model_data.get('exposure_log_stdized', model_data['exposure_log'])
            # 转换为NumPy浮点数组并做清洗，避免 NaN/Inf 传播到似然
            population_full = model_data['full_data']['Population'].to_numpy(dtype='float64')
            # 将缺失/无穷替换为合理的最小值，并裁剪到安全范围内
            population_full = np.nan_to_num(population_full, nan=1.0, posinf=1e12, neginf=1.0)
            population_full = np.clip(population_full, 1.0, 1e12)

            # --- 使用对数偏移项的稳定实现 ---
            # 对数偏移：log(population)，并对0做保护
            log_offset = pt.log(population_full)

            # 基础 log-rate（不含 offset）
            exposure_tensor = pt.as_tensor_variable(exposure)
            log_rate_base = intercept + beta_exposure * exposure_tensor

            # 协变量效应（若存在）：X · beta_covariates
            if n_covariates > 0:
                # 形状与对齐检查
                assert hasattr(X, 'shape'), "X 必须是NumPy数组或可转换为张量的矩阵"
                assert X.ndim == 2, f"X 应为二维矩阵，当前维度: {X.ndim}"
                assert X.shape[1] == n_covariates, (
                    f"X列数({X.shape[1]})应与协变量数({n_covariates})一致"
                )
                assert X.shape[0] == len(idx_all), (
                    f"X行数({X.shape[0]})应与full_data长度({len(idx_all)})一致"
                )
                # 确保 X 的行数与 full_data 对齐
                # 这里假设 prepare_model_data 已构造与 full_data 同长度的 X
                X_tensor = pt.as_tensor_variable(X)
                log_rate_base = log_rate_base + pt.dot(X_tensor, beta_covariates)

            # 加入时空随机效应
            log_rate_base = log_rate_base + phi_spatial[idx_all] + phi_temporal[idx_year_all]

            # 完整的 log(lambda) = log_offset + log_rate_base
            log_lambda_full = log_offset + log_rate_base

            # 约束以防数值溢出（在log尺度更稳定）
            stable_log_lambda_full = pt.clip(log_lambda_full, -25, 25)

            # 计算 lambda
            lambda_full = pt.exp(stable_log_lambda_full)

            # 分离观测和审查数据的索引
            is_censored = model_data['full_data']['is_censored'].values
            idx_obs = np.where(~is_censored)[0]
            idx_cens = np.where(is_censored)[0]

            # 选择对应的 lambda
            lambda_obs = lambda_full[idx_obs]
            lambda_cens = lambda_full[idx_cens]
            
            # 观测数据似然（使用长度判断，避免布尔歧义）
            if len(y_obs) > 0:
                y_obs_tensor = pt.as_tensor_variable(y_obs)
                pm.Potential('y_obs_likelihood', 
                    pt.sum(pm.logp(pm.Poisson.dist(mu=lambda_obs), y_obs_tensor))
                )

            # 审查数据似然：使用对区间内精确质量的 logsumexp 聚合，避免 CDF 差导致的灾难性抵消
            if len(y_cens_lower) > 0:
                # 全局K范围，并按行掩码保证仅累加各自区间内的质量
                k_min = int(np.min(y_cens_lower))
                k_max = int(np.max(y_cens_upper))
                # 保护：若k范围无效（极端情况），直接跳过
                if k_max >= k_min:
                    k = pt.arange(k_min, k_max + 1)  # 形如 [k_min, ..., k_max]

                    # 计算每个审查点在每个k上的对数概率（形状：n_cens x K）
                    k_row = k[None, :]
                    lambda_cens_mat = lambda_cens[:, None]
                    logp_matrix = pm.logp(pm.Poisson.dist(mu=lambda_cens_mat), k_row)

                    # 构造掩码，仅保留 [lower_i, upper_i] 内的项
                    lower_expanded = pt.as_tensor_variable(y_cens_lower)[:, None]
                    upper_expanded = pt.as_tensor_variable(y_cens_upper)[:, None]
                    in_range = (k_row >= lower_expanded) & (k_row <= upper_expanded)

                    # 将区间外的logp置为 -inf，确保对数求和不计入
                    logp_masked = pt.switch(in_range, logp_matrix, -np.inf)

                    # 对每一行进行logsumexp，得到 log P(lower_i <= Y <= upper_i)
                    log_interval_prob = pm.math.logsumexp(logp_masked, axis=1)

                    pm.Potential('y_cens_likelihood', pt.sum(log_interval_prob))
        
        # 存储模型数据用于后续分析
        model.model_data = model_data
        
        return model
    
    def fit_model(self, model: pm.Model) -> az.InferenceData:
        """
        拟合BYM2模型
        
        Parameters
        ----------
        model : pm.Model
            构建好的PyMC模型
            
        Returns
        -------
        az.InferenceData
            MCMC采样结果
        """
        import sys
        from datetime import datetime
        
        print(f"\n[{datetime.now().strftime('%H:%M:%S')}] ===== 开始MCMC采样 =====")
        print(f"[{datetime.now().strftime('%H:%M:%S')}] 采样配置: {self.sampling_config}")
        sys.stdout.flush()
        
        # 根据运行环境自动选择初始化策略与并行度
        try:
            cpu_env = int(os.environ.get('SLURM_CPUS_PER_TASK') or 0)
        except Exception:
            cpu_env = 0
        cpu_count = cpu_env if cpu_env > 0 else (os.cpu_count() or 1)
        # 初始化策略：CPU较多时使用 jitter+adapt_diag，提高初始点探索；CPU较少时使用 adapt_diag 更稳健
        init_strategy = 'jitter+adapt_diag' if cpu_count >= 8 else 'adapt_diag'
        # 限制并行度不过度超配
        req_chains = int(self.sampling_config.get('chains', 2))
        req_cores = int(self.sampling_config.get('cores', 2))
        sample_chains = max(1, min(req_chains, cpu_count))
        sample_cores = max(1, min(req_cores, cpu_count))
        
        print(f"[{datetime.now().strftime('%H:%M:%S')}] CPU检测: {cpu_count} 可用")
        print(f"[{datetime.now().strftime('%H:%M:%S')}] 采样策略: {init_strategy}")
        print(f"[{datetime.now().strftime('%H:%M:%S')}] 并行配置: {sample_chains} chains, {sample_cores} cores")
        print(f"[{datetime.now().strftime('%H:%M:%S')}] 开始初始化采样器...")
        sys.stdout.flush()
        
        with model:
            try:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] 启动NUTS采样器...")
                print(f"[{datetime.now().strftime('%H:%M:%S')}] 预计采样时间: ~{(self.sampling_config['draws'] + self.sampling_config['tune']) * sample_chains // 100:.1f} 分钟")
                sys.stdout.flush()
                
                # 运行MCMC采样
                trace = pm.sample(
                    draws=self.sampling_config['draws'],
                    tune=self.sampling_config['tune'],
                    chains=sample_chains,
                    cores=sample_cores,
                    target_accept=self.sampling_config['target_accept'],
                    init=init_strategy,
                    return_inferencedata=True,
                    progressbar=True,  # 确保显示进度条
                    random_seed=42  # 添加随机种子以便复现
                )
                
                print(f"[{datetime.now().strftime('%H:%M:%S')}] ✅ MCMC采样完成!")
                sys.stdout.flush()
                
            except Exception as e:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] ❌ MCMC采样失败: {e}")
                print(f"[{datetime.now().strftime('%H:%M:%S')}] 错误类型: {type(e).__name__}")
                if 'divergence' in str(e).lower():
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] 💡 提示: 发散问题，尝试降低target_accept或增加tune步数")
                elif 'memory' in str(e).lower():
                    print(f"[{datetime.now().strftime('%H:%M:%S')}] 💡 提示: 内存不足，尝试减少chains数量")
                sys.stdout.flush()
                raise
            
            # 添加后验预测检验
            try:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] 开始后验预测采样...")
                sys.stdout.flush()
                pm.sample_posterior_predictive(trace, extend_inferencedata=True)
                print(f"[{datetime.now().strftime('%H:%M:%S')}] ✅ 后验预测完成!")
                sys.stdout.flush()
            except Exception as e:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] ⚠️  后验预测采样失败: {e}")
                print(f"[{datetime.now().strftime('%H:%M:%S')}] 继续进行收敛诊断...")
                sys.stdout.flush()
        
        # 收敛诊断
        self._check_convergence(trace)
        
        return trace
    
    def _check_convergence(self, trace: az.InferenceData) -> None:
        """
        检查MCMC收敛性
        
        Parameters
        ----------
        trace : az.InferenceData
            MCMC采样结果
        """
        from datetime import datetime
        import sys
        
        print(f"\n[{datetime.now().strftime('%H:%M:%S')}] === 收敛诊断 ===")
        sys.stdout.flush()
        
        # 仅对核心参数计算 R-hat，避免在大型潜变量上造成内存压力
        try:
            core_vars = [
                'intercept', 'beta_exposure', 'beta_covariates',
                'sigma_spatial', 'rho', 'sigma_temporal'
            ]
            rhat_ds = az.rhat(trace, var_names=[v for v in core_vars if v in trace.posterior])
            # 将所有变量堆叠为单向量后取最大值
            if hasattr(rhat_ds, 'to_array'):
                max_rhat_val = rhat_ds.to_array().max().values
            else:
                max_rhat_val = rhat_ds.max().values if hasattr(rhat_ds, 'max') else float(rhat_ds)
            max_rhat = float(np.asarray(max_rhat_val))
            print(f"最大R̂: {max_rhat:.4f}")
            if np.isfinite(max_rhat) and max_rhat > 1.01:
                print("⚠️  警告: R̂ > 1.01，模型可能未收敛")
            elif np.isfinite(max_rhat):
                print("✅ R̂ < 1.01，收敛良好")
            else:
                print("ℹ️ R̂ 无法可靠估计（样本过少或变量不可用）")
        except Exception as e:
            print(f"ℹ️ R̂ 计算失败或不可靠: {e}")
        
        # 有效样本量（bulk/tail），对Dataset取最小值
        try:
            core_vars = [
                'intercept', 'beta_exposure', 'beta_covariates',
                'sigma_spatial', 'rho', 'sigma_temporal'
            ]
            ess_bulk_ds = az.ess(trace, method="bulk", var_names=[v for v in core_vars if v in trace.posterior])
            ess_tail_ds = az.ess(trace, method="tail", var_names=[v for v in core_vars if v in trace.posterior])
            if hasattr(ess_bulk_ds, 'to_array'):
                min_ess_bulk_val = ess_bulk_ds.to_array().min().values
            else:
                min_ess_bulk_val = ess_bulk_ds.min().values if hasattr(ess_bulk_ds, 'min') else ess_bulk_ds
            if hasattr(ess_tail_ds, 'to_array'):
                min_ess_tail_val = ess_tail_ds.to_array().min().values
            else:
                min_ess_tail_val = ess_tail_ds.min().values if hasattr(ess_tail_ds, 'min') else ess_tail_ds
            min_ess_bulk = float(np.asarray(min_ess_bulk_val))
            min_ess_tail = float(np.asarray(min_ess_tail_val))
            print(f"最小ESS (bulk): {min_ess_bulk:.0f}")
            print(f"最小ESS (tail): {min_ess_tail:.0f}")
            if np.isfinite(min_ess_bulk) and min_ess_bulk < 400:
                print("⚠️  警告: ESS (bulk) < 400，需要更多样本")
            if np.isfinite(min_ess_tail) and min_ess_tail < 400:
                print("⚠️  警告: ESS (tail) < 400，需要更多样本")
            if (np.isfinite(min_ess_bulk) and np.isfinite(min_ess_tail) and 
                min_ess_bulk >= 400 and min_ess_tail >= 400):
                print("✅ ESS充足")
        except Exception as e:
            print(f"ℹ️ ESS 计算失败或不可靠: {e}")
    
    def run_analysis(self, model_data: Dict) -> Tuple[pm.Model, az.InferenceData]:
        """
        运行完整的BYM2分析
        
        Parameters
        ---------- 
        model_data : dict
            模型数据字典
            
        Returns
        -------
        Tuple[pm.Model, az.InferenceData]
            模型对象和采样结果
        """
        print(f"\n=== BYM2模型分析 ===")
        print(f"疾病: {model_data['disease_code']}")
        print(f"化合物: {model_data['compound']}")
        print(f"模型: {model_data['model_type']}")
        
        # 构建模型
        model = self.build_bym2_model(model_data)
        
        # 拟合模型
        trace = self.fit_model(model)
        
        return model, trace


def quick_test_model():
    """快速测试模型构建"""
    from Utils_Data import WDPDataLoader
    
    print("快速测试BYM2模型...")
    
    try:
        # 加载数据
        loader = WDPDataLoader()
        model_data = loader.prepare_model_data(
            disease_code="C81-C96",
            compound="24D", 
            model_type="M0",  # 简单模型测试
            lag_years=5
        )
        
        # 测试模型构建
        fitter = BYM2ModelFitter(sampling_config={
            'draws': 100,
            'tune': 50,
            'chains': 1,
            'cores': 1,
            'target_accept': 0.8
        })
        
        model, trace = fitter.run_analysis(model_data)
        
        print("✅ 模型测试成功！")
        return model, trace
        
    except Exception as e:
        print(f"❌ 模型测试失败: {e}")
        import traceback
        traceback.print_exc()
        return None, None


if __name__ == "__main__":
    # 运行快速测试
    model, trace = quick_test_model()