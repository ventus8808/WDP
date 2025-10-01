#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM总控模块 - 统一的分析流程控制

主要功能:
1. 协调数据处理、模型拟合、结果输出三大模块
2. 支持标准分析和滞后分析两种模式
3. 批量处理多种癌症类型和EQI变量
4. 分层分析 (整体 + 城市化/RUCC分层)
5. 生成标准化结果和报告

分析流程:
数据加载 → 模型拟合 → 结果提取 → 文件输出

支持的分析类型:
- 标准分析: EQI vs 癌症死亡率
- 滞后分析: 多EQI变量 vs 癌症死亡率 (支持4种滞后场景)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import logging
import time
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
import traceback
import argparse

# 添加当前目录到Python路径
current_dir = Path(__file__).resolve().parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

# 导入自定义模块
from LMM_Data import LMMDataProcessor
from LMM_Model import LMMAnalyzer
from LMM_Result import LMMResultProcessor

# 设置日志
def setup_logging(log_level="INFO", log_file=None):
    """设置日志配置"""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    if log_file:
        logging.basicConfig(
            level=getattr(logging, log_level),
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(level=getattr(logging, log_level), format=log_format)

logger = logging.getLogger(__name__)

class LMMPipeline:
    """LMM分析主流程管理器 - 统一的分析流程控制"""
    
    def __init__(self, eqi_period: str = "0610", 
                 aamr_period: str = "2016_2020",
                 enable_stratification: bool = True,
                 run_all_combinations: bool = False):
        """
        初始化LMM分析流程 - 横断面分析
        
        参数:
            eqi_period: EQI时间跨度 ("0005"=2000-2005合并, "0610"=2006-2010合并)
            aamr_period: AAMR时间跨度 ("2006_2010", "2011_2015", "2016_2020"合并)
            enable_stratification: 是否进行分层分析
            run_all_combinations: 是否运行所有4种年份组合
            
        注: 所有数据都是时间跨度合并的横断面数据，不是从多年份中选择
        """
        self.eqi_period = eqi_period
        self.aamr_period = aamr_period
        self.enable_stratification = enable_stratification
        self.run_all_combinations = run_all_combinations
        
        # 设置项目根目录
        self.project_root = Path(__file__).resolve().parents[2]
        
        # 初始化组件 (横断面分析)
        self.data_processor = LMMDataProcessor(
            eqi_period=eqi_period,
            aamr_period=aamr_period
        )
        self.model_analyzer = LMMAnalyzer()
        self.result_processor = LMMResultProcessor(
            project_root=self.project_root
        )
        
        # 定义所有4种年份组合
        self.all_period_combinations = [
            ("0005", "2006_2010"),
            ("0005", "2016_2020"), 
            ("0610", "2006_2010"),
            ("0610", "2016_2020")
        ]
        
        # 分析配置
        self.cancer_types_config = {
            'primary': ['C00_C97', 'C34', 'C50', 'C81_C96'],
            'secondary': ['C15_C26', 'C18_C21', 'C25', 'C30_C39', 'C61', 'C64_C68'],
            'all': ['C00_C97', 'C15_C26', 'C18_C21', 'C25', 'C30_C39', 
                   'C34', 'C50', 'C51_C58', 'C60_C63', 'C61', 'C64_C68', 
                   'C76_C80', 'C81_C96']
        }
        
        # 分层配置 (横断面分析使用城市化类型分层)
        self.stratification_config = {
            'variable': 'Urbanization_Type', 
            'strata': ['Large Metro', 'Medium Metro', 'Small Metro', 'Noncore']
        }
        
    def run_complete_analysis(self, cancer_subset: str = "all") -> Dict[str, Any]:
        """
        执行完整的LMM分析流程
        
        参数:
            cancer_subset: 癌症类型子集 ("primary", "secondary", "all")
            
        返回:
            包含所有分析结果的字典
        """
        start_time = time.time()
        print(f"开始LMM横断面分析流程")
        print(f"EQI时间跨度: {self.eqi_period} (合并数据)")
        print(f"AAMR时间跨度: {self.aamr_period} (合并数据)")
        print(f"癌症子集: {cancer_subset}")
        print("=" * 60)
        
        try:
            # 步骤1: 数据准备
            print("\n步骤1: 数据加载和预处理...")
            
            # 先执行完整的数据处理流程
            if not self.data_processor.process_all():
                print("数据处理失败")
                return {}
            
            analysis_data = self.data_processor.get_analysis_ready_data()
            if analysis_data is None:
                print("获取分析数据失败")
                return {}
                
            print(f"数据加载完成: {analysis_data.shape[0]} 条记录, {analysis_data.shape[1]} 个特征")
            
            # 步骤2: 主要分析
            print("\n步骤2: 执行主要分析...")
            cancer_types = self.cancer_types_config[cancer_subset]
            main_results = self._run_main_analysis(analysis_data, cancer_types)
            
            # 步骤3: 分层分析 (可选)
            if self.enable_stratification:
                print("\n步骤3: 执行分层分析...")
                stratified_results = self._run_stratified_analysis(analysis_data, cancer_types)
                main_results['stratified'] = stratified_results
            
            # 步骤4: 结果保存和输出
            print("\n步骤4: 保存分析结果...")
            self._save_all_results(main_results)
            
            # 完成分析
            elapsed_time = time.time() - start_time
            print(f"\n分析完成! 总耗时: {elapsed_time:.2f} 秒")
            print("=" * 60)
            
            return main_results
            
        except Exception as e:
            print(f"分析过程中出现错误: {str(e)}")
            traceback.print_exc()
            return {}
    
    def run_all_period_combinations(self, cancer_subset: str = "all") -> Dict[str, Any]:
        """
        运行所有4种年份组合的分析
        
        参数:
            cancer_subset: 癌症类型子集
            
        返回:
            包含所有年份组合结果的字典
        """
        start_time = time.time()
        print(f"开始运行所有4种年份组合分析")
        print(f"癌症子集: {cancer_subset}")
        print("=" * 60)
        
        all_combination_results = {}
        
        for i, (eqi_period, aamr_period) in enumerate(self.all_period_combinations, 1):
            print(f"\n[{i}/4] 运行年份组合: EQI {eqi_period} vs AAMR {aamr_period}")
            print("-" * 50)
            
            try:
                # 重新初始化数据处理器用于当前组合
                self.data_processor = LMMDataProcessor(
                    eqi_period=eqi_period,
                    aamr_period=aamr_period
                )
                
                # 处理数据
                if not self.data_processor.process_all():
                    print(f"数据处理失败: EQI {eqi_period} vs AAMR {aamr_period}")
                    continue
                
                analysis_data = self.data_processor.get_analysis_ready_data()
                if analysis_data is None:
                    print(f"获取分析数据失败: EQI {eqi_period} vs AAMR {aamr_period}")
                    continue
                
                print(f"数据加载完成: {analysis_data.shape[0]} 条记录, {analysis_data.shape[1]} 个特征")
                
                # 运行分析
                cancer_types = self.cancer_types_config[cancer_subset]
                results = self._run_main_analysis(analysis_data, cancer_types)
                
                # 存储结果
                result_key = f"EQI{eqi_period}_AAMR{aamr_period}"
                all_combination_results[result_key] = results.get('overall', {})
                
                print(f"✅ 完成年份组合: EQI {eqi_period} vs AAMR {aamr_period}")
                
            except Exception as e:
                print(f"❌ 年份组合失败: EQI {eqi_period} vs AAMR {aamr_period} - {str(e)}")
                continue
        
        # 生成简洁表格
        print(f"\n生成简洁结果表格...")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        saved_files = self.result_processor.save_all_compact_tables(all_combination_results, timestamp)
        
        elapsed_time = time.time() - start_time
        print(f"\n所有年份组合分析完成!")
        print(f"总耗时: {elapsed_time:.2f} 秒")
        print(f"生成了 {len(saved_files)} 个简洁表格文件")
        print("=" * 60)
        
        return all_combination_results
    
    def _run_main_analysis(self, data: pd.DataFrame, cancer_types: List[str]) -> Dict[str, Any]:
        """执行主要分析（整体，不分层）"""
        results = {'overall': {}}
        
        print(f"    正在分析 {len(cancer_types)} 种癌症类型...")
        
        # 使用统一的模型分析器
        overall_results = self.model_analyzer.fit_multiple_models(
            data=data, 
            cancer_types=cancer_types
        )
        
        if overall_results:
            results['overall'] = overall_results
            print(f"    整体分析完成: {len(overall_results)} 个模型")
        else:
            print("    整体分析失败")
        
        return results
    
    def _run_stratified_analysis(self, data: pd.DataFrame, cancer_types: List[str]) -> Dict[str, Any]:
        """执行分层分析"""
        print(f"    按 {self.stratification_config['variable']} 进行分层分析...")
        
        stratified_results = {}
        
        for stratum in self.stratification_config['strata']:
            print(f"      分析分层: {stratum}")
            
            # 筛选该分层的数据
            stratum_data = data[data[self.stratification_config['variable']] == stratum].copy()
            
            if len(stratum_data) > 0:
                stratum_results = self.model_analyzer.fit_multiple_models(
                    data=stratum_data,
                    cancer_types=cancer_types
                )
                if stratum_results:
                    stratified_results[stratum] = stratum_results
                    print(f"        完成: {len(stratum_results)} 个模型")
                else:
                    print(f"        失败: {stratum}")
            else:
                print(f"        跳过: {stratum} (无数据)")
        
        return stratified_results
    
    def _save_all_results(self, results: Dict[str, Any]) -> None:
        """保存所有分析结果"""
        try:
            # 生成时间戳
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # 保存详细结果
            if 'overall' in results:
                self.result_processor.save_results(
                    results['overall'], 
                    analysis_type="Overall",
                    timestamp=timestamp
                )
            
            # 保存分层结果
            if 'stratified' in results:
                for stratum, stratum_results in results['stratified'].items():
                    self.result_processor.save_results(
                        stratum_results,
                        analysis_type=f"Stratified_{stratum}",
                        timestamp=timestamp
                    )
            
            print(f"    所有结果已保存到: Result/EQI_GLMM/ (时间戳: {timestamp})")
            
        except Exception as e:
            print(f"    保存结果时出错: {e}")
            traceback.print_exc()
    

def main():
    """主函数 - 执行LMM横断面分析"""
    parser = argparse.ArgumentParser(description='LMM横断面分析流程 - 两种核心模型')
    parser.add_argument('--eqi-period', default='0610',
                       choices=['0005', '0610'],
                       help='EQI时间跨度: 0005(2000-2005合并), 0610(2006-2010合并)')
    parser.add_argument('--aamr-period', default='2016_2020',
                       choices=['2006_2010', '2011_2015', '2016_2020'],
                       help='AAMR时间跨度: 合并的癌症死亡率数据')
    parser.add_argument('--cancer-subset', default='all',
                       choices=['primary', 'secondary', 'all'],
                       help='癌症类型子集')
    parser.add_argument('--no-stratification', action='store_true',
                       help='禁用分层分析')
    parser.add_argument('--all-combinations', action='store_true',
                       help='运行所有4种年份组合分析并生成简洁表格')
    
    args = parser.parse_args()
    
    print("LMM分析说明:")
    print("1. 主效应模型: 总体EQI vs 癌症死亡率")
    print("2. 领域探索模型: 五大环境领域 vs 癌症死亡率")
    print("注: 所有数据都是时间跨度合并的横断面数据\n")
    
    # 创建分析流程
    pipeline = LMMPipeline(
        eqi_period=args.eqi_period,
        aamr_period=args.aamr_period,
        enable_stratification=not args.no_stratification,
        run_all_combinations=args.all_combinations
    )
    
    # 执行分析
    if args.all_combinations:
        print("运行模式: 所有4种年份组合分析")
        results = pipeline.run_all_period_combinations(cancer_subset=args.cancer_subset)
    else:
        print(f"运行模式: 单一年份组合 (EQI {args.eqi_period} vs AAMR {args.aamr_period})")
        results = pipeline.run_complete_analysis(cancer_subset=args.cancer_subset)
    
    if results:
        print("🎉 分析成功完成!")
        return True
    else:
        print("❌ 分析失败!")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)



    
    def _run_stratified_analysis(self, cancer_types: List[str]) -> bool:
        """运行分层分析"""
        logger.info("步骤3: 分层分析")
        
        stratify_var = self.stratification_config['variable']
        target_strata = self.stratification_config['strata']
        
        total_success = 0
        
        for cancer_type in cancer_types:
            logger.info(f"  分层分析癌症类型: {cancer_type}")
            
            try:
                # 获取该癌症类型的完整数据
                full_data = self.data_processor.get_analysis_ready_data(
                    cancer_type=cancer_type
                )
                
                if full_data is None or len(full_data) == 0:
                    logger.warning(f"    癌症类型 {cancer_type} 无数据用于分层分析")
                    continue
                
                # 检查可用分层
                available_strata = full_data[stratify_var].dropna().unique()
                logger.info(f"    可用分层: {list(available_strata)}")
                
                cancer_success = 0
                
                for stratum in target_strata:
                    if stratum not in available_strata:
                        logger.warning(f"      分层 {stratum} 不存在，跳过")
                        continue
                    
                    logger.info(f"    处理分层: {stratum}")
                    
                    # 获取分层数据
                    stratum_data = self.data_processor.get_analysis_ready_data(
                        cancer_type=cancer_type,
                        stratum=stratum
                    )
                    
                    if stratum_data is None or len(stratum_data) < 100:
                        logger.warning(f"      分层 {stratum} 数据不足 ({len(stratum_data) if stratum_data is not None else 0})，跳过")
                        continue
                    
                    # 检查州数量
                    state_count = stratum_data['State_FIPS'].nunique()
                    if state_count < 3:
                        logger.warning(f"      分层 {stratum} 州数量不足 ({state_count})，跳过")
                        continue
                    
                    try:
                        # 拟合模型
                        result = self.model_analyzer.fit_lmm_model(
                            stratum_data,
                            model_name=f"stratified_{cancer_type}_{stratum}"
                        )
                        
                        if result is None:
                            logger.warning(f"      分层 {stratum} 模型拟合失败")
                            continue
                        
                        # 提取系数
                        coefficients = self.model_analyzer.extract_eqi_coefficients(
                            result, f"stratified_{cancer_type}_{stratum}"
                        )
                        
                        # 运行诊断
                        diagnostics = self.model_analyzer.run_model_diagnostics(
                            result, stratum_data, f"stratified_{cancer_type}_{stratum}"
                        )
                        
                        # 存储结果
                        self.result_processor.add_analysis_result(
                            cancer_type=cancer_type,
                            stratum=stratum,
                            coefficients=coefficients,
                            diagnostics=diagnostics,
                            model_info={
                                'type': 'stratified',
                                'stratify_var': stratify_var,
                                'n_obs': len(stratum_data),
                                'n_states': state_count
                            }
                        )
                        
                        cancer_success += 1
                        total_success += 1
                        logger.info(f"      ✅ 分层 {stratum} 完成")
                        
                    except Exception as e:
                        logger.error(f"      ❌ 分层 {stratum} 失败: {e}")
                        continue
                
                logger.info(f"    癌症 {cancer_type} 分层分析: {cancer_success} 个分层成功")
                
            except Exception as e:
                logger.error(f"    癌症类型 {cancer_type} 分层分析失败: {e}")
                continue
        
        logger.info(f"分层分析完成: 总计 {total_success} 个成功分析")
        return total_success > 0
    
    def _run_sensitivity_analysis(self, cancer_types: List[str]) -> bool:
        """运行敏感性分析"""
        logger.info("步骤4: 敏感性分析")
        
        # 选择一个代表性癌症类型进行敏感性分析
        representative_cancer = cancer_types[0] if cancer_types else 'C00_C97'
        logger.info(f"  使用 {representative_cancer} 进行敏感性分析")
        
        try:
            # 获取分析数据
            analysis_data = self.data_processor.get_analysis_ready_data(
                cancer_type=representative_cancer
            )
            
            if analysis_data is None or len(analysis_data) == 0:
                logger.warning("敏感性分析无可用数据")
                return False
            
            # 运行敏感性模型
            sensitivity_results = self.model_analyzer.fit_sensitivity_models(analysis_data)
            
            if not sensitivity_results:
                logger.warning("敏感性分析未产生有效结果")
                return False
            
            # 创建模型比较表
            comparison_df = self.model_analyzer.compare_models(sensitivity_results)
            
            if not comparison_df.empty:
                # 保存敏感性分析结果
                sensitivity_path = self.output_dir / f"sensitivity_analysis_{representative_cancer}_{self.time_period}.csv"
                comparison_df.to_csv(sensitivity_path, index=False)
                logger.info(f"敏感性分析结果已保存: {sensitivity_path}")
            
            logger.info("敏感性分析完成")
            return True
            
        except Exception as e:
            logger.error(f"敏感性分析失败: {e}")
            return False
    
    def _generate_final_results(self) -> bool:
        """生成最终结果"""
        logger.info("步骤5: 生成最终结果")
        
        try:
            # 创建宽格式摘要表
            summary_table = self.result_processor.create_wide_format_table()
            
            if summary_table.empty:
                logger.error("无法生成摘要表 - 没有有效结果")
                return False
            
            # 生成时间戳
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # 保存结果文件
            saved_files = self.result_processor.save_results(timestamp)
            
            # 生成可视化
            try:
                # 森林图 (Overall分层)
                forest_plot_path = self.result_processor.create_forest_plot(
                    stratum='Overall'
                )
                if forest_plot_path:
                    logger.info(f"森林图已生成: {forest_plot_path}")
                
                # 系数热图 (EQI Q5 vs Q1)
                heatmap_path = self.result_processor.create_coefficient_heatmap(
                    eqi_term='EQI_Q5_vs_Q1'
                )
                if heatmap_path:
                    logger.info(f"系数热图已生成: {heatmap_path}")
                
            except Exception as e:
                logger.warning(f"可视化生成失败: {e}")
            
            # 生成摘要报告
            summary_report = self.result_processor.generate_summary_report()
            
            # 保存摘要报告
            report_path = self.output_dir / f"LMM_Analysis_Report_{timestamp}.txt"
            with report_path.open('w', encoding='utf-8') as f:
                f.write(summary_report)
            
            # 输出关键信息
            logger.info("=" * 60)
            logger.info("分析完成! 主要输出文件:")
            for file_type, file_path in saved_files.items():
                logger.info(f"  {file_type}: {file_path}")
            logger.info(f"  分析报告: {report_path}")
            logger.info("=" * 60)
            
            # 输出摘要报告到控制台
            print("\n" + summary_report)
            
            return True
            
        except Exception as e:
            logger.error(f"生成最终结果失败: {e}")
            return False
    
    def get_results_summary(self) -> Dict:
        """获取结果摘要"""
        return {
            'timestamp': datetime.now().isoformat(),
            'time_period': self.time_period,
            'output_directory': str(self.output_dir),
            'results_data': self.result_processor.results_data
        }


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='LMM死亡率差分析')
    parser.add_argument('--time-period', default='2016_2020',
                       help='分析时间周期 (默认: 2016_2020)')
    parser.add_argument('--cancer-subset', default='primary',
                       choices=['primary', 'secondary', 'all'],
                       help='癌症类型子集 (默认: primary)')
    parser.add_argument('--no-stratification', action='store_true',
                       help='禁用分层分析')
    parser.add_argument('--no-sensitivity', action='store_true',
                       help='禁用敏感性分析')
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='日志级别 (默认: INFO)')
    parser.add_argument('--output-dir', type=str,
                       help='输出目录路径')
    
    args = parser.parse_args()
    
    # 设置日志
    setup_logging(log_level=args.log_level)
    
    # 创建分析管道
    pipeline = LMMPipeline(
        time_period=args.time_period,
        output_dir=Path(args.output_dir) if args.output_dir else None,
        enable_stratification=not args.no_stratification,
        enable_sensitivity=not args.no_sensitivity
    )
    
    # 运行分析
    success = pipeline.run_full_analysis(cancer_subset=args.cancer_subset)
    
    if success:
        logger.info("🎉 分析成功完成!")
        sys.exit(0)
    else:
        logger.error("❌ 分析失败!")
        sys.exit(1)


if __name__ == "__main__":
    main()