#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WDP PyMC主脚本
主程序，操控整个流程，支持传入参数：疾病、模型、滞后等等
Author: WDP Analysis Team
Date: 2025-09-26
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict, Optional
import traceback
import os

# 添加当前目录到Python路径
sys.path.insert(0, str(Path(__file__).parent))

# 导入自定义模块
from Utils_Data import WDPDataLoader
from Utils_Model import BYM2ModelFitter
from Utils_Result import BYM2ResultExtractor
from Utils_Others import (
    get_pymc_config, get_sampling_config, validate_disease_code, 
    validate_model_type, parse_compound_list, parse_model_list,
    parse_lag_years, check_data_availability, print_analysis_summary
)


class WDPPyMCAnalysis:
    """WDP PyMC分析主控制器"""
    
    def __init__(self, config_path: Optional[Path] = None, 
                 output_dir: Optional[Path] = None,
                 sampling_mode: str = 'test',
                 cores: Optional[str | int] = None,
                 chains: Optional[str | int] = None,
                 draws: Optional[int] = None,
                 tune: Optional[int] = None,
                 target_accept: Optional[float] = None):
        """
        初始化分析控制器
        
        Parameters
        ----------
        config_path : Path, optional
            配置文件路径
        output_dir : Path, optional
            输出目录
        sampling_mode : str
            采样模式 ('test' 或 'production')
        """
        self.config_path = config_path
        self.output_dir = output_dir
        self.sampling_mode = sampling_mode
        
        # 初始化组件
        self.data_loader = WDPDataLoader(config_path)
        
        # 获取采样配置
        sampling_config = get_sampling_config(sampling_mode, config_path)

        # Auto-detect cores and apply overrides
        cpu_count = os.cpu_count() or 1

        def _to_int_or_auto(v):
            if v is None:
                return None
            if isinstance(v, int):
                return v
            v_str = str(v).strip().lower()
            if v_str == 'auto':
                return 'auto'
            if v_str.isdigit():
                return int(v_str)
            return None

        cores = _to_int_or_auto(cores)
        chains = _to_int_or_auto(chains)

        # Default to max cores; chains align to cores unless overridden
        if cores in (None, 'auto'):
            sampling_config['cores'] = cpu_count
        else:
            sampling_config['cores'] = max(1, int(cores))

        if chains in (None, 'auto'):
            sampling_config['chains'] = sampling_config['cores']
        else:
            sampling_config['chains'] = max(1, int(chains))

        if draws is not None:
            sampling_config['draws'] = int(draws)
        if tune is not None:
            sampling_config['tune'] = int(tune)
        if target_accept is not None:
            sampling_config['target_accept'] = float(target_accept)

        self.model_fitter = BYM2ModelFitter(sampling_config)
        
        self.result_extractor = BYM2ResultExtractor(output_dir)
        
        # 初始化信息
        print(f"WDP PyMC分析系统初始化完成")
        print(f"采样模式: {self.sampling_mode}")
        print(f"输出目录: {self.result_extractor.output_dir}")
        print(
            "采样配置: "
            f"draws={self.model_fitter.sampling_config['draws']}, "
            f"tune={self.model_fitter.sampling_config['tune']}, "
            f"chains={self.model_fitter.sampling_config['chains']}, "
            f"cores={self.model_fitter.sampling_config['cores']}, "
            f"target_accept={self.model_fitter.sampling_config['target_accept']}"
        )
    
    def run_single_analysis(self, disease_code: str, compound: str, 
                          model_type: str, lag_years: int,
                          measure_type: str = 'Weight',
                          estimate_type: str = 'avg') -> Optional[Path]: # <--- 修改
        """
        运行单个分析
        
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
            测量类型
            
        Returns
        -------
        Optional[Path]
            输出文件路径，失败时返回None
        """
        try:
            print(f"\n{'='*80}")
            print(f"开始分析: {disease_code} - {compound} - {model_type}")
            print(f"{'='*80}")
            
            # 1. 检查数据可用性
            print("检查数据可用性...")
            data_report = check_data_availability(disease_code, compound, self.config_path)
            
            if not data_report['data_available']:
                print(f"❌ 数据不可用:")
                for issue in data_report['issues']:
                    print(f"   - {issue}")
                return None
            
            print("✅ 数据检查通过")
            
            # 2. 准备模型数据
            print("准备模型数据...")
            model_data = self.data_loader.prepare_model_data(
                disease_code=disease_code,
                compound=compound,
                model_type=model_type,
                lag_years=lag_years,
                measure_type=measure_type,
                estimate_type=estimate_type  # <--- 修改
            )
            
            # 3. 拟合模型
            print("拟合BYM2模型...")
            model, trace = self.model_fitter.run_analysis(model_data)
            
            # 4. 提取结果
            print("提取分析结果...")
            output_file = self.result_extractor.process_single_analysis(
                trace, model, model_data
            )
            
            # 5. 打印分析摘要
            print_analysis_summary(model_data, trace)
            
            print(f"✅ 分析完成！结果保存到: {output_file}")
            return output_file
            
        except Exception as e:
            print(f"❌ 分析失败: {e}")
            traceback.print_exc()
            return None
    
    def run_batch_analysis(self, disease_codes: List[str], 
                          compounds: List[str], model_types: List[str],
                          lag_years_list: List[int],
                          measure_types: List[str] = ['Weight'],
                          estimate_types: List[str] = ['avg']) -> Dict: # <--- 修改
        """
        运行批量分析
        
        Parameters
        ----------
        disease_codes : List[str]
            疾病编码列表
        compounds : List[str]
            化合物列表
        model_types : List[str]
            模型类型列表
        lag_years_list : List[int]
            滞后年份列表
        measure_types : List[str]
            测量类型列表
            
        Returns
        -------
        Dict
            批量分析结果报告
        """
        print(f"\n{'='*80}")
        print(f"WDP PyMC批量分析")
        print(f"{'='*80}")
        print(f"疾病: {', '.join(disease_codes)}")
        print(f"化合物: {', '.join(compounds[:5])}{'...' if len(compounds) > 5 else ''}")
        print(f"模型: {', '.join(model_types)}")
        print(f"滞后: {', '.join(map(str, lag_years_list))}")
        print(f"测量: {', '.join(measure_types)}")
        print(f"估算: {', '.join(estimate_types)}")
        
        # 计算总分析数
        total_analyses = (len(disease_codes) * len(compounds) * 
                         len(model_types) * len(lag_years_list) * 
                         len(measure_types) * len(estimate_types)) # <--- 修改
        print(f"总分析数: {total_analyses}")
        
        # 结果统计
        results = {
            'total': 0,
            'successful': 0,
            'failed': 0,
            'output_files': [],
            'errors': []
        }
        
        # 执行批量分析
        from itertools import product
        
        # 使用product生成所有组合
        combinations = list(product(
            disease_codes,
            compounds,
            model_types,
            lag_years_list,
            measure_types,
            estimate_types # <--- 修改
        ))
        
        for i, (disease, compound, model_type, lag_years, measure_type, estimate_type) in enumerate(combinations): # <--- 修改
            results['total'] += 1
            
            # --- START: 添加的日志 ---
            print(f"\n{'='*40}")
            print(f"===> 开始分析 [{i+1}/{len(combinations)}]: {disease}, {compound}, {model_type}, lag={lag_years}, measure={measure_type}, estimate={estimate_type}")
            print(f"{'='*40}")
            # --- END: 添加的日志 ---
            
            output_file = self.run_single_analysis(
                disease, compound, model_type, 
                lag_years, measure_type, estimate_type # <--- 修改
            )
            
            if output_file:
                # --- START: 添加的日志 ---
                print(f"===> 分析 [{i+1}/{len(combinations)}] 成功完成")
                # --- END: 添加的日志 ---
                results['successful'] += 1
                results['output_files'].append(output_file)
            else:
                # --- START: 添加的日志 ---
                print(f"===> 分析 [{i+1}/{len(combinations)}] 失败")
                # --- END: 添加的日志 ---
                results['failed'] += 1
                results['errors'].append({
                    'disease': disease,
                    'compound': compound,
                    'model': model_type,
                    'lag': lag_years,
                    'measure': measure_type,
                    'estimate': estimate_type # <--- 新增
                })
        
        # 打印批量分析摘要
        print(f"\n{'='*80}")
        print(f"批量分析完成")
        print(f"{'='*80}")
        print(f"总分析数: {results['total']}")
        print(f"成功: {results['successful']}")
        print(f"失败: {results['failed']}")
        print(f"成功率: {results['successful']/results['total']*100:.1f}%")
        
        if results['failed'] > 0:
            print(f"\n失败的分析:")
            for error in results['errors']:
                print(f"  - {error['disease']}-{error['compound']}-{error['model']}")
        
        # 创建汇总表
        if results['successful'] > 0:
            print("\n创建汇总表...")
            summary_df = self.result_extractor.create_summary_table()
            print(f"汇总表包含 {len(summary_df)} 条结果")
        
        return results


def create_parser() -> argparse.ArgumentParser:
    """
    创建命令行参数解析器
    
    Returns
    -------
    argparse.ArgumentParser
        参数解析器
    """
    parser = argparse.ArgumentParser(
        description="WDP PyMC贝叶斯空间时间分析",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 单个分析
  python Main.py --disease C81-C96 --compound 24D --model M0 --lag 5
  
  # 多模型分析
  python Main.py --disease C81-C96 --compound 24D --model M0,M1,M2,M3 --lag 5,10
  
  # 批量分析
  python Main.py --disease C81-C96,C50 --compound 24D,Atrazine --model M0,M3 --lag 5
  
  # 生产模式
  python Main.py --disease C81-C96 --compound 24D --model M3 --lag 5 --sampling-mode production
        """
    )
    
    # Required arguments
    parser.add_argument('--disease', '--disease-code', 
                       type=str, required=True,
                       help='Disease code (e.g., C81-C96, C50, C34)')
    
    parser.add_argument('--compound', '--exposure',
                       type=str, required=True, 
                       help='Compound name (e.g., 24D, Atrazine, Glyphosate)')
    
    # Optional arguments
    parser.add_argument('--model', '--model-type',
                       type=str, default='M0',
                       help='Model type, comma-separated (default: M0)')
    
    parser.add_argument('--lag', '--lag-years',
                       type=str, default='5',
                       help='Lag years, comma-separated (default: 5)')
    
    parser.add_argument('--measure', '--measure-type',
                       type=str, default='Weight',
                       choices=['Weight', 'Density', 'Weight,Density'],
                       help='Measure type (default: Weight)')
    
    parser.add_argument('--estimate', '--estimate-types',
                       type=str, default='avg',
                       help='Exposure estimate types, comma-separated (e.g., min,avg,max)')
    
    parser.add_argument('--sampling-mode', 
                       type=str, default='test',
                       choices=['test', 'production'],
                       help='Sampling mode (default: test)')
    
    parser.add_argument('--cores',
                        type=str, default=None,
                        help='Number of cores, integer or auto (default: auto=all available)')
    
    parser.add_argument('--chains',
                        type=str, default=None,
                        help='Number of chains, integer or auto (default: auto=same as cores)')
    
    parser.add_argument('--draws',
                        type=int, default=None,
                        help='Number of draws per chain')
    
    parser.add_argument('--tune',
                        type=int, default=None,
                        help='Number of tuning steps per chain')
    
    parser.add_argument('--target-accept',
                        type=float, default=None,
                        help='Target acceptance rate (e.g., 0.9)')
    
    parser.add_argument('--output-dir',
                       type=str, default=None,
                       help='Output directory (default: Result/PyMC_Results)')
    
    parser.add_argument('--config-path',
                       type=str, default=None,
                       help='Config file path (default: project_root/config.yaml)')
    
    parser.add_argument('--verbose', '-v',
                       action='store_true',
                       help='Verbose output')
    
    parser.add_argument('--dry-run',
                       action='store_true',
                       help='Check data availability only, do not run analysis')
    
    return parser


def main():
    """主函数"""
    # 解析命令行参数
    parser = create_parser()
    args = parser.parse_args()
    
    # 设置日志级别
    if args.verbose:
        import logging
        logging.basicConfig(level=logging.DEBUG)
    
    try:
        # 处理路径参数
        config_path = Path(args.config_path) if args.config_path else None
        output_dir = Path(args.output_dir) if args.output_dir else None
        
        # 验证疾病编码
        disease_codes = [d.strip() for d in args.disease.split(',')]
        for disease in disease_codes:
            if not validate_disease_code(disease, config_path):
                print(f"❌ 无效的疾病编码: {disease}")
                sys.exit(1)
        
        # 解析参数
        compounds = parse_compound_list(args.compound)
        model_types = parse_model_list(args.model)
        lag_years_list = parse_lag_years(args.lag)
        measure_types = [m.strip() for m in args.measure.split(',')]
        
        # <--- 新增：解析 estimate 参数 --->
        estimate_types = [e.strip().lower() for e in args.estimate.split(',')]
        
        # 验证模型类型
        for model_type in model_types:
            if not validate_model_type(model_type, config_path):
                print(f"❌ 无效的模型类型: {model_type}")
                sys.exit(1)
        
        print(f"WDP PyMC分析系统启动")
        print(f"疾病: {disease_codes}")
        print(f"化合物: {compounds}")
        print(f"模型: {model_types}")
        print(f"滞后: {lag_years_list}")
        print(f"测量: {measure_types}")
        
        # Dry run 模式
        if args.dry_run:
            print("\n=== 数据可用性与暴露列名解析 (Dry Run) ===")
            loader = WDPDataLoader(
                config_path=config_path
            )
            try:
                pesticide_df_weight = loader.load_pesticide_data("Weight")
                pesticide_df_density = loader.load_pesticide_data("Density")
            except Exception as e:
                print(f"❌ 无法加载农药数据: {e}")
                pesticide_df_weight = None
                pesticide_df_density = None
            
            # 使用第一个lag年份进行测试
            first_lag_year = lag_years_list[0]
            print(f"使用滞后年份 {first_lag_year} 进行列名解析测试")
            
            from itertools import product
            combinations = list(product(disease_codes, compounds, measure_types, estimate_types))

            for disease, compound, measure, estimate in combinations:
                print(f"\n--- 检查: Disease={disease}, Compound={compound}, Measure={measure}, Estimate={estimate} ---")
                report = check_data_availability(disease, compound, config_path)
                status = "✅" if report.get('data_available', False) else "❌"
                print(f"  数据可用性: {status}")
                if report.get('issues'):
                    for issue in report['issues']:
                        print(f"    - 问题: {issue}")

                pesticide_df = pesticide_df_weight if measure == 'Weight' else pesticide_df_density
                if pesticide_df is not None and report.get('pesticide_available', True):
                    try:
                        # 使用第一个lag年份进行测试
                        _lagged_df, found_col, _est_type = loader.calculate_lagged_exposure(
                            pesticide_df, compound, lag_years=first_lag_year, estimate_type=estimate
                        )
                        print(f"  解析结果: 输入 '{compound}' (estimate='{estimate}') 将使用数据列 '{found_col}'")
                    except ValueError as e:
                        print(f"  解析结果: ❌ 无法解析 '{compound}' (estimate='{estimate}'): {e}")
            return
        
        # 初始化分析系统
        analysis = WDPPyMCAnalysis(
            config_path=config_path,
            output_dir=output_dir,
            sampling_mode=args.sampling_mode,
            cores=args.cores,
            chains=args.chains,
            draws=args.draws,
            tune=args.tune,
            target_accept=args.target_accept
        )
        
        # 运行分析
        if (len(disease_codes) == 1 and len(compounds) == 1 and 
            len(model_types) == 1 and len(lag_years_list) == 1 and 
            len(measure_types) == 1 and len(estimate_types) == 1): # <--- 修改此条件
            
            # 单个分析
            output_file = analysis.run_single_analysis(
                disease_codes[0], compounds[0], model_types[0],
                lag_years_list[0], measure_types[0], estimate_types[0] # <--- 传递新参数
            )
            
            if output_file:
                print(f"\n🎉 分析成功完成！")
                print(f"结果文件: {output_file}")
            else:
                print(f"\n💥 分析失败!")
                sys.exit(1)
        
        else:
            # 批量分析
            results = analysis.run_batch_analysis(
                disease_codes, compounds, model_types,
                lag_years_list, measure_types, estimate_types # <--- 传递新参数
            )
            
            if results['successful'] > 0:
                print(f"\n🎉 批量分析完成！")
                print(f"成功分析: {results['successful']}/{results['total']}")
            else:
                print(f"\n💥 所有分析都失败了!")
                sys.exit(1)
    
    except KeyboardInterrupt:
        print(f"\n⏹  用户中断分析")
        sys.exit(1)
    
    except Exception as e:
        print(f"\n💥 系统错误: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()