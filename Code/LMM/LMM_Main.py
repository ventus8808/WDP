#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LMM主控制脚本 - 完整分析流程协调

主要功能:
1. 协调数据整合、模型分析和结果输出的完整流程
2. 支持命令行参数控制分析范围
3. 输出4个结果表（2时间跨度×2滞后期）
4. 提供进度监控和错误处理

使用方法:
python LMM_Main.py [--cancer-types C00_C97,C34,C50] [--output-dir /path/to/results]

输出文件:
- LMM_5_year_lag_Results_YYYYMMDD_HHMMSS.csv
- LMM_10_year_lag_Results_YYYYMMDD_HHMMSS.csv
- LMM_EQI0005_Results_YYYYMMDD_HHMMSS.csv
- LMM_EQI0610_Results_YYYYMMDD_HHMMSS.csv
"""

import sys
import argparse
from pathlib import Path
import logging
from datetime import datetime
from typing import List, Optional

# 导入自定义模块
try:
    from LMM_Data import LMMDataIntegrator
    from LMM_Model import LMMAnalyzer
    from LMM_Result import LMMResultFormatter
except ImportError as e:
    print(f"导入模块失败: {e}")
    print("请确保在正确的目录中运行脚本")
    sys.exit(1)

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('lmm_analysis.log', encoding='utf-8')
    ]
)
logger = logging.getLogger(__name__)

class LMMPipeline:
    """LMM分析完整流程管理器"""
    
    def __init__(self, data_file: Optional[str] = None, output_dir: Optional[str] = None, use_mice: bool = False):
        """初始化流程管理器"""
        self.use_mice = use_mice
        self.project_root = Path(__file__).resolve().parents[2]
        
        # 支持自定义数据文件
        if data_file:
            self.data_file = Path(data_file)
        else:
            self.data_file = self.project_root / "Data" / "Processed" / "df_EQI_AAMR" / "EQI_AAMR_Point.csv"
        
        # 支持自定义输出目录
        if output_dir:
            self.default_output_dir = Path(output_dir)
        else:
            self.default_output_dir = self.project_root / "Result" / "EQI_LMM"
        
        # 分析配置
        self.all_scenarios = [
            'EQI0005_AAMR2006_2010',
            'EQI0005_AAMR2011_2015',
            'EQI0610_AAMR2011_2015', 
            'EQI0610_AAMR2016_2020'
        ]
        
        self.all_cancer_types = [
            'C00_C97',     # All Cancers
            'C15_C26',     # Digestive Organs
            'C18_C21',     # Colorectal
            'C25',         # Pancreas
            'C30_C39',     # Respiratory
            'C34',         # Lung Cancer
            'C50',         # Breast Cancer
            'C51_C58',     # Female Genital
            'C60_C63',     # Male Genital
            'C61',         # Prostate
            'C64_C68',     # Urinary Tract
            'C76_C80',     # Ill-defined
            'C81_C96'      # Lymphoid/Hematologic
        ]
        
        # 组件实例
        self.data_integrator = None
        self.model_analyzer = None
        self.result_formatter = None
        
    def setup_components(self):
        """设置各组件"""
        logger.info("初始化分析组件...")
        
        # 初始化结果格式化器
        self.result_formatter = LMMResultFormatter(use_mice=self.use_mice)
        logger.info("结果格式化器初始化完成")
        
    def check_data_availability(self) -> bool:
        """检查数据文件可用性"""
        logger.info("检查数据文件...")
        
        if not self.data_file.exists():
            logger.error(f"数据文件不存在: {self.data_file}")
            logger.info("尝试重新生成数据...")
            
            # 尝试重新生成数据
            return self.regenerate_data()
        else:
            logger.info(f"数据文件存在: {self.data_file}")
            return True
    
    def regenerate_data(self) -> bool:
        """重新生成数据文件"""
        logger.info("开始重新生成数据...")
        
        try:
            self.data_integrator = LMMDataIntegrator(use_mice=self.use_mice)
            success = self.data_integrator.process_all()
            
            if success and self.data_file.exists():
                logger.info("数据重新生成成功")
                return True
            else:
                logger.error("数据重新生成失败")
                return False
                
        except Exception as e:
            logger.error(f"数据生成过程出错: {e}")
            return False
    
    def run_model_analysis(self, cancer_types: List[str]) -> Optional[dict]:
        """运行模型分析"""
        logger.info("=== 开始模型分析阶段 ===")
        
        try:
            # 初始化模型分析器
            self.model_analyzer = LMMAnalyzer(str(self.data_file))
            
            # 加载数据
            if not self.model_analyzer.load_data():
                logger.error("模型分析器数据加载失败")
                return None
            
            # 运行分析
            logger.info(f"开始分析 {len(self.all_scenarios)} 个场景, {len(cancer_types)} 种癌症类型")
            results = self.model_analyzer.run_full_analysis(self.all_scenarios, cancer_types)
            
            if results:
                logger.info("模型分析完成")
                return results
            else:
                logger.error("模型分析失败")
                return None
                
        except Exception as e:
            logger.error(f"模型分析过程出错: {e}")
            return None
    
    def run_model_analysis_with_incremental_save(self, cancer_types: List[str], output_dir: str, apply_correction: bool = True) -> bool:
        """运行模型分析阶段，边分析边保存"""
        logger.info("=== 开始增量模型分析阶段 ===")
        
        try:
            # 初始化模型分析器
            self.model_analyzer = LMMAnalyzer(str(self.data_file))
            
            # 加载数据
            if not self.model_analyzer.load_data():
                logger.error("模型分析器数据加载失败")
                return False
            
            logger.info(f"开始分析 {len(self.all_scenarios)} 个场景, {len(cancer_types)} 种癌症类型")
            
            # 为每种癌症类型分别运行并保存
            for cancer_type in cancer_types:
                logger.info(f"=== 开始分析癌症类型: {cancer_type} ===")
                
                # 运行单个癌症类型的分析
                single_cancer_results = self.model_analyzer.run_full_analysis(self.all_scenarios, [cancer_type])
                
                if single_cancer_results:
                    # 立即保存结果
                    saved_files = self.result_formatter.save_results(
                        single_cancer_results, 
                        [cancer_type], 
                        output_dir,
                        apply_correction=apply_correction
                    )
                    
                    if saved_files:
                        logger.info(f"癌症类型 {cancer_type} 分析完成并保存")
                        for file_key, file_path in saved_files.items():
                            logger.info(f"  保存文件: {file_path}")
                    else:
                        logger.warning(f"癌症类型 {cancer_type} 保存失败")
                else:
                    logger.error(f"癌症类型 {cancer_type} 分析失败")
            
            logger.info("=== 增量分析完成 ===")
            return True
                
        except Exception as e:
            logger.error(f"增量分析过程出错: {e}")
            return False

    def save_results(self, model_results: dict, cancer_types: List[str], output_dir: str, apply_correction: bool = True) -> dict:
        """保存分析结果"""
        logger.info("=== 开始结果输出阶段 ===")
        
        if apply_correction:
            logger.info("将应用多重校正（FDR）")
        else:
            logger.info("跳过多重校正，仅输出原始p值")
        
        try:
            saved_files = self.result_formatter.save_results(
                model_results, 
                cancer_types, 
                output_dir,
                apply_correction=apply_correction
            )
            
            logger.info(f"结果已保存到 {len(saved_files)} 个文件")
            return saved_files
            
        except Exception as e:
            logger.error(f"结果保存过程出错: {e}")
            return {}
    
    def run_full_pipeline(self, cancer_types: List[str], output_dir: str, apply_correction: bool = True) -> bool:
        """运行完整分析流程"""
        logger.info("=== 开始LMM完整分析流程 ===")
        logger.info(f"分析癌症类型: {', '.join(cancer_types)}")
        logger.info(f"输出目录: {output_dir}")
        
        start_time = datetime.now()
        
        try:
            # 1. 设置组件
            self.setup_components()
            
            # 2. 检查数据可用性
            if not self.check_data_availability():
                logger.error("数据准备失败")
                return False
            
            # 3. 运行模型分析
            model_results = self.run_model_analysis(cancer_types)
            if not model_results:
                logger.error("模型分析失败")
                return False
            
            # 4. 保存结果
            saved_files = self.save_results(model_results, cancer_types, output_dir, apply_correction)
            if not saved_files:
                logger.error("结果保存失败")
                return False
            
            # 5. 总结
            end_time = datetime.now()
            duration = end_time - start_time
            
            logger.info("=== 分析流程完成 ===")
            logger.info(f"总耗时: {duration}")
            logger.info(f"输出文件数: {len(saved_files)}")
            
            # 显示输出文件列表
            logger.info("生成的结果文件:")
            for file_type, file_path in saved_files.items():
                logger.info(f"  {file_type}: {Path(file_path).name}")
            
            return True
            
        except Exception as e:
            logger.error(f"分析流程出错: {e}")
            return False
    
    def run_incremental_pipeline(self, cancer_types: List[str], output_dir: str, apply_correction: bool = True) -> bool:
        """运行增量分析流程（边分析边保存）"""
        logger.info("=== 开始LMM增量分析流程 ===")
        logger.info(f"分析癌症类型: {', '.join(cancer_types)}")
        logger.info(f"输出目录: {output_dir}")
        
        start_time = datetime.now()
        
        try:
            # 1. 设置组件
            self.setup_components()
            
            # 2. 检查数据可用性
            if not self.check_data_availability():
                logger.error("数据准备失败")
                return False
            
            # 3. 运行增量模型分析（边分析边保存）
            success = self.run_model_analysis_with_incremental_save(cancer_types, output_dir, apply_correction)
            if not success:
                logger.error("增量分析失败")
                return False
            
            # 4. 总结
            end_time = datetime.now()
            duration = end_time - start_time
            
            logger.info("=== 增量分析流程完成 ===")
            logger.info(f"总耗时: {duration}")
            logger.info(f"成功分析癌症类型数: {len(cancer_types)}")
            
            return True
            
        except Exception as e:
            logger.error(f"增量分析流程出错: {e}")
            return False
    
    def run_quick_test(self) -> bool:
        """运行快速测试（仅分析总癌症）"""
        logger.info("=== 运行快速测试 ===")
        
        test_cancer_types = ['C00_C97']  # 仅分析总癌症
        test_output_dir = self.default_output_dir
        
        return self.run_full_pipeline(test_cancer_types, str(test_output_dir))


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="LMM多层回归分析主程序",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python LMM_Main.py                                    # 运行所有癌症类型分析（原始数据）
  python LMM_Main.py --mice                            # 使用MICE插补数据分析
  python LMM_Main.py --cancer-types C00_C97            # 仅分析总癌症
  python LMM_Main.py --cancer-types C00_C97,C34,C50    # 分析特定癌症类型
  python LMM_Main.py --mice --cancer-types C00_C97     # 用MICE数据分析总癌症
  python LMM_Main.py --test                             # 运行快速测试
  python LMM_Main.py --data-file /path/to/data.csv     # 指定数据文件
  python LMM_Main.py --output-dir /custom/path         # 指定输出目录
        """
    )
    
    parser.add_argument(
        '--cancer-types',
        type=str,
        help='要分析的癌症类型，用逗号分隔 (如: C00_C97,C34,C50)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        help='结果输出目录路径'
    )
    
    parser.add_argument(
        '--test',
        action='store_true',
        help='运行快速测试（仅分析总癌症）'
    )
    
    parser.add_argument(
        '--list-cancer-types',
        action='store_true',
        help='列出所有可用的癌症类型'
    )
    
    parser.add_argument(
        '--data-file',
        type=str,
        help='指定输入数据文件路径 (默认: EQI_LMM_Delete_df.csv)'
    )
    
    parser.add_argument(
        '--mice',
        action='store_true',
        help='使用MICE插补数据进行分析'
    )
    
    parser.add_argument(
        '--no-correction',
        action='store_true',
        help='跳过多重校正，仅输出原始p值结果'
    )
    
    parser.add_argument(
        '--correction-method',
        type=str,
        default='fdr_bh',
        choices=['fdr_bh', 'bonferroni', 'holm', 'sidak'],
        help='多重校正方法 (默认: fdr_bh)'
    )
    
    return parser.parse_args()


def main():
    """主函数"""
    print("=" * 60)
    print("LMM多层回归分析 - 县级环境质量与癌症死亡率关联研究")
    print("=" * 60)
    
    # 解析命令行参数
    args = parse_arguments()
    
    # 确定数据文件
    data_file = None
    if args.mice:
        # 使用MICE插补数据
        project_root = Path(__file__).resolve().parents[2]
        data_file = str(project_root / "Data" / "Processed" / "df_EQI_AAMR" / "EQI_AAMR_Point_MICE.csv")
        print(f"📊 使用MICE插补数据: EQI_AAMR_Point_MICE.csv")
    elif args.data_file:
        data_file = args.data_file
        print(f"📊 使用指定数据文件: {data_file}")
    
    # 确定输出目录
    output_dir = None
    if args.mice and not args.output_dir:
        # 使用MICE数据时，默认输出到MICE结果目录
        project_root = Path(__file__).resolve().parents[2]
        output_dir = str(project_root / "Result" / "EQI_LMM_MICE")
        print(f"📁 输出目录: Result/EQI_LMM_MICE")
    elif args.output_dir:
        output_dir = args.output_dir
        print(f"📁 输出目录: {output_dir}")
    
    # 创建流程管理器
    pipeline = LMMPipeline(data_file=data_file, output_dir=output_dir, use_mice=args.mice)
    
    # 处理列出癌症类型的请求
    if args.list_cancer_types:
        print("可用的癌症类型:")
        for i, cancer_type in enumerate(pipeline.all_cancer_types, 1):
            print(f"{i:2d}. {cancer_type}")
        return
    
    # 处理快速测试
    if args.test:
        success = pipeline.run_quick_test()
        if success:
            print("✅ 快速测试完成!")
        else:
            print("❌ 快速测试失败!")
        return
    
    # 确定要分析的癌症类型
    if args.cancer_types:
        cancer_types = [ct.strip() for ct in args.cancer_types.split(',')]
        # 验证癌症类型有效性
        invalid_types = [ct for ct in cancer_types if ct not in pipeline.all_cancer_types]
        if invalid_types:
            print(f"❌ 无效的癌症类型: {', '.join(invalid_types)}")
            print("使用 --list-cancer-types 查看所有可用类型")
            return
    else:
        cancer_types = pipeline.all_cancer_types
    
    # 确定最终输出目录
    if not output_dir:
        output_dir = str(pipeline.default_output_dir)
    
    # 确定多重校正设置
    apply_correction = not args.no_correction
    
    # 运行增量分析（边分析边保存）
    print(f"开始分析 {len(cancer_types)} 种癌症类型...")
    if apply_correction:
        print(f"将应用多重校正: {args.correction_method}")
    else:
        print("跳过多重校正")
    
    success = pipeline.run_incremental_pipeline(cancer_types, output_dir, apply_correction)
    
    if success:
        print("✅ 分析完成!")
    else:
        print("❌ 分析失败!")
        print("请检查日志文件 lmm_analysis.log 获取详细信息")


if __name__ == "__main__":
    main()