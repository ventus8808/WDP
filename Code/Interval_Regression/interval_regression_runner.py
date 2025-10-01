#!/usr/bin/env python3
"""
Interval Regression Analysis Runner
==================================

区间回归分析的主控制器

功能：
- 协调Python数据准备和R模型拟合
- 管理不同分析场景
- 整合结果和报告生成
"""

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict
import yaml
import json
from datetime import datetime

from interval_data_loader import IntervalRegressionDataLoader


class IntervalRegressionAnalyzer:
    """区间回归分析器"""
    
    def __init__(self):
        """初始化分析器"""
        self.project_root = Path(__file__).resolve().parents[2]
        self.code_dir = Path(__file__).parent
        self.data_dir = self.code_dir / "data"
        self.results_dir = self.code_dir / "results"
        
        # 创建必要目录
        self.data_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)
        
        self.config = self._load_config()
        self.data_loader = IntervalRegressionDataLoader()
        
    def _load_config(self) -> Dict:
        """加载配置文件"""
        config_path = self.project_root / 'config.yaml'
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    
    def run_analysis_scenario(self,
                            scenario_name: str,
                            cancer_types: List[str],
                            analysis_type: str = "total_eqi",
                            rucc_filter: Optional[List[int]] = None,
                            eqi_domains: Optional[List[str]] = None) -> Dict:
        """
        运行单个分析场景
        
        Parameters:
        -----------
        scenario_name : str
            场景名称，用于文件命名
        cancer_types : List[str]
            要分析的癌症类型
        analysis_type : str
            分析类型: "total_eqi", "domain_specific", "rucc_stratified"
        rucc_filter : List[int], optional
            城乡分类筛选
        eqi_domains : List[str], optional
            EQI领域列表（用于领域特异性分析）
            
        Returns:
        --------
        Dict
            分析结果摘要
        """
        print(f"\n🎯 运行分析场景: {scenario_name}")
        print("=" * 50)
        
        results = {
            'scenario_name': scenario_name,
            'analysis_type': analysis_type,
            'cancer_types': cancer_types,
            'timestamp': datetime.now().isoformat(),
            'success': False,
            'results_files': [],
            'errors': []
        }
        
        try:
            # 1. 准备数据
            print("📊 准备分析数据...")
            analysis_data = self.data_loader.prepare_analysis_data(
                cancer_types=cancer_types,
                analysis_scenario='EQI0610_AAMR2016_2020',
                rucc_filter=rucc_filter
            )
            
            if len(analysis_data) == 0:
                raise ValueError("分析数据为空")
            
            # 2. 导出R分析数据
            r_data_file = self.data_dir / f"{scenario_name}_data.csv"
            self.data_loader.export_for_r_analysis(r_data_file)
            
            # 3. 准备R分析参数
            r_params = self._prepare_r_parameters(
                scenario_name, analysis_type, cancer_types, eqi_domains
            )
            
            # 4. 运行R分析
            r_results = self._run_r_analysis(r_data_file, r_params, scenario_name)
            
            # 5. 处理结果
            results.update(r_results)
            results['success'] = True
            
            print(f"✅ 场景 {scenario_name} 分析完成")
            
        except Exception as e:
            error_msg = f"场景 {scenario_name} 分析失败: {str(e)}"
            print(f"❌ {error_msg}")
            results['errors'].append(error_msg)
        
        return results
    
    def _prepare_r_parameters(self, 
                            scenario_name: str,
                            analysis_type: str,
                            cancer_types: List[str],
                            eqi_domains: Optional[List[str]] = None) -> Dict:
        """准备R分析参数"""
        
        params = {
            'scenario_name': scenario_name,
            'analysis_type': analysis_type,
            'cancer_types': cancer_types,
            'data_file': str(self.data_dir / f"{scenario_name}_data.csv"),
            'output_dir': str(self.results_dir),
            'chains': 4,
            'iter': 2000,
            'warmup': 1000,
            'cores': 4
        }
        
        if eqi_domains:
            params['eqi_domains'] = eqi_domains
        
        return params
    
    def _run_r_analysis(self, 
                      data_file: Path,
                      parameters: Dict,
                      scenario_name: str) -> Dict:
        """运行R分析"""
        
        print(f"🔄 运行R区间回归分析...")
        
        # 保存参数到JSON文件
        params_file = self.data_dir / f"{scenario_name}_params.json"
        with open(params_file, 'w') as f:
            json.dump(parameters, f, indent=2)
        
        # R脚本路径
        r_script = self.code_dir / "interval_regression_analysis.R"
        
        if not r_script.exists():
            # 如果R脚本不存在，创建一个基本版本
            self._create_r_analysis_script(r_script)
        
        try:
            # 运行R脚本
            cmd = [
                "Rscript", 
                str(r_script),
                str(params_file)
            ]
            
            print(f"  执行命令: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                cwd=str(self.project_root),
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时
            )
            
            if result.returncode != 0:
                raise RuntimeError(f"R脚本执行失败:\n{result.stderr}")
            
            print(f"✅ R分析完成")
            
            # 处理R输出
            r_results = {
                'r_output': result.stdout,
                'r_errors': result.stderr if result.stderr else None,
                'results_files': self._find_result_files(scenario_name)
            }
            
            return r_results
            
        except subprocess.TimeoutExpired:
            raise RuntimeError("R分析超时")
        except Exception as e:
            raise RuntimeError(f"R分析执行错误: {str(e)}")
    
    def _create_r_analysis_script(self, script_path: Path):
        """创建R分析脚本"""
        
        r_code = '''#!/usr/bin/env Rscript
# Interval Regression Analysis Script
# ==================================

# 加载必要的包
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(ggplot2)
})

# 读取参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("请提供参数文件路径")
}

params_file <- args[1]
params <- fromJSON(params_file)

cat("🔄 区间回归分析\\n")
cat("================================\\n")
cat("场景名称:", params$scenario_name, "\\n")
cat("分析类型:", params$analysis_type, "\\n")

# 读取数据
cat("📁 读取数据:", params$data_file, "\\n")
data <- read_csv(params$data_file, show_col_types = FALSE)

cat("📊 数据维度:", nrow(data), "×", ncol(data), "\\n")

# 数据预处理
data <- data %>%
  mutate(
    # 确保分类变量是factor
    EQI_quintile = as.factor(EQI_quintile),
    State = as.factor(State),
    Cancer_Type = as.factor(Cancer_Type),
    # 创建区间变量
    cens = ifelse(AAMR_lower == AAMR_upper, "none", "interval"),
    AAMR_response = ifelse(cens == "none", AAMR_lower, NA)
  )

# 运行不同癌症类型的分析
results_list <- list()

for (cancer in params$cancer_types) {
  cat("\\n🎯 分析癌症类型:", cancer, "\\n")
  
  # 筛选数据
  cancer_data <- data %>% filter(Cancer_Type == cancer)
  
  if (nrow(cancer_data) < 100) {
    cat("⚠️  数据量不足，跳过\\n")
    next
  }
  
  cat("  数据量:", nrow(cancer_data), "\\n")
  
  # 构建模型公式
  if (params$analysis_type == "total_eqi") {
    formula <- bf(
      AAMR_response | cens(cens, AAMR_lower, AAMR_upper) ~ 
        EQI_quintile + Smoking_Rate_std + (1 | State),
      family = gaussian()
    )
  } else {
    # 其他分析类型的公式可以在这里扩展
    formula <- bf(
      AAMR_response | cens(cens, AAMR_lower, AAMR_upper) ~ 
        EQI_quintile + Smoking_Rate_std + (1 | State),
      family = gaussian()
    )
  }
  
  cat("📐 模型公式:", deparse(formula$formula), "\\n")
  
  # 拟合模型
  cat("⏳ 拟合模型...\\n")
  
  tryCatch({
    fit <- brm(
      formula = formula,
      data = cancer_data,
      chains = params$chains,
      iter = params$iter,
      warmup = params$warmup,
      cores = params$cores,
      seed = 12345,
      control = list(adapt_delta = 0.95),
      silent = TRUE,
      refresh = 0
    )
    
    # 提取结果
    fixed_effects <- fixef(fit, summary = TRUE)
    
    # 保存结果
    result_df <- data.frame(
      Scenario = params$scenario_name,
      Cancer_Type = cancer,
      Parameter = rownames(fixed_effects),
      Estimate = fixed_effects[, "Estimate"],
      Lower_CI = fixed_effects[, "Q2.5"],
      Upper_CI = fixed_effects[, "Q97.5"],
      Rhat = fixed_effects[, "Rhat"]
    )
    
    results_list[[cancer]] <- result_df
    
    # 保存详细结果
    output_file <- file.path(params$output_dir, 
                           paste0(params$scenario_name, "_", cancer, "_results.csv"))
    write_csv(result_df, output_file)
    
    cat("✅ 完成，结果已保存\\n")
    
  }, error = function(e) {
    cat("❌ 模型拟合失败:", conditionMessage(e), "\\n")
  })
}

# 合并所有结果
if (length(results_list) > 0) {
  combined_results <- do.call(rbind, results_list)
  
  # 保存合并结果
  combined_file <- file.path(params$output_dir, 
                           paste0(params$scenario_name, "_combined_results.csv"))
  write_csv(combined_results, combined_file)
  
  cat("\\n💾 合并结果已保存:", combined_file, "\\n")
  
  # 显示关键结果
  cat("\\n📈 关键结果摘要:\\n")
  eqi_results <- combined_results %>% 
    filter(grepl("EQI_quintile", Parameter)) %>%
    select(Cancer_Type, Parameter, Estimate, Lower_CI, Upper_CI)
  
  print(eqi_results)
}

cat("\\n🎉 分析完成!\\n")
'''
        
        with open(script_path, 'w') as f:
            f.write(r_code)
        
        print(f"📝 创建R分析脚本: {script_path}")
    
    def _find_result_files(self, scenario_name: str) -> List[str]:
        """查找结果文件"""
        result_files = []
        
        # 查找该场景的所有结果文件
        pattern = f"{scenario_name}*results.csv"
        for file_path in self.results_dir.glob(pattern):
            result_files.append(str(file_path))
        
        return result_files
    
    def run_comprehensive_analysis(self) -> Dict:
        """运行综合分析"""
        print("🚀 开始综合区间回归分析")
        print("=" * 60)
        
        # 定义分析场景
        scenarios = [
            {
                'name': 'primary_cancers_total_eqi',
                'cancer_types': ['C00_C97', 'C34', 'C50', 'C61'],
                'analysis_type': 'total_eqi',
                'description': '主要癌症类型 - 总EQI分析'
            },
            {
                'name': 'digestive_cancers_total_eqi', 
                'cancer_types': ['C15_C26', 'C18_C21', 'C25'],
                'analysis_type': 'total_eqi',
                'description': '消化系统癌症 - 总EQI分析'
            },
            {
                'name': 'all_cancers_urban_rural',
                'cancer_types': ['C00_C97'],
                'analysis_type': 'total_eqi',
                'rucc_filter': [1, 2, 3],  # 城市地区
                'description': '全部癌症 - 城市地区分析'
            }
        ]
        
        comprehensive_results = {
            'start_time': datetime.now().isoformat(),
            'scenarios': [],
            'summary': {},
            'errors': []
        }
        
        # 运行每个场景
        for scenario in scenarios:
            print(f"\n📋 {scenario['description']}")
            
            try:
                result = self.run_analysis_scenario(
                    scenario_name=scenario['name'],
                    cancer_types=scenario['cancer_types'],
                    analysis_type=scenario['analysis_type'],
                    rucc_filter=scenario.get('rucc_filter')
                )
                
                comprehensive_results['scenarios'].append(result)
                
            except Exception as e:
                error_msg = f"场景 {scenario['name']} 失败: {str(e)}"
                print(f"❌ {error_msg}")
                comprehensive_results['errors'].append(error_msg)
        
        # 生成摘要
        comprehensive_results['end_time'] = datetime.now().isoformat()
        comprehensive_results['summary'] = self._generate_analysis_summary(
            comprehensive_results['scenarios']
        )
        
        # 保存综合结果
        summary_file = self.results_dir / "comprehensive_analysis_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(comprehensive_results, f, indent=2, ensure_ascii=False)
        
        print(f"\n💾 综合分析摘要已保存: {summary_file}")
        
        return comprehensive_results
    
    def _generate_analysis_summary(self, scenario_results: List[Dict]) -> Dict:
        """生成分析摘要"""
        
        summary = {
            'total_scenarios': len(scenario_results),
            'successful_scenarios': sum(1 for r in scenario_results if r.get('success', False)),
            'failed_scenarios': sum(1 for r in scenario_results if not r.get('success', False)),
            'total_result_files': sum(len(r.get('results_files', [])) for r in scenario_results)
        }
        
        return summary


def main():
    """主函数 - 运行区间回归分析"""
    print("🔄 区间回归分析主程序")
    print("=" * 50)
    
    try:
        # 创建分析器
        analyzer = IntervalRegressionAnalyzer()
        
        # 运行综合分析
        results = analyzer.run_comprehensive_analysis()
        
        # 显示摘要
        summary = results['summary']
        print(f"\n📊 分析摘要:")
        print(f"  总场景数: {summary['total_scenarios']}")
        print(f"  成功场景: {summary['successful_scenarios']}")
        print(f"  失败场景: {summary['failed_scenarios']}")
        print(f"  结果文件: {summary['total_result_files']}")
        
        if results['errors']:
            print(f"\n⚠️  错误:")
            for error in results['errors']:
                print(f"    {error}")
        
        print(f"\n🎉 区间回归分析完成!")
        
    except Exception as e:
        print(f"❌ 分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()