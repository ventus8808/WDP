#!/usr/bin/env python3
"""
Simple orchestrator for brms scenarios defined in config.yaml.
- Runs 01_prepare_data.R once (optional, safe to re-run).
- Iterates active scenarios and calls 02_run_brms_model.R with scenario name only.

职责反转后的设计原则:
- Python只做调度：读取配置，遍历场景，调用R脚本
- Python不再关心数据内容或生成复杂参数组合
- 每个场景的完整定义都在config.yaml中
- R脚本成为智能引擎，根据场景名称自我配置

用法:
python Code/brms/03_main_runner.py [--cancer-types C00_C97,C34,C50] [--verbose]
"""
from __future__ import annotations
import subprocess
import sys
import argparse
from pathlib import Path
import yaml

try:
    from tqdm import tqdm
except ImportError:
    # Fallback no-op progress wrapper
    def tqdm(iterable, total=None, desc=None, ncols=None):
        return iterable

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'config.yaml'
CODE_DIR = ROOT / 'Code' / 'brms'


def run_command(cmd: list[str], description: str = "", verbose: bool = False) -> bool:
    """Run a command and return success status."""
    try:
        if verbose or description:
            print(f"\n[brms] Running: {' '.join(cmd)}")
            if description:
                print(f"[brms] {description}")
        
        # Using pipe to show R output in real-time (only in verbose mode)
        if verbose:
            with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1) as proc:
                for line in proc.stdout:
                    print(line, end='')
        else:
            # Silent mode - capture output but don't show it unless there's an error
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"\n[brms] ERROR in {cmd[0]}: {result.stderr}")
                return False
            return True
        
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd)
            
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n[brms] ERROR: Command failed with return code {e.returncode}: {' '.join(cmd)}")
        return False
    except FileNotFoundError:
        print(f"\n[brms] ERROR: Command not found. Is Rscript in your PATH?")
        return False

def load_config() -> dict:
    """Load configuration."""
    with open(CONFIG, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def generate_scenarios(config, selected_cancer_types=None):
    """根据分析参数自动生成所有场景组合 - 与LMM系统完全对齐"""
    brms_config = config.get('brms_analysis', {})
    params = brms_config.get('analysis_parameters', {})
    
    # 如果指定了特定癌症类型，使用指定的；否则使用配置中的全部
    all_cancer_types = params.get('cancer_types', ['C00_C97'])
    if selected_cancer_types:
        # 验证指定的癌症类型是否在配置中
        valid_types = []
        for cancer in selected_cancer_types:
            if cancer in all_cancer_types:
                valid_types.append(cancer)
            else:
                print(f"[brms] 警告: 癌症类型 {cancer} 不在配置的癌症类型列表中")
        cancer_types = valid_types
        if not cancer_types:
            print(f"[brms] 错误: 没有有效的癌症类型")
            return []
    else:
        cancer_types = all_cancer_types
    
    scenario_combos = params.get('scenario_combinations', [])
    model_types = params.get('model_types', [{'name': 'EQI', 'formula_type': 'total_eqi'}])
    rucc_strata = params.get('rucc_strata', [{'name': None, 'rucc_codes': None}])
    
    scenarios = []
    
    print(f"[brms] 生成场景 (正确的模型结构):")
    print(f"[brms]   - 癌症类型: {cancer_types}")
    print(f"[brms]   - 时期组合: {len(scenario_combos)}种")
    print(f"[brms]   - 模型类型: {len(model_types)}种 (总EQI + EQI细分联合)") 
    print(f"[brms]   - RUCC分层: {len(rucc_strata)}种")
    
    total_scenarios = len(cancer_types) * len(scenario_combos) * len(model_types) * len(rucc_strata)
    print(f"[brms]   - 预计总场景: {total_scenarios}")
    
    # 新逻辑：每个癌症类型+时期组合+RUCC分层下运行2种模型
    for cancer in cancer_types:
        for combo in scenario_combos:
            for rucc in rucc_strata:
                for model_type in model_types:
                    # 生成场景名称
                    name_parts = [cancer.replace('_', '')]
                    
                    # 添加RUCC分层信息（如果有）
                    if rucc['name']:
                        name_parts.append(rucc['name'])
                    
                    # 添加模型类型
                    name_parts.append(model_type['name'])
                    
                    # 添加时期和滞后信息
                    name_parts.append(f"Lag{combo['lag_years']}")
                    
                    # 如果不是默认时期组合，添加时期信息
                    if combo['eqi_period'] != '0005' or combo['aamr_period'] != '2006-2010':
                        eqi_formatted = combo['eqi_period']
                        aamr_formatted = combo['aamr_period'].replace('-', '')
                        name_parts.append(f"{eqi_formatted}_{aamr_formatted}")
                    
                    scenario_name = '_'.join(name_parts)
                    
                    scenario = {
                        'name': scenario_name,
                        'cancer_type': cancer,
                        'eqi_period': combo['eqi_period'],
                        'aamr_period': combo['aamr_period'],
                        'lag_years': combo['lag_years'],
                        'rucc_filter': rucc['rucc_codes'][0] if rucc['rucc_codes'] else None,
                        'model_type': model_type['name'],
                        'formula_type': model_type['formula_type']
                    }
                    scenarios.append(scenario)
    
    print(f"[brms] ✅ 实际生成了 {len(scenarios)} 个分析场景")
    return scenarios

def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='BRMS贝叶斯分析流水线 - 与LMM系统对齐',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  python Code/brms/03_main_runner.py                           # 运行所有癌症类型
  python Code/brms/03_main_runner.py --cancer-types C00_C97   # 只运行全癌症
  python Code/brms/03_main_runner.py --cancer-types C34,C50   # 运行肺癌和乳腺癌
  python Code/brms/03_main_runner.py --verbose                # 详细输出模式
  python Code/brms/03_main_runner.py --skip-prep              # 跳过数据准备步骤
        """
    )
    
    parser.add_argument(
        '--cancer-types',
        type=str,
        help='指定要分析的癌症类型，用逗号分隔 (例如: C00_C97,C34,C50)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='详细输出模式，显示更多调试信息'
    )
    
    parser.add_argument(
        '--skip-prep',
        action='store_true', 
        help='跳过数据准备步骤（如果数据已经准备好）'
    )
    
    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()
    cfg = load_config()
    
    # 解析癌症类型参数
    selected_cancer_types = None
    if args.cancer_types:
        selected_cancer_types = [t.strip() for t in args.cancer_types.split(',')]
        print(f"[brms] 🎯 指定癌症类型: {selected_cancer_types}")
    
    # 1. Prepare data (optional, safe to re-run)
    if not args.skip_prep:
        prepare_script = CODE_DIR / '01_prepare_data.R'
        if prepare_script.exists():
            print("[brms] --- Step 1: Preparing Data ---")
            if not run_command(['Rscript', str(prepare_script)], "Running data preparation script", True):
                print("[brms] FATAL: Data preparation failed. Aborting.")
                sys.exit(1)
    else:
        print("[brms] --- Step 1: Skipped Data Preparation (--skip-prep) ---")
            
    # 2. Find and execute active scenarios
    print("\n[brms] --- Step 2: Running Model Scenarios ---")
    active_scenarios = generate_scenarios(cfg, selected_cancer_types)
    
    if not active_scenarios:
        print("[brms] 没有生成任何分析场景，退出")
        return
    
    if args.verbose:
        print(f"[brms] 📋 生成的分析场景详情:")
        for i, scenario in enumerate(active_scenarios, 1):
            print(f"[brms]   {i:3d}. {scenario['name']} (Cancer: {scenario['cancer_type']}, Domain: {scenario['domain']})")
    
    failures = 0
    successes = 0
    
    # 使用tqdm显示进度条，并估算时间
    import time
    start_time = time.time()
    
    progress_bar = tqdm(
        active_scenarios, 
        desc="🧬 BRMS Analysis", 
        ncols=100,
        bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]"
    )
    
    for i, scenario in enumerate(progress_bar):
        scenario_name = scenario['name'] if isinstance(scenario, dict) else scenario
        
        # 更新进度条描述，添加时间估算
        short_name = scenario_name[:25] + "..." if len(scenario_name) > 25 else scenario_name
        
        # 计算平均时间和预估剩余时间
        if i > 0:
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            remaining_scenarios = len(active_scenarios) - i
            est_remaining = avg_time * remaining_scenarios
            est_remaining_str = f"{int(est_remaining//60):02d}:{int(est_remaining%60):02d}"
            progress_bar.set_description(f"🧬 {short_name} (~{est_remaining_str})")
        else:
            progress_bar.set_description(f"🧬 {short_name}")
        
        cmd = ['Rscript', str(CODE_DIR / '02_run_brms_model.R'), '--scenario', scenario_name]
        
        success = run_command(cmd, f"Scenario: {scenario_name}" if args.verbose else "", args.verbose)
        if success:
            successes += 1
        else:
            failures += 1
            print(f"\n[brms] ❌ ERROR: Scenario failed: {scenario_name}")
    
    # 3. Results are already saved in LMM-compatible format
    print("\n[brms] --- Step 3: Results Summary ---")
    print("[brms] All results have been saved in LMM-compatible format to Result/brms/")
    
    # 4. Final Report
    print("\n[brms] --- Analysis Complete ---")
    total_scenarios = len(active_scenarios)
    if failures == 0:
        print(f"[brms] 🎉 All {total_scenarios} scenarios completed successfully!")
    else:
        print(f"[brms] ⚠️  Completed with {failures} failures and {successes} successes out of {total_scenarios} scenarios")
        print(f"[brms] 📊 Success rate: {successes/total_scenarios*100:.1f}%")


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f"\n[brms] An unexpected error occurred in the main runner: {e}")
        sys.exit(1)
