#!/usr/bin/env python3
"""
Simple orchestrator for brms scenarios defined in config.yaml.
- Ensures output directories exist.
- Runs 01_prepare_data.R once (optional, safe to re-run).
- Iterates active scenarios and calls 02_run_brms_model.R with scenario name only.

职责反转后的设计原则:
- Python只做调度：读取配置，遍历场景，调用R脚本
- Python不再关心数据内容或生成复杂参数组合
- 每个场景的完整定义都在config.yaml中
- R脚本成为智能引擎，根据场景名称自我配置
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path
import yaml

try:
    from tqdm import tqdm
except ImportError:
    # Fallback no-op progress wrapper
    def tqdm(iterable, total=None, desc=None):
        return iterable

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / 'config.yaml'
CODE_DIR = ROOT / 'Code' / 'brms'


def run_command(cmd: list[str], description: str = "") -> bool:
    """Run a command and return success status."""
    try:
        print(f"[brms] Running: {' '.join(cmd)}")
        if description:
            print(f"[brms] {description}")
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[brms] Command failed with return code {e.returncode}: {' '.join(cmd)}")
        return False


def load_config() -> dict:
    """Load configuration from config.yaml."""
    with open(CONFIG, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def ensure_directories(cfg: dict) -> None:
    """Create output directories if they don't exist."""
    brms_config = cfg.get('brms_analysis', {})
    results = brms_config.get('results', {})
    
    for dir_key, dir_path in results.items():
        full_path = ROOT / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"[brms] Ensured directory: {full_path}")


def get_active_scenarios(cfg: dict) -> list[dict]:
    """Extract active scenarios from configuration."""
    brms_config = cfg.get('brms_analysis', {})
    scenarios = brms_config.get('scenarios', [])
    
    active_scenarios = [s for s in scenarios if s.get('active', False)]
    
    print(f"[brms] Found {len(active_scenarios)} active scenarios out of {len(scenarios)} total")
    for scenario in active_scenarios:
        print(f"[brms]   - {scenario['name']}")
    
    return active_scenarios


def main():
    """
    Main orchestration function - simplified to pure scheduling.
    No data reading, no complex parameter generation - just config-driven scenario execution.
    """
    print("[brms] Starting brms analysis orchestrator (refactored for simplicity)")
    
    # Load configuration
    cfg = load_config()
    
    # Ensure output directories exist
    ensure_directories(cfg)
    
    # Optional: Run data preparation script once
    prep_script = CODE_DIR / '01_prepare_data.R'
    if prep_script.exists():
        print(f"[brms] Running preparation script: {prep_script}")
        if not run_command(['Rscript', str(prep_script)], "Data preparation"):
            print("[brms] Warning: Preparation script failed, but continuing...")
    
    # Get active scenarios
    active_scenarios = get_active_scenarios(cfg)
    
    if not active_scenarios:
        print("[brms] No active scenarios found. Nothing to do.")
        return
    
    # Execute scenarios
    print(f"[brms] Executing {len(active_scenarios)} active scenarios...")
    
    failures = 0
    for scenario in tqdm(active_scenarios, desc="brms scenarios"):
        scenario_name = scenario['name']
        
        # Simple command: just pass the scenario name
        cmd = ['Rscript', str(CODE_DIR / '02_run_brms_model.R'), '--scenario', scenario_name]
        
        success = run_command(cmd, f"Scenario: {scenario_name}")
        if not success:
            failures += 1
            print(f"[brms] Scenario failed: {scenario_name}")
    
    # Report results
    if failures == 0:
        print(f"[brms] All {len(active_scenarios)} scenarios completed successfully!")
    else:
        print(f"[brms] Completed with {failures} failures out of {len(active_scenarios)} scenarios")
    
    # Optional: post-process results
    post_script = CODE_DIR / '04_process_results.R'
    if post_script.exists():
        print(f"[brms] Running post-processing script: {post_script}")
        if not run_command(['Rscript', str(post_script)], "Post-processing results"):
            print("[brms] Warning: Post-processing script failed")


if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        sys.exit(e.returncode)
