#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MICE插补进度监控脚本
"""

import os
import time
from pathlib import Path

def monitor_mice_progress():
    """监控MICE插补进度"""
    project_root = Path("/Users/ventus/Repository/WDP")
    mi_datasets_dir = project_root / "Data/df/MI_Datasets"
    diagnostics_dir = project_root / "Result/EQI_LMM_MI/MI_Diagnostics"
    
    print("🔍 MICE+PMM 插补进度监控")
    print("=" * 50)
    
    start_time = time.time()
    last_count = 0
    
    while True:
        try:
            # 检查插补数据集
            if mi_datasets_dir.exists():
                csv_files = list(mi_datasets_dir.glob("MI_dataset_*.csv"))
                current_count = len(csv_files)
                
                if current_count != last_count:
                    elapsed = time.time() - start_time
                    print(f"⏰ {time.strftime('%H:%M:%S')} - 已完成: {current_count}/20 个数据集 (用时: {elapsed:.1f}秒)")
                    
                    if current_count > 0:
                        avg_time = elapsed / current_count
                        remaining = (20 - current_count) * avg_time
                        print(f"📊 平均每个数据集: {avg_time:.1f}秒, 预计还需: {remaining:.0f}秒")
                    
                    last_count = current_count
                
                if current_count >= 20:
                    print("✅ 所有插补数据集已完成!")
                    break
            
            # 检查诊断文件
            if diagnostics_dir.exists():
                diag_files = list(diagnostics_dir.glob("*.png"))
                summary_file = diagnostics_dir / "imputation_summary.json"
                
                if summary_file.exists():
                    print(f"📈 诊断文件: {len(diag_files)} 个图表")
                    print("🎉 MICE+PMM 插补流程完全完成!")
                    break
            
            time.sleep(10)  # 每10秒检查一次
            
        except KeyboardInterrupt:
            print("\n⏸️ 监控中断")
            break
        except Exception as e:
            print(f"❌ 监控错误: {e}")
            time.sleep(5)

if __name__ == "__main__":
    monitor_mice_progress()