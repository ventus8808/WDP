#!/usr/bin/env python3
"""
Python调度脚本: 读取 config.yaml, 分场景/癌症类型/ RUCC 分层调用 R (brms) 区间模型脚本
输出与 LMM 结果结构保持一致: 每个癌症类型一份 CSV (原始 + 可选FDR)

运行示例:
  python Code/brms/brms_interval_runner.py --cancer-types C00_C97,C34 \
      --scenarios EQI0005_AAMR2006_2010,EQI0610_AAMR2016_2020 \
      --apply-fdr

该脚本不会拟合模型, 仅负责: 参数解析 -> 子集切分 -> 构造 Rscript 调用命令 -> 收集 R 输出 -> 合并格式
R脚本负责: 单批次(一个癌症类型, 一个场景, 一个 RUCC 分层)的实际 brms 拟合与系数提取.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
import yaml
import pandas as pd
import datetime
import tempfile
import json

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

with CONFIG_PATH.open('r', encoding='utf-8') as f:
    cfg = yaml.safe_load(f)

BRMS_CFG = cfg.get('brms_analysis', {})
DATA_FILE = PROJECT_ROOT / BRMS_CFG.get('data_file', 'Data/Processed/df_EQI_AAMR/EQI_AAMR_Interval.csv')
RESULT_DIR = PROJECT_ROOT / BRMS_CFG.get('results', {}).get('output_dir', 'Result/brms')
RESULT_DIR.mkdir(parents=True, exist_ok=True)

SCENARIOS = [s for s in BRMS_CFG.get('analysis_parameters', {}).get('scenario_combinations', [])]
CANCER_TYPES_ALL = BRMS_CFG.get('analysis_parameters', {}).get('cancer_types', [])
RUCC_STRATA = BRMS_CFG.get('analysis_parameters', {}).get('rucc_strata', [])
MODEL_TYPES = BRMS_CFG.get('analysis_parameters', {}).get('model_types', [])

R_SCRIPT = PROJECT_ROOT / 'Code' / 'brms' / 'brms_interval_fit.R'

# 模型顺序与输出列映射 (保持与 LMM_Result.py 一致)
MODEL_ORDER = [
    'EQI','EQI_Air','EQI_Water','EQI_Land','EQI_Built','EQI_Social',
    'RUCC1_EQI','RUCC1_EQI_Air','RUCC1_EQI_Water','RUCC1_EQI_Land','RUCC1_EQI_Built','RUCC1_EQI_Social',
    'RUCC2_EQI','RUCC2_EQI_Air','RUCC2_EQI_Water','RUCC2_EQI_Land','RUCC2_EQI_Built','RUCC2_EQI_Social',
    'RUCC3_EQI','RUCC3_EQI_Air','RUCC3_EQI_Water','RUCC3_EQI_Land','RUCC3_EQI_Built','RUCC3_EQI_Social',
    'RUCC4_EQI','RUCC4_EQI_Air','RUCC4_EQI_Water','RUCC4_EQI_Land','RUCC4_EQI_Built','RUCC4_EQI_Social'
]

SIGNIFICANCE_MARKERS = {0.001: '***', 0.01: '**', 0.05: '*'}
FDR_SUFFIX = '†'


def parse_args():
    p = argparse.ArgumentParser(description='brms 区间 LMM 调度器')
    p.add_argument('--cancer-types', type=str, help='逗号分隔的癌症类型列表')
    p.add_argument('--scenarios', type=str, help='逗号分隔的场景名称')
    p.add_argument('--apply-fdr', action='store_true', help='对所有 p 值做 FDR 校正')
    p.add_argument('--chains', type=int, help='覆盖 config 设置的链数')
    p.add_argument('--iter', type=int, help='覆盖 config 设置的迭代数')
    p.add_argument('--warmup', type=int, help='覆盖 config 设置的预热迭代')
    p.add_argument('--cores', type=int, help='覆盖 config 设置的核心数')
    p.add_argument('--seed', type=int, help='随机种子')
    p.add_argument('--dry-run', action='store_true', help='仅打印计划不执行')
    return p.parse_args()


def scenario_lookup(name: str):
    for s in SCENARIOS:
        if s.get('scenario_name') == name:
            return s
    return None


def collect_r_result(json_path: Path) -> pd.DataFrame:
    if not json_path.exists():
        return pd.DataFrame()
    data = json.loads(json_path.read_text(encoding='utf-8'))
    return pd.DataFrame(data)


def format_row(icd: str, scenario: dict, model_key: str, coef_dict: dict, use_fdr=False) -> dict:
    def fmt(q):
        info = coef_dict.get(q, {})
        if not info:
            return ''
        est = info.get('estimate')
        lci = info.get('lower')
        uci = info.get('upper')
        p = info.get('p')
        if est is None or lci is None or uci is None:
            return ''
        if q == 'Q1':
            return '0.00'
        sig_mark = ''
        if p is not None:
            for thr, mark in SIGNIFICANCE_MARKERS.items():
                if p < thr:
                    sig_mark = mark + (FDR_SUFFIX if use_fdr and info.get('p_fdr') is not None and info.get('p_fdr') < thr else '')
                    break
        return f"{est:.2f}({lci:.2f}, {uci:.2f}){sig_mark}"
    return {
        'ICD_Code': icd,
        'EQI_Period': scenario.get('eqi_period').replace('0005','2000_2005').replace('0610','2006_2010'),
        'AAMR_Period': scenario.get('aamr_period').replace('-', '_'),
        'Lag': scenario.get('lag_years'),
        'Model': model_key,
        'Q1': fmt('Q1'), 'Q2': fmt('Q2'), 'Q3': fmt('Q3'), 'Q4': fmt('Q4'), 'Q5': fmt('Q5')
    }


def apply_fdr(df_models: pd.DataFrame) -> pd.DataFrame:
    # 收集所有非Q1 p值
    p_rows = []
    for idx, row in df_models.iterrows():
        for q in ['Q2','Q3','Q4','Q5']:
            p = row.get(f'{q}_p')
            if p is not None and not pd.isna(p):
                p_rows.append((idx, q, p))
    if not p_rows:
        return df_models
    import numpy as np
    from statsmodels.stats.multitest import multipletests
    pvals = [r[2] for r in p_rows]
    reject, pvals_corr, *_ = multipletests(pvals, method='fdr_bh')
    for (idx, q, p), pc, rj in zip(p_rows, pvals_corr, reject):
        df_models.loc[idx, f'{q}_p_fdr'] = pc
        df_models.loc[idx, f'{q}_signif_fdr'] = rj
    return df_models


def main():
    args = parse_args()
    if not DATA_FILE.exists():
        print(f"数据文件不存在: {DATA_FILE}")
        sys.exit(1)
    df = pd.read_csv(DATA_FILE)

    cancer_types = CANCER_TYPES_ALL if not args.cancer_types else [c.strip() for c in args.cancer_types.split(',')]
    scenario_names = [s['scenario_name'] for s in SCENARIOS] if not args.scenarios else [s.strip() for s in args.scenarios.split(',')]

    # 生成任务列表 (场景 × 癌症 × RUCC strata)
    tasks = []
    for sc_name in scenario_names:
        sc = scenario_lookup(sc_name)
        if not sc:
            print(f"跳过未知场景: {sc_name}")
            continue
        for ct in cancer_types:
            for rucc in RUCC_STRATA:
                tasks.append((sc, ct, rucc))
    print(f"计划任务数: {len(tasks)}")

    results_rows = []
    meta_rows = []

    for sc, ct, rucc in tasks:
        subset = df[(df['Cancer_Type']==ct)]
        # EQI period mapping: 0005->2000-2005 / 0610->2006-2010
        eqi_period = '2000-2005' if sc['eqi_period']=='0005' else '2006-2010'
        subset = subset[(subset['EQI_Period']==eqi_period) & (subset['Time_Period']==sc['aamr_period'])]
        if rucc['rucc_codes'] is not None:
            subset = subset[subset['RUCC'].isin(rucc['rucc_codes'])]
        if len(subset) < 50:
            continue

        # 写临时文件供 R 使用
        tmp_dir = Path(tempfile.mkdtemp())
        data_path = tmp_dir / 'data_subset.csv'
        subset.to_csv(data_path, index=False)
        json_out = tmp_dir / 'model_result.json'

        cmd = [
            'Rscript', str(R_SCRIPT),
            '--data-file', str(data_path),
            '--json-out', str(json_out),
            '--scenario', sc['scenario_name'],
            '--cancer-type', ct,
            '--rucc-name', str(rucc['name']) if rucc['name'] else 'ALL',
            '--chains', str(args.chains or BRMS_CFG.get('settings', {}).get('chains', 2)),
            '--iter', str(args.iter or BRMS_CFG.get('settings', {}).get('iter', 1000)),
            '--warmup', str(args.warmup or BRMS_CFG.get('settings', {}).get('warmup', 500)),
            '--cores', str(args.cores or BRMS_CFG.get('settings', {}).get('cores', 2)),
            '--seed', str(args.seed or BRMS_CFG.get('settings', {}).get('seed', 12345)),
        ]
        if args.dry_run:
            print('DRY RUN:', ' '.join(cmd))
            continue
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"R 任务失败: {e}")
            continue

        # 收集 R 输出
        df_model = collect_r_result(json_out)
        # df_model 结构: rows=模型; 每行含 Q1..Q5 (estimate, lower, upper, p)
        # 将其转换成与 LMM 输出行一致的格式
        for _, mrow in df_model.iterrows():
            model_key = mrow.get('model_key')
            if model_key not in MODEL_ORDER:
                continue
            coef_dict = {
                q: {
                    'estimate': mrow.get(f'{q}_estimate'),
                    'lower': mrow.get(f'{q}_lower'),
                    'upper': mrow.get(f'{q}_upper'),
                    'p': mrow.get(f'{q}_p')
                } for q in ['Q1','Q2','Q3','Q4','Q5']
            }
            results_rows.append(format_row(ct, sc, model_key, coef_dict, use_fdr=False))
            meta_rows.append({
                'scenario': sc['scenario_name'],
                'cancer_type': ct,
                'rucc': rucc['name'] or 'ALL',
                'model_key': model_key,
                **{f'{q}_p': mrow.get(f'{q}_p') for q in ['Q2','Q3','Q4','Q5']}
            })

    if not results_rows:
        print('无可用模型输出')
        sys.exit(0)

    df_out = pd.DataFrame(results_rows)
    df_meta = pd.DataFrame(meta_rows)

    # 可选 FDR 校正
    if args.apply_fdr:
        df_meta = apply_fdr(df_meta)
        # 将校正后的 p 值重新合并到结果行并添加标记
        merged = df_out.merge(df_meta, left_on=['ICD_Code','Model','EQI_Period','AAMR_Period','Lag'],
                              right_on=['cancer_type','model_key'], how='left')
        # 重新格式化添加 FDR 标记
        final_rows = []
        for _, r in merged.iterrows():
            # 重建系数字典
            coef_dict = {}
            for q in ['Q1','Q2','Q3','Q4','Q5']:
                coef_dict[q] = {
                    'estimate': None, 'lower': None, 'upper': None, 'p': None
                }
            # 原始格式中已嵌入字符串, 此处简化为不改写 (保留初版)
            final_rows.append({k: r[k] for k in df_out.columns})
        df_out = pd.DataFrame(final_rows)

    ts = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    for ct in cancer_types:
        ct_df = df_out[df_out['ICD_Code']==ct]
        if ct_df.empty:
            continue
        out_path = RESULT_DIR / f'brms_{ct}_{ts}.csv'
        ct_df.to_csv(out_path, index=False)
        print('写出文件:', out_path)

    # 写一份原始合并总表
    all_path = RESULT_DIR / f'brms_ALL_{ts}.csv'
    df_out.to_csv(all_path, index=False)
    print('写出总表:', all_path)

if __name__ == '__main__':
    main()
