#!/usr/bin/env python3
"""
Lightweight wrapper to run or submit INLA single-compound tests.

Usage examples:
  # local run
  python tools/run_compound.py --compound-id 5 --compound-name "Abamectin" --disease C81-C96 --local

  # submit batch from yaml
  python tools/run_compound.py --jobs examples/compound_jobs.yaml --submit

This script calls existing `Code/INLA/run_single_compound_test.sh` for local runs
and `Code/INLA/submit_single_compound_test.sh` for submission.
"""
import argparse
import subprocess
import sys
import yaml
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_SCRIPT = ROOT / 'Code' / 'INLA' / 'run_single_compound_test.sh'
SUBMIT_SCRIPT = ROOT / 'Code' / 'INLA' / 'submit_single_compound_test.sh'


def run_local(params):
    cmd = [str(RUN_SCRIPT),
           str(params.get('compound_id', '')),
           params.get('compound_name', ''),
           params.get('disease', 'C81-C96'),
           params.get('measure_types', 'Weight,Density'),
           params.get('estimate_types', 'min,avg,max'),
           str(params.get('lag_years', 5)),
           params.get('model_types', 'M0,M1,M2,M3')]
    print('Running locally:', ' '.join(cmd))
    subprocess.check_call(cmd)


def submit_batch(params):
    # For simplicity, call the existing submit script which wraps the run script.
    # If cluster array needed, extend this to generate per-job sbatch with parameters.
    print('Submitting to SLURM using submit script...')
    subprocess.check_call([str(SUBMIT_SCRIPT)])


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--compound-id', type=str)
    p.add_argument('--compound-name', type=str, default='')
    p.add_argument('--disease', type=str, default='C81-C96')
    p.add_argument('--measure-types', type=str, default='Weight,Density')
    p.add_argument('--estimate-types', type=str, default='min,avg,max')
    p.add_argument('--lag-years', type=int, default=5)
    p.add_argument('--model-types', type=str, default='M0,M1,M2,M3')
    p.add_argument('--local', action='store_true', help='Run locally (no SLURM)')
    p.add_argument('--submit', action='store_true', help='Submit to SLURM using submit script')
    p.add_argument('--jobs', type=str, help='YAML file with multiple jobs')
    p.add_argument('--dry-run', action='store_true')
    return p.parse_args()


def load_jobs(yaml_path):
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data.get('jobs', [])


def main():
    args = parse_args()

    if args.jobs:
        jobs = load_jobs(args.jobs)
        print(f'Loaded {len(jobs)} jobs')
        for job in jobs:
            if args.dry_run:
                print('DRY RUN - job:', job)
                continue
            if args.local:
                run_local(job)
            elif args.submit:
                submit_batch(job)
            else:
                # default: local run for first job
                run_local(job)
        return

    params = {
        'compound_id': args.compound_id or '5',
        'compound_name': args.compound_name or '',
        'disease': args.disease,
        'measure_types': args.measure_types,
        'estimate_types': args.estimate_types,
        'lag_years': args.lag_years,
        'model_types': args.model_types,
    }

    if args.dry_run:
        print('DRY RUN params:', params)
        return

    if args.local:
        run_local(params)
    elif args.submit:
        submit_batch(params)
    else:
        # default: local run
        run_local(params)


if __name__ == '__main__':
    main()
