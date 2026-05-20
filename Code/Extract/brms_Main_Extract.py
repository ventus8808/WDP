#!/usr/bin/env python3
"""Extract brms_Main results (EQI 2000-2005 only) with full disease names from config."""

import yaml
import pandas as pd
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = PROJECT_ROOT / "Result" / "brms_Main"
OUTPUT_BASE = PROJECT_ROOT / "Result" / "Extraction"
CONFIG_FILE = PROJECT_ROOT / "config.yaml"

MODEL_RENAME = {
    "EQI+SM+PA+OB": "EQI",
    "EQI_Air+SM+PA+OB": "Air",
    "EQI_Water+SM+PA+OB": "Water",
    "EQI_Land+SM+PA+OB": "Land",
    "EQI_Built+SM+PA+OB": "Built",
    "EQI_Social+SM+PA+OB": "Social",
}


def collect_entries(config: dict):
    """Yield (group_shortname, icd_code, fullname, shortname) for every disease entry."""
    for group_key, group in config.get("diseases", {}).items():
        overall = group.get("overall", {})
        group_shortname = overall.get("shortname", group_key)

        if "icd_code" in overall:
            yield group_shortname, overall["icd_code"], overall.get("fullname", ""), overall.get("shortname", "")

        for entry in group.get("subtypes", {}).values():
            if "icd_code" in entry:
                yield group_shortname, entry["icd_code"], entry.get("fullname", ""), entry.get("shortname", "")


def main():
    with open(CONFIG_FILE) as f:
        config = yaml.safe_load(f)

    groups: dict[str, list] = defaultdict(list)
    for group_shortname, icd_code, fullname, shortname in collect_entries(config):
        groups[group_shortname].append((icd_code, fullname, shortname))

    for group_shortname, entries in groups.items():
        all_dfs = []

        for icd_code, fullname, shortname in entries:
            csv_path = INPUT_DIR / f"{icd_code}_Main.csv"
            if not csv_path.exists():
                continue

            df = pd.read_csv(csv_path)
            df = df[df["EQI_Period"] == "2000_2005"].copy()
            if df.empty:
                continue

            df["Outcome"] = shortname if shortname else fullname
            df["Model"] = df["Model"].map(MODEL_RENAME).fillna(df["Model"])
            df = df[["ICD_Code", "Outcome", "Lag", "Model", "Q1", "Q2", "Q3", "Q4", "Q5"]]
            all_dfs.append(df)

        if not all_dfs:
            continue

        result = pd.concat(all_dfs, ignore_index=True)
        OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
        out_file = OUTPUT_BASE / f"{group_shortname}_Main.csv"
        result.to_csv(out_file, index=False)
        print(f"[OK] {out_file.relative_to(PROJECT_ROOT)}  ({len(result)} rows)")


if __name__ == "__main__":
    main()
