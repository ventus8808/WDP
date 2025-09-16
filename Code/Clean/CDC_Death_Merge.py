#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# 导入配置文件
import sys
from pathlib import Path

# 读取 YAML 配置（必须存在）
try:
    import yaml  # type: ignore
except Exception:
    print("ERROR: 需要 PyYAML。请安装: pip install pyyaml", file=sys.stderr)
    sys.exit(1)

sys.path.append(str(Path(__file__).parent.parent))

# ============ 手动指定疾病子目录（相对输入根目录，例："C00-C97"） ============
MANUAL_ICD_GROUP = "C00-C97"
# =====================================================================

# 解析配置，计算输入输出路径
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = PROJECT_ROOT / "config.yaml"

if not CONFIG_PATH.exists():
    print(f"ERROR: 未找到配置文件: {CONFIG_PATH}", file=sys.stderr)
    sys.exit(1)

with CONFIG_PATH.open("r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

ds = (cfg.get("data_sources") or {}).get("cdc_wonder") or {}
base_rel = ds.get("integrity_base_dir")
processed_rel = ds.get("processed")

if not base_rel or not processed_rel:
    print("ERROR: config.yaml 缺少 data_sources.cdc_wonder.integrity_base_dir 或 processed", file=sys.stderr)
    sys.exit(1)

_base_dir = (PROJECT_ROOT / base_rel).resolve()
_output_dir = (PROJECT_ROOT / processed_rel).resolve()

if not MANUAL_ICD_GROUP.strip():
    print("ERROR: 请在脚本顶部设置 MANUAL_ICD_GROUP，例如 'C00-C97'", file=sys.stderr)
    sys.exit(1)

_input_dir = (_base_dir / MANUAL_ICD_GROUP.strip()).resolve()

if not _base_dir.exists():
    print(f"ERROR: 输入根目录不存在: {_base_dir}", file=sys.stderr)
    sys.exit(1)

if not _input_dir.exists():
    print(f"ERROR: 疾病目录不存在: {_input_dir}", file=sys.stderr)
    sys.exit(1)

# ==================== 在这里修改要检查的文件夹路径 ====================
FOLDER_PATH = _input_dir
# =====================================================================

# ==================== 数据完整性检查开关 ====================
SKIP_INTEGRITY_CHECK = True  # 设为True跳过数据完整性检查
# =====================================================================

import pandas as pd


def check_data_integrity(folder_path):
    """调用现有的数据完整性检查函数（保持兼容，默认跳过）"""
    try:
        from CDC_Data_Integrity_Checker import check_data_integrity as original_check  # type: ignore
        return original_check(folder_path)
    except Exception:
        print("[WARN] 完整性检查脚本接口可能已变更，跳过调用")
        return None


def merge_cancer_data(folder_path):
    """合并癌症数据"""
    print(f"开始合并文件夹: {folder_path}")

    # 需要的列
    needed_columns = ['Year', 'County', 'County Code', 'Sex Code', 'Race',
                     'Ten-Year Age Groups', 'Deaths', 'Population', 'Crude Rate Standard Error']

    all_data = []

    # 处理所有CSV文件
    for file_path in sorted(folder_path.glob("*.csv")):
        filename = file_path.name
        print(f"处理文件: {filename}")

        try:
            # 尝试不同的编码格式
            encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252', 'windows-1252']
            df = None

            for encoding in encodings:
                try:
                    df = pd.read_csv(file_path, encoding=encoding)
                    print(f"  ✅ 使用编码: {encoding}")
                    break
                except UnicodeDecodeError:
                    continue

            if df is None:
                print(f"  ❌ 无法读取文件，尝试了所有编码格式")
                continue

            df = df.dropna(how='all')

            # 检查必需的列是否存在
            missing_cols = [col for col in needed_columns if col not in df.columns]
            if missing_cols:
                print(f"  ❌ 缺少必需列: {missing_cols}")
                continue

            # 只保留需要的列
            df = df[needed_columns].copy()

            # 移除Year列为空或非数字的行（这些通常是标题行或无效数据）
            df = df.dropna(subset=['Year'])
            df = df[df['Year'] != 'Year']  # 移除标题行

            # 确保县代码是5位数字符串
            df['County Code'] = df['County Code'].astype(str).str.replace('.0', '').str.zfill(5)

            # 确保年份是整数
            df['Year'] = pd.to_numeric(df['Year'], errors='coerce').fillna(0).astype(int)

            # 只保留1999-2019年的数据
            df = df[(df['Year'] >= 1999) & (df['Year'] <= 2019)].copy()

            # 确保Deaths和Population是整数
            df['Deaths'] = pd.to_numeric(df['Deaths'], errors='coerce').fillna(0).astype(int)
            df['Population'] = pd.to_numeric(df['Population'], errors='coerce').fillna(0).astype(int)

            all_data.append(df)
            print(f"  ✅ 成功处理，{len(df)} 行")

        except Exception as e:
            print(f"  ❌ 错误: {e}")
            continue

    if not all_data:
        print("没有成功处理任何文件")
        return None

    # 合并所有数据
    print("\n合并所有数据...")
    merged_df = pd.concat(all_data, ignore_index=True)

    # 重命名列
    merged_df = merged_df.rename(columns={
        'County Code': 'COUNTY_FIPS',
        'Year': 'Year',
        'County': 'County',
        'Sex Code': 'Sex',
        'Race': 'Race',
        'Ten-Year Age Groups': 'Age',
        'Deaths': 'Deaths',
        'Population': 'Population',
        'Crude Rate Standard Error': 'SD'
    })

    print(f"合并完成，总共 {len(merged_df)} 行")

    # 保存合并后的数据
    output_dir = _output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"{MANUAL_ICD_GROUP}.csv"
    merged_df.to_csv(output_file, index=False)
    print(f"数据已保存到: {output_file}")

    return merged_df


def main():
    """主函数"""
    if not SKIP_INTEGRITY_CHECK:
        print("进行数据完整性检查...")
        check_data_integrity(FOLDER_PATH)
        print("\n" + "="*50 + "\n")

    # 合并数据
    merged_data = merge_cancer_data(FOLDER_PATH)

    if merged_data is not None:
        print("\n数据合并成功!")
        print(f"数据形状: {merged_data.shape}")
        print(f"年份范围: {merged_data['Year'].min()} - {merged_data['Year'].max()}")
        print(f"县数量: {merged_data['COUNTY_FIPS'].nunique()}")
    else:
        print("\n数据合并失败!")


if __name__ == "__main__":
    main()
