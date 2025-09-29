#!/bin/bash
# 调试专用：极简最小测试，用于排查集群环境问题
#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_Debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --time=30:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -eo pipefail

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# 增强的项目根目录检测
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""

if [ -f "${SCRIPT_DIR}/../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd -P)"
elif command -v git >/dev/null 2>&1; then
  GIT_ROOT="$(cd "${SCRIPT_DIR}" && git rev-parse --show-toplevel 2>/dev/null || true)"
  if [ -n "$GIT_ROOT" ] && [ -f "$GIT_ROOT/config.yaml" ]; then
    PROJECT_ROOT="$GIT_ROOT"
  fi
fi

if [ -z "$PROJECT_ROOT" ] && [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
fi

if [ -z "$PROJECT_ROOT" ] && [ -f "$PWD/config.yaml" ]; then
  PROJECT_ROOT="$PWD"
fi

if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "无法确定项目根目录"
  log ERROR "脚本目录: $SCRIPT_DIR"
  log ERROR "SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR-未设置}"
  log ERROR "当前工作目录: $PWD"
  exit 1
fi

cd "$PROJECT_ROOT" || { log ERROR "无法切换到: $PROJECT_ROOT"; exit 1; }
log INFO "项目根目录: $PROJECT_ROOT"

# 环境检查
log INFO "=== 环境检查 ==="
log INFO "Python路径: $(which python || echo '未找到')"
log INFO "PyMC安装: $(python -c 'import pymc; print(pymc.__version__)' 2>/dev/null || echo '未安装')"
log INFO "NumPy版本: $(python -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo '未安装')"
log INFO "可用内存: $(free -h | grep '^Mem:' | awk '{print $7}' || echo '未知')"

# 激活conda
set +u || true
export ADDR2LINE="${ADDR2LINE-}"
if [ "${CONDA_DEFAULT_ENV-}" != "pymc" ]; then
  if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
  else
    log ERROR "找不到conda"; exit 1
  fi
  conda activate pymc || { log ERROR "激活pymc失败"; exit 1; }
fi

log INFO "激活后Python: $(which python)"
log INFO "PyMC版本: $(python -c 'import pymc; print(pymc.__version__)')"

# 数据结构检查
log INFO "=== 数据结构检查 ==="
export PYTHONUNBUFFERED=1

python -u -c "
import sys
sys.path.insert(0, 'Code/PYMC')
import pandas as pd
from pathlib import Path

print('=== 检查CDC数据文件结构 ===')
try:
    cdc_file = Path('Data/Processed/CDC/AAMR_C81-C96.csv')
    if cdc_file.exists():
        df = pd.read_csv(cdc_file, nrows=5)  # 只读前5行检查结构
        print(f'文件存在: {cdc_file}')
        print(f'数据形状: {df.shape}')
        print(f'列名: {list(df.columns)}')
        print('前几行数据:')
        print(df.head())
        
        # 检查是否有Deaths相关列
        deaths_cols = [col for col in df.columns if 'death' in col.lower()]
        print(f'Deaths相关列: {deaths_cols}')
        
        # 检查是否有Type相关列
        type_cols = [col for col in df.columns if 'type' in col.lower()]
        print(f'Type相关列: {type_cols}')
        
    else:
        print(f'❌ 文件不存在: {cdc_file}')
        
except Exception as e:
    print(f'❌ 数据文件检查失败: {e}')
    import traceback
    traceback.print_exc()

print('\\n=== 详细数据加载测试 ===')
try:
    from Utils_Data import WDPDataLoader
    
    loader = WDPDataLoader()
    print('✅ 数据加载器初始化成功')
    
    # 直接尝试加载CDC数据
    print('直接加载CDC数据...')
    cdc_df = loader.load_mortality_data('C81-C96')
    print(f'✅ CDC数据加载成功，形状: {cdc_df.shape}')
    print(f'列名: {list(cdc_df.columns)}')
    print(f'Deaths_Type值: {cdc_df[\"Deaths_Type\"].value_counts()}')
    
    # 加载其他数据
    print('\\n加载协变量数据...')
    cov_df = loader.load_covariate_data()
    print(f'✅ 协变量数据加载成功，形状: {cov_df.shape}')
    
    print('\\n加载农药数据...')
    pest_df = loader.load_pesticide_data()
    print(f'✅ 农药数据加载成功，形状: {pest_df.shape}')
    
    print('\\n加载空间邻接数据...')
    adj_matrix, county_fips = loader.load_spatial_adjacency()
    print(f'✅ 空间数据加载成功，县数: {len(county_fips)}')
    
    print('\\n=== 综合数据可用性检查 ===')
    from Utils_Others import check_data_availability
    report = check_data_availability('C81-C96', '2')
    print(f'综合检查结果: {report[\"data_available\"]}')
    if not report['data_available']:
        for issue in report['issues']:
            print(f'  问题: {issue}')
    else:
        print('✅ 所有数据检查通过')
    
except Exception as e:
    print(f'❌ 数据加载失败: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
"

# 简单模型测试（仅M0，极少样本）
log INFO "=== 极简模型测试 ==="
python -u Code/PYMC/main.py \
  --disease "C81-C96" \
  --compound "2" \
  --model "M0" \
  --lag "5" \
  --measure "Weight" \
  --estimate "avg" \
  --sampling-mode "test" \
  --draws 100 --tune 50 --chains 1 --cores 1 --target-accept 0.8 \
  --config-path "config.yaml" --verbose \
  2>&1 | tee "debug_${SLURM_JOB_ID:-$$}.log"

RESULT=$?
if [ $RESULT -eq 0 ]; then
  log INFO "✅ 调试测试成功"
else
  log ERROR "❌ 调试测试失败，退出码: $RESULT"
  log ERROR "详细日志: debug_${SLURM_JOB_ID:-$$}.log"
  exit $RESULT
fi