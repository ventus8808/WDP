#!/bin/bash

# PyMC Test Shell Script
# Usage: ./pymc_test.sh --disease <disease> --compound <compound> --model <model> --lag <lag> [--sampling-mode <mode>]
# Example: ./pymc_test.sh --disease C81-C96 --compound 24D --model M0 --lag 5 --sampling-mode test

# ========================
# WONDER PyMC分析 单任务CPU测试脚本 (本地环境适配版)
# ========================

# ========================
# 项目目录设置
# ========================

# 项目根目录
PROJECT_ROOT="/Users/ventus/Repository/WDP"
PYMC_DIR="$PROJECT_ROOT/Code/PYMC"
echo "保持工作目录在项目根目录: ${PROJECT_ROOT}"
cd ${PROJECT_ROOT}

# 检查PyMC目录是否存在
if [ ! -d "$PYMC_DIR" ]; then
    echo "❌ Error: PyMC directory not found at $PYMC_DIR"
    exit 1
fi

# 检查main.py是否存在
if [ ! -f "$PYMC_DIR/main.py" ]; then
    echo "❌ Error: main.py not found at $PYMC_DIR/main.py"
    exit 1
fi

# 检查config.yaml是否存在
if [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
    echo "❌ Error: config.yaml not found at $PROJECT_ROOT/config.yaml"
    exit 1
fi

# 设置Python解释器路径（优先使用conda，否则使用默认）
PYTHON_CMD="python"

# 尝试检测并激活conda环境
if [ -f "~/miniconda3/etc/profile.d/conda.sh" ]; then
    source ~/miniconda3/etc/profile.d/conda.sh
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source /opt/anaconda3/etc/profile.d/conda.sh
elif [ -f "/Users/ventus/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source /Users/ventus/opt/anaconda3/etc/profile.d/conda.sh
fi

# 尝试激活pymc环境
if conda info --envs | grep -q "pymc"; then
    echo "激活Conda环境..."
    conda activate pymc
    echo "Conda环境 'pymc' 已激活。"
else
    echo "⚠️  Conda环境 'pymc' 未找到，使用默认Python环境"
fi

echo "Python路径: $(which python)"
echo "[诊断] Python环境与PyMC包可用性："
python -c "import pymc; print('PyMC available:', pymc.__version__)" || { echo "❌ PyMC导入失败"; exit 1; }
python -c "import arviz; print('ArviZ available:', arviz.__version__)" || { echo "❌ ArviZ导入失败"; exit 1; }

# ========================
# 运行PyMC分析命令
# ========================

# Parse command line arguments
DISEASE=""
COMPOUND=""
MODEL=""
LAG=""
SAMPLING_MODE="production"
DRY_RUN=false
MEASURE="Weight"
OUTPUT_DIR=""
CONFIG_PATH=""
VERBOSE=false
CORES=""
CHAINS=""
DRAWS=""
TUNE=""
TARGET_ACCEPT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --disease)
            DISEASE="$2"
            shift 2
            ;;
        --compound)
            COMPOUND="$2"
            shift 2
            ;;
        --model)
            MODEL="$2"
            shift 2
            ;;
        --lag)
            LAG="$2"
            shift 2
            ;;
        --sampling-mode)
            SAMPLING_MODE="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --measure)
            MEASURE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --config-path)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --verbose|-v)
            VERBOSE=true
            shift
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --chains)
            CHAINS="$2"
            shift 2
            ;;
        --draws)
            DRAWS="$2"
            shift 2
            ;;
        --tune)
            TUNE="$2"
            shift 2
            ;;
        --target-accept)
            TARGET_ACCEPT="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 --disease <disease> --compound <compound> --model <model> --lag <lag> [options]"
            echo ""
            echo "Required arguments:"
            echo "  --disease      Disease/ICD group (e.g., C81-C96)"
            echo "  --compound     Pesticide compound (e.g., 24D)"
            echo "  --model        Model type (M0, M1, M2, M3)"
            echo "  --lag          Lag period in years (e.g., 5)"
            echo ""
            echo "Optional arguments:"
            echo "  --sampling-mode  Sampling mode: production, test (default: production)"
            echo "  --measure        Measure type: Weight, Density (default: Weight)"
            echo "  --output-dir     Output directory (default: Result/PyMC_Results)"
            echo "  --config-path    Path to config.yaml (default: project root config.yaml)"
            echo "  --verbose, -v    Verbose output"
            echo "  --cores          Number of CPU cores or 'auto' (default: auto)"
            echo "  --chains         Number of chains or 'auto' (default: auto=cores)"
            echo "  --draws          Posterior draws per chain"
            echo "  --tune           Tuning steps per chain"
            echo "  --target-accept  Target accept rate, e.g., 0.9"
            echo "  --dry-run       Data validation only, no model fitting"
            echo "  --help, -h      Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 --disease C81-C96 --compound 24D --model M0 --lag 5"
            echo "  $0 --disease C81-C96 --compound 24D --model M0 --lag 5 --dry-run"
            echo "  $0 --disease C81-C96 --compound 24D --model M0 --lag 5 --sampling-mode test --measure Density"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$DISEASE" || -z "$COMPOUND" || -z "$MODEL" || -z "$LAG" ]]; then
    echo "❌ Error: Missing required arguments"
    echo "Required: --disease, --compound, --model, --lag"
    echo "Use --help for usage information"
    exit 1
fi

# Print configuration
echo "---"
echo "分析参数："
echo "  Disease/ICD: $DISEASE"
echo "  Compound: $COMPOUND"
echo "  Model: $MODEL"
echo "  Lag: $LAG years"
echo "  Sampling Mode: $SAMPLING_MODE"
echo "  Dry Run: $DRY_RUN"
echo "  Measure: $MEASURE"
if [ -n "$OUTPUT_DIR" ]; then echo "  Output Dir: $OUTPUT_DIR"; fi
if [ -n "$CONFIG_PATH" ]; then echo "  Config Path: $CONFIG_PATH"; fi
echo "  Verbose: $VERBOSE"
if [ -n "$CORES" ]; then echo "  Cores: $CORES"; fi
if [ -n "$CHAINS" ]; then echo "  Chains: $CHAINS"; fi
if [ -n "$DRAWS" ]; then echo "  Draws: $DRAWS"; fi
if [ -n "$TUNE" ]; then echo "  Tune: $TUNE"; fi
if [ -n "$TARGET_ACCEPT" ]; then echo "  Target Accept: $TARGET_ACCEPT"; fi
echo ""

# --- 配置文件路径 ---
CONFIG_PATH="${PROJECT_ROOT}/config.yaml"
echo "使用配置文件: $(realpath ${CONFIG_PATH})"

echo "开始执行PyMC分析脚本..."
echo "-------------------------------------"

# Print configuration
echo "---"
echo "分析参数："
echo "  Disease/ICD: $DISEASE"
echo "  Compound: $COMPOUND"
echo "  Model: $MODEL"
echo "  Lag: $LAG years"
echo "  Sampling Mode: $SAMPLING_MODE"
echo "  Dry Run: $DRY_RUN"
echo "  Measure: $MEASURE"
if [ -n "$OUTPUT_DIR" ]; then echo "  Output Dir: $OUTPUT_DIR"; fi
if [ -n "$CONFIG_PATH" ]; then echo "  Config Path: $CONFIG_PATH"; fi
echo "  Verbose: $VERBOSE"
if [ -n "$CORES" ]; then echo "  Cores: $CORES"; fi
if [ -n "$CHAINS" ]; then echo "  Chains: $CHAINS"; fi
if [ -n "$DRAWS" ]; then echo "  Draws: $DRAWS"; fi
if [ -n "$TUNE" ]; then echo "  Tune: $TUNE"; fi
if [ -n "$TARGET_ACCEPT" ]; then echo "  Target Accept: $TARGET_ACCEPT"; fi
echo ""

# Build Python command
PYTHON_CMD="python ${PYMC_DIR}/main.py \
  --disease \"${DISEASE}\" \
  --compound \"${COMPOUND}\" \
  --model \"${MODEL}\" \
  --lag \"${LAG}\" \
  --measure \"${MEASURE}\" \
  --sampling-mode \"${SAMPLING_MODE}\" \
  --config-path \"${CONFIG_PATH}\" \
  --verbose"

if [ "$DRY_RUN" = true ]; then
    PYTHON_CMD="$PYTHON_CMD --dry-run"
fi

if [ -n "$OUTPUT_DIR" ]; then
    PYTHON_CMD="$PYTHON_CMD --output-dir \"${OUTPUT_DIR}\""
fi

if [ -n "$CORES" ]; then
    PYTHON_CMD="$PYTHON_CMD --cores \"${CORES}\""
fi

if [ -n "$CHAINS" ]; then
    PYTHON_CMD="$PYTHON_CMD --chains \"${CHAINS}\""
fi

if [ -n "$DRAWS" ]; then
    PYTHON_CMD="$PYTHON_CMD --draws \"${DRAWS}\""
fi

if [ -n "$TUNE" ]; then
    PYTHON_CMD="$PYTHON_CMD --tune \"${TUNE}\""
fi

if [ -n "$TARGET_ACCEPT" ]; then
    PYTHON_CMD="$PYTHON_CMD --target-accept \"${TARGET_ACCEPT}\""
fi

echo "执行命令: $PYTHON_CMD"
echo "-------------------------------------"

# Execute Python script
eval $PYTHON_CMD

status=$?
if [ $status -ne 0 ]; then
    echo "！！！Python脚本执行失败，请检查错误日志"
else
    echo "--- Python脚本执行成功 ---"
    echo "请检查项目根目录下的日志文件"
fi
