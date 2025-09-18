#!/bin/bash
#
# 自动化脚本：创建并配置用于 WDP BYM INLA 生产分析的 Conda 环境
#
# 功能:
# 1. 创建一个名为 'INLA' 的新 Conda 环境。
# 2. 指定安装 R v4.3，以确保与 INLA 的良好兼容性。
# 3. 使用 Conda 安装所有已知的 R 依赖包 (dplyr, readr, yaml, etc.)。
# 4. 在 R 环境内部，使用 remotes 包安装特定版本的 INLA。
# 5. 主动添加备用仓库，以解决 fmesher 等依赖可能出现的二进制包下载问题。
#

# --- 配置 ---
# R 版本 (选择 4.3 是因为它有明确兼容的 INLA 版本)
R_VERSION="4.3"
# INLA 版本 (24.03.29 是与 R 4.3 兼容的版本)
INLA_VERSION="24.03.29"
# Conda 环境名称
ENV_NAME="INLA"

# --- 脚本开始 ---
set -e # 如果任何命令失败，立即退出脚本

echo "### 步骤 1/4: 创建 Conda 环境 '$ENV_NAME' 并安装 R 及其依赖包 ###"

# 从 conda-forge 频道创建环境并安装所有依赖
# 这包括 R 本体，您脚本所需的包，以及 remotes 和 Bioconductor 的基础包
conda create --name "$ENV_NAME" -c conda-forge -y \
  r-base="$R_VERSION" \
  r-dplyr \
  r-readr \
  r-yaml \
  r-argparse \
  r-progress \
  r-remotes \
  bioconductor-graph \
  bioconductor-rgraphviz

echo "✓ Conda 环境创建成功。"
echo

# --- 步骤 2: 激活环境 ---
echo "### 步骤 2/4: 激活 Conda 环境 '$ENV_NAME' ###"
# 为了在脚本中可靠地激活 conda，需要 source hook
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# 验证 R 版本
echo "✓ 环境已激活。当前 R 版本信息:"
R --version
echo

# --- 步骤 3: 在 R 中安装 INLA ---
echo "### 步骤 3/4: 使用 remotes 在 R 中安装 INLA v$INLA_VERSION ###"
echo "这将需要一些时间来下载和编译，请耐心等待..."

# 使用 Rscript 执行 R 命令
# - `options(repos=...)`: 我们主动添加了 inlabru-org 仓库，这是为了防止 fmesher 包在 CRAN 上找不到
#   适用于您服务器系统的二进制版本而导致安装失败，这是一个常见的坑。
# - `remotes::install_version`: 精确安装指定版本的 INLA，并处理其所有 R 依赖 (`dep=TRUE`)。
Rscript -e '
  # 设置仓库，优先使用 r-universe 来获取 fmesher，并添加 INLA 官方仓库
  options(repos = c(
    inlabruorg = "https://inlabru-org.r-universe.dev",
    INLA = "https://inla.r-inla-download.org/R/testing",
    CRAN = "https://cran.rstudio.com"
  ))

  # 安装指定版本的 INLA
  remotes::install_version(
    "INLA",
    version = "'"$INLA_VERSION"'",
    repos = getOption("repos"),
    dep = TRUE,
    upgrade = "never"
  )
'

echo "✓ INLA 安装命令已执行。"
echo

# --- 步骤 4: 最终验证 ---
echo "### 步骤 4/4: 验证所有核心包是否成功安装 ###"

Rscript -e '
  packages <- c("INLA", "fmesher", "dplyr", "readr", "yaml", "argparse", "progress")
  cat("正在检查包...\n")
  installed <- sapply(packages, requireNamespace, quietly = TRUE)

  if (all(installed)) {
    cat("\n✅✅✅ 恭喜！所有包都已成功安装！\n")
    cat("INLA 版本: ", as.character(packageVersion("INLA")), "\n")
    cat("fmesher 版本: ", as.character(packageVersion("fmesher")), "\n")
  } else {
    cat("\n❌❌❌ 安装失败！以下包未能加载:\n")
    cat(paste(names(installed[!installed]), collapse = "\n"), "\n")
    quit(status = 1)
  }
'

echo
echo "🚀 环境 '$ENV_NAME' 已准备就绪！"
echo "现在，您可以随时使用 'conda activate $ENV_NAME' 命令进入此环境，"
echo "然后运行您的 R 脚本，例如: Rscript BYM_INLA_Production.R --help"
echo

# --- 脚本结束 ---