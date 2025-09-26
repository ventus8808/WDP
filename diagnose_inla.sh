#!/bin/bash
# ========================
# INLA 环境诊断和修复脚本
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=INLA_Diagnosis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA

echo "=== INLA 环境诊断 ==="

# 系统信息
echo "1. 系统信息:"
uname -a
cat /etc/os-release | head -5

# INLA 版本和路径
echo -e "\n2. INLA 信息:"
Rscript -e "
library(INLA)
INLA:::inla.version()
cat('INLA binary path:', INLA:::inla.getOption('inla.call'), '\n')
cat('INLA binary exists:', file.exists(INLA:::inla.getOption('inla.call')), '\n')
"

# 库依赖检查
echo -e "\n3. 库依赖检查:"
INLA_BIN=$(Rscript -e "library(INLA); cat(INLA:::inla.getOption('inla.call'))")
echo "INLA binary: $INLA_BIN"
if [ -f "$INLA_BIN" ]; then
    echo "检查动态库依赖:"
    ldd "$INLA_BIN" | head -10
else
    echo "INLA binary 不存在!"
fi

# 尝试不同的INLA设置
echo -e "\n4. 测试不同的INLA配置:"

# 测试1: 禁用MKL
echo "测试1: 禁用MKL"
export INLA_DISABLE_MKL=1
Rscript -e "
library(INLA)
set.seed(123)
data <- data.frame(y = rpois(10, 5), x = rnorm(10))
tryCatch({
    model <- inla(y ~ x, family = 'poisson', data = data, verbose = FALSE)
    cat('✓ 禁用MKL成功\n')
}, error = function(e) cat('❌ 禁用MKL失败:', e\$message, '\n'))
"

# 测试2: 使用不同的数值库
echo -e "\n测试2: 使用基础数值库"
unset INLA_DISABLE_MKL
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
Rscript -e "
library(INLA)
inla.setOption(num.threads = '1:1')
set.seed(123)
data <- data.frame(y = rpois(10, 5), x = rnorm(10))
tryCatch({
    model <- inla(y ~ x, family = 'poisson', data = data, verbose = FALSE)
    cat('✓ 单线程成功\n')
}, error = function(e) cat('❌ 单线程失败:', e\$message, '\n'))
"

# 测试3: 最小化设置
echo -e "\n测试3: 最小化设置"
Rscript -e "
library(INLA)
# 清除所有环境设置
inla.setOption(inla.mode = 'classic')
inla.setOption(num.threads = '1:1')
inla.setOption(smtp = 'pardiso')
set.seed(123)
data <- data.frame(y = rpois(10, 5), x = rnorm(10))
tryCatch({
    model <- inla(y ~ x, family = 'poisson', data = data,
                  control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
                  verbose = FALSE)
    cat('✓ 最小化设置成功\n')
    print(summary(model))
}, error = function(e) cat('❌ 最小化设置失败:', e\$message, '\n'))
"

# 测试4: 重新安装INLA
echo -e "\n测试4: 检查INLA安装"
Rscript -e "
# 检查INLA安装位置和版本
cat('INLA installed at:', find.package('INLA'), '\n')
cat('R version:', R.version.string, '\n')

# 尝试重新安装INLA (稳定版本)
cat('尝试安装stable版本的INLA...\n')
tryCatch({
    install.packages('INLA', repos = c(getOption('repos'), 
                     INLA = 'https://inla.r-inla-download.org/R/stable'), 
                     dep = TRUE)
    cat('✓ INLA重新安装成功\n')
}, error = function(e) cat('❌ INLA重新安装失败:', e\$message, '\n'))
"

echo -e "\n=== 诊断完成 ==="