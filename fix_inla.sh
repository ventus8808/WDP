#!/bin/bash
# ========================
# INLA 修复策略脚本
# 尝试多种方法让INLA正常工作
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=INLA_Fix
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

source ~/miniconda3/etc/profile.d/conda.sh

echo "=== INLA 修复策略 ==="

# 策略1: 重新创建conda环境
echo "策略1: 重新创建INLA环境"
conda deactivate
conda remove -n INLA --all -y
conda create -n INLA -c conda-forge r-base=4.3 r-devtools -y
conda activate INLA

# 安装R包
Rscript -e "
# 设置CRAN镜像
options(repos = c(CRAN = 'https://cloud.r-project.org/'))

# 安装基础包
install.packages(c('here', 'dplyr', 'readr', 'yaml', 'argparse'), quiet = TRUE)

# 安装stable版本的INLA
install.packages('INLA', repos = c(getOption('repos'), 
                 INLA = 'https://inla.r-inla-download.org/R/stable'), 
                 dep = TRUE, quiet = TRUE)
"

echo "测试新环境:"
Rscript -e "
library(INLA)
INLA:::inla.version()

# 设置保守的INLA选项
inla.setOption(inla.mode = 'classic')
inla.setOption(num.threads = '1:1')
inla.setOption(smtp = 'pardiso')

# 测试最简单的模型
set.seed(123)
n <- 20
data <- data.frame(
    y = rpois(n, 3),
    x = rnorm(n)
)

tryCatch({
    model <- inla(y ~ x, 
                  family = 'poisson', 
                  data = data,
                  control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
                  verbose = FALSE)
    cat('✅ 策略1成功 - 新环境可以运行INLA\n')
    
    # 保存成功的配置
    cat('成功的配置:\n')
    cat('- inla.mode: classic\n')
    cat('- num.threads: 1:1\n')  
    cat('- smtp: pardiso\n')
    cat('- 所有compute选项关闭\n')
    
}, error = function(e) {
    cat('❌ 策略1失败:', e\$message, '\n')
})
"

# 如果策略1失败，尝试策略2：使用Docker/Singularity
if [ $? -ne 0 ]; then
    echo -e "\n策略2: 检查是否可以使用容器"
    
    # 检查singularity
    if command -v singularity &> /dev/null; then
        echo "✓ Singularity 可用"
        echo "建议使用容器运行INLA："
        echo "singularity pull docker://rocker/r-ver:4.3"
        echo "singularity exec r-ver_4.3.sif R --slave -e \"install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'))\""
    else
        echo "❌ 没有容器环境"
    fi
    
    # 策略3：使用系统R
    echo -e "\n策略3: 尝试系统R"
    if command -v R &> /dev/null; then
        echo "找到系统R:"
        R --version | head -1
        
        R --slave -e "
        if (!require('INLA', quietly = TRUE)) {
            install.packages('INLA', repos = c(getOption('repos'), 
                           INLA = 'https://inla.r-inla-download.org/R/stable'))
        }
        
        library(INLA)
        inla.setOption(inla.mode = 'classic')
        inla.setOption(num.threads = '1:1')
        
        set.seed(123)
        data <- data.frame(y = rpois(10, 3), x = rnorm(10))
        tryCatch({
            model <- inla(y ~ x, family = 'poisson', data = data, verbose = FALSE)
            cat('✅ 系统R成功\n')
        }, error = function(e) cat('❌ 系统R失败:', e\$message, '\n'))
        "
    else
        echo "❌ 没有找到系统R"
    fi
fi

echo -e "\n=== 修复尝试完成 ==="