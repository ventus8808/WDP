#!/bin/bash
# 简单的R脚本测试 - 用于验证路径修复
# 在HPC集群上运行此脚本来测试路径是否正确

PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
cd ${PROJECT_ROOT}

echo "测试项目根目录: $(pwd)"
echo "测试.here文件存在: $(test -f .here && echo 'YES' || echo 'NO')"

# 测试R中的here包是否工作正常
Rscript -e "
if (!require(here, quietly = TRUE)) {
  install.packages('here', repos='http://cran.us.r-project.org')
  library(here, quietly = TRUE)
}

cat('Project root detected by here():', here(), '\n')

# 检查工具文件是否存在
utils_dir <- here('Code', 'INLA', 'INLA_Utils')
cat('Utils directory:', utils_dir, '\n')
cat('Utils directory exists:', dir.exists(utils_dir), '\n')

if (dir.exists(utils_dir)) {
  util_files <- c('INLA_Utils_Data.R', 'INLA_Utils_Model.R', 'INLA_Utils_Results.R')
  for (f in util_files) {
    path <- file.path(utils_dir, f)
    cat('  ', f, ':', file.exists(path), '\n')
  }
}

# 检查配置文件
config_path <- here('Code', 'INLA', 'INLA_Config', 'analysis_config.yaml')
cat('Config file exists:', file.exists(config_path), '\n')
"