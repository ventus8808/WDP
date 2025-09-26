#!/bin/bash
# ========================
# WONDER R分析 调试版本 - 简化模型测试
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=WONDER_R_Debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=30:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

# ========================
# Conda环境设置
# ========================
echo "激活Conda环境..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA

# ========================
# 运行简化测试
# ========================

PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
cd ${PROJECT_ROOT}

# 创建简化的测试R脚本
cat > test_simple_inla.R << 'EOF'
#!/usr/bin/env Rscript

# 加载必要的包
library(INLA)
library(here)
library(dplyr)
library(readr)

# 设置项目根目录
here::set_here()

# 简单测试数据
set.seed(123)
n <- 100
test_data <- data.frame(
  y = rpois(n, lambda = 10),
  x = rnorm(n),
  region = rep(1:10, each = 10),
  expected = rep(10, n)
)

test_data$log_expected <- log(test_data$expected)

# 最简单的泊松模型
cat("测试1: 基础泊松模型\n")
formula1 <- y ~ 1 + x
model1 <- tryCatch({
  inla(formula1, 
       family = "poisson", 
       data = test_data,
       offset = test_data$log_expected,
       control.compute = list(dic = TRUE),
       verbose = TRUE)
}, error = function(e) {
  cat("模型1失败:", e$message, "\n")
  NULL
})

if (!is.null(model1)) {
  cat("✓ 基础模型成功\n")
  print(summary(model1))
} else {
  cat("❌ 基础模型失败\n")
}

# 测试简单的IID随机效应
cat("\n测试2: IID随机效应模型\n")
formula2 <- y ~ 1 + x + f(region, model = "iid")
model2 <- tryCatch({
  inla(formula2, 
       family = "poisson", 
       data = test_data,
       offset = test_data$log_expected,
       control.compute = list(dic = TRUE),
       verbose = TRUE)
}, error = function(e) {
  cat("模型2失败:", e$message, "\n")
  NULL
})

if (!is.null(model2)) {
  cat("✓ IID模型成功\n")
} else {
  cat("❌ IID模型失败\n")
}

cat("INLA基础测试完成\n")
EOF

echo "运行INLA基础测试..."
Rscript test_simple_inla.R

status=$?
if [ $status -eq 0 ]; then
    echo "✓ INLA基础测试成功"
else
    echo "❌ INLA基础测试失败"
fi