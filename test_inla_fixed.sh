#!/bin/bash
# ========================
# 最简INLA测试 - 修复版
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=INLA_Minimal_Test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=15:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA

PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
cd ${PROJECT_ROOT}

# 设置环境变量禁用问题库
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export INLA_DISABLE_MKL=1

echo "=== 最简INLA测试 (修复版) ==="

Rscript -e "
# 加载库
library(INLA)

# 显示版本信息
INLA:::inla.version()

# 设置最保守的INLA选项
cat('设置INLA选项...\n')
inla.setOption(inla.mode = 'classic')
inla.setOption(num.threads = '1:1')
inla.setOption(smtp = 'pardiso')
inla.setOption(safe = TRUE)
inla.setOption(keep = TRUE)

# 创建最简单的测试数据
cat('创建测试数据...\n')
set.seed(42)
n <- 50
test_data <- data.frame(
  y = rpois(n, lambda = 5),
  x1 = rnorm(n, 0, 1),
  expected = rep(5, n)
)
test_data\$log_expected <- log(test_data\$expected)

cat('数据摘要:\n')
print(summary(test_data))

# 测试1: 最简单的截距模型
cat('\n测试1: 截距模型\n')
tryCatch({
  model1 <- inla(y ~ 1, 
                 family = 'poisson', 
                 data = test_data,
                 offset = test_data\$log_expected,
                 control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
                 control.predictor = list(compute = FALSE),
                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE)
  
  cat('✅ 截距模型成功!\n')
  cat('固定效应:\n')
  print(model1\$summary.fixed)
  
}, error = function(e) {
  cat('❌ 截距模型失败:', e\$message, '\n')
})

# 测试2: 简单回归模型  
cat('\n测试2: 简单回归模型\n')
tryCatch({
  model2 <- inla(y ~ 1 + x1, 
                 family = 'poisson', 
                 data = test_data,
                 offset = test_data\$log_expected,
                 control.compute = list(dic = FALSE, waic = FALSE, cpo = FALSE),
                 control.predictor = list(compute = FALSE),
                 control.inla = list(strategy = 'gaussian'),
                 verbose = FALSE)
  
  cat('✅ 回归模型成功!\n')
  cat('固定效应:\n')
  print(model2\$summary.fixed)
  
  # 计算基础结果
  coef_x1 <- model2\$summary.fixed['x1', ]
  rr <- exp(coef_x1\$mean)
  rr_lower <- exp(coef_x1\$'0.025quant')
  rr_upper <- exp(coef_x1\$'0.975quant')
  
  cat('x1的相对风险: RR =', round(rr, 4), 
      '[', round(rr_lower, 4), ',', round(rr_upper, 4), ']\n')
      
}, error = function(e) {
  cat('❌ 回归模型失败:', e\$message, '\n')
})

cat('\n=== 测试完成 ===\n')
"

status=$?
if [ $status -eq 0 ]; then
    echo "✅ 修复测试成功！可以继续主分析"
else
    echo "❌ 修复测试失败，需要进一步诊断"
fi