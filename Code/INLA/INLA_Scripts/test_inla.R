#!/usr/bin/env Rscript
# INLA Installation and Functionality Test
# 测试INLA是否正确安装和配置

cat("🔍 INLA诊断测试\n")
cat("================\n\n")

# 1. Check INLA package
cat("📦 [1] 检查INLA包安装...\n")
if (require("INLA", quietly = TRUE)) {
  cat("✅ INLA包已安装\n")
  cat(sprintf("   版本: %s\n", packageVersion("INLA")))
} else {
  cat("❌ INLA包未找到\n")
  quit(status = 1)
}

# 2. Check INLA binary
cat("\n🔧 [2] 检查INLA二进制文件...\n")
inla_binary <- inla.getOption("inla.call")
cat(sprintf("   INLA调用命令: %s\n", inla_binary))

if (file.exists(inla_binary) || Sys.which(inla_binary) != "") {
  cat("✅ INLA二进制文件可访问\n")
} else {
  cat("⚠️ INLA二进制文件可能有问题\n")
}

# 3. Test basic INLA functionality
cat("\n🧪 [3] 测试基础INLA功能...\n")
tryCatch({
  # Create simple test data
  n <- 100
  test_data <- data.frame(
    y = rpois(n, lambda = 2),
    x = rnorm(n),
    idx = 1:n
  )
  
  # Try a very simple model
  simple_result <- inla(y ~ x, 
                       data = test_data, 
                       family = "poisson",
                       verbose = FALSE)
  
  if (!is.null(simple_result)) {
    cat("✅ 基础INLA模型测试成功\n")
  } else {
    cat("❌ 基础INLA模型测试失败\n")
  }
  
}, error = function(e) {
  cat(sprintf("❌ INLA测试失败: %s\n", e$message))
})

# 4. Test working directory
cat("\n📁 [4] 检查工作目录...\n")
work_dir <- inla.getOption("working.directory")
temp_dir <- tempdir()

cat(sprintf("   INLA工作目录: %s\n", work_dir))
cat(sprintf("   系统临时目录: %s\n", temp_dir))

# Test write permissions
test_file <- file.path(temp_dir, "inla_test.txt")
tryCatch({
  writeLines("test", test_file)
  if (file.exists(test_file)) {
    file.remove(test_file)
    cat("✅ 临时目录写入权限正常\n")
  } else {
    cat("❌ 无法在临时目录创建文件\n")
  }
}, error = function(e) {
  cat(sprintf("❌ 临时目录权限测试失败: %s\n", e$message))
})

# 5. Check memory and system info
cat("\n💻 [5] 系统信息...\n")
cat(sprintf("   R版本: %s\n", R.version.string))
cat(sprintf("   平台: %s\n", R.version$platform))
cat(sprintf("   系统: %s\n", Sys.info()["sysname"]))

# Memory info (basic)
memory_info <- gc()
cat(sprintf("   内存使用: %.1f MB\n", sum(memory_info[, 2])))

cat("\n✅ 诊断完成\n")
cat("如果所有测试通过但模型仍失败，问题可能在于:\n")
cat("  • 数据复杂度过高\n")
cat("  • 空间模型配置问题\n") 
cat("  • 内存不足\n")
cat("  • 临时文件系统问题\n")