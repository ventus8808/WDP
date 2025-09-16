# A minimal R script to test the batch execution environment.

# 打印一条信息，证明脚本开始运行
cat("=========================\n")
cat("Hello, R World from Slurm!\n")
cat("=========================\n\n")

# 打印R的会话信息，这可以告诉我们脚本是在哪个环境下运行的
print(sessionInfo())

cat("\nScript finished successfully.\n")
