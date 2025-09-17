#!/bin/bash

#================================================================================
# SLURM 资源配置
#================================================================================
#SBATCH --partition=kshdtest                 # 指定分区
#SBATCH --job-name=WDP_Analysis              # 任务名 (会被提交脚本覆盖)
#SBATCH --nodes=1                            # 使用 1 个节点
#SBATCH --ntasks-per-node=1                  # 每个节点运行 1 个任务
#SBATCH --cpus-per-task=4                    # 为每个任务分配 4 个 CPU 核心 (匹配 R 代码中的 inla.setOption)
#SBATCH --mem=64G                            # 分配 64GB 内存
#SBATCH --time=12:00:00                      # 任务最长运行时间 12 小时
#SBATCH --output=slurm_logs/%x-%j.out        # 标准输出日志 (%x=任务名, %j=任务ID)
#SBATCH --error=slurm_logs/%x-%j.err         # 标准错误日志

#================================================================================
# 任务执行环境
#================================================================================

# 接收从 sbatch 或 submit_jobs.sh 传入的参数
if [ "$#" -ne 2 ]; then
    echo "错误: 需要提供两个参数: 疾病代码 和 化合物ID"
    exit 1
fi

DISEASE_CODE="$1"
COMPOUND_ID="$2"

echo "================================================="
echo "WONDER 研究计划 - Slurm 任务开始"
echo "开始时间: $(date)"
echo "运行节点: $(hostname)"
echo "任务 ID: $SLURM_JOB_ID"
echo "任务名称: $SLURM_JOB_NAME"
echo "-------------------------------------------------"
echo "分析参数:"
echo "  疾病代码 (Disease Code): ${DISEASE_CODE}"
echo "  化合物 ID (Compound ID): ${COMPOUND_ID}"
echo "================================================="

# 加载 R 模块 (重要: 请根据服务器实际情况修改模块名)
# 您可以使用 'module avail' 命令查看可用的 R 模块
echo "加载 R 模块..."
module load R/4.2.2-gfbr-2022b # 这是一个示例，请务必替换为您的服务器上正确的 R 模块

# 定义项目路径和脚本路径
# SLURM_SUBMIT_DIR 是提交任务时所在的目录，这里应该是 /public/home/acf4pijnzl/WDP
PROJECT_DIR=$SLURM_SUBMIT_DIR
SCRIPT_PATH="${PROJECT_DIR}/Code/BYM_INLA_Production.R"
CONFIG_PATH="${PROJECT_DIR}/config.yaml"

# 切换到项目根目录，以确保所有相对路径正确
cd $PROJECT_DIR

#================================================================================
# 执行 R 分析脚本
#================================================================================
echo "开始执行 R 分析脚本..."

Rscript ${SCRIPT_PATH} \
  --config ${CONFIG_PATH} \
  --disease-code "${DISEASE_CODE}" \
  --pesticide-category "compound:${COMPOUND_ID}" \
  --measure-type "Weight,Density" \
  --estimate-types "avg" \
  --lag-years "5,10" \
  --model-types "M0,M1,M2,M3" \
  --verbose

# 检查 R 脚本的退出状态
status=$?
if [ $status -ne 0 ]; then
    echo "错误: R 脚本执行失败，退出代码: $status"
    exit $status
fi

echo "================================================="
echo "R 脚本执行完成。"
echo "结束时间: $(date)"
echo "WONDER 研究计划 - Slurm 任务结束"
echo "================================================="
