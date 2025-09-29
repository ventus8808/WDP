# WDP 集群作业脚本使用指南# WDP 集群作业脚本使用指南# WD## 目录



> **集群配置**: 32核心/123GB节点，Hygon C86处理器，Slurm调度器- `debug_test.sh`：**调试专用**极简测试（环境检查+最小模型，30分钟）



本目录提供四个优化的作业脚本，**全部使用32核心配置**，充分利用服务器资源。> **集群配置**: 32核心/123GB节点，Hygon C86处理器，Slurm调度器- `test0_small.sh`：最小可行烟囱测试（单化合物、少量模型、短采样）



## 📋 脚本总览- `test1_completed.sh`：完整"单化合物×全部模型"测试（较长采样）



| 脚本 | 用途 | 时间 | 资源 | 采样配置 | 适用场景 |本目录提供四个优化的作业脚本，按复杂度递增设计，适配集群环境。- `test2_chemical.sh`：完整"多化合物×全部模型"生产跑

|------|------|------|------|----------|----------|

| `debug_test.sh` | 🔧 环境诊断 | 1小时 | 32核96G | 4链×32核 | 首次运行，故障排查 |

| `test0_small.sh` | 🧪 基础验证 | 2小时 | 32核96G | 4链×32核 | 双阶段模型验证 |

| `test1_completed.sh` | 📊 单化合物全测 | 8小时 | 32核96G | 4链×32核 | 完整模型套件 |## 📋 脚本总览> 以上脚本都会自动切换到项目根目录，激活 `pymc` conda 环境，并写日志。**首次使用建议先运行调试脚本排查问题**。脚本使用说明

| `test2_chemical.sh` | 🚀 生产批量 | 2天 | 32核98G | 4链×32核 | 多化合物生产分析 |



> **推荐顺序**: debug → test0 → test1 → test2  

> **性能提升**: 全部32核配置，比之前配置提升4-8倍并行度| 脚本 | 用途 | 时间 | 资源 | 适用场景 |本目录提供三种层次的作业脚本，覆盖从最小烟囱测试到完整化合物×模型的生产运行。



---|------|------|------|------|----------|



## debug_test.sh — 环境诊断测试 ⭐ **首次必运行**| `debug_test.sh` | 🔧 环境诊断 | 1小时 | 4核8G | 首次运行，故障排查 |## 目录



### 功能特性| `test0_small.sh` | 🧪 基础验证 | 2小时 | 8核24G | 双阶段模型验证 |- `test0_small.sh`：最小可行烟囱测试（单化合物、少量模型、短采样）

- **环境检查**: Python版本、PyMC安装、内存状态

- **数据验证**: 逐个检查CDC、协变量、农药、空间数据| `test1_completed.sh` | 📊 单化合物全测 | 8小时 | 16核48G | 完整模型套件 |- `test1_completed.sh`：完整“单化合物×全部模型”测试（较长采样）

- **最小模型**: M0基础模型，100/50采样

- **全核心**: 4链×32核心并行采样| `test2_chemical.sh` | 🚀 生产批量 | 2天 | 32核98G | 多化合物生产分析 |- `test2_chemical.sh`：完整“多化合物×全部模型”生产跑

- **快速完成**: 1小时内完成



### 提交命令

```bash> **推荐顺序**: debug → test0 → test1 → test2> 以上脚本都会自动切换到项目根目录，激活 `pymc` conda 环境，并写日志到 `logs/`。

cd /path/to/WDP

sbatch Shell/debug_test.sh

```

------

**日志文件**: `WDP_Debug-<jobid>.out/.err` + `debug_<jobid>.log`



---

## debug_test.sh — 环境诊断测试 ⭐ **首次必运行**## debug_test.sh — 调试专用测试 ✨ **推荐首次使用**

## test0_small.sh — 基础验证测试

- **适用场景**：集群环境首次运行或故障排查。

### 功能特性

- **双阶段验证**: M0基础模型 → M1社会脆弱性模型### 功能特性- **功能**：

- **中等采样**: 500/200采样，4链×32核心并行

- **快速反馈**: 2小时内完成两个模型- **环境检查**: Python版本、PyMC安装、内存状态  - 环境检查（Python、PyMC、内存状态）

- **错误隔离**: 分阶段执行，便于定位问题

- **数据验证**: 逐个检查CDC、协变量、农药、空间数据  - 详细数据结构检查（CDC、协变量、农药、空间数据）

### 环境变量配置

```bash- **最小模型**: M0基础模型，100/50采样  - 最小模型测试（M0，100/50样本）

# 可选覆盖默认值

export DISEASE="C81-C96"        # 疾病编码- **快速完成**: 1小时内完成- **时间**：30分钟内完成。

export COMPOUND="2"             # 化合物ID

```- **资源**：2核心，3G内存。



### 提交命令### 提交命令

```bash

cd /path/to/WDP```bash提交命令：

sbatch Shell/test0_small.sh

```cd /path/to/WDP```bash



**日志文件**: `WDP_Basic_Test-<jobid>.out/.err` + `smoke_test_<jobid>.log`sbatch Shell/debug_test.shcd /path/to/WDP



---```sbatch Shell/debug_test.sh



## test1_completed.sh — 单化合物完整测试```



### 功能特性**日志文件**: `WDP_Debug-<jobid>.out/.err` + `debug_<jobid>.log`日志：`WDP_Debug-<jobid>.out/.err` + `debug_<jobid>.log`

- **全模型测试**: M0, M1, M2, M3四种模型配置

- **多滞后分析**: 5年和10年滞后效应

- **中等采样**: 1500/750采样，4链×32核心并行

- **详细诊断**: 完整的收敛和效应分析------



### 环境变量配置

```bash

export DISEASE="C81-C96"        # 疾病编码  ## test0_small.sh — 基础验证测试## test0_small.sh — 小规模烟囱测试

export COMPOUND="2"             # 化合物ID

export MODELS="M0,M1,M2,M3"     # 模型类型- 适用场景：首次上集群验证能否跑通端到端。

```

### 功能特性- 默认：`disease=C81-C96`，`compound=2`，`models=M5_SVI,M6_ENV1`，`lag=10`，`estimate=avg`。

### 提交命令

```bash- **双阶段验证**: M0基础模型 → M1社会脆弱性模型- 资源：`--cpus-per-task=4`，`--mem-per-cpu=2G`，`--time=1:00:00`。

cd /path/to/WDP

sbatch Shell/test1_completed.sh- **适中采样**: 500/200采样，2链并行

```

- **快速反馈**: 2小时内完成两个模型提交命令（**请确保在仓库根目录 WDP 下提交**）：

**日志文件**: `WDP_Single_Compound-<jobid>.out/.err`

- **错误隔离**: 分阶段执行，便于定位问题```bash

---

# 确保在项目根目录

## test2_chemical.sh — 生产级批量分析

### 环境变量配置cd /path/to/WDP

### 功能特性

- **大规模批量**: 多化合物×多模型×多滞后×多测量类型```bashsbatch Shell/test0_small.sh

- **生产采样**: 4000/2000高质量采样，4链×32核心并行  

- **全资源利用**: 32核心，98GB内存，48小时时间# 可选覆盖默认值```

- **完整覆盖**: Weight/Density测量，avg/max估算

export DISEASE="C81-C96"        # 疾病编码日志：`WONDER_PyMC_Smoke_Sub-<jobid>.out/.err`（位于项目根目录）

### 环境变量配置

```bashexport COMPOUND="2"             # 化合物ID

export DISEASE="C81-C96"                # 疾病编码

export COMPOUNDS="2,9,cat21,cat33"      # 多个化合物```> **重要提示**: 脚本会自动检测项目根目录，但请务必在WDP仓库根目录下提交作业，以确保路径检测正确。

export MODELS="M0,M1,M2,M3"             # 全部基础模型

```



### 提交命令### 提交命令---

```bash

cd /path/to/WDP  ```bash

sbatch Shell/test2_chemical.sh

```cd /path/to/WDP## test1_completed.sh — 完整单模型测试



**日志文件**: `WDP_Production-<jobid>.out/.err`sbatch Shell/test0_small.sh- 适用场景：对一个化合物跑全套模型，采样更充足以观测诊断。



---```- 环境变量可覆盖：



## 🚀 32核心性能优势  - `DISEASE`（默认 `C81-C96`）



### 并行度对比**日志文件**: `WDP_Basic_Test-<jobid>.out/.err` + `smoke_test_<jobid>.log`  - `COMPOUND`（默认 `2`）

- **之前配置**: 1-16核心，单链或2链采样

- **现在配置**: 32核心，4链并行采样  - `MODELS`（默认 `M0,M1,M2,M3,M5_SVI,M6_ENV1`）

- **性能提升**: 

  - 调试测试: 4倍提升 (1核→32核)---- 资源：`--cpus-per-task=16`，`--mem-per-cpu=3G`，`--time=12:00:00`。

  - 基础验证: 4倍提升 (8核→32核)  

  - 单化合物: 2倍提升 (16核→32核)

  - 生产批量: 保持32核心满负荷

## test1_completed.sh — 单化合物完整测试提交命令：

### MCMC采样优化

- **chains=4**: 4条独立采样链，提高收敛可靠性```bash

- **cores=32**: 每条链充分利用多核心计算

- **内存优化**: 96-98GB内存，支持大数据集高并发### 功能特性# 确保在项目根目录



---- **全模型测试**: M0, M1, M2, M3四种模型配置cd /path/to/WDP



## 🎯 使用建议- **多滞后分析**: 5年和10年滞后效应sbatch Shell/test1_completed.sh



### 新用户工作流- **中等采样**: 1500/750采样，4链并行```

1. **环境验证** → 运行 `debug_test.sh`

2. **基础测试** → 运行 `test0_small.sh`  - **详细诊断**: 完整的收敛和效应分析日志：`WONDER_PyMC_OneChem_AllModels-<jobid>.out/.err`（位于项目根目录）

3. **模型验证** → 运行 `test1_completed.sh`

4. **批量生产** → 运行 `test2_chemical.sh`



### 资源利用策略### 环境变量配置---

- **全核心利用**: 所有脚本都使用32核心，无资源浪费

- **内存充足**: 96-98GB配置，支持大数据集```bash

- **时间合理**: 从1小时调试到2天生产，梯度递增

export DISEASE="C81-C96"        # 疾病编码  ## test2_chemical.sh — 完整化合物批量（生产）

### 采样质量建议

- **快速验证**: debug (100/50), test0 (500/200)export COMPOUND="2"             # 化合物ID- 适用场景：多个化合物×全模型配置的正式运行。

- **标准分析**: test1 (1500/750)

- **发表质量**: test2 (4000/2000)export MODELS="M0,M1,M2,M3"     # 模型类型- 环境变量可覆盖：



---```  - `DISEASE`（默认 `C81-C96`）



## 📊 结果输出  - `COMPOUNDS`（默认 `2,9,cat21,cat33`）



### 文件结构### 提交命令  - `MODELS`（默认 `M0,M1,M2,M3,M5_SVI,M6_ENV1`）

```

Result/PyMC_Results/```bash- 资源：单节点满配：`--cpus-per-task=32` + `--mem=98G` + `--time=1-00:00:00`。

├── <DISEASE>_<COMPOUND>_Results.csv    # 主要结果文件

├── All_Results_Summary.csv             # 汇总表cd /path/to/WDP

└── debug_<jobid>.log                   # 调试日志

```sbatch Shell/test1_completed.sh提交命令：



### 结果内容``````bash

- **效应量**: RR per SD, RR per IQR, 四分位对比

- **协变量效应**: 社会脆弱性、环境因子系数# 确保在项目根目录

- **诊断指标**: R-hat, ESS, WAIC收敛诊断

- **元数据**: 时间戳、参数配置、样本量**日志文件**: `WDP_Single_Compound-<jobid>.out/.err`cd /path/to/WDP



---sbatch Shell/test2_chemical.sh



## 🛠️ 常用SLURM命令---```



### 作业管理日志：`WONDER_PyMC_AllChem-<jobid>.out/.err`（位于项目根目录）

```bash

squeue -u $USER                          # 查看队列## test2_chemical.sh — 生产级批量分析

scontrol show job <jobid>                # 作业详情  

scancel <jobid>                          # 取消作业---

sstat -j <jobid> --format=JobID,MaxRSS   # 运行中资源

sacct -u $USER --format=JobID,State,MaxRSS,Elapsed  # 历史资源### 功能特性

```

- **大规模批量**: 多化合物×多模型×多滞后×多测量类型## 结果输出

### 实时监控

```bash- **生产采样**: 4000/2000高质量采样，4链并行  - 结果会追加写入 `Result/PyMC_Results/<DISEASE>_<COMPOUND>_Results.csv`。

# 查看作业输出

tail -f WDP_Debug-<jobid>.out- **全资源利用**: 32核心，98GB内存，48小时时间- 图与表（若生成）对应写入 `Result/` 子目录。



# 查看错误日志- **完整覆盖**: Weight/Density测量，avg/max估算

tail -f WDP_Debug-<jobid>.err

## 资源与稳定性建议

# 查看节点资源使用

ssh <node> htop### 环境变量配置- 小规模测试先行，再逐步放大至 Test（1000/1000）乃至 Production（2000~4000/2000）。

```

```bash- 样本不多时 R-hat/ESS 警告属正常；正式跑请使用更高 `draws/tune` 和 `target_accept>=0.9`。

---

export DISEASE="C81-C96"                # 疾病编码- 代码已根据 `SLURM_CPUS_PER_TASK` 自动选择初始化策略与并发度，无需手工改；如需强制，脚本内的 `--chains/--cores` 仍可覆盖。

## ⚠️ 故障排除

export COMPOUNDS="2,9,cat21,cat33"      # 多个化合物

### 常见错误1：数据列缺失（已修复 ✅）

**状态**: 已修复。数据加载逻辑已更新，能正确处理CDC数据文件结构。export MODELS="M0,M1,M2,M3"             # 全部基础模型## 常用SLURM命令



### 常见错误2：项目根目录检测失败``````bash

**错误信息**: `未找到 Code/PYMC/main.py`

squeue -u $USER                          # 查看队列

**解决方案**:

```bash### 提交命令scontrol show job <jobid>                # 作业详情

# 确保在正确位置提交

cd /path/to/WDP```bashsstat -j <jobid> --format=JobID,MaxRSS   # 运行中资源

ls config.yaml Code/PYMC/main.py  # 确认文件存在

sbatch Shell/debug_test.shcd /path/to/WDP  sacct -u $USER --format=JobID,State,...  # 历史资源

```

sbatch Shell/test2_chemical.sh```

### 常见错误3：Conda环境问题

**解决方案**:```

```bash

conda create -n pymc python=3.11 -y## 故障排除

conda activate pymc

conda install -c conda-forge pymc=5.23.0**日志文件**: `WDP_Production-<jobid>.out/.err`

```

### 常见错误1：数据列缺失（已修复 ✅）

### 常见错误4：作业运行几分钟后失败

**症状**: 作业提交后运行几分钟就失败，无实时输出---**错误信息**：`缺少必需的列: ['Deaths_Type']`



**解决方案**:

1. **首先运行调试脚本**: `sbatch Shell/debug_test.sh`

2. **查看详细日志**: `cat WDP_Debug-<jobid>.out`## 🎯 使用建议**状态**：已修复。数据加载逻辑已更新，能正确处理CDC数据文件的所有列结构。

3. **逐步增加复杂度**: debug → test0 → test1 → test2



### 常见错误5：内存不足（已优化 ✅）

**状态**: 已优化。所有脚本都配置96-98GB内存，充足应对大数据集。### 新用户工作流### 常见错误2：项目根目录检测失败



---1. **环境验证** → 运行 `debug_test.sh`**错误信息**：`未找到 Code/PYMC/main.py，请检查项目根目录是否正确`



## 📈 性能优化2. **基础测试** → 运行 `test0_small.sh`  



### 集群特定优化3. **模型验证** → 运行 `test1_completed.sh`**原因**：脚本在集群环境下无法正确检测项目根目录。

- **CPU**: **32核心全利用** - 基于Hygon C86架构优化

- **内存**: **96-98GB配置** - 使用123GB节点内存，为系统预留25GB4. **批量生产** → 运行 `test2_chemical.sh`

- **网络**: 利用InfiniBand高速互联（11.2.200.x网段）

- **存储**: 优化大数据集I/O（37TB /data分区）**解决方案**：



### 采样策略优化### 资源优化策略1. **确保在正确位置提交**：

- **调试**: 32核，小样本 (100/50) - 快速环境验证

- **验证**: 32核，中样本 (500-1500/200-750) - 模型可行性确认- **小数据集**: 使用test0或test1   ```bash

- **生产**: 32核，大样本 (4000/2000) - 高质量研究结果

- **单化合物深度分析**: 使用test1   cd /path/to/your/WDP  # 切换到WDP仓库根目录

### MCMC并行优化

```bash- **多化合物对比**: 使用test2   ls config.yaml        # 确认config.yaml存在

# 所有脚本统一配置

--chains 4              # 4条独立采样链- **调试问题**: 总是先运行debug   sbatch Shell/test0_small.sh

--cores 32               # 32核心并行计算  

--mem-per-cpu 3G         # 每核心3GB内存   ```

--cpus-per-task 32       # SLURM 32核心分配

```### 采样质量建议



---- **快速验证**: test0 (500/200)2. **检查目录结构**：



## 💡 最佳实践- **标准分析**: test1 (1500/750)   ```bash



### 首次使用- **发表质量**: test2 (4000/2000)   # 项目根目录应包含以下文件/目录：

1. 运行 `debug_test.sh` 确认环境正常

2. 检查输出日志确认32核心被正确使用   config.yaml

3. 观察内存使用情况，确保无OOM

---   Code/PYMC/main.py

### 生产使用

1. 使用 `test1_completed.sh` 验证单化合物效果   Shell/test0_small.sh

2. 确认收敛质量满意后，运行 `test2_chemical.sh`

3. 定期检查 `sstat` 确认资源使用率## 📊 结果输出   ```



### 性能监控

```bash

# 实时查看CPU使用率### 文件结构3. **手动指定项目路径**（如果自动检测仍失败）：

sstat -j <jobid> --format=JobID,AveCPU,MaxRSS

```   ```bash

# 查看内存峰值

sacct -j <jobid> --format=JobID,MaxRSS,AveRSSResult/PyMC_Results/   # 编辑脚本，在PROJECT_ROOT检测逻辑后添加：

```

├── <DISEASE>_<COMPOUND>_Results.csv    # 主要结果文件   PROJECT_ROOT="/your/actual/WDP/path"  # 手动指定

---

├── All_Results_Summary.csv             # 汇总表   ```

*文档版本: v4.0 - 全32核心优化版*  

*更新时间: 2025年9月29日*  └── debug_<jobid>.log                   # 调试日志

*配置状态: 32核心×4脚本，Hygon C86集群优化*
```### 常见错误3：Conda环境激活失败

**错误信息**：`激活pymc失败`

### 结果内容

- **效应量**: RR per SD, RR per IQR, 四分位对比**解决方案**：

- **协变量效应**: 社会脆弱性、环境因子系数1. 确保pymc环境已创建：

- **诊断指标**: R-hat, ESS, WAIC收敛诊断   ```bash

- **元数据**: 时间戳、参数配置、样本量   conda create -n pymc python=3.11 -y

   conda activate pymc

---   # 安装必需包...

   ```

## 🛠️ 常用SLURM命令

2. 检查conda初始化脚本路径是否正确（脚本会自动检测常见位置）

### 作业管理

```bash### 常见错误4：作业运行几分钟后失败且无实时输出

squeue -u $USER                          # 查看队列**症状**：作业提交后运行4-5分钟，然后失败，看不到MCMC采样进度。

scontrol show job <jobid>                # 作业详情  

scancel <jobid>                          # 取消作业**原因**：

sstat -j <jobid> --format=JobID,MaxRSS   # 运行中资源1. **输出缓冲**：集群环境下Python输出被缓冲，不会实时显示

sacct -u $USER --format=JobID,State,MaxRSS,Elapsed  # 历史资源2. **数值问题**：模型参数或数据导致采样失败

```3. **内存不足**：复杂模型消耗内存过大



### 实时监控**解决方案**：

```bash1. **首先运行调试脚本**：

# 查看作业输出   ```bash

tail -f WDP_Debug-<jobid>.out   sbatch Shell/debug_test.sh  # 30分钟内完成，提供详细诊断

   ```

# 查看错误日志

tail -f WDP_Debug-<jobid>.err2. **查看详细日志**：

   ```bash

# 查看节点资源使用   # 检查Slurm日志

ssh <node> htop   cat WDP_Debug-<jobid>.out

```   cat WDP_Debug-<jobid>.err

   

---   # 检查Python详细日志

   cat debug_<jobid>.log

## ⚠️ 故障排除   ```



### 常见错误1：数据列缺失（已修复 ✅）3. **如果调试脚本成功**，再尝试烟囱测试：

**状态**: 已修复。数据加载逻辑已更新，能正确处理CDC数据文件结构。   ```bash

   sbatch Shell/test0_small.sh

### 常见错误2：项目根目录检测失败   ```

**错误信息**: `未找到 Code/PYMC/main.py`

4. **逐步增加复杂度**：

**解决方案**:   - M0 模型 → M1,M2 → 交互模型

```bash   - 少样本 → 多样本

# 确保在正确位置提交   - 单核 → 多核

cd /path/to/WDP

ls config.yaml Code/PYMC/main.py  # 确认文件存在### 常见错误5：权限或路径问题

sbatch Shell/debug_test.sh**解决方案**：

```1. 确保对项目目录有读写权限

2. 检查SLURM_SUBMIT_DIR环境变量是否正确设置

### 常见错误3：Conda环境问题3. 使用绝对路径提交作业：

**解决方案**:   ```bash

```bash   sbatch /full/path/to/WDP/Shell/debug_test.sh

conda create -n pymc python=3.11 -y   ```

conda activate pymc
conda install -c conda-forge pymc=5.23.0
```

### 常见错误4：作业运行几分钟后失败
**症状**: 作业提交后运行几分钟就失败，无实时输出

**解决方案**:
1. **首先运行调试脚本**: `sbatch Shell/debug_test.sh`
2. **查看详细日志**: `cat WDP_Debug-<jobid>.out`
3. **逐步增加复杂度**: debug → test0 → test1 → test2

### 常见错误5：内存不足
**解决方案**:
```bash
# 减少并行度
--chains 2 --cores 8

# 增加内存分配  
--mem-per-cpu=4G
```

---

## 📈 性能优化

### 集群特定优化
- **CPU**: 基于32核心节点优化的核心分配
- **内存**: 使用123GB节点内存，为系统预留25GB
- **网络**: 利用InfiniBand高速互联（11.2.200.x网段）
- **存储**: 优化大数据集I/O（37TB /data分区）

### 采样策略
- **调试**: 1核，小样本 (100/50)
- **验证**: 多核，中样本 (500-1500/200-750)  
- **生产**: 满核，大样本 (4000/2000)

---

*文档版本: v3.0*  
*更新时间: 2025年9月29日*  
*集群状态: Hygon C86, 317节点, Slurm 22.05.8*