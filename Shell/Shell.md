# WD## 目录
- `debug_test.sh`：**调试专用**极简测试（环境检查+最小模型，30分钟）
- `test0_small.sh`：最小可行烟囱测试（单化合物、少量模型、短采样）
- `test1_completed.sh`：完整"单化合物×全部模型"测试（较长采样）
- `test2_chemical.sh`：完整"多化合物×全部模型"生产跑

> 以上脚本都会自动切换到项目根目录，激活 `pymc` conda 环境，并写日志。**首次使用建议先运行调试脚本排查问题**。脚本使用说明

本目录提供三种层次的作业脚本，覆盖从最小烟囱测试到完整化合物×模型的生产运行。

## 目录
- `test0_small.sh`：最小可行烟囱测试（单化合物、少量模型、短采样）
- `test1_completed.sh`：完整“单化合物×全部模型”测试（较长采样）
- `test2_chemical.sh`：完整“多化合物×全部模型”生产跑

> 以上脚本都会自动切换到项目根目录，激活 `pymc` conda 环境，并写日志到 `logs/`。

---

## debug_test.sh — 调试专用测试 ✨ **推荐首次使用**
- **适用场景**：集群环境首次运行或故障排查。
- **功能**：
  - 环境检查（Python、PyMC、内存状态）
  - 详细数据结构检查（CDC、协变量、农药、空间数据）
  - 最小模型测试（M0，100/50样本）
- **时间**：30分钟内完成。
- **资源**：2核心，3G内存。

提交命令：
```bash
cd /path/to/WDP
sbatch Shell/debug_test.sh
```
日志：`WDP_Debug-<jobid>.out/.err` + `debug_<jobid>.log`

---

## test0_small.sh — 小规模烟囱测试
- 适用场景：首次上集群验证能否跑通端到端。
- 默认：`disease=C81-C96`，`compound=2`，`models=M5_SVI,M6_ENV1`，`lag=10`，`estimate=avg`。
- 资源：`--cpus-per-task=4`，`--mem-per-cpu=2G`，`--time=1:00:00`。

提交命令（**请确保在仓库根目录 WDP 下提交**）：
```bash
# 确保在项目根目录
cd /path/to/WDP
sbatch Shell/test0_small.sh
```
日志：`WONDER_PyMC_Smoke_Sub-<jobid>.out/.err`（位于项目根目录）

> **重要提示**: 脚本会自动检测项目根目录，但请务必在WDP仓库根目录下提交作业，以确保路径检测正确。

---

## test1_completed.sh — 完整单模型测试
- 适用场景：对一个化合物跑全套模型，采样更充足以观测诊断。
- 环境变量可覆盖：
  - `DISEASE`（默认 `C81-C96`）
  - `COMPOUND`（默认 `2`）
  - `MODELS`（默认 `M0,M1,M2,M3,M5_SVI,M6_ENV1`）
- 资源：`--cpus-per-task=16`，`--mem-per-cpu=3G`，`--time=12:00:00`。

提交命令：
```bash
# 确保在项目根目录
cd /path/to/WDP
sbatch Shell/test1_completed.sh
```
日志：`WONDER_PyMC_OneChem_AllModels-<jobid>.out/.err`（位于项目根目录）

---

## test2_chemical.sh — 完整化合物批量（生产）
- 适用场景：多个化合物×全模型配置的正式运行。
- 环境变量可覆盖：
  - `DISEASE`（默认 `C81-C96`）
  - `COMPOUNDS`（默认 `2,9,cat21,cat33`）
  - `MODELS`（默认 `M0,M1,M2,M3,M5_SVI,M6_ENV1`）
- 资源：单节点满配：`--cpus-per-task=32` + `--mem=98G` + `--time=1-00:00:00`。

提交命令：
```bash
# 确保在项目根目录
cd /path/to/WDP
sbatch Shell/test2_chemical.sh
```
日志：`WONDER_PyMC_AllChem-<jobid>.out/.err`（位于项目根目录）

---

## 结果输出
- 结果会追加写入 `Result/PyMC_Results/<DISEASE>_<COMPOUND>_Results.csv`。
- 图与表（若生成）对应写入 `Result/` 子目录。

## 资源与稳定性建议
- 小规模测试先行，再逐步放大至 Test（1000/1000）乃至 Production（2000~4000/2000）。
- 样本不多时 R-hat/ESS 警告属正常；正式跑请使用更高 `draws/tune` 和 `target_accept>=0.9`。
- 代码已根据 `SLURM_CPUS_PER_TASK` 自动选择初始化策略与并发度，无需手工改；如需强制，脚本内的 `--chains/--cores` 仍可覆盖。

## 常用SLURM命令
```bash
squeue -u $USER                          # 查看队列
scontrol show job <jobid>                # 作业详情
sstat -j <jobid> --format=JobID,MaxRSS   # 运行中资源
sacct -u $USER --format=JobID,State,...  # 历史资源
```

## 故障排除

### 常见错误1：数据列缺失（已修复 ✅）
**错误信息**：`缺少必需的列: ['Deaths_Type']`

**状态**：已修复。数据加载逻辑已更新，能正确处理CDC数据文件的所有列结构。

### 常见错误2：项目根目录检测失败
**错误信息**：`未找到 Code/PYMC/main.py，请检查项目根目录是否正确`

**原因**：脚本在集群环境下无法正确检测项目根目录。

**解决方案**：
1. **确保在正确位置提交**：
   ```bash
   cd /path/to/your/WDP  # 切换到WDP仓库根目录
   ls config.yaml        # 确认config.yaml存在
   sbatch Shell/test0_small.sh
   ```

2. **检查目录结构**：
   ```bash
   # 项目根目录应包含以下文件/目录：
   config.yaml
   Code/PYMC/main.py
   Shell/test0_small.sh
   ```

3. **手动指定项目路径**（如果自动检测仍失败）：
   ```bash
   # 编辑脚本，在PROJECT_ROOT检测逻辑后添加：
   PROJECT_ROOT="/your/actual/WDP/path"  # 手动指定
   ```

### 常见错误3：Conda环境激活失败
**错误信息**：`激活pymc失败`

**解决方案**：
1. 确保pymc环境已创建：
   ```bash
   conda create -n pymc python=3.11 -y
   conda activate pymc
   # 安装必需包...
   ```

2. 检查conda初始化脚本路径是否正确（脚本会自动检测常见位置）

### 常见错误4：作业运行几分钟后失败且无实时输出
**症状**：作业提交后运行4-5分钟，然后失败，看不到MCMC采样进度。

**原因**：
1. **输出缓冲**：集群环境下Python输出被缓冲，不会实时显示
2. **数值问题**：模型参数或数据导致采样失败
3. **内存不足**：复杂模型消耗内存过大

**解决方案**：
1. **首先运行调试脚本**：
   ```bash
   sbatch Shell/debug_test.sh  # 30分钟内完成，提供详细诊断
   ```

2. **查看详细日志**：
   ```bash
   # 检查Slurm日志
   cat WDP_Debug-<jobid>.out
   cat WDP_Debug-<jobid>.err
   
   # 检查Python详细日志
   cat debug_<jobid>.log
   ```

3. **如果调试脚本成功**，再尝试烟囱测试：
   ```bash
   sbatch Shell/test0_small.sh
   ```

4. **逐步增加复杂度**：
   - M0 模型 → M1,M2 → 交互模型
   - 少样本 → 多样本
   - 单核 → 多核

### 常见错误5：权限或路径问题
**解决方案**：
1. 确保对项目目录有读写权限
2. 检查SLURM_SUBMIT_DIR环境变量是否正确设置
3. 使用绝对路径提交作业：
   ```bash
   sbatch /full/path/to/WDP/Shell/debug_test.sh
   ```
