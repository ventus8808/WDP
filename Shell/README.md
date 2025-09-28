# WDP 集群提交脚本使用说明

本目录提供三种层次的作业脚本，覆盖从最小烟囱测试到完整化合物×模型的生产运行。

## 目录
- `test0_small.sh`：最小可行烟囱测试（单化合物、少量模型、短采样）
- `test1_completed.sh`：完整“单化合物×全部模型”测试（较长采样）
- `test2_chemical.sh`：完整“多化合物×全部模型”生产跑

> 以上脚本都会自动切换到项目根目录，激活 `pymc` conda 环境，并写日志到 `logs/`。

---

## test0_small.sh — 小规模烟囱测试
- 适用场景：首次上集群验证能否跑通端到端。
- 默认：`disease=C81-C96`，`compound=2`，`models=M5_SVI,M6_ENV1`，`lag=10`，`estimate=avg`。
- 资源：`--cpus-per-task=4`，`--mem-per-cpu=2G`，`--time=1:00:00`。

提交命令（在仓库根目录 WDP 下）：
```bash
sbatch Shell/test0_small.sh
```
日志：`logs/WONDER_Smoke-<jobid>.out/.err`

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
sbatch Shell/test1_completed.sh
```
日志：`logs/WONDER_OneChem_<jobid>.out/.err`

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
sbatch Shell/test2_chemical.sh
```
日志：`logs/WONDER_AllChem_<jobid>.out/.err`

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
