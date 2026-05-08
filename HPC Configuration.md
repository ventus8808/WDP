# HPC 环境安装手册

> 适用于：登录节点有外网、计算节点无外网、CentOS 7 (glibc 2.17)、SLURM 集群
> 
> 关键约束：
> 
> - 登录节点不能跑重计算（会被 kill），只能用于下载和轻量安装
> - 计算节点无外网，只能用本地文件
> - 系统 glibc 2.17 较旧，conda-forge 预编译的新版 cmdstan 有 ABI 兼容性问题
> - 解决方案：用 cmdstanr 从源码编译 cmdstan 2.35.0，在计算节点上编译

---

## 1. 安装 micromamba

```bash
curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
mkdir -p ~/micromamba
mv bin/micromamba ~/micromamba/micromamba
~/micromamba/micromamba shell init --shell bash --root-prefix ~/micromamba
source ~/.bashrc
~/micromamba/micromamba --version
# 预期输出：2.6.0
```

---

## 2. 配置 Git 与 GitHub

```bash
# 安装 git
~/micromamba/micromamba install -n base -c conda-forge git --root-prefix ~/micromamba -y
~/micromamba/micromamba run -n base git --version
```

```bash
# 生成 SSH 密钥,回车三次
ssh-keygen -t ed25519 -C "Ventus8808@iCloud.com"

# 复制公钥，添加到 GitHub → Settings → SSH Keys
cat ~/.ssh/id_ed25519.pub

# 测试连接
ssh -T git@github.com

# 确认身份
git config --global user.email "ventus8808@icloud.com"

# 克隆项目
git clone git@github.com:ventus8808/WDP.git
```

---

## 3. 创建 brms 环境（约 10 分钟）

不安装 cmdstan，让 conda-forge 自选 gcc，避免版本冲突：

```bash
~/micromamba/micromamba create -n brms -c conda-forge \
  r-base \
  r-optparse \
  r-data.table \
  r-dplyr \
  r-stringr \
  r-tidyr \
  r-readr \
  r-purrr \
  r-posterior \
  r-jsonlite r-processx r-ps r-callr r-checkmate \
  --root-prefix ~/micromamba \
  -y
  
```

---
## 4. 本地下载并SFTP cmdstanr和cmdstan源码
```bash

wget https://github.com/stan-dev/cmdstanr/archive/refs/heads/master.tar.gz -O ~/cmdstanr.tar.gz

wget https://github.com/stan-dev/cmdstan/releases/download/v2.35.0/cmdstan-2.35.0.tar.gz -O ~/cmdstan-2.35.0.tar.gz
```
---
## 5. 安装 cmdstanr

```bash
micromamba run -n brms Rscript -e "install.packages('~/cmdstanr.tar.gz', repos=NULL, type='source')"
```

---
## 6. 在计算节点编译 cmdstan

登录节点内存不足，编译必须在计算节点上进行：

```bash
# 申请计算节点（交互式） 乌镇-通用计算
srun --partition=wzhctest --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=12G --time=01:00:00 --pty bash

# 申请计算节点（交互式） 昆山-通用计算
srun --partition=kshctest02 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=12G --time=01:00:00 --pty bash

# 申请计算节点（交互式） 核心分区-异构加速
srun --partition=hx1hdnormal --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=12G --time=01:00:00 --gres=dcu:1 --pty bash

# 进入计算节点后：
cd ~
mkdir -p ~/.cmdstan
tar -xzf ~/cmdstan-2.35.0.tar.gz -C ~/.cmdstan/
cd ~/.cmdstan/cmdstan-2.35.0

# 写入 TBB 编译配置
cat > make/local << EOF
TBB_CXX_TYPE=gcc
LDFLAGS += -L$(pwd)/stan/lib/stan_math/lib/tbb -Wl,-rpath,$(pwd)/stan/lib/stan_math/lib/tbb
EOF

# 编译（约 10-20 分钟）
make build -j4

# 编译成功后验证
micromamba run -n brms Rscript -e "
library(cmdstanr)
set_cmdstan_path('~/.cmdstan/cmdstan-2.35.0')
writeLines('parameters { real y; } model { y ~ normal(0,1); }', '/test.stan')
m <- cmdstan_model('/test.stan')
fit <- m\$sample(chains=1, iter_sampling=100, iter_warmup=100)
print(fit\$summary())
"
# 预期：MCMC 正常运行并输出结果表格

exit  # 退出计算节点
```

---

## 7. 安装 brms（非常慢，至少30min）

```bash
micromamba run -n brms Rscript -e \
  "install.packages('brms', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
```

---

## 8. 验证所有包

```bash
micromamba run -n brms Rscript -e "
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(cmdstanr)
  library(posterior)
  library(brms)
})
cat('所有包加载成功\n')
"
# 预期输出：所有包加载成功
```

---

## R 脚本可移植性配置

**不要硬编码路径**，在 R 脚本开头加以下代码，自动找 `~/.cmdstan/` 下的最新版本：

```r
suppressPackageStartupMessages(library(cmdstanr))
candidates <- list.dirs(path.expand("~/.cmdstan"), recursive = FALSE)
if (length(candidates) == 0) stop("未找到 cmdstan，请先编译安装 cmdstan")
cmdstanr::set_cmdstan_path(tail(sort(candidates), 1))
```

---

## SLURM 提交脚本说明

- 分区：`kshctest02`，时间上限：`2-00:00:00`
- 用 `micromamba run -n brms Rscript` 替代直接调用 `Rscript`
- 不需要 `module load devtoolset-8`，不需要手动指定 gcc
- 激活环境方式：

```bash
export MAMBA_EXE="$HOME/micromamba/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate brms
```

---

## 注意事项

- `~/.cmdstan/` 下只保留一个版本，避免自动选择出错（删除多余版本：`rm -rf ~/.cmdstan/cmdstan-2.38.0`）
- 计算节点无外网，所有下载必须在登录节点完成
- 登录节点内存有限，Stan 模型编译必须在计算节点上进行
- MCMC 运行中出现 `Informational Message` 是正常现象，不影响结果