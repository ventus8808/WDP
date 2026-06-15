#!/bin/bash
# Combined launcher: submit the SVI MRR pipeline, then the SVI RL pipeline.
# RL is chained with a Slurm dependency (afterok) on the MRR array, so it only
# starts once every MRR task has finished successfully.
#
# This is a controller only — run it on a login node:
#   bash Code/brms_submit/SVI_submit_Main_MRR_RL.sh
#
# It array-submits the existing per-pipeline worker scripts directly (with
# --array set, they run in worker mode), so there is no duplicated job logic.

set -eo pipefail
log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$1] - $2"; }

# --- Locate project root (expects config.yaml there) ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT=""
if [ -f "${SCRIPT_DIR}/../../config.yaml" ]; then
  PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
elif [ -n "${SLURM_SUBMIT_DIR-}" ] && [ -f "${SLURM_SUBMIT_DIR}/config.yaml" ]; then
  PROJECT_ROOT="$SLURM_SUBMIT_DIR"
else
  if [ -f "config.yaml" ]; then PROJECT_ROOT="$(pwd -P)"; fi
fi
if [ -z "$PROJECT_ROOT" ] || [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  log ERROR "无法确定项目根目录 (找不到 config.yaml)"; exit 1
fi
cd "$PROJECT_ROOT"
log INFO "项目根目录: $PROJECT_ROOT"

ENV_NAME="${ENV_NAME:-brms}"
MRR_SCRIPT="Code/brms_submit/SVI_submit_Main_MRR.sh"
RL_SCRIPT="Code/brms_submit/SVI_submit_Main_RL.sh"
MRR_LIST="outcome.list"            # MRR: all outcomes
RL_LIST="outcome_overall.list"     # RL : overall outcomes

for f in "$MRR_SCRIPT" "$RL_SCRIPT" "$MRR_LIST" "$RL_LIST"; do
  if [ ! -f "$f" ]; then log ERROR "找不到: $f"; exit 1; fi
done

count_lines() { grep -Ev '^[[:space:]]*$' "$1" | wc -l | tr -d ' '; }

# --- Phase 1: MRR array ---------------------------------------------------------
N_MRR=$(count_lines "$MRR_LIST")
if [ "$N_MRR" -le 0 ]; then log ERROR "$MRR_LIST is empty"; exit 1; fi
log INFO "提交 MRR 数组任务: 0-$((N_MRR-1)) (共 $N_MRR 个outcomes)"
MRR_JOBID=$(sbatch --parsable --array=0-$((N_MRR-1)) \
  --export=ALL,OUTCOME_FILE="$PROJECT_ROOT/$MRR_LIST",ENV_NAME="$ENV_NAME" \
  "$MRR_SCRIPT")
log INFO "MRR 已提交: JobID=$MRR_JOBID"

# --- Phase 2: RL array (runs after MRR completes successfully) ------------------
N_RL=$(count_lines "$RL_LIST")
if [ "$N_RL" -le 0 ]; then log ERROR "$RL_LIST is empty"; exit 1; fi
log INFO "提交 RL 数组任务: 0-$((N_RL-1)) (共 $N_RL 个outcomes), 依赖 afterok:$MRR_JOBID"
RL_JOBID=$(sbatch --parsable --dependency=afterok:"$MRR_JOBID" --kill-on-invalid-dep=yes \
  --array=0-$((N_RL-1)) \
  --export=ALL,OUTCOME_FILE="$PROJECT_ROOT/$RL_LIST",ENV_NAME="$ENV_NAME" \
  "$RL_SCRIPT")
log INFO "RL 已提交: JobID=$RL_JOBID (等待 MRR JobID=$MRR_JOBID 全部成功后开始)"

log INFO "提交完成。使用 squeue 查看进度。"
log INFO "  MRR: JobID=$MRR_JOBID   RL: JobID=$RL_JOBID"
