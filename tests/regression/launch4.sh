#!/bin/bash
# launch4.sh JOBS : run regression cases in parallel with run_case.py.
#   JOBS: one job per line, "CASEDIR RUNDIR [BIN]" (BIN defaults to $BIN). Existing RUNDIRs are skipped.
#   env: P      parallel jobs (default 4)
#        BIN    default binary
#        RC     run_case.py (default: this repo's tests/regression/run_case.py)
#        CONDA_ENV  conda env to activate for the runtime libraries (optional; on calib-elchin
#                   ~/miniconda3/envs/dvmdostem-build); REG_LD_LIBRARY_PATH defaults to its lib/.
#        DRY=1  print the jobs only
# Output: RUNDIR/{meta.json,spin000/,main/} and RUNDIR.out; "ALL_DONE <date>" at the end.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
J=${1:?usage: launch4.sh JOBS}
P=${P:-4}; RC=${RC:-$HERE/run_case.py}; B=${BIN:-}
if [ -n "${CONDA_ENV:-}" ]; then
  source "$(dirname "$(dirname "$CONDA_ENV")")/etc/profile.d/conda.sh" 2>/dev/null || source "$HOME/miniconda3/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
  export REG_LD_LIBRARY_PATH=${REG_LD_LIBRARY_PATH:-$CONDA_PREFIX/lib}
fi
if [ "${DRY:-0}" = 1 ]; then grep -v '^\s*#' "$J" | sed '/^\s*$/d; s/^/DRY /'; exit 0; fi
grep -v '^\s*#' "$J" | sed '/^\s*$/d' | xargs -P "$P" -L 1 bash -c 'b=${2:-'"$B"'}; [ -n "$b" ] || { echo "NO_BIN $1"; exit 1; }
  if [ -e "$1" ]; then echo SKIP_EXISTS $1; exit 0; fi; mkdir -p "$(dirname "$1")"; echo START $1
  nice -n 5 python3 '"$RC"' "$b" "$0" "$1" > "$1.out" 2>&1; echo END rc=$? $1'
echo ALL_DONE $(date -Is)
