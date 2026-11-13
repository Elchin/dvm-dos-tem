#!/bin/bash
# run_all.sh: rebuild the reference binaries and rerun the fire-balance identity and balance tests.
#
#   tests/regression/reburn/run_all.sh --work DIR [options]
#
# Case data: the 8 calibrated regression cases (leg/) are on calib-elchin only for now
# (~/tem-dev/srd-work/cs/leg, 144 MB; not in the repo). Everything else is generated from them:
# s4_f1 and all smoke variants by mk_variants.py (variants.json), the smoke bases F0g/F1c/F1p by
# mk_smoke_cases.py and E1 by mk_explicit_case.py.
#
# Steps (all under --work, nothing else is written):
#   1. build   git worktree (detached) of 19af750f, 06c41688, f94fceaf -> WORK/src/<sha>;
#              scripts/get_external.sh (Fireweed 31051b0, Burnup 6d56e29, needs GitHub);
#              make -> WORK/bin/dvmdostem_fb2|fb3|fb4   (fb2 = 19af750f, fb3 = 06c41688, fb4 = f94fceaf)
#   2. cases   WORK/cases/reg/s4_f1 (hard links to REGCASES/leg) and WORK/cases/smoke/* (bases + variants)
#   3. runs    identity: 16 regression cases (leg + s4_f1) and E1/F1c/F1p/F0g/F1cdbg with fb2, fb3, fb4,
#                        F1c_on/E1on with fb3, fb4;
#              balance:  the fire-balance, N-leak, overshoot and restart-split runs
#              (20-yr spin + 75-yr main each, run_case.py via ../launch4.sh, P jobs in parallel)
#   4. check   ident.py: fb2 vs fb3 (42 stage folders) and fb3 vs fb4 (46), all outputs and restarts,
#              expected 88/88 identical; analysis scripts -> WORK/summaries, compared with the committed
#              docs/figures/{fire_balance,reburn}/*.json by compare_summaries.py (expected all PASS).
# Cost: about 1.5 core-h on calib-elchin (3 builds ~0.3, ~90 runs of ~1 min each ~1.2), ~2 GB in WORK.
#
# Options:
#   --work DIR        scratch root (required)
#   --regcases DIR    root holding leg/ (default $TEM_REG_CASES or $TEM_WORK/srd-work/cs, TEM_WORK=~/tem-dev)
#   --sets LIST       build,cases,identity,balance,check (default all)
#   --reg-cases LIST  regression cases to run (default all 8), e.g. TVC_Warming_nfeed
#   --reg-sets LIST   leg,s4_f1 (default both)
#   --smoke LIST      smoke cases for the identity set (default E1,F1c,F1p,F0g,F1cdbg)
#   --fixes-on LIST   fixes-on cases for the fb3 vs fb4 identity (default F1c_on,E1on; "" for none)
#   --bins DIR        use prebuilt dvmdostem_fb2|fb3|fb4 from DIR instead of building
#   --conda ENV       conda env with the build toolchain/runtime (calib-elchin: ~/miniconda3/envs/dvmdostem-build)
#   -P N              parallel runs and make -j (default 4)
#   --dry-run         print what would be done, do nothing
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../../.." && pwd); REG=$REPO/tests/regression
TEM_WORK=${TEM_WORK:-$HOME/tem-dev}
WORK=""; REGCASES=${TEM_REG_CASES:-$TEM_WORK/srd-work/cs}; SETS=build,cases,identity,balance,check
REG_CASES="HPC_NoWarming_nfeed HPC_NoWarming_nonfeed HPC_Warming_nfeed HPC_Warming_nonfeed TVC_NoWarming_nfeed TVC_NoWarming_nonfeed TVC_Warming_nfeed TVC_Warming_nonfeed"
REG_SETS="leg s4_f1"; SMOKE="E1 F1c F1p F0g F1cdbg"; FIXON="F1c_on E1on"; BINS=""; CONDA_ENV=${CONDA_ENV:-}; P=4; DRY=0
while [ $# -gt 0 ]; do case "$1" in
  --work) WORK=$2; shift;; --regcases) REGCASES=$2; shift;; --sets) SETS=$2; shift;;
  --reg-cases) REG_CASES=${2//,/ }; shift;; --reg-sets) REG_SETS=${2//,/ }; shift;; --smoke) SMOKE=${2//,/ }; shift;; --fixes-on) FIXON=${2//,/ }; shift;;
  --bins) BINS=$2; shift;; --conda) CONDA_ENV=$2; shift;; -P) P=$2; shift;; --dry-run) DRY=1;;
  -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d; s/^# \{0,1\}//'; exit 0;;
  *) echo "unknown option $1" >&2; exit 2;; esac; shift; done
[ -n "$WORK" ] || { echo "--work DIR is required" >&2; exit 2; }
WORK=$(mkdir -p "$WORK" 2>/dev/null; cd "$WORK" 2>/dev/null && pwd || echo "$WORK")
has() { [[ ",$SETS," == *",$1,"* ]]; }
say() { echo "[run_all $(date +%H:%M:%S)] $*"; }
run() { if [ $DRY = 1 ]; then echo "DRY: $*" >&2; else "$@"; fi; }
declare -A SHA=([fb2]=19af750f [fb3]=06c41688 [fb4]=f94fceaf)
BIN=${BINS:-$WORK/bin}; CASES=$WORK/cases/smoke; RCASES=$WORK/cases/reg; RUNS=$WORK/runs; JOBS=$WORK/jobs; SUM=$WORK/summaries
if [ -n "$CONDA_ENV" ] && [ $DRY = 0 ]; then
  source "$(dirname "$(dirname "$CONDA_ENV")")/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"
  export REG_LD_LIBRARY_PATH=${REG_LD_LIBRARY_PATH:-$CONDA_PREFIX/lib}
fi
export CONDA_ENV
say "repo $REPO, work $WORK, regcases $REGCASES, sets $SETS, P=$P, dry=$DRY"
[ -d "$REGCASES/leg" ] || { echo "no regression cases at $REGCASES/leg (calib-elchin only for now)" >&2; [ $DRY = 1 ] || exit 2; }

# ---- 1. build ----
if has build && [ -z "$BINS" ]; then
  run mkdir -p "$WORK/src" "$BIN"
  for L in fb2 fb3 fb4; do s=${SHA[$L]}; wt=$WORK/src/$s
    say "build $L = $s"
    [ -d "$wt" ] || run git -C "$REPO" worktree add --detach "$wt" "$s"
    run bash "$wt/scripts/get_external.sh"
    if [ -n "$CONDA_ENV" ]; then
      run make -C "$wt" -j"$P" CC="${CXX:-g++}" SITE_SPECIFIC_INCLUDES="-I${CONDA_PREFIX:-}/include" \
          SITE_SPECIFIC_LIBS="-L${CONDA_PREFIX:-}/lib -Wl,-rpath,${CONDA_PREFIX:-}/lib"
    else run make -C "$wt" -j"$P"; fi
    run cp "$wt/dvmdostem" "$BIN/dvmdostem_$L"
  done
fi

# ---- 2. cases ----
if has cases; then
  run mkdir -p "$CASES" "$RCASES"
  say "cases -> $WORK/cases"
  for x in leg; do [ -e "$RCASES/$x" ] || run ln -s "$REGCASES/$x" "$RCASES/$x"; done
  run python3 "$HERE/mk_smoke_cases.py" --regcases "$REGCASES" --cases "$CASES" --figdir "$SUM/tmp"
  run python3 "$HERE/mk_explicit_case.py" --regcases "$REGCASES" --cases "$CASES" --figdir "$SUM/tmp"
  # s4_f1 is hard-linked from the real leg/ (same file system as WORK required for cp --link)
  run python3 "$HERE/mk_variants.py" --regcases "$REGCASES" --cases "$CASES" --out "$WORK/cases/_v" --figdir "$SUM/tmp" --only s4_f1
  [ -e "$RCASES/s4_f1" ] || run ln -s "$WORK/cases/_v/reg/s4_f1" "$RCASES/s4_f1"
  run python3 "$HERE/mk_variants.py" --regcases "$REGCASES" --cases "$CASES" --figdir "$SUM/tmp" \
      --only "$(python3 -c "import json;print(','.join(json.load(open('$HERE/variants.json'))['smoke']))")"
fi

# ---- 3. runs ----
mkjobs() { : > "$JOBS/$1.txt"; }
job() { echo "$1 $2 $BIN/dvmdostem_$3" >> "$JOBS/$4.txt"; }
run mkdir -p "$JOBS" "$RUNS/reg" "$RUNS/smoke"
if has identity; then
  [ $DRY = 1 ] && mkdir -p "$JOBS"; mkjobs identity
  for L in fb2 fb3 fb4; do
    for s in $REG_SETS; do for c in $REG_CASES; do job "$RCASES/$s/$c" "$RUNS/reg/${L^^}_$s/$c" $L identity; done; done
    for c in $SMOKE; do job "$CASES/$c" "$RUNS/smoke/${c}_$L" $L identity; done
  done
  for L in fb3 fb4; do for c in $FIXON; do job "$CASES/$c" "$RUNS/smoke/${c}_$L" $L identity; done; done
  say "identity: $(wc -l < "$JOBS/identity.txt") runs"
  env P=$P DRY=$DRY bash "$REG/launch4.sh" "$JOBS/identity.txt"
fi
if has balance; then
  [ $DRY = 1 ] && mkdir -p "$JOBS"; mkjobs balance
  # fire_balance_summary / smoke_summary run names (19af750f; with all new switches 0 it equals the older binaries)
  for c in F1c_N F1c_C F1c_on F1c_L F1c_onL F1p_on F1p_L F1p_onL; do job "$CASES/$c" "$RUNS/smoke/$c" fb2 balance; done
  for p in "F1c F1c" "F0g F0g" "F1p F1p_c" "E1on E1on_fb2"; do set -- $p; job "$CASES/$1" "$RUNS/smoke/$2" fb2 balance; done
  for c in F1c_on E1on E1; do job "$CASES/$c" "$RUNS/smoke/${c}_fb3" fb3 balance; done
  for c in F1c_onL F1c_onLK F1c_onLD F1c_onLR F1p_onL F1p_onLK E1on E1onL; do job "$CASES/$c" "$RUNS/smoke/${c}_fb4" fb4 balance; done
  # reuse identical identity runs instead of rerunning them (launch4 skips existing run dirs)
  if [ $DRY = 0 ]; then for p in "F1c_fb2 F1c" "F0g_fb2 F0g" "F1p_fb2 F1p_c"; do set -- $p
    [ -d "$RUNS/smoke/$1" ] && [ ! -e "$RUNS/smoke/$2" ] && ln -s "$1" "$RUNS/smoke/$2"; done; fi
  say "balance: $(wc -l < "$JOBS/balance.txt") runs (minus reused) + 3 restart splits"
  env P=$P DRY=$DRY bash "$REG/launch4.sh" "$JOBS/balance.txt"
  for c in E1on_fb3:fb3:E1on E1on_fb2:fb2:E1on E1_fb3:fb3:E1; do IFS=: read r L k <<< "$c"
    [ -e "$RUNS/smoke/${r}_split" ] || run python3 "$HERE/split_run.py" "$BIN/dvmdostem_$L" "$CASES/$k" "$RUNS/smoke/$r" "$RUNS/smoke/${r}_split" 35 75 &
    while [ "$(jobs -rp | wc -l)" -ge "$P" ]; do sleep 2; done
  done; wait
fi

# ---- 4. check ----
if has check; then
  run mkdir -p "$SUM/fire_balance" "$SUM/reburn"
  : > "$WORK/pairs.txt" 2>/dev/null || true
  if [ $DRY = 0 ]; then
    for s in $REG_SETS; do for c in $REG_CASES; do
      echo "$RUNS/reg/FB2_$s/$c $RUNS/reg/FB3_$s/$c"; echo "$RUNS/reg/FB3_$s/$c $RUNS/reg/FB4_$s/$c"; done; done > "$WORK/pairs.txt"
    for c in $SMOKE; do echo "$RUNS/smoke/${c}_fb2 $RUNS/smoke/${c}_fb3"; echo "$RUNS/smoke/${c}_fb3 $RUNS/smoke/${c}_fb4"; done >> "$WORK/pairs.txt"
    for c in $FIXON; do echo "$RUNS/smoke/${c}_fb3 $RUNS/smoke/${c}_fb4"; done >> "$WORK/pairs.txt"
  fi
  rc=0
  run python3 "$REG/ident.py" --pairs "$WORK/pairs.txt" --json "$SUM/identity.json" || rc=1
  A="--runs $RUNS/smoke --work $WORK/analysis"
  if has balance || [ -d "$RUNS/smoke/F1c_onLK_fb4" ]; then
    run python3 "$HERE/analyze_fire_balance.py" $A --figdir "$SUM/fire_balance" > "$WORK/analysis_fire_balance.log" || rc=1
    run python3 "$HERE/analyze_n_leaks.py" $A --figdir "$SUM/fire_balance" > "$WORK/analysis_n_leaks.log" || rc=1
    run python3 "$HERE/redep_overshoot.py" $A --figdir "$SUM/fire_balance" > "$WORK/analysis_overshoot.log" || rc=1
    run python3 "$HERE/compare_split.py" $A --figdir "$SUM/fire_balance" > "$WORK/analysis_split.log" || rc=1
    run python3 "$HERE/plot_restart_continuity.py" $A --figdir "$SUM/fire_balance" > /dev/null || rc=1
    run python3 "$HERE/analyze_smoke.py" $A --figdir "$SUM/reburn" > "$WORK/analysis_smoke.log" || rc=1
    F=$REPO/docs/figures
    run python3 "$HERE/compare_summaries.py" \
      "$SUM/fire_balance/fire_balance_summary.json" "$F/fire_balance/fire_balance_summary.json" \
      "$SUM/fire_balance/n_leak_summary.json" "$F/fire_balance/n_leak_summary.json" \
      "$SUM/fire_balance/overshoot.json" "$F/fire_balance/overshoot.json" \
      "$SUM/fire_balance/split_summary.json" "$F/fire_balance/split_summary.json" \
      "$SUM/reburn/smoke_summary.json" "$F/reburn/smoke_summary.json" || rc=1
  fi
  say "check done, rc=$rc (identity: $SUM/identity.json; summaries: $SUM)"
  exit $rc
fi
