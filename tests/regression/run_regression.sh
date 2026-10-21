#!/bin/bash
# Regression harness for daily-mode dvm-dos-tem: compare two binaries on a case set.
#
# usage: run_regression.sh [-j N] [--no-expect-identical] [--no-repeat] [--strict] BIN_A BIN_B CASESET OUTDIR
#
#   BIN_A, BIN_B  two dvmdostem binaries (A = reference, B = candidate)
#   CASESET       directory made by prepare_caseset.py (OUTDIR/<case>/case.json ...)
#   OUTDIR        new directory for the runs and the report (must not exist)
#
# For every case: 20-yr spin-up on the 1950-69 forcing from the case restart, then sc -n 75
# (run_case.py; --no-output-cleanup, OPENBLAS_NUM_THREADS=1, fixed environment). Runs:
#   A/<case>  with BIN_A,   B/<case>  with BIN_B,   R/<case>  BIN_B again (determinism).
# (A, B, R have equal-length paths on purpose.)
# Checks (summarize.py): A vs B bit-identity of spin-up and main outputs (required unless
# --no-expect-identical), B vs R bit-identity (determinism, always required), C/water/N closure of
# B (water |resid| < 1% P, C max |annual| < 0.5 gC/m2, N |mean| < 1 gN/m2/yr), run_status == 100
# and no fail_log for A and B. Reported as WARN (FAIL with --strict): B has more integrator
# safety-net / veg give-up events than A, or B's water/N closure is worse than A's.
# Writes OUTDIR/report.json and OUTDIR/report.txt; exit 0 = PASS.
#
# Environment: runs get ONLY the variables in run_case.ENV; LD_LIBRARY_PATH is
# $REG_LD_LIBRARY_PATH (default /workspace/wf4-runtime/lib). Bit-identity is only meaningful
# between runs made on the same machine with the same runtime libraries (on the GCP VM, output
# depended on the size of the inherited environment; use the box with the wf4 runtime).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
J=8; EXPECT=1; REPEAT=1; STRICT=0
while [[ $# -gt 0 && "$1" == -* ]]; do
  case "$1" in
    -j) J="$2"; shift 2;;
    --no-expect-identical) EXPECT=0; shift;;
    --no-repeat) REPEAT=0; shift;;
    --strict) STRICT=1; shift;;
    -h|--help) sed -n '2,25p' "$0"; exit 0;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
[[ $# -eq 4 ]] || { sed -n '2,8p' "$0"; exit 2; }
BIN_A="$(readlink -f "$1")"; BIN_B="$(readlink -f "$2")"; CS="$(readlink -f "$3")"; OUT="$4"
[[ -x "$BIN_A" && -x "$BIN_B" ]] || { echo "binaries not executable" >&2; exit 2; }
[[ -f "$CS/caseset.json" ]] || { echo "$CS is not a case set (run prepare_caseset.py)" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "$OUT exists" >&2; exit 2; }
mkdir -p "$OUT"; OUT="$(readlink -f "$OUT")"
CASES=$(python3 -c "import json,sys;print(' '.join(json.load(open(sys.argv[1]))['cases']))" "$CS/caseset.json")
echo "A=$BIN_A ($(md5sum < "$BIN_A" | cut -c1-32))" | tee "$OUT/binaries.txt"
echo "B=$BIN_B ($(md5sum < "$BIN_B" | cut -c1-32))" | tee -a "$OUT/binaries.txt"
JOBS="$OUT/jobs.txt"; : > "$JOBS"
for c in $CASES; do
  echo "$BIN_A $CS/$c $OUT/A/$c" >> "$JOBS"
  echo "$BIN_B $CS/$c $OUT/B/$c" >> "$JOBS"
  [[ $REPEAT -eq 1 ]] && echo "$BIN_B $CS/$c $OUT/R/$c" >> "$JOBS"
done
echo "running $(wc -l < "$JOBS") case runs with -j $J ..."
t0=$(date +%s)
# a failing model run is reported by summarize.py, so don't abort here
xargs -P "$J" -L 1 python3 "$HERE/run_case.py" < "$JOBS" > "$OUT/runs.log" 2>&1 || true
echo "runs done in $(( $(date +%s) - t0 )) s"
python3 "$HERE/summarize.py" "$OUT" --expect-identical "$EXPECT" --repeat "$REPEAT" --strict "$STRICT" $CASES
