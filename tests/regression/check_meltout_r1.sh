#!/bin/bash
# Regression check for the melt-out closure fix (design doc section 15g): re-run Calibration-expert's r1 case
# (TVC_Warming_nfeed, CMT51 snwdennew 75, nfactor(w) 0.3, thermal_solver 1) and apply the daily enthalpy gate.
# usage: check_meltout_r1.sh BIN [CASESET] [RUNDIR]    (calib-elchin paths by default; needs REG_LD_LIBRARY_PATH)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1
CS=${2:-$HOME/tem-dev/snw-work/stageE_tests/cs/r1_s0/TVC_Warming_nfeed}
RD=${3:-$(mktemp -d)/r1_s0}
python3 "$HERE/run_case.py" "$BIN" "$CS" "$RD" > "$RD.out" 2>&1 || { cat "$RD.out"; exit 2; }
python3 "$HERE/closure.py" "$RD" --etol-enth 1e-6 > "$RD.closure.json" || true
python3 - "$RD" "$RD.closure.json" <<'PY'
import csv, json, sys
rows = list(csv.DictReader(open(sys.argv[1] + "/main/output/energy_daily.csv")))
w = max(rows, key=lambda r: abs(float(r["resday"])))
d = [r for r in rows if r["day"] == "27162"]
j = json.load(open(sys.argv[2]))
print("r1 worst day", w["day"], "resday", w["resday"], "J/m2; 1 Jun 2024 (day 27162):", d[0]["resday"] if d else "n/a",
      "; Enth_resday_max_abs_Wm2", j.get("Enth_resday_max_abs_Wm2"))
ok = j.get("Enth_resday_max_abs_Wm2", 1.0) <= 1e-6 and j.get("Enth_res_max_abs_Wm2", 1.0) <= 1e-6 and j.get("Enth_fallback_failures", 1) == 0
print("PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
PY
