#!/usr/bin/env python3
"""Run one regression case with one binary (same protocol and launcher as vrun_seb.one, 1 segment):
  spin000: 20-yr sc spin-up on the 1950-69 forcing from the case's restart (restart-tr.nc),
  main:    sc -n 75 (1950-2024) from the spin-up's restart-sc.nc.
argv: -p 0 -e 0 -s 0 -t 0 -n N -f config/config.js -l disabled --no-output-cleanup
env : fixed dict (no inherited environment), OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1,
      LD_LIBRARY_PATH=$REG_LD_LIBRARY_PATH (default /workspace/wf4-runtime/lib), stdin=/dev/null.
usage: run_case.py BINARY CASEDIR RUNDIR      -> RUNDIR/{meta.json, spin000/, main/}"""
import hashlib, json, re, shutil, subprocess, sys, time
import os
from pathlib import Path
ENV = {"PATH": "/usr/bin:/bin", "HOME": "/tmp",
       "LD_LIBRARY_PATH": os.environ.get("REG_LD_LIBRARY_PATH", "/workspace/wf4-runtime/lib"),
       "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "LANG": "C"}

def md5(p):
    h = hashlib.md5(); h.update(Path(p).read_bytes()); return h.hexdigest()

def stage(case, cmeta, name, src, dst, restart):
    shutil.copytree(case / src, dst, ignore=shutil.ignore_patterns("output"))
    shutil.copytree(case / "parameters", dst / "parameters", dirs_exist_ok=True)
    (dst / "output").mkdir(); shutil.copy2(restart, dst / "output/restart-tr.nc")
    c = dst / "config/config.js"; t = c.read_text().replace(cmeta["stage_root"][src], str(dst))
    assert cmeta["stage_root"][src] not in t; c.write_text(t)

def run(d, n, b):
    t0 = time.time()
    with open(d / "run.log", "w") as f:
        f.write(f"BIN={b} N={n}\n"); f.flush()
        rc = subprocess.call([b, "-p", "0", "-e", "0", "-s", "0", "-t", "0", "-n", str(n), "-f", "config/config.js",
                              "-l", "disabled", "--no-output-cleanup"], cwd=d, stdout=f, stderr=subprocess.STDOUT,
                             env=ENV, stdin=subprocess.DEVNULL)
        f.write(f"EXIT_CODE={rc} SEC={time.time()-t0:.0f}\n")
    return rc

def check(d):
    import numpy as np, netCDF4 as nc
    o = d / "output"; r = {}
    r["fail_logs"] = sorted(p.name for p in o.glob("fail_log*"))
    try: r["run_status"] = [int(x) for x in np.unique(nc.Dataset(o / "run_status.nc")["run_status"][:])]
    except Exception as e: r["run_status"] = f"ERR {e}"
    log = (d / "run.log").read_text()
    m = re.findall(r"soil_safetynet_events=(\d+)", log); r["soil_safetynet_events"] = sum(int(x) for x in m) if m else 0
    g = re.findall(r"veg_integrator_giveups=(\d+)", log); r["veg_giveups"] = sum(int(x) for x in g) if g else 0
    csv = o / "integrator_safety_net.csv"
    r["safetynet_csv_rows"] = max(0, len(csv.read_text().strip().splitlines()) - 1) if csv.exists() else 0
    return r

def ok(ck): return ck["exit"] == 0 and not ck["fail_logs"] and ck["run_status"] == [100]

def main(b, case, rd):
    b = str(Path(b).resolve()); case = Path(case).resolve(); rd = Path(rd).resolve()
    cmeta = json.loads((case / "case.json").read_text())
    assert not rd.exists(), f"{rd} exists"
    rd.mkdir(parents=True)
    meta = dict(case=cmeta["case"], bin=b, bin_md5=md5(b), env=ENV, caseset=str(case),
                started=time.strftime("%Y-%m-%d %H:%M:%S %z"))
    sp = rd / "spin000"; stage(case, cmeta, "spin", "spinup", sp, case / "restart0.nc")
    rc = run(sp, cmeta["n_spin"], b); ck = check(sp); ck["exit"] = rc; meta["spin"] = ck
    if not ok(ck):
        meta["status"] = "SPIN_FAIL"
    else:
        mn = rd / "main"; stage(case, cmeta, "main", "main", mn, sp / "output/restart-sc.nc")
        rc = run(mn, cmeta["n_main"], b); ck = check(mn); ck["exit"] = rc; meta["main"] = ck
        meta["status"] = "OK" if ok(ck) else "MAIN_FAIL"
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S %z")
    (rd / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(cmeta["case"], meta["status"], flush=True)
    return 0 if meta["status"] == "OK" else 1

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
