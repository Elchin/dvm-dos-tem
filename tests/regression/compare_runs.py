#!/usr/bin/env python3
"""Bit-identity of two runs of the same case (logic of /workspace/v5i/scripts/bitid.py):
every variable of every NetCDF file in spin000/output and main/output must be identical
(NaN == NaN; run_status.nc:total_runtime ignored), the file lists must match, and
integrator_safety_net.csv must be byte-identical (or absent in both).
usage: compare_runs.py RUN_A RUN_B   (prints JSON; exit 0 if identical)"""
import json, sys
from pathlib import Path
import numpy as np, netCDF4 as nc

def cmp(a, b):
    out = []
    for st in ("spin000", "main"):
        da, db = Path(a) / st / "output", Path(b) / st / "output"
        fa = sorted(p.name for p in da.glob("*.nc")); fb = sorted(p.name for p in db.glob("*.nc"))
        nv = nd = 0; diffs = []
        if fa != fb: diffs.append(["FILELIST", sorted(set(fa) ^ set(fb))]); nd += 1
        for f in fa:
            if f not in fb: continue
            with nc.Dataset(da / f) as A, nc.Dataset(db / f) as B:
                if set(A.variables) != set(B.variables): diffs.append([f, "VARLIST"]); nd += 1
                for v in A.variables:
                    if v not in B.variables or (f == "run_status.nc" and v == "total_runtime"): continue
                    nv += 1
                    x = np.ma.getdata(A[v][:]); y = np.ma.getdata(B[v][:])
                    same = x.shape == y.shape and (np.array_equal(x, y, equal_nan=True) if x.dtype.kind == "f" else np.array_equal(x, y))
                    if not same:
                        nd += 1
                        md = float(np.nanmax(np.abs(x.astype(float) - y.astype(float)))) if x.shape == y.shape and x.dtype.kind in "fiu" else None
                        diffs.append([f, v, md])
        csv = [d / "integrator_safety_net.csv" for d in (da, db)]
        c_same = (not csv[0].exists() and not csv[1].exists()) or (csv[0].exists() and csv[1].exists() and csv[0].read_bytes() == csv[1].read_bytes())
        if not c_same: nd += 1
        out.append(dict(stage=st, files=len(fa), vars=nv, ndiff=nd, safetynet_csv_same=c_same, diffs=diffs[:20]))
    return out

if __name__ == "__main__":
    r = cmp(sys.argv[1], sys.argv[2]); print(json.dumps(r))
    sys.exit(0 if all(s["ndiff"] == 0 for s in r) else 1)
