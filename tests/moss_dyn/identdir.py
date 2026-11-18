#!/usr/bin/env python3
"""Bitwise compare every NetCDF variable (data only) of all *.nc files under two run dirs.
  identdir.py A B  -> prints JSON {n_files, n_vars, n_identical, differing:[...]}; exit 0 iff all identical"""
import sys, os, json, glob, numpy as np, netCDF4
A, B = sys.argv[1], sys.argv[2]
fa = sorted(os.path.relpath(p, A) for p in glob.glob(os.path.join(A, '**', '*.nc'), recursive=True))
fb = sorted(os.path.relpath(p, B) for p in glob.glob(os.path.join(B, '**', '*.nc'), recursive=True))
fa = [f for f in fa if not f.endswith('fri-fire.nc')]; fb = [f for f in fb if not f.endswith('fri-fire.nc')]
res = dict(n_files=0, n_vars=0, n_identical=0, differing=[], only_one=sorted(set(fa) ^ set(fb)))
for f in sorted(set(fa) & set(fb)):
    res['n_files'] += 1
    da, db = netCDF4.Dataset(os.path.join(A, f)), netCDF4.Dataset(os.path.join(B, f))
    for v in da.variables:
        if f.endswith('run_status.nc') and v == 'total_runtime': continue  # wall-clock seconds
        res['n_vars'] += 1
        if v not in db.variables: res['differing'].append(f + ':' + v + ' (missing)'); continue
        x, y = da[v][:], db[v][:]
        xa = np.ma.getdata(x); ya = np.ma.getdata(y)
        same = xa.shape == ya.shape and xa.tobytes() == ya.tobytes() and \
               np.array_equal(np.ma.getmaskarray(x), np.ma.getmaskarray(y))
        if same: res['n_identical'] += 1
        else:
            try: md = float(np.nanmax(np.abs(xa.astype(float) - ya.astype(float))))
            except Exception: md = None
            res['differing'].append('%s:%s maxabs=%s' % (f, v, md))
print(json.dumps(res, indent=1))
sys.exit(0 if not res['differing'] and not res['only_one'] else 1)
