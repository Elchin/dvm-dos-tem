#!/usr/bin/env python3
"""Per-cell RTS area blend for analysis: X = (1 - A) * X_baseline + A * X_rts.

X_baseline / X_rts are monthly outputs of the stage-2b `--mode baseline` and
`--mode rts` runs of scripts/run_rts_two_stage.py (identical ice init / spin-up,
CH4 on or off). A = RTS-affected area fraction of the cell (0..1). Pools and
fluxes blend linearly, so the blended C balance closes whenever both runs close
(dC = GPP - RA - TKRHCO2 - RHDWD - TKCH4EMIS - TKLATC - TKRESPC, C incl. TKCH4POOL).

usage: rts_area_blend.py <run_workdir> <A> [VAR ...]   (prints annual sums/Dec values as CSV)
"""
import os, sys
import numpy as np
import netCDF4 as nc

POOLS = {'SOC', 'VEGC', 'DEADC', 'DWDC', 'TKCH4POOL', 'SHLWC', 'DEEPC', 'TKEXICE', 'TKSUBSID', 'TKPONDSAT'}

def read(stage_dir, var, stage='sc'):
    f = os.path.join(stage_dir, '%s_monthly_%s.nc' % (var, stage))
    d = nc.Dataset(f); a = np.ma.filled(np.array(d[var][:], dtype=float), np.nan); d.close()
    return np.nan_to_num(a.reshape(a.shape[0], -1)[:, 0])

def blend(base, rts, A):
    """Linear area blend of two arrays (or dicts of arrays). A=0 returns `base` exactly."""
    if not 0.0 <= A <= 1.0: raise ValueError('A must be in [0, 1]')
    if isinstance(base, dict): return {k: blend(base[k], rts[k], A) for k in base}
    base, rts = np.asarray(base, dtype=float), np.asarray(rts, dtype=float)
    if A == 0.0: return base.copy()
    if A == 1.0: return rts.copy()
    return (1.0 - A) * base + A * rts

def blend_workdir(workdir, A, variables, base_dir=None):
    b = base_dir or os.path.join(workdir, 'stage2b_baseline')
    r = os.path.join(workdir, 'stage2b_rts')
    return {v: blend(read(b, v), read(r, v), A) for v in variables}

if __name__ == '__main__':
    wd, A = sys.argv[1], float(sys.argv[2])
    vs = sys.argv[3:] or ['TKCH4EMIS', 'TKRHCO2', 'TKLATC', 'TKLATCH4', 'SOC']
    out = blend_workdir(wd, A, vs)
    print('year_index,' + ','.join(vs))
    n = len(out[vs[0]]) // 12
    for y in range(n):
        row = [out[v][12 * y + 11] if v in POOLS else out[v][12 * y:12 * y + 12].sum() for v in vs]
        print('%d,' % y + ','.join('%.6g' % x for x in row))
