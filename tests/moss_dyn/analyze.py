#!/usr/bin/env python3
"""moss_dyn run summary: yearly MOSSDZ / ALD / YSD, moss PFT C, and the monthly C balance
residual r = d(VEGC+SOC+DEADC+DWDC) + TKNETC  (stage 2a: TKNETC = RECO - GPP when RTS is off;
fire months add BURNVEG2AIRC + BURNSOIL2AIRC). Prints JSON.
  analyze.py RUNDIR [--moss-pfts 5,6,7]"""
import sys, os, json, numpy as np, netCDF4
R = sys.argv[1]
mp = [int(x) for x in (sys.argv[sys.argv.index('--moss-pfts') + 1].split(',') if '--moss-pfts' in sys.argv else [])]

def rd(stage, v, res):
    sub, tag = {'tr': ('stage2a_hist', 'tr'), 'sc': ('stage2b_rts', 'sc')}[stage]
    p = os.path.join(R, sub, '%s_%s_%s.nc' % (v, res, tag))
    if not os.path.exists(p): return None
    x = np.ma.filled(netCDF4.Dataset(p)[v][:].astype(float), np.nan)
    return x.reshape(x.shape[0], -1) if x.ndim > 3 else x.reshape(x.shape[0])

def series(v, res):
    a, b = rd('tr', v, res), rd('sc', v, res)
    return None if a is None else (np.concatenate([a, b]) if b is not None else a)

out = {}
ntr = len(rd('tr', 'MOSSDZ', 'yearly')); ny = len(series('MOSSDZ', 'yearly'))
start = {'G': 1950, 'L': 1950}.get(json.load(open(os.path.join(R, 'meta.json')))['site'], 1901)
years = np.arange(start, start + ny)
out['years'] = [int(years[0]), int(years[-1])]
for v in ('MOSSDZ', 'ALD', 'YSD', 'SHLWDZ'):
    s = series(v, 'yearly'); out[v] = [round(float(x), 5) for x in s]
veg = series('VEGC', 'monthly')  # (t, pft)
out['moss_c_dec'] = [round(float(x), 3) for x in veg[11::12][:, mp].sum(axis=1)] if mp else None
stock = veg.sum(axis=1) + series('SOC', 'monthly') + series('DEADC', 'monthly') + series('DWDC', 'monthly')
flux = series('TKNETC', 'monthly')
g = series('GPP', 'monthly'); g = g.sum(axis=1) if g.ndim > 1 else g
nee = series('RECO', 'monthly') - g
burn = np.nan_to_num(series('BURNVEG2AIRC', 'monthly')) + np.nan_to_num(series('BURNSOIL2AIRC', 'monthly'))
r = np.diff(stock) + flux[1:] + burn[1:]
r_nee = np.diff(stock) + nee[1:] + burn[1:]
m = np.arange(1, len(stock))
out['cbal'] = dict(max_abs_resid=float(np.nanmax(np.abs(r))), mean_resid=float(np.nanmean(r)),
                   max_abs_resid_month=int(m[np.nanargmax(np.abs(r))]), max_abs_resid_nee_form=float(np.nanmax(np.abs(r_nee))),
                   n_months=int(len(r)), total_stock_change=float(stock[-1] - stock[0]))
out['ntr_years'] = ntr
json.dump(out, sys.stdout); print()
