#!/usr/bin/env python3
"""demo10 analysis: per-cell identity on vs off_auto (tr+sc outputs and restarts minus tk_*),
C closure per cell, RTS effect in the slump cell, restart tk_* presence."""
import os, sys, json, numpy as np, netCDF4 as nc
W = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.getcwd(); O = lambda n: os.path.join(W, n, 'output')
def arr(f, v):
    with nc.Dataset(f) as x: return np.ma.getdata(x[v][:]).astype(float)
res = {'cells_differ': {}, 'restart': {}}
files = sorted(f for f in os.listdir(O('on')) if f.endswith('.nc'))
for f in files:
    a, b = os.path.join(O('on'), f), os.path.join(O('off_auto'), f)
    if not os.path.exists(b): res.setdefault('only_in_on', []).append(f); continue
    with nc.Dataset(a) as x, nc.Dataset(b) as y:
        for v in x.variables:
            if v.startswith('tk_') or v not in y.variables or x[v].ndim < 2: continue
            A, B = np.ma.getdata(x[v][:]), np.ma.getdata(y[v][:])
            if A.shape != B.shape: res['cells_differ'].setdefault('shape', []).append(f + ':' + v); continue
            # cell axes are the last two for outputs (time,...,y,x) and the first two for restarts (y,x,...)
            if f.startswith('restart'):
                d = (A != B).reshape(A.shape[0], A.shape[1], -1).any(-1)
            else:
                d = (A != B).reshape(-1, A.shape[-2], A.shape[-1]).any(0)
            for (i, j) in np.argwhere(d):
                if A.shape[-2] == 10 or f.startswith('restart'):
                    res['cells_differ'].setdefault('%d,%d' % (i, j), set()).add(f.split('_')[0] if not f.startswith('restart') else f)
res['cells_differ'] = {k: sorted(v) if isinstance(v, set) else v for k, v in res['cells_differ'].items()}
for n in ('on', 'off_auto', 'off'):
    f = os.path.join(O(n), 'restart-tr.nc')
    if os.path.exists(f):
        with nc.Dataset(f) as x: res['restart'][n] = sorted(v for v in x.variables if v.startswith('tk_'))
with nc.Dataset(os.path.join(O('on'), 'restart-sc.nc')) as x:
    res['restart_on_sc_tk_init_row0'] = np.ma.getdata(x['tk_init'][0, :]).tolist()
    res['restart_on_sc_tk_drainage_row0'] = np.ma.getdata(x['tk_drainage_type'][0, :]).tolist()
def m(n, v, st='tr'):
    f = os.path.join(O(n), '%s_monthly_%s.nc' % (v, st))
    a = arr(f, v); a = a.reshape(a.shape[0], -1, 10, 10) if a.ndim > 3 else a[:, None]
    return np.nansum(np.where(np.abs(a) > 1e20, 0, a), axis=1)  # (t, y, x)
clos = {}
for n in ('on', 'off_auto'):
    st = m(n, 'VEGC') + m(n, 'SOC') + m(n, 'DEADC') + m(n, 'DWDC')
    nee = m(n, 'RECO') - m(n, 'GPP')
    lat = m(n, 'TKLATC') + m(n, 'TKRESPC') - m(n, 'TKREGROWC')
    rr = np.diff(st, axis=0) + (nee + lat)[1:]
    clos[n] = dict(max_abs_any_cell=float(np.max(np.abs(rr[:, 0, :]))), cell00_max=float(np.max(np.abs(rr[:, 0, 0]))))
res['C_closure_tr_row0'] = clos
tk = m('on', 'TKREMTHICK'); i = np.argwhere(tk[:, 0, 0] > 0)
res['slump_month_cell00'] = int(i[0, 0]) if len(i) else None
if len(i):
    k = int(i[0, 0])
    for v in ('SOC', 'VEGC'):
        e = m('on', v)[:, 0, 0] - m('off_auto', v)[:, 0, 0]
        res['cell00_slump_drop_' + v] = float(e[k] - e[k - 1]); res['cell00_end_effect_' + v] = float(e[-1])
    res['cell00_TKLATC_slump'] = float(m('on', 'TKLATC')[k, 0, 0]); res['cell00_TKRESPC_slump'] = float(m('on', 'TKRESPC')[k, 0, 0])
    g, gb = m('on', 'GPP')[:, 0, 0], m('off_auto', 'GPP')[:, 0, 0]
    res['cell00_gpp_last10yr_frac_of_off'] = float(g[-120:].sum() / gb[-120:].sum())
json.dump(res, open(os.path.join(W, 'demo10_summary.json'), 'w'), indent=1); print(json.dumps(res, indent=1))
