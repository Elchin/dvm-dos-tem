#!/usr/bin/env python3
"""First-look legacy vs enthalpy comparison: cmp_ts.py LEGACY_RUNROOT ENTH_RUNROOT  (dirs with <case>/main/output)"""
import sys, os, json, numpy as np, netCDF4 as nc
L, E = sys.argv[1], sys.argv[2]
def rd(d, v):
    f = nc.Dataset(os.path.join(d, v)); var = [k for k in f.variables if k not in ('time','x','y','albers_conical_equal_area','lat','lon')][0]
    a = np.array(f.variables[var][:], dtype=float); f.close()
    a = a[..., 0, 0] if a.ndim >= 3 else a
    return a
def tdepth(d, depths):
    T = rd(d, 'TLAYER_monthly_sc.nc'); Z = rd(d, 'LAYERDEPTH_monthly_sc.nc'); DZ = rd(d, 'LAYERDZ_monthly_sc.nc')
    out = np.full((T.shape[0], len(depths)), np.nan)
    for t in range(T.shape[0]):
        ok = DZ[t] > 0
        zc = (Z[t] + 0.5 * DZ[t])[ok]; tt = T[t][ok]
        out[t] = np.interp(depths, zc, tt)
    return out
def ann(m):  # monthly -> annual sums
    n = m.shape[0] // 12; return m[:n*12].reshape(n, 12, *m.shape[1:]).sum(1)
cases = sorted(c for c in os.listdir(E) if os.path.isdir(os.path.join(E, c, 'main', 'output')))
dep = [0.2, 0.5, 1.0, 2.0]
rows = []
for c in cases:
    dl, de = os.path.join(L, c, 'main', 'output'), os.path.join(E, c, 'main', 'output')
    r = {'case': c}
    for tag, d in (('L', dl), ('E', de)):
        ald = rd(d, 'ALD_yearly_sc.nc'); r['ALD_'+tag] = ald[-30:].mean(); r['ALDall_'+tag] = ald.mean()
        g = rd(d, 'GPP_monthly_sc.nc'); g = g.reshape(g.shape[0], -1).sum(1) if g.ndim > 1 else g
        re = rd(d, 'RECO_monthly_sc.nc'); re = re.reshape(re.shape[0], -1).sum(1) if re.ndim > 1 else re
        r['GPP_'+tag] = ann(g)[-30:].mean(); r['NEE_'+tag] = (ann(re) - ann(g))[-30:].mean()
        td = tdepth(d, dep); r['T_'+tag] = td[-360:].mean(0).tolist()
        tdj = td[-360:].reshape(30, 12, -1); r['Tjan_'+tag] = tdj[:, 0].mean(0).tolist(); r['Tjul_'+tag] = tdj[:, 6].mean(0).tolist()
    rows.append(r)
json.dump(rows, open(E.rstrip('/') + '_cmp.json', 'w'), indent=1)
print('%-22s %7s %7s %6s | %7s %7s | %7s %7s | %s' % ('case (last 30 yr)', 'ALD_L', 'ALD_E', 'dALD', 'GPP_L', 'GPP_E', 'NEE_L', 'NEE_E', 'dT annual @0.2/0.5/1/2 m ; dT Jan ; dT Jul'))
for r in rows:
    dT = np.array(r['T_E']) - np.array(r['T_L']); dJ = np.array(r['Tjan_E']) - np.array(r['Tjan_L']); dU = np.array(r['Tjul_E']) - np.array(r['Tjul_L'])
    f = lambda a: '/'.join('%+.2f' % x for x in a)
    print('%-22s %7.3f %7.3f %+6.3f | %7.1f %7.1f | %+7.1f %+7.1f | %s ; %s ; %s' % (r['case'], r['ALD_L'], r['ALD_E'], r['ALD_E']-r['ALD_L'], r['GPP_L'], r['GPP_E'], r['NEE_L'], r['NEE_E'], f(dT), f(dJ), f(dU)))
