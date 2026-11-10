#!/usr/bin/env python3
"""Stage D: model vs observations (climatological; forcing is bias-corrected CanRCM4, so years are
not the observed weather). usage: stageD_eval.py OBS_JSON OUT_JSON LABEL=RUNROOT [LABEL=RUNROOT ...]
RUNROOT/<case>/main/output. Model window: last 15 years of the 1950-2024 main stage (2010-2024)."""
import sys, os, json, numpy as np, netCDF4 as nc
obs = json.load(open(sys.argv[1])); out = sys.argv[2]; runs = [a.split('=', 1) for a in sys.argv[3:]]
NY = 15
def rd(d, v):
    with nc.Dataset(os.path.join(d, v)) as f:
        k = [x for x in f.variables if x not in ('time', 'x', 'y', 'albers_conical_equal_area', 'lat', 'lon')][0]
        a = np.array(f.variables[k][:], dtype=float)
    return a[..., 0, 0] if a.ndim >= 3 else a
def tprof(d, depths):
    T = rd(d, 'TLAYER_monthly_sc.nc'); Z = rd(d, 'LAYERDEPTH_monthly_sc.nc'); DZ = rd(d, 'LAYERDZ_monthly_sc.nc')
    o = np.full((T.shape[0], len(depths)), np.nan)
    for t in range(T.shape[0]):
        ok = (DZ[t] > 0) & (DZ[t] < 100) & (T[t] > -100)
        zc = (Z[t] + 0.5 * DZ[t])[ok]; o[t] = np.where(np.array(depths) <= zc.max(), np.interp(depths, zc, T[t][ok]), np.nan)
    return o[-12 * NY:].reshape(NY, 12, -1)
def snow(d):
    swe = rd(d, 'SWE_daily_sc.nc').ravel(); n = len(swe) // 365; swe = swe[:n * 365].reshape(n, 365)[-NY:]
    # hydrological-year peak (Sep..Aug): peak of Sep(y-1)..Aug(y)
    flat = swe.ravel(); peaks = [flat[365 * i - 122: 365 * i + 243].max() for i in range(1, NY)]
    off = []
    for y in range(NY):
        s = swe[y]; c = 0; v = None
        for k in range(90, 196):
            c = c + 1 if s[k] <= 0.0 else 0
            if c == 5: v = k - 4 + 1; break
        if v: off.append(v)
    st = rd(d, 'SNOWTHICK_monthly_sc.nc').ravel()[-12 * NY:].reshape(NY, 12)
    return float(np.mean(peaks)), float(np.mean(off)) if off else None, float(st[:, 0:4].mean())
SEAS = {'DJF': [11, 0, 1], 'MAM': [2, 3, 4], 'JJA': [5, 6, 7], 'SON': [8, 9, 10]}
res = []
for lab, root in runs:
    for c in sorted(os.listdir(root)):
        d = os.path.join(root, c, 'main', 'output')
        if not os.path.isdir(d): continue
        site = c[:3]; ob = obs[site]; r = {'run': lab, 'case': c, 'site': site}
        ald = rd(d, 'ALD_yearly_sc.nc')[-NY:]; r['ALD'] = float(ald.mean()); r['ALD_bias'] = r['ALD'] - ob['ALD']['primary' if site == 'TVC' else 'thaw_tube']
        pk, off, sd = snow(d); r['peakSWE'] = pk; r['peakSWE_bias'] = pk - ob['peak_SWE_mm'][0]
        r['snowoff'] = off; r['snowoff_bias'] = (off - ob['snowoff_doy']['mean']) if off else None; r['snowdepth_JFMA'] = sd
        if site == 'TVC':
            dm = ob['soilT_monthly']['depth_m']; names = list(dm); tp = tprof(d, [dm[k] for k in names]).mean(0)   # 12 x nd
            for j, k in enumerate(names):
                o = np.array([ob['soilT_monthly']['values'][k][str(m + 1)][0] for m in range(12)])
                e = tp[:, j] - o; r[f'T{k}_bias'] = float(e.mean()); r[f'T{k}_rmse'] = float(np.sqrt((e ** 2).mean()))
                for s, ms in SEAS.items(): r[f'T{k}_{s}_bias'] = float(e[ms].mean())
        else:
            o = np.array([ob['surfT_monthly_3_5cm'][str(m + 1)] for m in range(12)])
            tp = tprof(d, [0.04]).mean(0)[:, 0]; e = tp - o
            r['T4cm_bias'] = float(e.mean()); r['T4cm_rmse'] = float(np.sqrt((e ** 2).mean()))
            for s, ms in SEAS.items(): r[f'T4cm_{s}_bias'] = float(e[ms].mean())
            env = ob['envelope_01TC02']; dz = {'0.5 m': 0.5, '1 m': 1.0, '1.5 m': 1.5, '2 m': 2.0, '3 m': 3.0}
            tpa = tprof(d, list(dz.values()))   # NY x 12 x nd (monthly means: damped extremes)
            for j, k in enumerate(dz):
                omax = np.mean([float(x['max']) for x in env[k]]); omin = np.mean([float(x['min']) for x in env[k]])
                r[f'env{k}_max_bias'] = float(np.nanmean(tpa[:, :, j].max(1)) - omax); r[f'env{k}_min_bias'] = float(np.nanmean(tpa[:, :, j].min(1)) - omin)
        res.append(r)
json.dump(res, open(out, 'w'), indent=1); print('wrote', out, len(res))
