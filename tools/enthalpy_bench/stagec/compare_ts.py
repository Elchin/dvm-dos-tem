#!/usr/bin/env python3
"""First-look comparison legacy (A) vs thermal_solver = 1 (B), main stage outputs.
usage: compare_ts.py ADIR BDIR [--last 20]   (ADIR/BDIR contain <case>/main/output)"""
import sys, json, argparse
from pathlib import Path
import numpy as np, netCDF4 as nc
ap = argparse.ArgumentParser(); ap.add_argument("A"); ap.add_argument("B"); ap.add_argument("--last", type=int, default=20)
a = ap.parse_args()
def rd(d, v, f): return np.ma.filled(np.ma.masked_invalid(nc.Dataset(Path(d) / "main/output" / f"{v}_{f}_sc.nc")[v][:].astype(float)), np.nan)
def flat(x): return x.reshape(x.shape[0], -1)
def tot(x): x = flat(x); return np.nansum(np.where(np.abs(x) < 1e20, x, 0), 1)
DEP = (0.1, 0.3, 0.5, 1.0, 2.0, 5.0)
def tprof(d):
    T = flat(rd(d, "TLAYER", "monthly")); z = flat(rd(d, "LAYERDEPTH", "monthly")); dz = flat(rd(d, "LAYERDZ", "monthly"))
    out = np.full((T.shape[0], len(DEP)), np.nan)
    for m in range(T.shape[0]):
        ok = (dz[m] > 0) & (dz[m] < 100) & (np.abs(T[m]) < 100)
        zt, zb, t = z[m][ok], z[m][ok] + dz[m][ok], T[m][ok]
        for k, D in enumerate(DEP):
            i = np.where((zt <= D) & (zb > D))[0]
            if len(i): out[m, k] = t[i[0]]
    return out
def metrics(d, L):
    r = {}
    ald = flat(rd(d, "ALD", "yearly"))[:, 0]; r["ALD_m"] = float(np.nanmean(ald[-L:])); r["ALD_max_m"] = float(np.nanmax(ald[-L:]))
    g = tot(rd(d, "GPP", "monthly")); re = tot(rd(d, "RECO", "monthly")); ny = len(g) // 12
    G = g[: ny * 12].reshape(ny, 12).sum(1); N = (re - g)[: ny * 12].reshape(ny, 12).sum(1)
    r["GPP_gm2yr"] = float(G[-L:].mean()); r["NEE_gm2yr"] = float(N[-L:].mean()); r["NEE_cum_gm2"] = float(N.sum())
    tp = tprof(d)[-12 * L:]
    r["T_mean"] = {str(D): float(np.nanmean(tp[:, k])) for k, D in enumerate(DEP)}
    mon = tp.reshape(L, 12, -1)
    r["T_amp"] = {str(D): float(np.nanmean(np.nanmax(mon[:, :, k], 1) - np.nanmin(mon[:, :, k], 1))) for k, D in enumerate(DEP)}
    r["SOC_end"] = float(tot(rd(d, "SOC", "monthly"))[-1])
    return r
res = {}
cases = sorted(p.name for p in Path(a.B).iterdir() if (p / "main/output/TLAYER_monthly_sc.nc").exists())
hdr = f"{'case':24s} {'ALD A':>6s} {'ALD B':>6s} {'GPP A':>7s} {'GPP B':>7s} {'NEE A':>7s} {'NEE B':>7s} " + " ".join(f"T{D}A/B".rjust(13) for D in (0.1, 0.5, 1.0, 2.0)) + f" {'SOC dB-A':>9s}"
print(f"last {a.last} years of the main stage (1950-2024): ALD m, GPP/NEE gC/m2/yr, T degC mean, SOC gC/m2 at end"); print(hdr)
for c in cases:
    if not (Path(a.A) / c / "main/output").exists(): continue
    ma, mb = metrics(Path(a.A) / c, a.last), metrics(Path(a.B) / c, a.last); res[c] = dict(A=ma, B=mb)
    print(f"{c:24s} {ma['ALD_m']:6.2f} {mb['ALD_m']:6.2f} {ma['GPP_gm2yr']:7.1f} {mb['GPP_gm2yr']:7.1f} {ma['NEE_gm2yr']:7.2f} {mb['NEE_gm2yr']:7.2f} "
          + " ".join(f"{ma['T_mean'][str(D)]:6.2f}/{mb['T_mean'][str(D)]:6.2f}" for D in (0.1, 0.5, 1.0, 2.0)) + f" {mb['SOC_end'] - ma['SOC_end']:9.1f}")
json.dump(res, open(Path(a.B) / "compare_vs_legacy.json", "w"), indent=1)
