#!/usr/bin/env python3
"""Mass-closure and run-health checks for one regression run (main stage, 1950-2024).
Formulas as /workspace/v5h/seb_explore/scripts/metrics_tvc.py:
  water (1951-2024): P - EET - QRUNOFF - QDRAINAGE - dS, dS = d(sum (LWC+IWC)*dz*1000) + dSWE  [mm]
                     PASS if |residual| < 1 % of P
  carbon (yearly):   d(SOC + VEGC + DEADC + DWDC) + NEE  [gC/m2/yr], NEE = RECO - GPP
  nitrogen (yearly): d(VEGNTOT + ORGN + AVLN + DEADN + DWDN) - (NINPUT - NLOST)  [gN/m2/yr]
                     PASS if |mean residual| < ntol (default 1.0 gN/m2/yr: a coarse absolute
                     bound; summarize.py also requires B to be no worse than A).
                     NLOST must be the AVLN leaching loss only (v5i); the old NLOST, which added the
                     never-removed baseline 'orgnlost' every month, gave residuals of order
                     10-1000s gN/m2/yr. NOTE: the N budget does not close exactly from the
                     outputs in the nfeed calibration cases: mean residual about -0.22 (HPC) and
                     +0.01..+0.05 (TVC) gN/m2/yr, single years up to ~0.9, identical in v5h and v5i
                     (unbooked N flux; not diagnosed here). Hence the coarse absolute bound.
usage: closure.py RUNDIR [--wtol 1.0] [--ctol 0.5] [--ntol 1.0]   (prints JSON; exit 0 if PASS)"""
import argparse, json, sys
from pathlib import Path
import numpy as np, netCDF4 as nc

def ma(x): return np.ma.filled(np.ma.masked_invalid(np.ma.asarray(x, dtype=float)), np.nan)
def tot(a): a = a.reshape(a.shape[0], -1); return np.nansum(np.where(np.abs(a) < 1e20, a, 0), 1)

def closure(rundir):
    o = Path(rundir) / "main/output"
    g = lambda v, f="monthly": ma(nc.Dataset(o / f"{v}_{f}_sc.nc")[v][:])
    m = lambda v: tot(g(v))
    yv = lambda v: tot(g(v, "yearly"))
    have = lambda v, f="yearly": (o / f"{v}_{f}_sc.nc").exists()
    nm = len(m("GPP")); ny = nm // 12; r = dict(n_years=ny)
    T = g("TLAYER").reshape(nm, -1); dz = g("LAYERDZ").reshape(nm, -1)
    lw = g("LWCLAYER").reshape(nm, -1); iw = g("IWCLAYER").reshape(nm, -1)
    ok = (dz > 0) & (dz < 100) & (T > -100)
    sw = g("SWE", "daily").ravel()[: ny * 365].reshape(ny, 365)
    stor = np.nansum(np.where(ok, (lw + iw) * dz * 1000., 0.), 1)
    dS = (stor[nm - 1] - stor[11]) + (sw[-1, -1] - sw[0, -1]); sl = slice(12, nm)
    P = m("RAINFALL")[sl].sum() + m("SNOWFALL")[sl].sum()
    res = P - m("EET")[sl].sum() - m("QRUNOFF")[sl].sum() - m("QDRAINAGE")[sl].sum() - dS
    r["P_mm_yr"] = P / (ny - 1); r["W_resid_mm_yr"] = res / (ny - 1); r["W_resid_pctP"] = 100 * res / P
    socy = m("SOC")[11::12]; st = socy + yv("VEGC")
    r["C_terms"] = ["SOC", "VEGC"] + [v for v in ("DEADC", "DWDC") if have(v)]
    for v in ("DEADC", "DWDC"):
        if have(v): st = st + yv(v)
    nee_y = (m("RECO") - m("GPP")).reshape(ny, 12).sum(1)
    rc = np.diff(st) + nee_y[1:]
    r["C_resid_mean_gm2yr"] = float(rc.mean()); r["C_resid_max_abs_annual"] = float(np.abs(rc).max())
    if all(have(v) for v in ("VEGNTOT", "ORGN", "AVLN", "NINPUT", "NLOST")):
        Nst = yv("VEGNTOT") + yv("ORGN") + yv("AVLN"); r["N_terms"] = ["VEGNTOT", "ORGN", "AVLN"]
        for v in ("DEADN", "DWDN"):
            if have(v): Nst = Nst + yv(v); r["N_terms"].append(v)
        nin, nl = yv("NINPUT"), yv("NLOST")
        rn = np.diff(Nst) - (nin[1:] - nl[1:])
        r["N_resid_mean_gm2yr"] = float(rn.mean()); r["N_resid_max_abs_annual"] = float(np.abs(rn).max())
        r["NLOST_mean_gm2yr"] = float(nl[1:].mean()); r["NINPUT_mean_gm2yr"] = float(nin[1:].mean())
    else:
        r["N_resid_mean_gm2yr"] = None
    return r

def evaluate(rundir, wtol=1.0, ctol=0.5, ntol=1.0):
    meta = json.loads((Path(rundir) / "meta.json").read_text())
    r = dict(run=str(rundir), case=meta["case"], status=meta["status"])
    for s in ("spin", "main"):
        if s in meta:
            ck = meta[s]; r[f"{s}_run_status"] = ck["run_status"]; r[f"{s}_fail_logs"] = len(ck["fail_logs"])
            r[f"{s}_safetynet"] = max(ck["soil_safetynet_events"], ck["safetynet_csv_rows"]); r[f"{s}_veg_giveups"] = ck["veg_giveups"]
    fails = []
    if meta["status"] != "OK": fails.append("status " + meta["status"])
    else:
        r.update(closure(rundir))
        if abs(r["W_resid_pctP"]) >= wtol: fails.append(f"water |resid| {r['W_resid_pctP']:.3f}% P >= {wtol}%")
        if r["C_resid_max_abs_annual"] >= ctol: fails.append(f"C max |annual resid| {r['C_resid_max_abs_annual']:.3g} >= {ctol}")
        if r["N_resid_mean_gm2yr"] is None: fails.append("N closure outputs missing")
        elif abs(r["N_resid_mean_gm2yr"]) >= ntol: fails.append(f"N |mean resid| {r['N_resid_mean_gm2yr']:.3g} >= {ntol} gN/m2/yr")
    r["closure_fail"] = fails; r["closure_pass"] = not fails
    return r

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--wtol", type=float, default=1.0)
    ap.add_argument("--ctol", type=float, default=0.5); ap.add_argument("--ntol", type=float, default=1.0)
    a = ap.parse_args(); r = evaluate(a.rundir, a.wtol, a.ctol, a.ntol); print(json.dumps(r))
    sys.exit(0 if r["closure_pass"] else 1)
