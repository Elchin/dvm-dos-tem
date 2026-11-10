#!/usr/bin/env python3
"""Mass-closure and run-health checks for one regression run (main stage, 1950-2024).
  energy (phase 2b, only when main/output/energy_daily.csv exists, i.e. energy_diag = 1):
                     daily thermal-step residual (h1-h0) - (gtop+gbot), see energy(); FAIL if
                     |run mean| >= etol (default 3 W/m2, phase 2c); summarize.py WARNs if B is worse than A.
Formulas as /workspace/v5h/seb_explore/scripts/metrics_tvc.py:
  water (1951-2024): P - EET - QRUNOFF - QDRAINAGE - dS, dS = d(sum (LWC+IWC)*dz*1000) + dSWE  [mm]
                     PASS if |residual| < 1 % of P
  carbon (yearly):   d(SOC + VEGC + DEADC + DWDC) + NEE  [gC/m2/yr], NEE = RECO - GPP
  nitrogen (yearly): d(VEGNTOT + ORGN + AVLN + DEADN + DWDN) - (NINPUT [+ BNF_VEG] - NLOST)  [gN/m2/yr]
                     v5q: NINPUT includes the soil part of biological N fixation (BNF output);
                     BNF_VEG (bnf_to_plant = 1: BNF given directly to plant labile N) is an
                     extra input, added when its yearly (else monthly, summed) file exists.
                     PASS if |mean residual| < ntol (default 1.0 gN/m2/yr: a coarse absolute
                     bound; summarize.py also requires B to be no worse than A).
                     v5o: the N tolerance is enforced (pass/fail) ONLY when all three N
                     conservation fixes are on in the run (n_rootless_fix 1|2, n_som_move_fix 1,
                     n_nonvasc_extract 1; read from main/config/config.js, missing keys take the
                     binary default: v5o = 2/1/1, or legacy 0/0/0 with --n-defaults legacy or
                     REG_N_DEFAULTS=legacy for pre-v5o binaries). Otherwise the N residual is
                     only reported as a warning (closure_warn) naming the three bugs.
                     NLOST must be the AVLN leaching loss only (v5i); the old NLOST, which added the
                     never-removed baseline 'orgnlost' every month, gave residuals of order
                     10-1000s gN/m2/yr. NOTE: the N budget does not close exactly from the
                     outputs in the nfeed calibration cases with the legacy N code: mean residual about -0.22
                     (HPC) and +0.01..+0.05 (TVC) gN/m2/yr, single years up to ~0.9 (6 in the
                     2025-2100 part). Diagnosed (v5n NFEED_RESIDUAL_REPORT.md): a real N mass error
                     from the known deltastate/uptake bugs N1 (rootless-layer netnmin lost), N2 (SOM
                     N move overwritten, sink) and N3 (moss N uptake not extracted, source); with
                     all three fixes on it closes to <= 2e-4 gN/m2/yr (float output round-off).
  enthalpy (Stage C, only when energy_daily.csv has the 'resday' column, i.e. thermal_solver = 1):
                     full-day energy balance of the enthalpy solver. Per day resday = res + equil + unattr
                     [J/m2] where res is the solve residual, equil the H-conserving equilibration of
                     off-curve cells (0 up to round-off) and unattr = (H at the start of the thermal step -
                     H at the end of the previous one) - (booked enthalpy changes of the snow, soil water,
                     layer/BGC and fire windows). PASS iff max_day |res|/86400 <= etol_enth AND
                     max_day |resday|/86400 <= etol_enth (default 1e-6 W/m2) AND no step needed the
                     dt-halving fallback beyond its limit (fallbacks column: number of failed steps).
usage: closure.py RUNDIR [--wtol 1.0] [--ctol 0.5] [--ntol 1.0] [--etol-enth 1e-6] [--n-defaults v5o|legacy]   (prints JSON; exit 0 if PASS)"""
import argparse, json, os, re, sys
from pathlib import Path
import numpy as np, netCDF4 as nc

def ma(x): return np.ma.filled(np.ma.masked_invalid(np.ma.asarray(x, dtype=float)), np.nan)
def tot(a): a = a.reshape(a.shape[0], -1); return np.nansum(np.where(np.abs(a) < 1e20, a, 0), 1)

N_FIX_DEFAULTS = {"v5o": dict(n_rootless_fix=2, n_som_move_fix=1, n_nonvasc_extract=1),
                  "legacy": dict(n_rootless_fix=0, n_som_move_fix=0, n_nonvasc_extract=0)}
N_BUGS = ("N1 n_rootless_fix (rootless-layer net mineralisation N lost), N2 n_som_move_fix (SOM N moved "
          "shallow->deep overwritten: N sink), N3 n_nonvasc_extract (non-vascular N uptake not extracted "
          "from the soil: N source)")

def n_fixes(rundir, defaults=None):
    """Effective N1-N3 switch values of a run (main/config/config.js; missing keys -> binary defaults)."""
    defaults = defaults or os.environ.get("REG_N_DEFAULTS", "v5o")
    eff = dict(N_FIX_DEFAULTS[defaults]); src = {k: "default:" + defaults for k in eff}
    c = Path(rundir) / "main/config/config.js"
    if c.exists():
        txt = re.sub(r"(?m)//.*$", "", c.read_text())
        ms = json.loads(txt).get("model_settings", {})
        for k in eff:
            if k in ms: eff[k] = int(ms[k]); src[k] = "config"
    on = eff["n_rootless_fix"] in (1, 2) and eff["n_som_move_fix"] == 1 and eff["n_nonvasc_extract"] == 1
    return eff, src, on

def n_budget(nst, nin, nl, bnfveg=None):
    """Yearly N residual d(N stocks) - (NINPUT [+ BNF_VEG] - NLOST), years 2..n (v5q: BNF_VEG optional)."""
    if bnfveg is not None: nin = nin + bnfveg
    return np.diff(nst) - (nin[1:] - nl[1:])

def energy(csvpath):
    """Phase 2b energy closure of the daily snow/soil thermal step from energy_daily.csv
    (written when model_settings.energy_diag = 1). Per day: res = (h1 - h0) - (gtop + gbot) - snowmk [J/m2]
    (phase 2d: snowmk = snow-melt marker set inside Stefan, energy supplied by Snow_Env; older files
    without the column: res as written),
    h = column enthalpy relative to liquid water at 0 C (solid and water sensible heat, ice carries
    -L_f), gtop/gbot = heat through the solver's virtual top/bottom boundary nodes (+ = into the
    column). Reported in W/m2 (/86400). Not included: heat advected with water (infiltration,
    drainage, ET), snowfall/snowmelt mass changes and layer restructuring between thermal steps
    (reported separately as E_other_mean_Wm2, not a closure term)."""
    import csv
    rows = list(csv.DictReader(open(csvpath)))
    if not rows: return dict(E_days=0)
    g = lambda k: np.array([float(x[k]) for x in rows])
    res, gt, gb, nf, sn, h0, h1 = g("res") / 86400., g("gtop") / 86400., g("gbot") / 86400., g("nfront"), g("snow"), g("h0"), g("h1")
    ny = len(res) // 365; r = dict(E_days=len(res))
    r["E_resid_mean_Wm2"] = float(res.mean())
    r["E_resid_rms_daily_Wm2"] = float(np.sqrt((res ** 2).mean()))
    r["E_resid_max_abs_annual_Wm2"] = float(np.abs(res[: ny * 365].reshape(ny, 365).mean(1)).max()) if ny else float("nan")
    r["E_G_mean_Wm2"] = float((gt + gb).mean()); r["E_G_mean_abs_Wm2"] = float(np.abs(gt).mean())
    r["E_gbot_mean_Wm2"] = float(gb.mean())
    r["E_resid_rel"] = float(np.abs(res).sum() / max(np.abs(gt).sum(), 1e-30))
    for lab, m in (("nofront", nf == 0), ("front", nf > 0), ("snow", sn == 1), ("nosnow", sn == 0)):
        r[f"E_resid_mean_{lab}_Wm2"] = float(res[m].mean()) if m.any() else None; r[f"E_days_{lab}"] = int(m.sum())
    r["E_other_mean_Wm2"] = float(((h0[1:] - h1[:-1]) / 86400.).mean()) if len(res) > 1 else 0.0
    r["E_dH_run_Wm2"] = float((h1[-1] - h0[0]) / 86400. / len(res))
    if "snowmk" in rows[0]: r["E_snowmk_mean_Wm2"] = float((g("snowmk") / 86400.).mean())
    if "resday" in rows[0]:   # Stage C enthalpy solver (thermal_solver = 1)
        rd = g("resday") / 86400.
        r["Enth_check"] = "on"
        r["Enth_res_max_abs_Wm2"] = float(np.abs(res).max())
        r["Enth_resday_max_abs_Wm2"] = float(np.abs(rd).max())
        r["Enth_resday_mean_Wm2"] = float(rd.mean())
        r["Enth_unattr_max_abs_Wm2"] = float(np.abs(g("unattr") / 86400.).max())
        r["Enth_fallback_failures"] = int(g("fallbacks").sum())
        r["Enth_solves"] = int(g("solves").sum())
        for k in ("src", "clip", "init", "equil", "bsnow", "bwater", "blayer", "bfire"):
            r[f"Enth_{k}_mean_Wm2"] = float((g(k) / 86400.).mean())
        if "rsnow" in rows[0]:   # Stage C2 energy-conserving remap
            for k in ("rsnow", "rwater"):
                r[f"Enth_{k}_mean_Wm2"] = float((g(k) / 86400.).mean())
                r[f"Enth_{k}_max_abs_Wm2"] = float(np.abs(g(k) / 86400.).max())
            r["Enth_wlayer_max_abs_kg"] = float(np.abs(g("wlayer")).max())
            r["Enth_remap_fallbacks"] = int(g("remapfb").max())
    else:
        r["Enth_check"] = "off"
    return r

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
        r["N_inputs"] = ["NINPUT"]; bv = None
        if have("BNF_VEG") or have("BNF_VEG", "monthly"):   # v5q Phase 2: BNF to plant labile N bypasses NINPUT
            bv = yv("BNF_VEG") if have("BNF_VEG") else m("BNF_VEG").reshape(ny, 12).sum(1)
            r["N_inputs"].append("BNF_VEG")
            r["BNF_VEG_mean_gm2yr"] = float(bv[1:].mean())
        rn = n_budget(Nst, nin, nl, bv)
        r["N_resid_mean_gm2yr"] = float(rn.mean()); r["N_resid_max_abs_annual"] = float(np.abs(rn).max())
        r["NLOST_mean_gm2yr"] = float(nl[1:].mean()); r["NINPUT_mean_gm2yr"] = float(nin[1:].mean())
    else:
        r["N_resid_mean_gm2yr"] = None
    ef = o / "energy_daily.csv"
    r["E_check"] = "on" if ef.exists() else "off"
    if ef.exists(): r.update(energy(ef))
    return r

def evaluate(rundir, wtol=1.0, ctol=0.5, ntol=1.0, n_defaults=None, etol=3.0, etol_enth=1e-6):
    meta = json.loads((Path(rundir) / "meta.json").read_text())
    r = dict(run=str(rundir), case=meta["case"], status=meta["status"])
    for s in ("spin", "main"):
        if s in meta:
            ck = meta[s]; r[f"{s}_run_status"] = ck["run_status"]; r[f"{s}_fail_logs"] = len(ck["fail_logs"])
            r[f"{s}_safetynet"] = max(ck["soil_safetynet_events"], ck["safetynet_csv_rows"]); r[f"{s}_veg_giveups"] = ck["veg_giveups"]
    fails = []; warns = []
    eff, src, n_on = n_fixes(rundir, n_defaults); r["N_fixes"] = eff; r["N_fixes_source"] = src
    r["N_check"] = "enforced" if n_on else "warn"
    if meta["status"] != "OK": fails.append("status " + meta["status"])
    else:
        r.update(closure(rundir))
        if abs(r["W_resid_pctP"]) >= wtol: fails.append(f"water |resid| {r['W_resid_pctP']:.3f}% P >= {wtol}%")
        if r["C_resid_max_abs_annual"] >= ctol: fails.append(f"C max |annual resid| {r['C_resid_max_abs_annual']:.3g} >= {ctol}")
        if r.get("Enth_check") == "on":   # Stage C: strict full-day balance replaces the legacy mean check
            for k in ("Enth_res_max_abs_Wm2", "Enth_resday_max_abs_Wm2"):
                if not np.isfinite(r[k]): r[k] = float("inf")   # NaN never passes
            if r["Enth_res_max_abs_Wm2"] > etol_enth:
                fails.append(f"enthalpy max |solve resid| {r['Enth_res_max_abs_Wm2']:.3g} W/m2 > {etol_enth}")
            if r["Enth_resday_max_abs_Wm2"] > etol_enth:
                fails.append(f"enthalpy max |full-day resid| {r['Enth_resday_max_abs_Wm2']:.3g} W/m2 > {etol_enth}")
            if r["Enth_fallback_failures"] > 0:
                fails.append(f"enthalpy solver failed to converge on {r['Enth_fallback_failures']} steps")
        elif r.get("E_resid_mean_Wm2") is not None and abs(r["E_resid_mean_Wm2"]) >= etol:
            fails.append(f"energy |mean thermal-step resid| {r['E_resid_mean_Wm2']:.3f} W/m2 >= {etol}")
        if r["N_resid_mean_gm2yr"] is None: fails.append("N closure outputs missing")
        elif n_on:
            if abs(r["N_resid_mean_gm2yr"]) >= ntol: fails.append(f"N |mean resid| {r['N_resid_mean_gm2yr']:.3g} >= {ntol} gN/m2/yr")
        else:   # legacy N code: known mass errors, report only
            warns.append(f"N closure not enforced (N1-N3 fixes not all on: {eff}); N mean resid "
                         f"{r['N_resid_mean_gm2yr']:.3g}, max |annual| {r['N_resid_max_abs_annual']:.3g} gN/m2/yr; "
                         f"known N mass-conservation bugs: {N_BUGS}")
    r["closure_fail"] = fails; r["closure_warn"] = warns; r["closure_pass"] = not fails
    return r

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--wtol", type=float, default=1.0)
    ap.add_argument("--ctol", type=float, default=0.5); ap.add_argument("--ntol", type=float, default=1.0)
    ap.add_argument("--etol", type=float, default=3.0)
    ap.add_argument("--etol-enth", type=float, default=1e-6)
    ap.add_argument("--n-defaults", choices=sorted(N_FIX_DEFAULTS), default=None,
                    help="defaults of the binary for missing N1-N3 keys (default $REG_N_DEFAULTS or v5o)")
    a = ap.parse_args(); r = evaluate(a.rundir, a.wtol, a.ctol, a.ntol, a.n_defaults, a.etol, a.etol_enth); print(json.dumps(r))
    sys.exit(0 if r["closure_pass"] else 1)
