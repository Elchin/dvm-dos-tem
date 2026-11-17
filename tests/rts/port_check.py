#!/usr/bin/env python3
"""RTS port (model_settings.thermokarst) checks on the 1-cell synthetic case (testing-data/rts-synthetic).

  python3 tests/rts/port_check.py BINARY WORKDIR [--cases off,noinit,on,ice,on_fx,off_fx,ice_fx] [-P 4]

Runs dvmdostem (RUN_ARGS of rtslib: pr 10, eq 300, sp 30, tr 30, sc 0; tr outputs) for the cases
  off     thermokarst 0
  noinit  thermokarst 1 + RTS_PARAMS, thermokarst_initiation false (module on, no slump: must equal off)
  on      thermokarst 1 + RTS_PARAMS (v0.3 defaults: floor + regrow on)            ~ reference/v03/rts
  ice     on + ICE_PARAMS (excess ground ice)                                     ~ reference/v03/rts_ice
  *_fx    the same with this fork's layer/N switches on (layer_bookkeep_fix, layer_n_remap_fix,
          soil_negpool_fix, soil_giveup_retry = 1)
with NBAL_DIAG=1 (monthly N balance from model state, thermokarst = 1 only) and writes
WORKDIR/summary.json: identity off vs noinit, monthly C / N / water closure residuals, and the RTS
effect (on - off) per variable next to rts-regrow v0.3 (reference/v03/X - reference/baseline).
"""
import argparse, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np, netCDF4 as nc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rtslib as r

FX = {'layer_bookkeep_fix': 1, 'layer_n_remap_fix': 1, 'soil_negpool_fix': 1, 'soil_giveup_retry': 1}
WATER = {v: 'm' for v in ('RAINFALL', 'SNOWFALL', 'EET', 'QRUNOFF', 'QDRAINAGE', 'SWE', 'IWCLAYER', 'TKWDRAIN', 'TKWLAT')}
NOUT = {v: 'm' for v in ('VEGNTOT', 'AVLN', 'DEADN', 'DWDN', 'NINPUT', 'NLOST')}
OUTS = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.REGROW_OUTPUTS, **r.ICE_OUTPUTS, **WATER, **NOUT}
CASES = {
  'off':    (False, {}),
  'noinit': (True,  {'thermokarst_initiation': False}),
  'on':     (True,  {}),
  'ice':    (True,  dict(r.ICE_PARAMS)),
  'off_fx': (False, dict(FX)),
  'on_fx':  (True,  dict(FX)),
  'ice_fx': (True,  {**r.ICE_PARAMS, **FX}),
  # thermokarst = 1 with the three auto-enabled layer/N switches forced off (shows why they are on)
  # thermokarst = 0 with the three switches that thermokarst = 1 auto-enables (base for the RTS effect)
  'off_auto': (False, {'soil_negpool_fix': 1, 'layer_n_remap_fix': 1, 'layer_bookkeep_fix': 1}),
  # safety-net N conservation (n_safetynet_conserve) and/or soil_giveup_retry on top of the defaults
  'on_snc': (True, {'n_safetynet_conserve': 1}),
  'ice_snc': (True, {**r.ICE_PARAMS, 'n_safetynet_conserve': 1}),
  'off_snc': (False, {'soil_negpool_fix': 1, 'layer_n_remap_fix': 1, 'layer_bookkeep_fix': 1, 'n_safetynet_conserve': 1}),
  # global n_safetynet_conserve 2 (L1 + exact soil accumulation, base fix)
  'on_snc2': (True, {'n_safetynet_conserve': 2}),
  'ice_snc2': (True, {**r.ICE_PARAMS, 'n_safetynet_conserve': 2}),
  'off_snc2': (False, {'soil_negpool_fix': 1, 'layer_n_remap_fix': 1, 'layer_bookkeep_fix': 1, 'n_safetynet_conserve': 2}),
  'on_ta09': (True, {'rts_regrow_t_active': 0.9}),   # stabilization ~1 month earlier (tested, not adopted)
  # thermokarst = 1 without the post-initiation safety-net N conservation (thermokarst_safetynet_conserve)
  'on_nosnc': (True, {'thermokarst_safetynet_conserve': False}),
  'ice_nosnc': (True, {**r.ICE_PARAMS, 'thermokarst_safetynet_conserve': False}),
  'on_noauto': (True, {'soil_negpool_fix': 0, 'layer_n_remap_fix': 0, 'layer_bookkeep_fix': 0}),
}
REF = os.path.join(r.REPO, 'tests', 'rts', 'reference')

def run_case(binary, work, name):
    on, extra = CASES[name]
    wd = os.path.join(work, name)
    os.makedirs(wd, exist_ok=True)
    cfg = r.make_config(wd, on, outputs=OUTS, extra=extra or None)
    env = dict(os.environ, NBAL_DIAG='1')
    with open(os.path.join(wd, 'stderr.txt'), 'w') as fe:
        p = subprocess.run([binary, '-f', cfg, '-l', 'fatal'] + r.RUN_ARGS, cwd=r.REPO, env=env,
                           stdout=subprocess.DEVNULL, stderr=fe, timeout=3600)
    return name, p.returncode

def rd(d, v, ts='monthly'):
    f = os.path.join(d, '%s_%s_tr.nc' % (v, ts))
    if not os.path.exists(f): return None
    x = nc.Dataset(f); a = np.array(x[v][:], dtype=float); x.close()
    a = np.where(np.abs(a) > 1e20, np.nan, a); a = np.where(a <= -9999, np.nan, a)
    return a

def cell(a):  # (time, [layer|pft...], y, x) -> (time, ...) for the single cell
    return a.reshape(a.shape[0], -1) if a.ndim > 1 else a

def colsum(a): return np.nansum(cell(a), axis=1)

def c_closure(d):
    st = colsum(rd(d, 'VEGC')) + colsum(rd(d, 'SOC')) + colsum(rd(d, 'DEADC')) + colsum(rd(d, 'DWDC'))
    nee = colsum(rd(d, 'RECO')) - colsum(rd(d, 'GPP'))
    lat = colsum(rd(d, 'TKLATC')) if rd(d, 'TKLATC') is not None else 0 * nee
    rsp = colsum(rd(d, 'TKRESPC')) if rd(d, 'TKRESPC') is not None else 0 * nee
    rg = colsum(rd(d, 'TKREGROWC')) if rd(d, 'TKREGROWC') is not None else 0 * nee
    res = np.diff(st) + (nee + lat + rsp - rg)[1:]
    return res

def w_closure_wb(wd):
    """Exact water check of the RTS step from the NBAL_DIAG 'WB' lines (tr stage): (W0 - W1) - TKWDRAIN -
    TKWLAT per month, kg/m2; excess-ice initialization (a source, first RTS month) excluded."""
    rows = [l.split()[1:] for l in open(os.path.join(wd, 'stderr.txt')) if l.startswith('WB ')]
    if not rows: return None
    yr, m, w0, wp, wr, w1, dr = np.array(rows[-360:], dtype=float).T
    res = (w0 - w1) - dr - (wp - wr)
    return res[1:]

def n_closure_nb(wd):
    """Monthly N residual from the NBAL_DIAG lines (tr stage = last 360 lines)."""
    rows = [l.split()[1:] for l in open(os.path.join(wd, 'stderr.txt')) if l.startswith('NB ')]
    if not rows: return None
    a = np.array(rows[-360:], dtype=float)
    yr, m, n0, n1, ain, oin, alost, olost, tlat, tresp, vegin, rgn = a.T
    return (n1 - n0) - (ain + oin - alost - tlat - tresp + vegin + rgn)

def w_closure(d, drain=True):
    dz = cell(rd(d, 'LAYERDZ')); lwc = cell(rd(d, 'LWCLAYER')); iwc = cell(rd(d, 'IWCLAYER'))
    stor = np.nansum((np.nan_to_num(lwc) + np.nan_to_num(iwc)) * np.nan_to_num(dz) * 1000., axis=1) + colsum(rd(d, 'SWE'))
    P = colsum(rd(d, 'RAINFALL')) + colsum(rd(d, 'SNOWFALL'))
    out = colsum(rd(d, 'EET')) + colsum(rd(d, 'QRUNOFF')) + colsum(rd(d, 'QDRAINAGE'))
    for v in (('TKWDRAIN', 'TKWLAT', 'TKEXMELT') if drain else ()):   # RTS water terms (TKEXMELT: source)
        x = rd(d, v)
        if x is not None: out = out + (-1 if v == 'TKEXMELT' else 1) * colsum(x)
    return np.diff(stor) - (P - out)[1:]

def ident(a, b, skip=('TK', 'DRAINTYPE')):
    bad, n = [], 0
    for f in sorted(os.listdir(a)):
        if not f.endswith('_tr.nc') or f.startswith(skip): continue
        n += 1
        x, y = nc.Dataset(os.path.join(a, f)), nc.Dataset(os.path.join(b, f))
        v = f.split('_')[0]
        if not np.array_equal(np.array(x[v][:]), np.array(y[v][:])): bad.append(f)
        x.close(); y.close()
    return n, bad

EFFECT_VARS = ['VEGC', 'SOC', 'SHLWC', 'DEEPC', 'ORGN', 'GPP', 'NPP', 'RHSOM', 'RECO', 'DEADC', 'DWDC', 'ALD',
               'MOSSDZ', 'SHLWDZ', 'DEEPDZ', 'TKLATC', 'TKRESPC', 'TKLATN', 'TKREMTHICK', 'TKICEMELT', 'DRAINTYPE',
               'TKSUBSID', 'TKEXICE', 'TKWDRAIN', 'TKWLAT']

def series(d, v):
    ts = r.BASE_OUTPUTS.get(v, 'm')
    a = rd(d, v, 'yearly' if ts == 'y' else 'monthly')
    if a is None: return None
    s = colsum(a) if a.ndim > 1 and cell(a).shape[1] > 1 and v not in ('ALD',) else cell(a)[:, 0]
    if ts == 'm' and v not in ('TKLATC', 'TKRESPC', 'TKLATN', 'TKREMTHICK', 'TKICEMELT', 'TKWDRAIN', 'TKWLAT'):
        s = s.reshape(-1, 12).mean(1)      # stocks/rates: yearly mean of the monthly values
    elif ts == 'm':
        s = s.reshape(-1, 12).sum(1)       # RTS pulses: yearly sums
    return s

def effect(on, off, ref_on, ref_off):
    out = {}
    for v in EFFECT_VARS:
        a, b, ra, rb = series(on, v), series(off, v), series(ref_on, v), series(ref_off, v)
        if a is None or ra is None: continue
        if b is None: b = 0 * a
        if rb is None: rb = 0 * ra
        e, re = a - b, ra - rb
        out[v] = dict(ours_end=float(e[-1]), v03_end=float(re[-1]), diff_end=float(e[-1] - re[-1]),
                      ours_mean=float(np.nanmean(e)), v03_mean=float(np.nanmean(re)),
                      ours_max_abs=float(np.nanmax(np.abs(e))), v03_max_abs=float(np.nanmax(np.abs(re))),
                      ours_off_end=float(b[-1]), v03_off_end=float(rb[-1]))
        if v == 'GPP':
            out[v]['ours_ratio_yearly'] = (a / b).tolist(); out[v]['v03_ratio_yearly'] = (ra / rb).tolist()
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('binary'); ap.add_argument('work')
    ap.add_argument('--cases', default=','.join(CASES)); ap.add_argument('-P', type=int, default=4)
    ap.add_argument('--no-run', action='store_true')
    a = ap.parse_args()
    names = a.cases.split(',')
    if not a.no_run:
        with ThreadPoolExecutor(a.P) as ex:
            for n, rc in ex.map(lambda n: run_case(os.path.abspath(a.binary), os.path.abspath(a.work), n), names):
                print('run', n, 'rc', rc); assert rc == 0, n
    O = lambda n: os.path.join(a.work, n, 'output')
    S = {'cases': {n: CASES[n][1] for n in names}, 'rts_params': r.RTS_PARAMS, 'run_args': r.RUN_ARGS}
    if 'off' in names and 'noinit' in names:
        n, bad = ident(O('off'), O('noinit')); S['ident_off_vs_noinit'] = dict(files=n, differ=bad)
    if 'off_auto' in names and 'noinit' in names:
        n, bad = ident(O('off_auto'), O('noinit')); S['ident_off_auto_vs_noinit'] = dict(files=n, differ=bad)
    S['closure'] = {}
    for n in names:
        c = c_closure(O(n)); w = w_closure(O(n)); nb = n_closure_nb(os.path.join(a.work, n))
        P = colsum(rd(O(n), 'RAINFALL')) + colsum(rd(O(n), 'SNOWFALL'))
        ent = dict(C_max_abs_month=float(np.nanmax(np.abs(c))), C_sum=float(np.nansum(c)),
                   W_max_abs_month=float(np.nanmax(np.abs(w))), W_sum=float(np.nansum(w)),
                   W_sum_pctP=float(100 * np.nansum(w) / np.nansum(P[1:])))
        w0 = w_closure(O(n), drain=False); wb = w_closure_wb(os.path.join(a.work, n))
        if wb is not None: ent.update(WB_rts_step_max_abs_month=float(np.max(np.abs(wb))))
        ent.update(W_nodrain_max_abs_month=float(np.nanmax(np.abs(w0))), W_nodrain_sum=float(np.nansum(w0)))
        if nb is not None: ent.update(N_max_abs_month=float(np.max(np.abs(nb))), N_sum=float(np.sum(nb)))
        ent['init_month'] = None
        tk = rd(O(n), 'TKREMTHICK')
        if tk is not None and np.nanmax(cell(tk)) > 0: ent['init_month'] = int(np.nanargmax(cell(tk)[:, 0] > 0))
        S['closure'][n] = ent
        np.savez(os.path.join(a.work, n, 'closure.npz'), C=c, W=w, N=nb if nb is not None else np.zeros(0))
    S['effect'] = {}
    for on, off, ref in (('on', 'off', 'v03/rts'), ('ice', 'off', 'v03/rts_ice'), ('on_fx', 'off_fx', 'v03/rts'),
                         ('ice_fx', 'off_fx', 'v03/rts_ice'), ('on_noauto', 'off', 'v03/rts'),
                         ('on', 'off_auto', 'v03/rts'), ('ice', 'off_auto', 'v03/rts_ice'),
                         ('on_ta09', 'off_auto', 'v03/rts'), ('on_nosnc', 'off_auto', 'v03/rts')):
        if on in names and off in names:
            S['effect'][on if off != 'off_auto' else on + '-off_auto'] = effect(O(on), O(off), os.path.join(REF, ref), os.path.join(REF, 'baseline'))
    json.dump(S, open(os.path.join(a.work, 'summary.json'), 'w'), indent=1)
    print(json.dumps({k: S[k] for k in ('ident_off_vs_noinit', 'ident_off_auto_vs_noinit', 'closure') if k in S}, indent=1))

if __name__ == '__main__':
    main()
