#!/usr/bin/env python3
"""Regression / closure checks for the phase-2 cohort runs (run_phase2.py tree).

Per run (ROOT/cNN/<run>/stage2b_{rts|baseline}/*_sc.nc):
  * NaN/inf scan of every output and restart-sc.nc
  * TKNETC identity vs components (NEE_month + TKLATC + TKRESPC - TKREGROWC, CH4
    soil term TKRHCO2+TKCH4EMIS when present) < 1e-9 -- same as test_rts_tknetc.py
  * carbon closure: monthly |dECOC + TKNETC| < 0.01 gC/m2/month and cumulative
    |sum TKNETC + (ECOC_end - ECOC_start)| < 1.0 gC/m2 (suite tolerances),
    ECOC = VEGC + SOC + DEADC + DWDC (+ TKCH4POOL)
  * N closure from nbal.txt (NBAL_DIAG=1): max |N_res| < 1e-3 gN/m2/month (check_stage1
    formula incl. the regrow N input column)
Per cell / folr (A=0 regressions):
  * base_f08 == base_f05 bit-for-bit (folr only enters at initiation)
  * notrig_fNN == base_fNN bit-for-bit (rts mode, all-zero mask)
  * yYYYY_fNN == base_fNN bit-for-bit for every month before Aug YYYY
    (yearly files: every year before YYYY), and differs from August YYYY on
Writes ROOT/checks.json (per-run rows + summary). Exit 1 on any failure.
"""
import argparse, glob, json, os, sys
import numpy as np
import netCDF4

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from rts_sus2case.check_stage1 import parse_nbal, n_residual, nan_vars  # noqa: E402

FIRST = 2016; AUG = 7
TOL_ID, TOL_M, TOL_CUM = 1e-9, 0.01, 1.0
TOL_N = 1e-3   # gN/m2/month, as the stage-1 gate (nbal.txt prints 6 decimals)


def read(d, v, ts='monthly'):
    p = os.path.join(d, f'{v}_{ts}_sc.nc')
    if not os.path.exists(p): return None
    with netCDF4.Dataset(p) as ds:
        return np.array(ds[v][:], dtype=float)


def outdir(run):
    g = glob.glob(os.path.join(run, 'stage2b_*'))
    return g[0] if g else None


def carbon(d):
    g = lambda v: (lambda a: None if a is None else a.reshape(a.shape[0], -1).sum(axis=1) if a.ndim > 1 else a)(read(d, v))
    ecoc = g('VEGC') + g('SOC') + g('DEADC') + g('DWDC')
    pool = g('TKCH4POOL')
    if pool is not None: ecoc = ecoc + pool
    soil = g('TKRHCO2') + g('TKCH4EMIS') if g('TKRHCO2') is not None else g('RHSOM')
    nee = g('RG') + g('RM') + soil + g('RHDWD') - g('GPP')
    expect = nee + g('TKLATC') + g('TKRESPC')
    if g('TKREGROWC') is not None: expect = expect - g('TKREGROWC')
    T = g('TKNETC')
    return closure(ecoc, T, expect)


def closure(ecoc, tknetc, expect):
    m = np.diff(ecoc) + tknetc[1:]
    return dict(id_max=float(np.max(np.abs(tknetc - expect))), c_month_max=float(np.max(np.abs(m))),
                c_cum=float(tknetc[1:].sum() + (ecoc[-1] - ecoc[0])))


def prefix_equal(a, b, n):
    """True if a[:n] and b[:n] are bit-identical (NaN-safe)."""
    return bool(np.array_equal(a[:n], b[:n], equal_nan=True))


def cutoff(year, ts):
    return (year - FIRST) * 12 + AUG if ts == 'monthly' else year - FIRST


def compare(a, b, year=None):
    """Compare all *_sc.nc outputs of two stage-2b dirs. year=None: full bit-for-bit.
    Else: bit-for-bit before the initiation cutoff. Returns (n_files, mismatches, any_diff_after)."""
    bad, after, n = [], False, 0
    for p in sorted(glob.glob(os.path.join(a, '*_sc.nc'))):
        f = os.path.basename(p)
        if f.startswith('restart') or f.startswith('run_status'): continue
        v, ts = f.split('_')[0], f.split('_')[1]
        x, y = read(a, v, ts), read(b, v, ts)
        if y is None: bad.append(f + ' missing'); continue
        n += 1
        k = len(x) if year is None else cutoff(year, ts)
        if not prefix_equal(x, y, k): bad.append(f)
        if year is not None and not np.array_equal(x[k:], y[k:], equal_nan=True): after = True
    return n, bad, after


def check_run(run):
    d = outdir(run); r = dict(run=os.path.relpath(run, os.path.dirname(os.path.dirname(run))))
    if d is None or not os.path.exists(os.path.join(run, 'timings.json')):
        r['error'] = 'not finished'; return r
    r.update(carbon(d))
    nb = os.path.join(run, 'nbal.txt')
    a, _ = parse_nbal(nb)
    if len(a): r['N_res_max'] = float(np.max(np.abs(n_residual(a)[0])))
    nan = {}
    for p in glob.glob(os.path.join(d, '*.nc')):
        b, _ = nan_vars(p)
        if b: nan[os.path.basename(p)] = b
    r['nan'] = nan
    r['ok'] = (r['id_max'] < TOL_ID and r['c_month_max'] < TOL_M and abs(r['c_cum']) < TOL_CUM and not nan
               and r.get('N_res_max', 0.0) < TOL_N)
    return r


def check_cell(cdir):
    out = dict(cell=os.path.basename(cdir), runs=[], ab=[])
    runs = sorted(x for x in os.listdir(cdir) if x.startswith(('base_', 'notrig_', 'y')) and os.path.isdir(os.path.join(cdir, x)))
    for x in runs: out['runs'].append(check_run(os.path.join(cdir, x)))
    od = lambda x: outdir(os.path.join(cdir, x))
    folrs = sorted({x.split('_')[-1] for x in runs})
    bases = {f: od('base_' + f) for f in folrs if 'base_' + f in runs}
    if len(bases) == 2:
        (f1, a), (f2, b) = sorted(bases.items())
        n, bad, _ = compare(a, b); out['ab'].append(dict(test=f'base_{f1}==base_{f2}', n=n, bad=bad, ok=not bad))
    for x in runs:
        f = x.split('_')[-1]
        if f not in bases or x.startswith('base_'): continue
        if x.startswith('notrig_'):
            n, bad, _ = compare(od(x), bases[f])
            out['ab'].append(dict(test=f'{x}==base_{f}', n=n, bad=bad, ok=not bad))
        else:
            y = int(x[1:5]); n, bad, after = compare(od(x), bases[f], y)
            out['ab'].append(dict(test=f'{x}==base_{f} before Aug {y}', n=n, bad=bad, differs_after=after,
                                  ok=not bad and after))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root'); ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--cells', default=None)
    a = ap.parse_args()
    cells = a.cells.split(',') if a.cells else sorted(x for x in os.listdir(a.root) if x.startswith('c') and os.path.isdir(os.path.join(a.root, x)))
    C = [check_cell(os.path.join(a.root, c)) for c in cells]
    R = [r for c in C for r in c['runs']]; AB = [t for c in C for t in c['ab']]
    fin = [r for r in R if 'error' not in r]
    S = dict(cells=len(C), runs=len(R), finished=len(fin), closure_fail=[r['run'] for r in fin if not r['ok']],
             unfinished=[r['run'] for r in R if 'error' in r],
             id_max=max((r['id_max'] for r in fin), default=None),
             c_month_max=max((r['c_month_max'] for r in fin), default=None),
             c_cum_absmax=max((abs(r['c_cum']) for r in fin), default=None),
             N_res_max=max((r.get('N_res_max', 0) for r in fin), default=None),
             runs_with_nan=[r['run'] for r in fin if r['nan']],
             ab_tests=len(AB), ab_fail=[c['cell'] + ':' + t['test'] for c in C for t in c['ab'] if not t['ok']])
    S['ok'] = not (S['closure_fail'] or S['unfinished'] or S['ab_fail'] or S['runs_with_nan'])
    json.dump(dict(summary=S, cells=C), open(a.out or os.path.join(a.root, 'checks.json'), 'w'), indent=1)
    print(json.dumps(S, indent=1)[:4000])
    sys.exit(0 if S['ok'] else 1)


if __name__ == '__main__':
    main()
