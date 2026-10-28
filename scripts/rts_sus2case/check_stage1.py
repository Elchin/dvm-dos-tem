#!/usr/bin/env python3
"""Checks for the stage-1 (pr/eq/sp, RTS off) runs written by run_stage1.py.

Per cell: run finished (timings.json, run.log), restart-{pr,eq,sp}.nc present,
no NaN/inf in any numeric restart variable, run_status, and the monthly N budget
from NBAL_DIAG=1 (nbal.txt; columns yr m N0 N1 avlninput orgninput avlnlost
orgnlost tk_lat tk_resp vegin regrow_n). Residual as in spatial10x10/analyze.py:
  N_res = (N1-N0) - (avlnin + orgnin + vegin - avlnlost - tk_lat - tk_resp)
(orgnlost reported separately; it is not part of analyze.py's residual).
The year counter restarts at 0 in each stage, which splits pr/eq/sp.
Give-up counters are only printed with soil_giveup_diag=true (kept off), so
they are not available here.
"""
import argparse, glob, json, os
import numpy as np

STAGES = ('pr', 'eq', 'sp')


def parse_nbal(path):
    rows = [l.split()[1:] for l in open(path) if l.startswith('NB ')]
    a = np.array(rows, float) if rows else np.zeros((0, 12))
    seg = np.concatenate([[0], np.cumsum((np.diff(a[:, 0]) < 0) & True)]) if len(a) else np.zeros(0, int)
    return a, seg


def n_residual(a):
    _, _, n0, n1, ai, oi, al, ol, tl, tv, vi = a.T[:11]
    return (n1 - n0) - (ai + oi + vi - al - tl - tv), ol


def nan_vars(path):
    import netCDF4
    bad = {}
    with netCDF4.Dataset(path) as ds:
        for k, v in ds.variables.items():
            if v.dtype.kind not in 'fiu': continue
            x = np.ma.filled(v[:].astype(float), 0.0)
            n = int((~np.isfinite(x)).sum())
            if n: bad[k] = n
        nvar = len(ds.variables)
    return bad, nvar


def check_cell(d):
    c = os.path.basename(d); r = dict(cell=c)
    t = os.path.join(d, 'timings.json')
    r['finished'] = os.path.exists(t)
    r['stage1_s'] = json.load(open(t)).get('stage1') if r['finished'] else None
    for s in STAGES:
        r[f'restart_{s}'] = os.path.exists(os.path.join(d, 'stage1', f'restart-{s}.nc'))
    sp = os.path.join(d, 'stage1', 'restart-sp.nc')
    if r['restart_sp']:
        r['nan_vars'], r['n_restart_vars'] = nan_vars(sp)
    nb = os.path.join(d, 'nbal.txt')
    if os.path.exists(nb):
        a, seg = parse_nbal(nb)
        r['nbal_months'] = int(len(a))
        res, ol = n_residual(a)
        for k, s in enumerate(STAGES):
            m = seg == k
            if m.any():
                r[f'N_res_max_{s}'] = float(np.abs(res[m]).max())
                r[f'N_drift_{s}'] = float(a[m][-1, 3] - a[m][0, 2])
                r[f'orgnlost_sum_{s}'] = float(ol[m].sum())
        r['nbal_nonNB_lines'] = sum(1 for l in open(nb) if not l.startswith('NB '))
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root'); ap.add_argument('-o', '--out')
    a = ap.parse_args()
    cells = sorted(d for d in glob.glob(os.path.join(a.root, 'c*')) if os.path.isdir(d))
    R = [check_cell(d) for d in cells]
    T = np.array([r['stage1_s'] for r in R if r['stage1_s']], float)
    summ = dict(n_cells=len(R), finished=sum(r['finished'] for r in R),
                restarts_sp=sum(r['restart_sp'] for r in R),
                restarts_all=sum(all(r[f'restart_{s}'] for s in STAGES) for r in R),
                cells_with_nan=[r['cell'] for r in R if r.get('nan_vars')],
                stage1_s=dict(min=float(T.min()), median=float(np.median(T)), max=float(T.max()), sum=float(T.sum())) if len(T) else None)
    for s in STAGES:
        v = [r.get(f'N_res_max_{s}', np.nan) for r in R]; dr = [r.get(f'N_drift_{s}', np.nan) for r in R]
        summ[f'N_res_max_{s}'] = float(np.nanmax(v)); summ[f'N_drift_absmax_{s}'] = float(np.nanmax(np.abs(dr)))
    out = dict(summary=summ, cells=R)
    if a.out: json.dump(out, open(a.out, 'w'), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
