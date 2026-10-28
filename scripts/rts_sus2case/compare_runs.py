#!/usr/bin/env python3
"""Compare stage-2b outputs of two run trees (e.g. e1000 vs e2000 equilibrium).

For each cell and run name present in both trees (run_phase2.py layout
ROOT/cNN/<run>/stage2b_*/), summarise ALT (ALD), ECOC (VEGC+SOC+DEADC+DWDC
[+TKCH4POOL]), soil C pools (SOC, SHLWC, DEEPC, mineral = SOC-SHLWC-DEEPC) and
NEE (loss-positive, RG+RM+soil+RHDWD-GPP; TKNETC adds the RTS terms):
  * values at the start (2016 mean / Jan 2016 pool) and end (2091-2100 mean / Dec 2100)
  * cumulative 2016-2100 NEE and TKNETC (gC/m2)
  * the RTS effect (yYYYY - base) for each tree, and its difference between trees
Writes JSON (and a markdown table with --md).
"""
import argparse, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from rts_sus2case.check_phase2 import read, outdir  # noqa: E402

POOLS = ('ECOC', 'SOC', 'SHLWC', 'DEEPC', 'MINEC')


def series(d):
    g = lambda v, ts='monthly': (lambda a: None if a is None else a.reshape(a.shape[0], -1).sum(axis=1))(read(d, v, ts))
    s = {k: g(k) for k in ('VEGC', 'SOC', 'SHLWC', 'DEEPC', 'DEADC', 'DWDC', 'TKNETC')}
    pool = g('TKCH4POOL')
    s['ECOC'] = s['VEGC'] + s['SOC'] + s['DEADC'] + s['DWDC'] + (0 if pool is None else pool)
    s['MINEC'] = s['SOC'] - s['SHLWC'] - s['DEEPC']
    soil = g('TKRHCO2') + g('TKCH4EMIS') if g('TKRHCO2') is not None else g('RHSOM')
    s['NEE'] = g('RG') + g('RM') + soil + g('RHDWD') - g('GPP')
    a = read(d, 'ALD', 'yearly'); s['ALD'] = a.reshape(a.shape[0], -1)[:, 0]
    return s


def summary(s):
    o = {'ALT_2016_2025': float(s['ALD'][:10].mean()), 'ALT_2091_2100': float(s['ALD'][-10:].mean())}
    for k in POOLS:
        o[f'{k}_start'] = float(s[k][0]); o[f'{k}_end'] = float(s[k][-1]); o[f'd{k}'] = float(s[k][-1] - s[k][0])
    o['NEE_cum'] = float(s['NEE'].sum()); o['TKNETC_cum'] = float(s['TKNETC'].sum())
    o['NEE_2016_2025_yr'] = float(s['NEE'][:120].sum() / 10); o['NEE_2091_2100_yr'] = float(s['NEE'][-120:].sum() / 10)
    return o


def diff(a, b):
    return {k: dict(a=a[k], b=b[k], d=b[k] - a[k], rel=(b[k] - a[k]) / abs(a[k]) if a[k] else None) for k in a}


def compare(root_a, root_b, cells=None):
    out = {}
    cells = cells or sorted(c for c in os.listdir(root_a) if os.path.isdir(os.path.join(root_b, c)))
    for c in cells:
        ra = {x for x in os.listdir(os.path.join(root_a, c)) if outdir(os.path.join(root_a, c, x))}
        rb = {x for x in os.listdir(os.path.join(root_b, c)) if outdir(os.path.join(root_b, c, x))}
        S = {}
        for x in sorted(ra & rb):
            S[x] = (summary(series(outdir(os.path.join(root_a, c, x)))), summary(series(outdir(os.path.join(root_b, c, x)))))
        res = {x: diff(*v) for x, v in S.items()}
        for x in S:
            if x.startswith('y') and 'base_' + x.split('_')[-1] in S:
                bx = 'base_' + x.split('_')[-1]
                ea = {k: S[x][0][k] - S[bx][0][k] for k in S[x][0]}
                eb = {k: S[x][1][k] - S[bx][1][k] for k in S[x][1]}
                res[f'effect_{x}'] = diff(ea, eb)
        out[c] = res
    return out


def to_md(res, la, lb, keys=('ALT_2091_2100', 'ECOC_start', 'ECOC_end', 'dECOC', 'SOC_end', 'SHLWC_end',
                            'DEEPC_end', 'MINEC_end', 'NEE_2016_2025_yr', 'NEE_2091_2100_yr', 'NEE_cum', 'TKNETC_cum')):
    L = [f'| cell | run | quantity | {la} | {lb} | {lb}-{la} | rel |', '|---|---|---|---|---|---|---|']
    for c, R in res.items():
        for x, D in R.items():
            for k in keys:
                if k not in D: continue
                v = D[k]; rel = '' if v['rel'] is None else f"{100 * v['rel']:+.1f}%"
                L.append(f"| {c} | {x} | {k} | {v['a']:.4g} | {v['b']:.4g} | {v['d']:+.4g} | {rel} |")
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root_a'); ap.add_argument('root_b'); ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--cells', default=None); ap.add_argument('--md', default=None)
    ap.add_argument('--labels', default='a,b')
    a = ap.parse_args()
    res = compare(a.root_a, a.root_b, a.cells.split(',') if a.cells else None)
    json.dump(res, open(a.out, 'w'), indent=1)
    if a.md: open(a.md, 'w').write(to_md(res, *a.labels.split(',')))


if __name__ == '__main__':
    main()
