#!/usr/bin/env python3
"""Effect of n_safetynet_conserve 2 on the regression cases: PAIRS file lines 'A B' (run dirs, A = reference,
B = level 2). Per case and stage (spin000, main): bit identity, C closure (closure.py), monthly N residual
d(VEGNTOT+ORGN+AVLN+DEADN+DWDN) - (NINPUT - NLOST) max |.|, and last-year VEGC / GPP / SOC (B - A, %).
usage: port_snc2_reg.py PAIRS OUT.json   (run from the reg dir)"""
import os, sys, json, subprocess, glob, numpy as np, netCDF4 as nc
HERE = os.path.dirname(os.path.abspath(__file__)); CL = os.path.join(HERE, '..', 'regression', 'closure.py')
def rd(o, v):
    f = glob.glob(os.path.join(o, v + '_monthly_*.nc')) or glob.glob(os.path.join(o, v + '_yearly_*.nc'))
    if not f: return None, None
    f = sorted(f)[0]
    with nc.Dataset(f) as x: a = np.ma.getdata(x[v][:]).astype(float)
    a = np.where(np.abs(a) > 1e20, 0.0, a)
    return a.reshape(a.shape[0], -1).sum(1), ('yearly' if '_yearly_' in f else 'monthly')
def last_year(o, v):
    a, ts = rd(o, v)
    if a is None: return None
    return float(a[-12:].mean() if ts == 'monthly' and v != 'GPP' else (a[-12:].sum() if ts == 'monthly' else a[-1]))
def n_resid(o):
    st = 0
    for v in ('VEGNTOT', 'ORGN', 'AVLN', 'DEADN', 'DWDN'):
        a, ts = rd(o, v)
        if a is None or ts != 'monthly': return None
        st = st + a
    i, _ = rd(o, 'NINPUT'); l, _ = rd(o, 'NLOST')
    if i is None or l is None: return None
    r = np.diff(st) - (i - l)[1:]
    return dict(max_abs_month=float(np.max(np.abs(r))), sum=float(r.sum()))
def ident(a, b):
    d = []
    for f in sorted(os.listdir(a)):
        if not f.endswith('.nc') or not os.path.exists(os.path.join(b, f)): continue
        with nc.Dataset(os.path.join(a, f)) as x, nc.Dataset(os.path.join(b, f)) as y:
            for v in x.variables:
                if v in y.variables and not np.array_equal(np.ma.getdata(x[v][:]), np.ma.getdata(y[v][:])): d.append(f); break
    return d
out = {}
for line in open(sys.argv[1]):
    if not line.strip(): continue
    A, B = line.split(); case = '/'.join(A.split('/')[-2:]); out[case] = {}
    cl = {k: json.loads(subprocess.run([sys.executable, CL, k], capture_output=True, text=True).stdout or '{}') for k in (A, B)}
    out[case]['closure'] = {k: {x: cl[r].get(x) for x in ('C_resid_mean_gm2yr', 'C_resid_max_abs_annual', 'main_safetynet', 'spin_safetynet', 'status')} for k, r in (('A', A), ('B', B))}
    for st in ('spin000', 'main'):
        oa, ob = os.path.join(A, st, 'output'), os.path.join(B, st, 'output'); e = {}
        e['differ_files'] = len(ident(oa, ob)); e['N_A'] = n_resid(oa); e['N_B'] = n_resid(ob)
        for v in ('VEGC', 'GPP', 'SOC'):
            a, b = last_year(oa, v), last_year(ob, v)
            if a is not None and b is not None: e[v] = dict(A=a, B=b, d=b - a, pct=100 * (b - a) / a if a else None)
        out[case][st] = e
    print(case, json.dumps(out[case]['main'].get('VEGC')), flush=True)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
