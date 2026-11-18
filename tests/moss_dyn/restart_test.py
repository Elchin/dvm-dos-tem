#!/usr/bin/env python3
"""Restart round-trip test: sp 60 yr in one process (A) vs sp 30 -> restart-sp.nc -> sp 30 from it (B2).
sp climate is year index % 30 and fire is off, so A and B2 must match bitwise if the restart carries
the whole state (incl. moss thickness). Compares every variable of A/restart-sp.nc and B2/restart-sp.nc.
  restart_test.py NAME BINARY TEMPLATE_CONFIG 'MS_JSON'"""
import json, os, subprocess, sys, numpy as np, netCDF4
W = os.path.expanduser('~/tem-dev/moss-work/rtest')
name, binary, tmpl, ms = sys.argv[1], os.path.expanduser(sys.argv[2]), sys.argv[3], json.loads(sys.argv[4])
base = json.load(open(tmpl)); base['model_settings'].update(ms)
root = os.path.join(W, name); os.makedirs(root)
def run(tag, years, restart=''):
    d = os.path.join(root, tag); os.makedirs(d)
    c = json.loads(json.dumps(base)); c['IO']['output_dir'] = d + '/'; c['IO']['restart_from'] = restart
    c['IO']['output_nc_eq'] = 0; c['IO']['output_nc_sp'] = 0
    p = os.path.join(d, 'config.js'); json.dump(c, open(p, 'w'), indent=1)
    r = subprocess.run(['nice', '-n', '10', binary, '-f', p, '-l', 'fatal'] + years, capture_output=True, text=True,
                       cwd=os.path.expanduser('~/tem-dev/rts-work/v4/src'))
    open(os.path.join(d, 'stdout.log'), 'w').write(r.stdout + r.stderr)
    if r.returncode: sys.exit('%s failed rc=%d' % (tag, r.returncode))
    return d
A = run('A', ['-p', '10', '-e', '300', '-s', '60', '-t', '0', '-n', '0'])
B1 = run('B1', ['-p', '10', '-e', '300', '-s', '30', '-t', '0', '-n', '0'])
B2 = run('B2', ['-p', '0', '-e', '0', '-s', '30', '-t', '0', '-n', '0'], os.path.join(B1, 'restart-sp.nc'))
da, db = netCDF4.Dataset(os.path.join(A, 'restart-sp.nc')), netCDF4.Dataset(os.path.join(B2, 'restart-sp.nc'))
n = same = 0; diff = []
for v in da.variables:
    n += 1; x, y = np.ma.getdata(da[v][:]), np.ma.getdata(db[v][:])
    if x.tobytes() == y.tobytes(): same += 1
    else: diff.append((v, float(np.nanmax(np.abs(x.astype(float) - y.astype(float))))))
b1 = netCDF4.Dataset(os.path.join(B1, 'restart-sp.nc'))
ty = np.ravel(da['TYPEsoil'][:]); dz = np.ravel(da['DZsoil'][:])
res = dict(name=name, ms=ms, n_vars=n, n_identical=same, differing=diff,
           moss_dz_A=float(dz[0]) if ty[0] == 0 else None,
           moss_dz_B1=float(np.ravel(b1['DZsoil'][:])[0]), yrsdist_A=int(np.ravel(da['yrsdist'][:])[0]))
print(json.dumps(res)); json.dump(res, open(os.path.join(root, 'result.json'), 'w'), indent=1)
