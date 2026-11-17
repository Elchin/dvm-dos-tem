#!/usr/bin/env python3
"""tr->sc split-run check of the RTS restart state (scratch). For each case: A = continuous
pr/eq/sp/tr(12 yr)/sc(10 yr) in one process; B = sc only, new process, from a COPY of A's restart-tr.nc.
usage: port_split.py BIN_A BIN_B WORK"""
import os, sys, json, shutil, subprocess
import numpy as np, netCDF4 as nc
from concurrent.futures import ThreadPoolExecutor
REPO = os.path.expanduser('~/tem-dev/rts'); sys.path.insert(0, os.path.join(REPO, 'tests', 'rts'))
import rtslib as r
binA, binB, W = [os.path.abspath(x) for x in sys.argv[1:4]]
OUTS = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.REGROW_OUTPUTS, 'NINPUT': 'm', 'AVLN': 'm'}
AUTO = {'soil_negpool_fix': 1, 'layer_n_remap_fix': 1, 'layer_bookkeep_fix': 1}
CASES = {'off_auto': (False, AUTO), 'on': (True, {})}
def cfg(wd, on, extra, restart=''):
    p = r.make_config(wd, on, outputs=OUTS, extra=extra)
    c = json.load(open(p)); c['IO']['output_nc_sc'] = 1; c['IO']['output_nc_tr'] = 1; c['IO']['restart_from'] = restart
    json.dump(c, open(p, 'w'), indent=2); return p
def run(b, p, args, wd):
    with open(os.path.join(wd, 'stderr.txt'), 'w') as fe:
        return subprocess.run([b, '-f', p, '-l', 'fatal'] + args, cwd=REPO, stdout=subprocess.DEVNULL, stderr=fe).returncode
def A(n):
    on, ex = CASES[n]; wd = os.path.join(W, n + '_A'); os.makedirs(wd, exist_ok=True)
    return run(binA, cfg(wd, on, ex), ['-p', '10', '-e', '300', '-s', '30', '-t', '12', '-n', '10'], wd)
def B(n, b, tag):
    on, ex = CASES[n]; wd = os.path.join(W, n + '_' + tag); os.makedirs(wd, exist_ok=True)
    rs = os.path.join(wd, 'restart-tr-copy.nc'); shutil.copy(os.path.join(W, n + '_A', 'output', 'restart-tr.nc'), rs)
    return run(b, cfg(wd, on, ex, restart=rs), ['-p', '0', '-e', '0', '-s', '0', '-t', '0', '-n', '10', '--no-output-cleanup'], wd)
with ThreadPoolExecutor(2) as e: print('A', list(e.map(A, CASES)))
jobs = [(n, binA, 'B') for n in CASES] + [('on', binB, 'Bold')]
with ThreadPoolExecutor(3) as e: print('B', list(e.map(lambda j: B(*j), jobs)))
for n, tag in (('off_auto', 'B'), ('on', 'B'), ('on', 'Bold')):
    a, b = os.path.join(W, n + '_A', 'output'), os.path.join(W, n + '_' + tag, 'output')
    diff, maxd = [], {}
    for f in sorted(os.listdir(a)):
        if not f.endswith('_sc.nc') or not os.path.exists(os.path.join(b, f)): continue
        v = f.split('_')[0]; x = np.array(nc.Dataset(os.path.join(a, f))[v][:], float); y = np.array(nc.Dataset(os.path.join(b, f))[v][:], float)
        if not np.array_equal(x, y): diff.append(f); maxd[v] = float(np.nanmax(np.abs(x - y)))
    print(n, tag, 'sc files differing:', len(diff), {k: round(v, 4) for k, v in sorted(maxd.items())})
