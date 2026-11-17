#!/usr/bin/env python3
"""RTS port: 10-cell demo-data case (toolik 10x10, row 0; RTS in cell (0,0) per the repo's
historic-explicit-thermokarst.nc). Scratch only: data dir = symlinks to demo-data + scratch run-mask.
usage: demo10.py BINARY WORK [-P 2]   runs off/on, -p 10 -e 100 -s 10 -t 115 -n 10, tr+sc outputs."""
import os, sys, json, subprocess, shutil
import numpy as np, netCDF4 as nc
REPO = os.path.expanduser('~/tem-dev/rts'); sys.path.insert(0, os.path.join(REPO, 'tests', 'rts'))
import rtslib as r
from concurrent.futures import ThreadPoolExecutor
DEMO = os.path.join(REPO, 'demo-data', 'cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10')
ARGS = ['-p', '10', '-e', '100', '-s', '10', '-t', '115', '-n', '10']
binary, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
P = int(sys.argv[sys.argv.index('-P') + 1]) if '-P' in sys.argv else 2
data = os.path.join(work, 'data'); os.makedirs(data, exist_ok=True)
for f in os.listdir(DEMO):
    if f.endswith('.nc') and f != 'run-mask.nc' and not os.path.exists(os.path.join(data, f)):
        os.symlink(os.path.join(DEMO, f), os.path.join(data, f))
rm = os.path.join(data, 'run-mask.nc')
if not os.path.exists(rm):
    shutil.copy(os.path.join(DEMO, 'run-mask.nc'), rm)
    d = nc.Dataset(rm, 'a'); m = np.zeros(d['run'].shape, dtype=d['run'].dtype); m[0, :] = 1; d['run'][:] = m; d.close()
OUTS = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.REGROW_OUTPUTS, 'RAINFALL': 'm', 'SNOWFALL': 'm', 'EET': 'm',
        'QRUNOFF': 'm', 'QDRAINAGE': 'm', 'SWE': 'm', 'IWCLAYER': 'm', 'VEGNTOT': 'm', 'AVLN': 'm', 'DEADN': 'm',
        'DWDN': 'm', 'NINPUT': 'm', 'NLOST': 'm'}
def run(name, on):
    wd = os.path.join(work, name); os.makedirs(wd, exist_ok=True)
    cfgp = r.make_config(wd, on, outputs=OUTS, data_dir=data)
    cfg = json.load(open(cfgp)); cfg['IO']['output_nc_sc'] = 1; json.dump(cfg, open(cfgp, 'w'), indent=2)
    with open(os.path.join(wd, 'stderr.txt'), 'w') as fe:
        p = subprocess.run([binary, '-f', cfgp, '-l', 'fatal'] + ARGS, cwd=REPO, stdout=subprocess.DEVNULL, stderr=fe, timeout=7200)
    return name, p.returncode
AUTO = {'soil_negpool_fix': 1, 'layer_n_remap_fix': 1, 'layer_bookkeep_fix': 1}
CASES = {'off': (False, None), 'on': (True, None), 'off_auto': (False, AUTO),
         'off_snc': (False, {**AUTO, 'n_safetynet_conserve': 1}),   # option A: global safety-net N conservation
         'off_snc2': (False, {**AUTO, 'n_safetynet_conserve': 2}),  # base fix: L1 + exact soil accumulation
         'on_snc2': (True, {'n_safetynet_conserve': 2})}
names = sys.argv[sys.argv.index('--cases') + 1].split(',') if '--cases' in sys.argv else ['off', 'on']
def run2(name):
    on, extra = CASES[name]
    wd = os.path.join(work, name); os.makedirs(wd, exist_ok=True)
    cfgp = r.make_config(wd, on, outputs=OUTS, data_dir=data, extra=extra)
    cfg = json.load(open(cfgp)); cfg['IO']['output_nc_sc'] = 1; json.dump(cfg, open(cfgp, 'w'), indent=2)
    with open(os.path.join(wd, 'stderr.txt'), 'w') as fe:
        p = subprocess.run([binary, '-f', cfgp, '-l', 'fatal'] + ARGS, cwd=REPO, stdout=subprocess.DEVNULL, stderr=fe, timeout=7200)
    return name, p.returncode
with ThreadPoolExecutor(P) as ex:
    for n, rc in ex.map(run2, names): print('run', n, 'rc', rc)
