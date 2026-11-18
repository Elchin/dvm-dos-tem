#!/usr/bin/env python3
"""Restart carry test for moss_dyn runs (driver layout runs/<name>):
 (1) write : in-memory moss dz at Dec 2015 (MOSSDZ_yearly_tr last) == restart-tr.nc DZsoil[0]
 (2) read  : a 2-yr sc run from that restart with moss_dyn 0 (no thickness update) shows the same
             moss dz in Jan of year 1 (LAYERDZ monthly layer 0) -> the restart rebuilt the layer as saved
 (3) state : the first sc-year January thickness of the original run == the mode formula evaluated on
             the restored restart state (yrsdist for mode 1, moss PFT vegc for mode 2), bit for bit
  rcarry.py NAME BINARY"""
import json, os, subprocess, sys, numpy as np, netCDF4
W = os.path.expanduser('~/tem-dev/moss-work')
name, binary = sys.argv[1], os.path.expanduser(sys.argv[2])
stg = sys.argv[3] if len(sys.argv) > 3 else 'tr'   # tr: restart-tr (Dec 2015); sc: restart-sc (end of scenario, post-RTS regrowth)
R = os.path.join(W, 'runs', name)
rst = os.path.join(R, 'stage2a_hist', 'restart-tr.nc') if stg == 'tr' else os.path.join(R, 'stage2b_rts', 'restart-sc.nc')
d = netCDF4.Dataset(rst)
ty = np.ravel(d['TYPEsoil'][:]); dz0 = float(np.ravel(d['DZsoil'][:])[0]); ysd = int(np.ravel(d['yrsdist'][:])[0])
vegc = np.ma.getdata(d['vegc'][:])[0, 0]   # (pftpart, pft)
mem = float(np.ravel(netCDF4.Dataset(os.path.join(R, 'stage2a_hist', 'MOSSDZ_yearly_tr.nc') if stg == 'tr' else
                     os.path.join(R, 'stage2b_rts', 'MOSSDZ_yearly_sc.nc'))['MOSSDZ'][:])[-1])
cfg = json.load(open(os.path.join(R, 'stage2a_hist', 'config.js')))
ms = cfg['model_settings']; mode = ms.get('moss_dyn', 0)
dmin = ms.get('moss_min_thick', 0.001); th = ms.get('moss_t_half', 5.0)
dmaxp = ms.get('moss_dmax', -1.0)
cmt_dmax = {'T': 0.30, 'G': 0.035}[json.load(open(os.path.join(R, 'meta.json')))['site']]
dmax = cmt_dmax if dmaxp < 0 else dmaxp
if mode == 1:
    t = float(ysd); f = dmax * t / (t + th) if (t >= 0 and th > 0) else dmax
    f = min(dmax, max(dmin, f))
else:
    nv = [5, 6, 7]  # CMT05 nonvascular PFTs (cmt_dimvegetation nonvasculartype 3, 2, 1)
    c = 0.0
    for ip in nv:
        for i in range(3): c += float(vegc[i, ip])
    f = min(dmax, max(dmin, max(c, 0.0) / (ms.get('moss_bulkden', 25.0) * 1000.0 * ms.get('moss_cfrac', 0.45))))
# (2) carry runs from the restart: moss_dyn 0 (read check) and the run's own mode (state check)
def carry_run(tag, m):
    out = os.path.join(W, 'rcarry', name + '_' + stg + '_' + tag); os.makedirs(out)
    c2 = json.loads(json.dumps(cfg)); c2['IO']['output_dir'] = out + '/'; c2['IO']['restart_from'] = rst
    c2['model_settings']['moss_dyn'] = m; c2['model_settings']['thermokarst'] = 0
    p = os.path.join(out, 'config.js'); json.dump(c2, open(p, 'w'), indent=1)
    r = subprocess.run(['nice', '-n', '10', binary, '-f', p, '-l', 'fatal', '-p', '0', '-e', '0', '-s', '0', '-t', '0', '-n', '2'],
                       capture_output=True, text=True, cwd=os.path.expanduser('~/tem-dev/rts-work/v4/src'))
    if r.returncode: sys.exit('carry run failed: ' + r.stderr[-500:])
    return np.ma.getdata(netCDF4.Dataset(os.path.join(out, 'LAYERDZ_monthly_sc.nc'))['LAYERDZ'][:])[:, 0, 0, 0]
carry = carry_run('m0', 0)
jan = float(carry_run('m%d' % mode, mode)[0])
res = dict(name=name, restart=stg, mode=mode, type0=int(ty[0]), mem_dec2015=mem, restart_dz0=dz0,
           write_ok=(mem == dz0), carry_jan=float(carry[0]), carry_all_months_equal=bool(np.all(carry == dz0)),
           read_ok=bool(carry[0] == dz0), yrsdist=ysd, jan_after_restart=jan, formula=f, state_ok=(jan == f),
           formula_absdiff=abs(jan - f))
print(json.dumps(res))
