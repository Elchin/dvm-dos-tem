#!/usr/bin/env python3
"""moss_dyn test harness (scratch, ~/tem-dev/moss-work). One case = one workdir under runs/<name>,
run through moss-wt/scripts/run_rts_two_stage.py (stage1 pr/eq/sp, 2a hist, 2b scenario with RTS 2021).
  run_case.py NAME SITE BINARY 'MS_JSON' [--pr 10 --eq 300 --sp 30 --hist N --sc 30 --mode rts]
              [--fri N]  (scratch copy of fri-fire.nc with FRI N, fire on in eq only)
Base config = 7b6a629a config/config.js (no moss keys); parameters = scratch copy of 7b6a629a parameters/.
Model settings = bench v4b (thermal_solver 1, n_safetynet_conserve 2, thermokarst_regrow, rmf {"70":0}) + MS_JSON."""
import argparse, csv, json, os, shutil, subprocess, sys, hashlib
W = os.path.expanduser('~/tem-dev/moss-work'); B = os.path.expanduser('~/tem-dev/rts-work/v4b/bench')
WT = os.path.expanduser('~/tem-dev/moss-wt')
MONTHLY = ['GPP', 'NPP', 'RECO', 'RHSOM', 'VEGC', 'SOC', 'DEADC', 'DWDC', 'TKNETC', 'TKLATC', 'TKRESPC',
           'TKREGROWC', 'BURNVEG2AIRC', 'BURNSOIL2AIRC', 'LAYERDZ', 'TLAYER', 'LAYERTYPE', 'LWCLAYER']
YEARLY = ['ALD', 'MOSSDZ', 'SHLWDZ', 'DEEPDZ', 'SHLWC', 'DEEPC', 'MINEC', 'YSD', 'LAI']
PFTV = {'VEGC', 'NPP', 'GPP'}
LAYERED = {'LAYERDZ', 'TLAYER', 'LAYERTYPE', 'LWCLAYER'}

def outspec(dst):
    rows = list(csv.reader(open(os.path.join(WT, 'config', 'output_spec.csv'))))
    h = rows[0]; c = {k: i for i, k in enumerate(h)}
    for r in rows[1:]:
        for col in ('Yearly', 'Monthly', 'Daily', 'PFT', 'Compartments', 'Layers'):
            if r[c[col]] not in ('invalid', 'forced'): r[c[col]] = ''
    for r in rows[1:]:
        v = r[0]
        if v in MONTHLY + YEARLY:
            if r[c['Yearly']] != 'invalid': r[c['Yearly']] = 'y'
            if v in MONTHLY and r[c['Monthly']] != 'invalid': r[c['Monthly']] = 'm'
            if v in PFTV and r[c['PFT']] != 'invalid': r[c['PFT']] = 'p'
            if v in LAYERED and r[c['Layers']] not in ('invalid',): r[c['Layers']] = 'l'
    csv.writer(open(dst, 'w', newline='')).writerows(rows)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('name'); ap.add_argument('site'); ap.add_argument('binary')
    ap.add_argument('ms', nargs='?', default='{}')
    for k, d in (('pr', 10), ('eq', 300), ('sp', 30), ('sc', 30)): ap.add_argument('--' + k, type=int, default=d)
    ap.add_argument('--hist', type=int, default=None); ap.add_argument('--mode', default='rts')
    ap.add_argument('--params', default='params'); ap.add_argument('--fri', type=int, default=0); ap.add_argument('--stop-after', default='stage2b')
    a = ap.parse_args()
    wd = os.path.join(W, 'runs', a.name)
    if os.path.exists(wd): sys.exit('exists: ' + wd)
    os.makedirs(wd)
    inp = os.path.join(B, 'inputs', a.site)
    ov = json.load(open(os.path.join(inp, 'two_stage_overrides.json')))
    for k, v in list(ov['IO'].items()):
        if k.endswith('_file'): ov['IO'][k] = os.path.join(inp, v)
    ov['IO']['parameter_dir'] = os.path.join(W, a.params) + '/'
    spec = os.path.join(wd, 'output_spec.csv'); outspec(spec); ov['IO']['output_spec_file'] = spec
    ov['IO']['output_nc_tr'] = 1; ov['IO']['output_nc_sc'] = 1
    ov['model_settings'].update({'thermal_solver': 1, 'n_safetynet_conserve': 2, 'thermokarst_regrow': True,
                                 'rhmoistfrozen_override': {'70': 0.0}})
    ov['model_settings'].update(json.loads(a.ms))
    if a.fri:
        import netCDF4
        f = os.path.join(wd, 'fri-fire.nc'); shutil.copy(os.path.join(inp, 'fri-fire.nc'), f)
        d = netCDF4.Dataset(f, 'a'); d['fri'][:] = a.fri; d['fri_severity'][:] = 2
        d['fri_jday_of_burn'][:] = 200; d['fri_area_of_burn'][:] = 1; d.close()
        ov['IO']['fri_fire_file'] = f
        ov['model_settings']['fire'] = 1
        ov.setdefault('stage_settings', {})['eq'] = {'fire': True, 'dsb': True}
    ovp = os.path.join(wd, 'overrides.json'); json.dump(ov, open(ovp, 'w'), indent=1)
    hist = a.hist or {'G': 66, 'L': 66}.get(a.site, 115)
    cmd = [sys.executable, os.path.expanduser("~/tem-dev/rts-work/v4/src/scripts/run_rts_two_stage.py"), '--binary', a.binary,
           '--config', os.path.join(W, 'config_ref.js'), '--config-overrides', ovp, '--workdir', wd,
           '--pr', str(a.pr), '--eq', str(a.eq), '--sp', str(a.sp), '--hist', str(hist), '--sc', str(a.sc),
           '--mode', a.mode, '--stop-after', a.stop_after]
    meta = dict(name=a.name, site=a.site, binary=a.binary, binary_md5=hashlib.md5(open(a.binary, 'rb').read()).hexdigest(),
                ms=json.loads(a.ms), cmd=cmd, fri=a.fri)
    json.dump(meta, open(os.path.join(wd, 'meta.json'), 'w'), indent=1)
    rc = subprocess.call(['nice', '-n', '10'] + cmd, stdout=open(os.path.join(wd, 'driver.log'), 'w'), stderr=subprocess.STDOUT)
    print(a.name, 'rc', rc, open(os.path.join(wd, 'driver.log')).read().strip().splitlines()[-1:])
    sys.exit(rc)

if __name__ == '__main__': main()
