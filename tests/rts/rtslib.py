"""Helpers to run dvmdostem on the 1-cell synthetic RTS case (testing-data/rts-synthetic)."""
import json, os, re, subprocess, shutil
import numpy as np
import netCDF4 as nc

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA = os.path.join(REPO, 'testing-data', 'rts-synthetic')

# Run lengths (years). Kept short so the whole suite runs in well under a minute.
RUN_ARGS = ['-p', '10', '-e', '300', '-s', '30', '-t', '30', '-n', '0']

# outputs: name -> timestep column ('y' yearly, 'm' monthly)
BASE_OUTPUTS = {'SOC': 'm', 'SHLWC': 'm', 'DEEPC': 'm', 'MOSSDZ': 'm', 'SHLWDZ': 'm',
                'DEEPDZ': 'm', 'ALD': 'y', 'LWCLAYER': 'm', 'TLAYER': 'm',
                'LAYERDZ': 'm', 'GPP': 'm', 'RH': 'm', 'VEGC': 'm', 'ORGN': 'm'}
RTS_OUTPUTS = {'TKLATC': 'm', 'TKLATN': 'm', 'TKREMTHICK': 'm', 'TKICEMELT': 'm', 'DRAINTYPE': 'm'}
ICE_OUTPUTS = {'TKSUBSID': 'm', 'TKEXICE': 'm', 'TKEXMELT': 'm', 'WATERTAB': 'm'}
# PLACEHOLDER ice-rich profile for the synthetic cell (tunable)
ICE_PARAMS = {'thermokarst_excess_fraction': 0.3, 'thermokarst_excess_top': 0.3,
              'thermokarst_excess_bottom': 2.0}

# PLACEHOLDER RTS parameters for the test (tunable, not calibrated)
RTS_PARAMS = {'thermokarst_folr': 0.8, 'thermokarst_c_resp_frac': 0.1,
              'thermokarst_veg_remove_frac': 0.99, 'thermokarst_thaw_depth': 1.0,
              'thermokarst_post_drainage': 1}

def strip_js_comments(txt):
    out = []
    for line in txt.splitlines():
        # remove // comments that are not inside strings (config.js has no '//' in strings except none)
        i = line.find('//')
        out.append(line if i < 0 else line[:i])
    return re.sub(r',(\s*[}\]])', r'\1', '\n'.join(out))

def write_outspec(path, outputs):
    src = os.path.join(REPO, 'config', 'output_spec.csv')
    lines = open(src).read().splitlines()
    res = [lines[0]]
    for ln in lines[1:]:
        f = ln.split(',')
        if f[0] in outputs:
            ts = outputs[f[0]]
            f[3] = 'y' if ts == 'y' else ''
            f[4] = 'm' if ts == 'm' else (f[4] if f[4] == 'invalid' else '')
        else:
            f[3] = '' if f[3] != 'invalid' else f[3]
            f[4] = '' if f[4] != 'invalid' else f[4]
        if len(f) > 5 and f[5] != 'invalid': f[5] = ''
        res.append(','.join(f))
    open(path, 'w').write('\n'.join(res) + '\n')

def make_config(workdir, rts_on, tk_file='historic-explicit-thermokarst.nc', outputs=None, extra=None):
    cfg = json.loads(strip_js_comments(open(os.path.join(REPO, 'config', 'config.js')).read()))
    io = cfg['IO']
    files = {'hist_climate_file': 'historic-climate.nc', 'proj_climate_file': 'projected-climate.nc',
             'veg_class_file': 'vegetation.nc', 'drainage_file': 'drainage.nc',
             'soil_texture_file': 'soil-texture.nc', 'co2_file': 'co2.nc',
             'proj_co2_file': 'projected-co2.nc', 'runmask_file': 'run-mask.nc',
             'topo_file': 'topo.nc', 'fri_fire_file': 'fri-fire.nc',
             'hist_exp_fire_file': 'historic-explicit-fire.nc',
             'proj_exp_fire_file': 'projected-explicit-fire.nc',
             'hist_exp_thermokarst_file': tk_file,
             'proj_exp_thermokarst_file': 'projected-explicit-thermokarst.nc'}
    for k, v in files.items():
        io[k] = os.path.join(DATA, v)
    io['parameter_dir'] = os.path.join(REPO, 'parameters') + '/'
    io['output_dir'] = os.path.join(workdir, 'output') + '/'
    io['output_spec_file'] = os.path.join(workdir, 'output_spec.csv')
    io['output_nc_eq'] = 0; io['output_nc_sp'] = 0; io['output_nc_tr'] = 1; io['output_nc_sc'] = 0
    io['restart_from'] = ''
    ms = cfg['model_settings']
    ms['thermokarst'] = bool(rts_on)
    if rts_on: ms.update(RTS_PARAMS)
    if extra: ms.update(extra)
    os.makedirs(os.path.join(workdir, 'output'), exist_ok=True)
    write_outspec(io['output_spec_file'], outputs or BASE_OUTPUTS)
    path = os.path.join(workdir, 'config.js')
    json.dump(cfg, open(path, 'w'), indent=2)
    return path

def run(binary, workdir, rts_on, **kw):
    if os.path.exists(workdir): shutil.rmtree(workdir)
    os.makedirs(workdir)
    cfgp = make_config(workdir, rts_on, **kw)
    r = subprocess.run([binary, '-f', cfgp, '-l', 'fatal'] + RUN_ARGS, cwd=REPO,
                       capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        raise RuntimeError('dvmdostem failed:\n' + r.stdout[-2000:] + r.stderr[-2000:])
    return os.path.join(workdir, 'output')

def read(outdir, var, ts='monthly', stage='tr'):
    d = nc.Dataset(os.path.join(outdir, '%s_%s_%s.nc' % (var, ts, stage)))
    a = np.array(d[var][:], dtype=float); d.close()
    return a

def compare_dirs(a, b, rtol=0.0, atol=0.0, skip=()):
    """Compare all .nc outputs in dirs a and b. Returns list of mismatch strings."""
    bad = []
    fa = sorted(f for f in os.listdir(a) if f.endswith('.nc') and not f.startswith('restart') and f != 'run_status.nc')
    fb = sorted(f for f in os.listdir(b) if f.endswith('.nc') and not f.startswith('restart') and f != 'run_status.nc')
    if fa != fb: bad.append('file sets differ: %s' % (sorted(set(fa) ^ set(fb)),))
    for f in set(fa) & set(fb):
        if any(f.startswith(s) for s in skip): continue
        da, db = nc.Dataset(os.path.join(a, f)), nc.Dataset(os.path.join(b, f))
        for v in da.variables:
            if v not in db.variables: bad.append('%s:%s missing' % (f, v)); continue
            x, y = np.array(da[v][:], dtype=float), np.array(db[v][:], dtype=float)
            if x.shape != y.shape or not np.allclose(x, y, rtol=rtol, atol=atol, equal_nan=True):
                diff = np.nanmax(np.abs(x - y)) if x.shape == y.shape else 'shape'
                bad.append('%s:%s differs (max abs diff %s)' % (f, v, diff))
        da.close(); db.close()
    return bad
