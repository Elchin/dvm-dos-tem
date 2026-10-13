"""Regression tests for the two-stage RTS workflow (scripts/run_rts_two_stage.py).

Short years for CI. Set RTS_FULL_TIMING=1 to also run the full stage 1
(-p 100 -e 1000 -s 200) and check it completes in <= 300 s for one cell.
"""
import json, os, shutil, subprocess, sys, time
import numpy as np
import pytest
import rtslib as r

DRIVER = os.path.join(r.REPO, 'scripts', 'run_rts_two_stage.py')
BIN = os.environ.get('DVMDOSTEM_BIN', os.path.join(r.REPO, 'dvmdostem'))
SHORT = ['--pr', '10', '--eq', '100', '--sp', '20', '--hist', '10', '--sc', '15']
SC_RTS_YEAR, RTS_MONTH = 5, 7          # projected-explicit-thermokarst-rts.nc
MS = SC_RTS_YEAR * 12 + RTS_MONTH
HAS_CH4 = os.path.exists(os.path.join(r.REPO, 'include', 'RTSMethane.h'))
OUTS = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS, **(r.CH4_OUTPUTS if HAS_CH4 else {})}

def _overrides(w):
    spec = os.path.join(w, 'output_spec.csv'); r.write_outspec(spec, OUTS)
    ov = json.load(open(os.path.join(r.DATA, 'two_stage_overrides.json')))
    for k, v in ov['IO'].items():
        if isinstance(v, str) and v.endswith('.nc'): ov['IO'][k] = os.path.join(r.DATA, v)
    ov['IO']['output_spec_file'] = spec
    ov['IO']['parameter_dir'] = os.path.join(r.REPO, 'parameters') + '/'
    p = os.path.join(w, 'overrides.json'); json.dump(ov, open(p, 'w')); return p

def _drive(w, *extra):
    os.makedirs(w, exist_ok=True)
    p = subprocess.run([sys.executable, DRIVER, '--binary', BIN, '--workdir', w,
                        '--config-overrides', _overrides(w)] + SHORT + list(extra),
                       capture_output=True, text=True, timeout=900)
    assert p.returncode == 0, p.stdout + p.stderr
    return w

@pytest.fixture(scope='module')
def wf(tmp_path_factory):
    w = tmp_path_factory.mktemp('two_stage')
    rts = _drive(str(w / 'rts'), *(['--ch4'] if HAS_CH4 else []))
    base = _drive(str(w / 'base'), '--mode', 'baseline',
                  '--stage1-restart', os.path.join(rts, 'stage1', 'restart-sp.nc'),
                  *(['--ch4'] if HAS_CH4 else []))
    # single continuous run (RTS off) for restart-continuity check
    cont = str(w / 'cont'); os.makedirs(cont)
    cfg = r.make_config(cont, False, tk_file='thermokarst-none.nc', outputs=OUTS)
    subprocess.run([BIN, '-f', cfg, '-l', 'fatal', '-p', '10', '-e', '100', '-s', '20', '-t', '10', '-n', '0'],
                   cwd=r.REPO, check=True, capture_output=True)
    return {'rts': rts, 'base': base, 'cont': os.path.join(cont, 'output')}

def _sc(d, v, mode):
    return r.read(os.path.join(d, 'stage2b_' + mode), v, 'monthly', 'sc').ravel()

def test_restart_chain_files(wf):
    for f in ['stage1/restart-sp.nc', 'stage2a_hist/restart-tr.nc', 'stage2b_rts/restart-sc.nc', 'timings.json']:
        assert os.path.exists(os.path.join(wf['rts'], f)), f

def test_stage1_and_historical_match_continuous_run(wf):
    """Stage 1 + 2a (RTS off) via restarts == one continuous original-model run, bit for bit."""
    assert r.compare_dirs(wf['cont'], os.path.join(wf['rts'], 'stage2a_hist'), 0.0, 0.0) == []

def test_ice_init_spinup_and_rts_in_scenario(wf):
    d = wf['rts']
    x, lat, sub, dr = (_sc(d, v, 'rts') for v in ('TKEXICE', 'TKLATC', 'TKSUBSID', 'DRAINTYPE'))
    assert x[0] > 0                                     # ice initialized in first scenario month
    assert np.all(sub[:MS] == 0) and np.all(lat[:MS] == 0)   # short spin-up: no slump yet
    assert lat[MS] > 0 and sub[MS] > 0.1                # RTS initiation in the flagged year
    assert np.all(dr[:MS] == 0) and np.all(dr[MS:] == 1)
    assert _sc(d, 'SOC', 'rts')[MS] < 0.9 * _sc(d, 'SOC', 'rts')[MS - 1]

def test_baseline_mode_no_slump_same_spinup(wf):
    b, d = wf['base'], wf['rts']
    assert np.all(_sc(b, 'TKLATC', 'baseline') == 0) and np.all(_sc(b, 'TKSUBSID', 'baseline') == 0)
    assert _sc(b, 'TKEXICE', 'baseline')[0] == _sc(d, 'TKEXICE', 'rts')[0]
    for v in ['SOC', 'VEGC', 'GPP', 'RHSOM', 'TLAYER']:
        assert np.array_equal(_sc(b, v, 'baseline')[:MS], _sc(d, v, 'rts')[:MS]), v

@pytest.mark.skipif(not HAS_CH4, reason='rts-ch4 only')
def test_ch4_in_workflow(wf):
    assert _sc(wf['rts'], 'TKCH4EMIS', 'rts').sum() > 0
    assert _sc(wf['base'], 'TKCH4EMIS', 'baseline').sum() > 0

@pytest.mark.skipif(os.environ.get('RTS_FULL_TIMING') != '1', reason='set RTS_FULL_TIMING=1')
def test_full_stage1_timing(tmp_path):
    w = str(tmp_path / 'full')
    p = subprocess.run([sys.executable, DRIVER, '--binary', BIN, '--workdir', w,
                        '--config-overrides', _overrides_dir(tmp_path), '--hist', '1', '--sc', '1'],
                       capture_output=True, text=True, timeout=1200)
    assert p.returncode == 0, p.stdout + p.stderr
    t = json.load(open(os.path.join(w, 'timings.json')))
    print('full stage-1 timing: %.1f s' % t['stage1'])
    assert t['stage1'] <= 300.0

def _overrides_dir(tmp_path):
    d = str(tmp_path / 'ov'); os.makedirs(d, exist_ok=True); return _overrides(d)
