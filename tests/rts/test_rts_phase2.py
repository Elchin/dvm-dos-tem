"""Tests for scripts/rts_sus2case/run_phase2.py and check_phase2.py.

Unit tests use synthetic files; the end-to-end test runs one Gydan cell (stage 2a
+ base/notrig/one start year, short sc) when the inputs, a stage-1 restart and a
binary are available (RTS_SUS2CASE_INPUTS, RTS_SUS2CASE_STAGE1, DVMDOSTEM_BIN).
"""
import json, os, subprocess, sys
import numpy as np
import netCDF4
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts')); sys.path.insert(0, HERE)
from rts_sus2case import run_phase2 as rp, check_phase2 as cp  # noqa: E402


def test_plan_counts_and_names():
    cells = ['c%02d' % i for i in range(90)]
    P = rp.plan(cells, rp.YEARS, [0.8, 0.5], checks=False)
    assert rp.YEARS == (2021, 2031, 2041, 2051, 2061, 2071, 2081, 2091)
    assert len(P) == 1620 and sum(p[1] == 'base' for p in P) == 180
    assert len(rp.plan(['c55'], rp.YEARS, [0.8], checks=True)) == 10
    assert rp.run_name('base', None, 0.8) == 'base_f08' and rp.run_name('notrig', None, 0.5) == 'notrig_f05'
    assert rp.run_name('rts', 2031, 0.5) == 'y2031_f05'
    names = {rp.run_name(k, y, f) for _, k, y, f in P}
    assert len(names) == 18


def test_outputs_reduced_but_complete_for_closure():
    o = rp.outputs()
    for v in rp.DROP: assert v not in o
    need = ('VEGC', 'SOC', 'DEADC', 'DWDC', 'SHLWC', 'DEEPC', 'GPP', 'RG', 'RM', 'RHSOM', 'RHDWD', 'TKLATC',
            'TKRESPC', 'TKNETC', 'TKREGROWC', 'TKRHCO2', 'TKCH4EMIS', 'TKCH4POOL', 'ALD', 'ORGN', 'AVLN')
    for v in need: assert v in o, v


def test_overrides_keep_n_switches_default(tmp_path):
    d = tmp_path / 'cell'; d.mkdir()
    json.dump({'IO': {'proj_exp_thermokarst_file': 'x.nc', 'other': 'y'},
               'model_settings': {'thermokarst_folr': 0.8}}, open(d / 'two_stage_overrides.json', 'w'))
    p = rp._overrides(str(d), str(tmp_path), 'spec.csv', 0.5, '/t/trig.nc')
    ov = json.load(open(p))
    assert ov['model_settings']['thermokarst_folr'] == 0.5 and ov['model_settings']['thermokarst_regrow'] is True
    assert ov['IO']['proj_exp_thermokarst_file'] == '/t/trig.nc' and ov['IO']['other'] == 'y'
    for k in ('n_nonvasc_extract', 'n_extract_limiter', 'soil_giveup_diag'): assert k not in ov['model_settings']
    json.dump({'IO': {}, 'model_settings': {'n_extract_limiter': True}}, open(d / 'two_stage_overrides.json', 'w'))
    with pytest.raises(AssertionError): rp._overrides(str(d), str(tmp_path), 's', 0.8)


def test_cutoff_and_closure():
    assert cp.cutoff(2021, 'monthly') == 67 and cp.cutoff(2021, 'yearly') == 5
    ecoc = np.array([10.0, 9.0, 9.5]); t = np.array([0.0, 1.0, -0.5])
    c = cp.closure(ecoc, t, t)
    assert c['id_max'] == 0 and c['c_month_max'] == 0 and c['c_cum'] == 0
    c = cp.closure(ecoc, t + [0, 0.02, 0], t)
    assert abs(c['c_month_max'] - 0.02) < 1e-12 and abs(c['id_max'] - 0.02) < 1e-12


def _nc(path, v, a):
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('time', len(a)); ds.createVariable(v, 'f8', ('time',))[:] = a


def test_compare_prefix(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'; a.mkdir(); b.mkdir()
    x = np.arange(85 * 12, dtype=float); y = x.copy(); y[cp.cutoff(2031, 'monthly'):] += 1
    _nc(a / 'SOC_monthly_sc.nc', 'SOC', x); _nc(b / 'SOC_monthly_sc.nc', 'SOC', y)
    xy = np.arange(85.0); _nc(a / 'ALD_yearly_sc.nc', 'ALD', xy); _nc(b / 'ALD_yearly_sc.nc', 'ALD', xy)
    n, bad, after = cp.compare(str(a), str(b), 2031)
    assert n == 2 and not bad and after
    n, bad, _ = cp.compare(str(a), str(b), 2041)
    assert bad == ['SOC_monthly_sc.nc']
    assert cp.compare(str(a), str(b))[1] == ['SOC_monthly_sc.nc']


INP = os.environ.get('RTS_SUS2CASE_INPUTS', os.path.expanduser('~/box_archive/russia-inputs/spatial10x10/inputs'))
ST1 = os.environ.get('RTS_SUS2CASE_STAGE1', '')
BIN = os.environ.get('DVMDOSTEM_BIN', os.path.join(REPO, 'dvmdostem'))


@pytest.mark.skipif(not (ST1 and os.path.exists(os.path.join(ST1, 'c55', 'stage1', 'restart-sp.nc'))
                         and os.path.exists(os.path.join(INP, 'c55')) and os.path.exists(BIN)),
                    reason='Gydan inputs / stage-1 restart / binary not available')
def test_end_to_end_one_cell(tmp_path):
    out = tmp_path / 'p2'
    p = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'rts_sus2case', 'run_phase2.py'),
                        '--inputs', INP, '--stage1', ST1, '--out', str(out), '--binary', BIN, '--nproc', '4',
                        '--cells', 'c55', '--years', '2021', '--folr', '0.8,0.5', '--sc', '10', '--checks'],
                       capture_output=True, text=True, timeout=3600)
    assert p.returncode == 0, (p.stdout + p.stderr)[-3000:]
    T = json.load(open(out / 'timings.json'))['summary']
    assert T['failed'] == 0 and T['runs'] == 6
    q = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'rts_sus2case', 'check_phase2.py'), str(out)],
                       capture_output=True, text=True)
    assert q.returncode == 0, q.stdout[-3000:]
    S = json.load(open(out / 'checks.json'))['summary']
    assert S['finished'] == 6 and S['ab_tests'] == 5 and not S['ab_fail']
