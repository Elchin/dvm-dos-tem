"""TKNETC combined total-carbon-loss diagnostic: identity, CH4 handling and closure.

TKNETC = NEE_month + TKLATC + TKRESPC - TKREGROWC (gC/m2/month, positive =
carbon leaving the ecosystem), where NEE_month is the loss-positive monthly
NEE of the RTS closures: RG+RM+soil+RHDWD-GPP. The soil term is RHSOM with
thermokarst_ch4 off, and TKRHCO2+TKCH4EMIS with it on -- in code
TKRHCO2 = RHSOM-prod+oxid and TKCH4EMIS = prod-oxid, so that sum is RHSOM
with CH4 carbon counted exactly once (not RHSOM+TKCH4EMIS, which would
double-count it, and not TKRHCO2 alone, which would drop it). TKLATCH4 is a
diagnostic split of TKLATC and is never added. ECOC below is
VEGC+SOC+DEADC+DWDC+TKCH4POOL (TKCH4POOL is 0 in CH4 scheme a), so monthly
closure is dECOC = -TKNETC in the diff/F[1:] alignment used by the suite.
"""
import json, os, subprocess, sys
import numpy as np
import pytest
import rtslib as r

BIN = os.environ.get('DVMDOSTEM_BIN', os.path.join(r.REPO, 'dvmdostem'))
GYDAN70 = os.environ.get('RTS_GYDAN70_DIR', '/workspace/russia-inputs/inputs/gydan_cmt70')
OUTS = {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.CH4_OUTPUTS, **r.REGROW_OUTPUTS}

def _g(d, v, stage='tr'):
    return r.read(d, v, 'monthly', stage).ravel()

def _ecoc(d, stage='tr'):
    c = _g(d, 'VEGC', stage) + _g(d, 'SOC', stage) + _g(d, 'DEADC', stage) + _g(d, 'DWDC', stage)
    if os.path.exists(os.path.join(d, 'TKCH4POOL_monthly_%s.nc' % stage)):
        c = c + _g(d, 'TKCH4POOL', stage)
    return c

def _nee(d, ch4, stage='tr'):
    soil = (_g(d, 'TKRHCO2', stage) + _g(d, 'TKCH4EMIS', stage)) if ch4 else _g(d, 'RHSOM', stage)
    return _g(d, 'RG', stage) + _g(d, 'RM', stage) + soil + _g(d, 'RHDWD', stage) - _g(d, 'GPP', stage)

def _identity_resid(d, ch4, regrow, stage='tr'):
    expect = _nee(d, ch4, stage) + _g(d, 'TKLATC', stage) + _g(d, 'TKRESPC', stage)
    if regrow:
        expect = expect - _g(d, 'TKREGROWC', stage)
    return _g(d, 'TKNETC', stage) - expect

@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    w = tmp_path_factory.mktemp('tknetc')
    yearly = {**OUTS, 'TKNETC': 'y'}
    return {
        'base': r.run(BIN, str(w / 'base'), False, outputs=OUTS),
        'rts': r.run(BIN, str(w / 'rts'), True, outputs=OUTS, extra=r.REGROW_ON),
        'rts_ch4': r.run(BIN, str(w / 'rts_ch4'), True, outputs=OUTS,
                         extra={**r.REGROW_ON, 'thermokarst_ch4': True}),
        'rts_yearly': r.run(BIN, str(w / 'rts_yearly'), True, outputs=yearly, extra=r.REGROW_ON),
    }

def test_tknetc_identity_from_components(runs):
    """(a) TKNETC recomputed from the component outputs matches TKNETC itself."""
    assert np.max(np.abs(_identity_resid(runs['rts'], ch4=False, regrow=True))) < 1e-9
    assert _g(runs['rts'], 'TKRESPC').max() > 0 and _g(runs['rts'], 'TKREGROWC').sum() > 0

def test_tknetc_ch4_counted_once(runs):
    d = runs['rts_ch4']
    # TKRHCO2 + TKCH4EMIS == RHSOM: CH4 C is neither double-counted nor dropped
    assert np.max(np.abs(_g(d, 'TKRHCO2') + _g(d, 'TKCH4EMIS') - _g(d, 'RHSOM'))) < 1e-9
    assert _g(d, 'TKCH4EMIS').sum() > 0
    assert np.max(np.abs(_identity_resid(d, ch4=True, regrow=True))) < 1e-9
    # the double-counting formula RHSOM+TKCH4EMIS would differ by the emission itself
    wrong = _identity_resid(d, ch4=True, regrow=True) + _g(d, 'TKCH4EMIS')
    assert np.max(np.abs(wrong)) > 0.1 * _g(d, 'TKCH4EMIS').max()

def test_tknetc_baseline_equals_nee(runs):
    """(c) No-RTS baseline: TK terms are zero and TKNETC == NEE_month."""
    d = runs['base']
    for v in ('TKLATC', 'TKRESPC', 'TKREGROWC'):
        assert np.all(_g(d, v) == 0), v
    assert np.max(np.abs(_g(d, 'TKNETC') - _nee(d, ch4=False))) < 1e-9

def test_tknetc_monthly_and_cumulative_closure_synthetic(runs):
    """(b) Monthly dECOC = -TKNETC (<0.01 gC/m2/month, as the suite's other
    closures) and cumulative TKNETC = -(ECOC_end-ECOC_start). Achieved on this
    case: monthly max and cumulative diff are asserted in the test output via
    the tight bounds below (cumulative is the sum of the monthly residuals)."""
    for k, ch4 in (('rts', False), ('rts_ch4', True)):
        d = runs[k]; C, T = _ecoc(d), _g(d, 'TKNETC')
        monthly = np.diff(C) + T[1:]
        assert np.max(np.abs(monthly)) < 0.01, (k, np.max(np.abs(monthly)))
        cum = T[1:].sum() - (-(C[-1] - C[0]))
        assert abs(cum) < 1.0, (k, cum)   # same cumulative bound as test_carbon_balance_after_initiation (achieved: 0.0053)

def test_tknetc_yearly_is_sum_of_monthly(runs):
    y = r.read(runs['rts_yearly'], 'TKNETC', 'yearly').ravel()
    m = _g(runs['rts'], 'TKNETC').reshape(len(y), 12).sum(axis=1)
    assert np.max(np.abs(y - m)) < 1e-6

def _gydan_regrow(tmp_path):
    w = str(tmp_path / 'c70rg'); os.makedirs(w)
    sf = os.path.join(w, 'output_spec.csv')
    r.write_outspec(sf, OUTS)
    ov = json.load(open(os.path.join(GYDAN70, 'two_stage_overrides.json')))
    for k, v in ov['IO'].items():
        if isinstance(v, str) and v.endswith('.nc'): ov['IO'][k] = os.path.join(GYDAN70, v)
    ov['IO']['output_spec_file'] = sf; ov['IO']['parameter_dir'] = os.path.join(r.REPO, 'parameters/')
    ov['model_settings'].update({'thermokarst_ch4': True, 'thermokarst_regrow': True})
    json.dump(ov, open(os.path.join(w, 'ov.json'), 'w'))
    p = subprocess.run([sys.executable, os.path.join(r.REPO, 'scripts', 'run_rts_two_stage.py'),
                        '--binary', BIN, '--workdir', w, '--config-overrides', os.path.join(w, 'ov.json'),
                        '--hist', '76', '--sc', '30', '--mode', 'rts'],
                       capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]
    return os.path.join(w, 'stage2b_rts')

@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_tknetc_gydan_cmt70_regrow_closure(tmp_path):
    """(b) Gydan CMT70 single-cell two-stage regrow (+CH4) run: TKNETC identity,
    monthly closure <0.01 gC/m2/month and cumulative closure vs ECOC change."""
    d = _gydan_regrow(tmp_path)
    assert _g(d, 'TKRESPC', 'sc').max() > 0 and _g(d, 'TKREGROWC', 'sc').sum() > 0
    assert np.max(np.abs(_identity_resid(d, ch4=True, regrow=True, stage='sc'))) < 1e-9
    C, T = _ecoc(d, 'sc'), _g(d, 'TKNETC', 'sc')
    monthly = np.diff(C) + T[1:]
    assert np.max(np.abs(monthly)) < 0.01, np.max(np.abs(monthly))
    cum = T[1:].sum() - (-(C[-1] - C[0]))
    assert abs(cum) < 1.0, cum   # achieved: -0.0087 (monthly max 0.00104)
