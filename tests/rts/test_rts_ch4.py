"""RTS CH4 (rts-light v0.2 prototypes a/b): switch, area blend and C balance incl. CH4.
Run: make rts-test (Gydan CMT70 test skipped when the inputs are absent)."""
import json, os, subprocess, sys
import numpy as np
import pytest
import rtslib as r

BIN = os.environ.get('DVMDOSTEM_BIN', os.path.join(r.REPO, 'dvmdostem'))
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference')
sys.path.insert(0, os.path.join(r.REPO, 'scripts'))
import rts_area_blend as blend
GYDAN70 = os.environ.get('RTS_GYDAN70_DIR', '/workspace/russia-inputs/inputs/gydan_cmt70')
CH4 = {'thermokarst_ch4': True}

def test_ch4_switch_off_identical_to_v01(tmp_path):
    """thermokarst_ch4=false given explicitly (v0.2 floor off): RTS-on synthetic run == stored v0.1 reference, exact."""
    o = r.run(BIN, str(tmp_path / 'off'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS},
              extra={**r.V01_OFF, 'thermokarst_ch4': False})
    assert r.compare_dirs(os.path.join(REF, 'rts'), o, 0.0, 0.0) == []

def test_area_blend_zero_is_baseline():
    rng = np.random.default_rng(1)
    b, x = rng.normal(size=48), rng.normal(size=48)
    assert np.array_equal(blend.blend(b, x, 0.0), b)
    assert np.array_equal(blend.blend(b, x, 1.0), x)
    assert np.allclose(blend.blend(b, x, 0.25), 0.75 * b + 0.25 * x, rtol=0, atol=1e-15)
    d = blend.blend({'a': b}, {'a': x}, 0.0); assert np.array_equal(d['a'], b)
    with pytest.raises(ValueError): blend.blend(b, x, 1.5)

def _two_stage(tmp_path, mode, stage1=None):
    w = str(tmp_path / ('c70_' + mode)); os.makedirs(w)
    sf = os.path.join(w, 'output_spec.csv')
    r.write_outspec(sf, {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.CH4_OUTPUTS})
    ov = json.load(open(os.path.join(GYDAN70, 'two_stage_overrides.json')))
    for k, v in ov['IO'].items():
        if isinstance(v, str) and v.endswith('.nc'): ov['IO'][k] = os.path.join(GYDAN70, v)
    ov['IO']['output_spec_file'] = sf; ov['IO']['parameter_dir'] = os.path.join(r.REPO, 'parameters/')
    ov['model_settings'].update(CH4)
    json.dump(ov, open(os.path.join(w, 'ov.json'), 'w'))
    cmd = [sys.executable, os.path.join(r.REPO, 'scripts', 'run_rts_two_stage.py'), '--binary', BIN,
           '--workdir', w, '--config-overrides', os.path.join(w, 'ov.json'), '--hist', '76', '--sc', '30',
           '--mode', mode]
    if stage1: cmd += ['--stage1-restart', stage1]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]
    return w

def _g(d, v):
    return r.read(d, v, 'monthly', 'sc').ravel()

def _balance(d):
    C = _g(d, 'VEGC') + _g(d, 'SOC') + _g(d, 'DEADC') + _g(d, 'DWDC') + _g(d, 'TKCH4POOL')
    F = (_g(d, 'GPP') - _g(d, 'RG') - _g(d, 'RM') - _g(d, 'TKRHCO2') - _g(d, 'RHDWD')
         - _g(d, 'TKCH4EMIS') - _g(d, 'TKLATC') - _g(d, 'TKRESPC'))
    return np.diff(C) - F[1:]

@pytest.fixture(scope='module')
def c70(tmp_path_factory):
    if not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')):
        pytest.skip('Gydan CMT70 inputs not available')
    t = tmp_path_factory.mktemp('ch4c70')
    rts = _two_stage(t, 'rts')
    base = _two_stage(t, 'baseline', os.path.join(rts, 'stage1', 'restart-sp.nc'))
    return {'rts': os.path.join(rts, 'stage2b_rts'), 'base': os.path.join(base, 'stage2b_baseline')}

def test_cmt70_rts_on_c_balance_with_ch4(c70):
    """CMT70 two-stage, RTS + CH4 on: monthly C closure incl. CH4 < 0.01 g/m2/month."""
    for k in ('rts', 'base'):
        res = _balance(c70[k])
        assert np.all(np.isfinite(res)) and np.max(np.abs(res)) < 0.01, (k, np.max(np.abs(res)))
    assert _g(c70['rts'], 'TKRESPC').max() > 0                       # the slump happened
    assert _g(c70['rts'], 'TKCH4EMIS').sum() > 0 and _g(c70['base'], 'TKCH4EMIS').sum() > 0

def test_cmt70_ch4_partition_consistency(c70):
    for k in ('rts', 'base'):
        d = c70[k]
        prod, ox, em = _g(d, 'TKCH4PROD'), _g(d, 'TKCH4OXID'), _g(d, 'TKCH4EMIS')
        assert np.all(ox >= 0) and np.all(em >= 0) and np.all(_g(d, 'TKCH4PLANT') <= em + 1e-12)
        # monthly CH4 mass balance: d(pool) = produced - oxidized - emitted
        pool = _g(d, 'TKCH4POOL')
        assert np.max(np.abs(np.diff(pool) - (prod - ox - em)[1:])) < 1e-9
    lat, latch4 = _g(c70['rts'], 'TKLATC'), _g(c70['rts'], 'TKLATCH4')
    assert np.all(latch4 >= 0) and np.all(latch4 <= lat + 1e-12) and latch4.max() > 0

def test_cmt70_ponding_after_slump(c70):
    p, lat = _g(c70['rts'], 'TKPONDSAT'), _g(c70['rts'], 'TKLATC')
    m = int(np.argmax(lat > 0))
    assert np.all(p[:m + 1] == 0) and p[m + 1] > 0 and np.all(np.diff(p[m + 1:]) <= 0)
    assert np.all(_g(c70['base'], 'TKPONDSAT') == 0)

def test_cmt70_blend_zero_is_baseline(c70):
    for v in ('TKCH4EMIS', 'SOC', 'TKRHCO2'):
        b, x = _g(c70['base'], v), _g(c70['rts'], v)
        assert np.array_equal(blend.blend(b, x, 0.0), b)


def test_floor_switches_off_identical_to_v01(tmp_path):
    """RTS floor switches given explicitly false: RTS-on synthetic run == stored v0.1 reference, exact."""
    off = {'thermokarst_rh_cap': False, 'thermokarst_floor_export': False, 'thermokarst_floor_anox': False,
           'thermokarst_floor_rmin': False, 'thermokarst_ch4': False}
    o = r.run(BIN, str(tmp_path / 'floor_off'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS}, extra=off)
    assert r.compare_dirs(os.path.join(REF, 'rts'), o, 0.0, 0.0) == []

def test_floor_on_no_effect_without_slump(tmp_path):
    """Floor terms act only after initiation: RTS module on but no initiation (baseline) -> v0.1 baseline
    run with the same settings, exact."""
    on = {'thermokarst_rh_cap': True, 'thermokarst_floor_export': True, 'thermokarst_floor_anox': True,
          'thermokarst_floor_rmin': True}
    a = r.run(BIN, str(tmp_path / 'b_off'), False, outputs=r.BASE_OUTPUTS)
    b = r.run(BIN, str(tmp_path / 'b_on'), False, outputs=r.BASE_OUTPUTS, extra=on)
    assert r.compare_dirs(a, b, 0.0, 0.0) == []
