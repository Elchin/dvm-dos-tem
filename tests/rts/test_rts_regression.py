"""RTS unit + regression tests on the 1-cell synthetic case.  Run: make rts-test

Tolerances for comparisons against stored references are set by env
RTS_RTOL / RTS_ATOL (default exact-ish 1e-9; CI may loosen for other compilers).
"""
import os, subprocess
import numpy as np
import pytest
import rtslib as r

BIN = os.environ.get('DVMDOSTEM_BIN', os.path.join(r.REPO, 'dvmdostem'))
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference')
RTOL = float(os.environ.get('RTS_RTOL', '1e-9'))
ATOL = float(os.environ.get('RTS_ATOL', '1e-9'))
RTS_YEAR, RTS_MONTH = 10, 7          # synthetic forcing: tr year 10, jday 213 -> Aug
M = RTS_YEAR * 12 + RTS_MONTH        # monthly index of initiation

@pytest.fixture(scope='session')
def runs(tmp_path_factory):
    w = tmp_path_factory.mktemp('rts')
    return {
        'base': r.run(BIN, str(w / 'base'), False),
        'noop': r.run(BIN, str(w / 'noop'), True, tk_file='thermokarst-none.nc'),
        'rts':  r.run(BIN, str(w / 'rts'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS}),
        # ice-rich cell, module on, no RTS: excess ice stays frozen -> baseline physics
        'ice_noop': r.run(BIN, str(w / 'ice_noop'), True, tk_file='thermokarst-none.nc',
                          outputs={**r.BASE_OUTPUTS, **r.ICE_OUTPUTS}, extra=r.ICE_PARAMS),
        'rts_ice': r.run(BIN, str(w / 'rts_ice'), True,
                         outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS}, extra=r.ICE_PARAMS),
    }

def test_unit_binary():
    exe = os.path.join(r.REPO, 'tests', 'rts', 'test_rts_units')
    if not os.path.exists(exe):
        pytest.skip('unit test binary not built (make rts-unit-test)')
    p = subprocess.run([exe], capture_output=True, text=True)
    print(p.stdout)
    assert p.returncode == 0

def test_baseline_matches_reference(runs):
    """RTS off: identical to stored baseline produced by the unmodified RTS-branch binary."""
    assert r.compare_dirs(os.path.join(REF, 'baseline'), runs['base'], RTOL, ATOL) == []

def test_noop_when_no_initiation(runs):
    """Module on, but no initiation year (rts_init_year=-1 equivalent): identical to baseline."""
    assert r.compare_dirs(runs['base'], runs['noop'], 0.0, 0.0) == []

def test_rts_matches_reference(runs):
    assert r.compare_dirs(os.path.join(REF, 'rts'), runs['rts'], RTOL, ATOL) == []

def test_unchanged_before_initiation(runs):
    for v in ['SOC', 'VEGC', 'GPP', 'TLAYER']:
        a, b = r.read(runs['base'], v), r.read(runs['rts'], v)
        assert np.array_equal(a[:M], b[:M]), v

def test_organic_c_drop_and_lateral_export(runs):
    o, f = runs['rts'], r.RTS_PARAMS['thermokarst_c_resp_frac']
    soc, vegc = r.read(o, 'SOC').ravel(), r.read(o, 'VEGC').ravel()
    latc = r.read(o, 'TKLATC').ravel()
    shlw, deep = r.read(o, 'SHLWC').ravel(), r.read(o, 'DEEPC').ravel()
    assert soc[M] < 0.9 * soc[M - 1]                 # big organic C drop
    assert shlw[M] < 0.05 * shlw[M - 1]               # fibric layer mostly removed
    assert deep[M] < deep[M - 1]
    assert latc[M] > 0 and np.all(latc[:M] == 0) and np.all(latc[M + 1:] == 0)
    # removed (soil + veg) ~= lateral / (1 - resp_frac); monthly fluxes add small noise
    # NOTE: VEGC output for month M is captured before the disturbance step,
    # so the veg drop appears between M and M+1 (SOC drops at M).
    removed = (soc[M - 1] - soc[M]) + (vegc[M] - vegc[M + 1])
    assert latc[M] / (1 - f) == pytest.approx(removed, rel=0.02)
    assert r.read(o, 'TKLATN').ravel()[M] > 0
    assert r.read(o, 'TKREMTHICK').ravel()[M] > 0

def test_thaw(runs):
    ice = r.read(runs['rts'], 'TKICEMELT').ravel()
    assert ice[M] > 0 and np.all(ice[:M] == 0)
    ald_b, ald_r = r.read(runs['base'], 'ALD', 'yearly').ravel(), r.read(runs['rts'], 'ALD', 'yearly').ravel()
    assert ald_r[RTS_YEAR] > ald_b[RTS_YEAR] + 0.2   # deeper thaw in initiation year

def test_drainage_switch(runs):
    d = r.read(runs['rts'], 'DRAINTYPE').ravel()
    assert np.all(d[:M] == 0)                         # synthetic cell starts well-drained
    assert np.all(d[M:] == r.RTS_PARAMS['thermokarst_post_drainage'])


# ---------------- excess ice / subsidence (ported from thermokarst-prototype) -------------
def test_rts_ice_matches_reference(runs):
    assert r.compare_dirs(os.path.join(REF, 'rts_ice'), runs['rts_ice'], RTOL, ATOL) == []

def test_excess_ice_noop_without_melt(runs):
    """Ice-rich profile below the ALD, no RTS: no melt, no subsidence, and all
    baseline variables identical to the RTS-off baseline (excess ice is passive)."""
    o = runs['ice_noop']
    assert np.all(r.read(o, 'TKSUBSID') == 0) and np.all(r.read(o, 'TKEXMELT') == 0)
    x = r.read(o, 'TKEXICE').ravel()
    assert x[0] > 0 and np.all(x == x[0])
    bad = r.compare_dirs(runs['base'], o, 0.0, 0.0, skip=('TK', 'WATERTAB'))
    assert [b for b in bad if not b.startswith('file sets differ')] == []

def test_initialized_excess_ice_amount(runs):
    p = r.ICE_PARAMS; f = p['thermokarst_excess_fraction']
    x0 = r.read(runs['rts_ice'], 'TKEXICE').ravel()[0]
    # all mineral layers in [top,bottom] are frozen in January of tr year 0 at this site
    assert x0 == pytest.approx(917.0 * (p['thermokarst_excess_bottom'] - p['thermokarst_excess_top']) * f / (1 - f), rel=1e-6)

def test_subsidence_equals_excess_melt(runs):
    o = runs['rts_ice']
    sub, melt, x = (r.read(o, v).ravel() for v in ('TKSUBSID', 'TKEXMELT', 'TKEXICE'))
    assert np.all(sub[:M] == 0)                                  # nothing before initiation
    assert sub[M] > 0.1                                          # RTS thaw releases excess ice
    assert np.all(np.diff(sub) >= 0)                             # subsidence is monotonic
    assert sub[-1] == pytest.approx(np.cumsum(melt)[-1] / 917.0, rel=1e-9)
    assert x[0] - x[-1] == pytest.approx(np.sum(melt), rel=1e-9) # excess-ice mass balance

def test_rts_ice_carbon_and_ald(runs):
    o, b = runs['rts_ice'], runs['base']
    soc = r.read(o, 'SOC').ravel()
    assert soc[M] < 0.9 * soc[M - 1]
    ald_b, ald_r = r.read(b, 'ALD', 'yearly').ravel(), r.read(o, 'ALD', 'yearly').ravel()
    assert ald_r[RTS_YEAR] > ald_b[RTS_YEAR] + 0.2

def test_water_table_shallows_after_subsidence(runs):
    # coupling: on the step after collapse the WTD is reduced by the subsidence increment
    wt_r = r.read(runs['rts_ice'], 'WATERTAB').ravel()
    wt_n = r.read(runs['ice_noop'], 'WATERTAB').ravel()
    assert np.array_equal(wt_r[:M], wt_n[:M])
    assert wt_r[M + 1] < wt_n[M + 1]
