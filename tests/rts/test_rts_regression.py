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
# reference/{baseline,rts,rts_ice,legacy_m2d}: rts-light-v0.1 outputs (floor off);
# reference/v02/{rts,rts_ice}: v0.2 defaults (floor on, mg0 = 1), see make_reference.py
REF02 = os.path.join(REF, 'v02')
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
        # v0.2 floor switches off -> must be bitwise rts-light-v0.1
        'rts_v01': r.run(BIN, str(w / 'rts_v01'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS},
                         extra=r.V01_OFF),
        'rts_ice_v01': r.run(BIN, str(w / 'rts_ice_v01'), True,
                             outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS},
                             extra={**r.ICE_PARAMS, **r.V01_OFF}),
        # ice-rich cell, module on, no RTS: excess ice stays frozen -> baseline physics
        'ice_noop': r.run(BIN, str(w / 'ice_noop'), True, tk_file='thermokarst-none.nc',
                          outputs={**r.BASE_OUTPUTS, **r.ICE_OUTPUTS}, extra=r.ICE_PARAMS),
        'rts_ice': r.run(BIN, str(w / 'rts_ice'), True,
                         outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS}, extra=r.ICE_PARAMS),
        # legacy monthly->daily (flag off) must reproduce the pre-m2d references exactly
        'leg_base': r.run(BIN, str(w / 'leg_base'), False, extra=r.LEGACY_M2D),
        'leg_rts': r.run(BIN, str(w / 'leg_rts'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS},
                         extra={**r.LEGACY_M2D, **r.V01_OFF}),
        'leg_rts_ice': r.run(BIN, str(w / 'leg_rts_ice'), True,
                             outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS},
                             extra={**r.LEGACY_ICE_PARAMS, **r.LEGACY_M2D, **r.V01_OFF}),
        # post-slump regrowth: explicit off must be bitwise v0.2; on: C/N closure with tracked inputs
        'rts_rgoff': r.run(BIN, str(w / 'rts_rgoff'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS},
                           extra={'thermokarst_regrow': False}),
        'rts_ice_rgoff': r.run(BIN, str(w / 'rts_ice_rgoff'), True,
                               outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.ICE_OUTPUTS},
                               extra={**r.ICE_PARAMS, 'thermokarst_regrow': False}),
        'rts_rg': r.run(BIN, str(w / 'rts_rg'), True, outputs={**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **r.REGROW_OUTPUTS},
                        extra=r.REGROW_ON),
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
    """v0.2 defaults (floor on, mg0 = 1) vs the v0.2 reference."""
    assert r.compare_dirs(os.path.join(REF02, 'rts'), runs['rts'], RTOL, ATOL) == []

@pytest.mark.parametrize('case', ['rts', 'rts_ice'])
def test_floor_off_bit_identical_to_v01(runs, case):
    """All v0.2 floor switches off + mg0 = 0.05: bitwise identical to the rts-light-v0.1 references."""
    assert r.compare_dirs(os.path.join(REF, case), runs[case + '_v01'], 0.0, 0.0) == []

def test_unchanged_before_initiation(runs):
    for v in ['SOC', 'VEGC', 'GPP', 'TLAYER']:
        a, b = r.read(runs['base'], v), r.read(runs['rts'], v)
        assert np.array_equal(a[:M], b[:M]), v

@pytest.mark.parametrize('run', ['rts_v01', 'rts'])
def test_organic_c_drop_and_lateral_export(runs, run):
    o, f = runs[run], r.RTS_PARAMS['thermokarst_c_resp_frac']
    soc, vegc = r.read(o, 'SOC').ravel(), r.read(o, 'VEGC').ravel()
    latc = r.read(o, 'TKLATC').ravel()
    shlw, deep = r.read(o, 'SHLWC').ravel(), r.read(o, 'DEEPC').ravel()
    assert soc[M] < 0.9 * soc[M - 1]                 # big organic C drop
    assert shlw[M] < 0.05 * shlw[M - 1]               # fibric layer mostly removed
    assert deep[M] < deep[M - 1]
    assert latc[M] > 0 and np.all(latc[:M] == 0)
    if run == 'rts_v01':
        assert np.all(latc[M + 1:] == 0)
    else:                                              # v0.2 floor export after the slump month
        assert np.all(latc[M + 1:] > 0) and np.all(latc[M + 1:] < latc[M])
    # removed (soil + veg) ~= lateral / (1 - resp_frac); monthly fluxes add small noise
    # VEGC and SOC both drop in month M (VEGC used to lag one month)
    removed = (soc[M - 1] - soc[M]) + (vegc[M - 1] - vegc[M])
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
    assert r.compare_dirs(os.path.join(REF02, 'rts_ice'), runs['rts_ice'], RTOL, ATOL) == []

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


# ---------------- carbon fluxes: RH / RA / GPP / NPP / balance -------------
FLUXES = ['GPP', 'NPP', 'RG', 'RM', 'RHSOM', 'RHDWD', 'RECO']

def _g(o, v):
    return r.read(o, v).ravel()

def _ra(o):
    return _g(o, 'RG') + _g(o, 'RM')

def test_fluxes_present_and_baseline_identical(runs):
    # RTS off (and module on without initiation) bit-identical for C fluxes
    for v in FLUXES:
        assert np.array_equal(_g(runs['base'], v), _g(runs['noop'], v)), v
    bad = r.compare_dirs(os.path.join(REF, 'baseline'), runs['base'], RTOL, ATOL)
    assert bad == []

@pytest.mark.parametrize('run', ['rts_v01', 'rts'])
def test_rh_drops_after_initiation(runs, run):
    rh_b, rh_r = _g(runs['base'], 'RHSOM'), _g(runs[run], 'RHSOM')
    # same calendar months one year after initiation: Aug-Sep of year 11
    win = slice(M + 12, M + 14)
    assert rh_r[win].sum() < 0.5 * rh_b[win].sum()
    if run == 'rts_v01':
        assert rh_r[M + 1] < 0.2 * rh_b[M + 1]
    else:   # v0.2 floor minimum respiration: >= rmin0 (0.3) x 5-yr pre-slump mean of the same month
        pre = np.mean([rh_r[M + 1 - 12 * k] for k in range(1, 6)])
        assert 0.3 * pre - 1e-9 <= rh_r[M + 1] < 0.5 * rh_b[M + 1]

def test_ra_drops_when_vegetation_removed(runs):
    ra_b, ra_r = _ra(runs['base']), _ra(runs['rts'])
    assert ra_r[M + 1] < 0.2 * ra_b[M + 1]

def test_gpp_minus_ra_equals_npp(runs):
    for k in ('base', 'rts'):
        o = runs[k]
        assert np.allclose(_g(o, 'GPP') - _ra(o), _g(o, 'NPP'), atol=1e-6), k

def test_reco_components(runs):
    o = runs['rts']
    assert np.allclose(_g(o, 'RECO'), _g(o, 'RHSOM') + _g(o, 'RHDWD') + _ra(o), atol=1e-6)

def test_direct_respired_pulse(runs):
    o, f = runs['rts'], r.RTS_PARAMS['thermokarst_c_resp_frac']
    resp, lat = _g(o, 'TKRESPC'), _g(o, 'TKLATC')
    assert resp[M] > 0 and np.all(resp[:M] == 0) and np.all(resp[M + 1:] == 0)
    # removed C = respired + lateral; respired = c_resp_frac * removed
    assert resp[M] == pytest.approx(f * (resp[M] + lat[M]), rel=1e-9)

def _total_c(o):
    return _g(o, 'VEGC') + _g(o, 'SOC') + _g(o, 'DEADC') + _g(o, 'DWDC')

def _net_flux(o):
    f = _g(o, 'GPP') - _ra(o) - _g(o, 'RHSOM') - _g(o, 'RHDWD')
    if os.path.exists(os.path.join(o, 'TKLATC_monthly_tr.nc')):
        f = f - _g(o, 'TKLATC') - _g(o, 'TKRESPC')
    return f

def test_carbon_balance_before_initiation(runs):
    for k in ('base', 'rts'):
        C, F = _total_c(runs[k]), _net_flux(runs[k])
        assert C[M - 1] - C[12] == pytest.approx(F[13:M].sum(), abs=0.05), k

def test_carbon_balance_at_initiation(runs):
    # window spanning the event
    o = runs['rts']; C, F = _total_c(o), _net_flux(o)
    removed = _g(o, 'TKLATC')[M] + _g(o, 'TKRESPC')[M]
    assert C[M + 1] - C[M - 1] == pytest.approx(F[M:M + 2].sum(), abs=0.01 * removed)


def test_carbon_balance_slump_month(runs):
    """Pool removal and TKLATC/TKRESPC fall in the same month (VEGC used to lag
    one month; dead-root N went to somcr; the re-created fibric layer got
    litter C out of nothing)."""
    o = runs['rts']; C, F = _total_c(o), _net_flux(o)
    assert abs((C[M] - C[M - 1]) - F[M]) < 0.01
    assert abs((C[M + 1] - C[M]) - F[M + 1]) < 0.01


def test_carbon_balance_after_initiation(runs):
    """Post-RTS balance closes (was a ~25 gC/m2/yr Jul-Sep leak: soil integrator
    'blackhole' fallback on negative mineral AVLN discarded litter input / RH)."""
    o = runs['rts']; C, F = _total_c(o), _net_flux(o)
    res = np.diff(C)[M + 1:] - F[M + 2:]
    assert np.max(np.abs(res)) < 0.01
    assert C[-1] - C[M + 1] == pytest.approx(F[M + 2:].sum(), abs=1.0)


# ---------------- monthly -> daily: legacy path + new path sanity -------------
LEG = os.path.join(REF, 'legacy_m2d')
# bit-identical locally; CI on another compiler may loosen via RTS_RTOL/RTS_ATOL
LRTOL = float(os.environ.get('RTS_RTOL', '0'))
LATOL = float(os.environ.get('RTS_ATOL', '0'))

@pytest.mark.parametrize('case,run', [('baseline', 'leg_base'), ('rts', 'leg_rts'), ('rts_ice', 'leg_rts_ice')])
def test_m2d_flag_off_bit_identical_to_legacy(runs, case, run):
    """monthly2daily_mean_preserving=false reproduces the pre-m2d references."""
    assert r.compare_dirs(os.path.join(LEG, case), runs[run], LRTOL, LATOL) == []

def test_m2d_new_path_changes_and_deepens_ald(runs):
    """Mean-preserving forcing restores lost thaw degree-days -> ALD not shallower."""
    a_new = r.read(runs['base'], 'ALD', 'yearly').ravel()
    a_old = r.read(runs['leg_base'], 'ALD', 'yearly').ravel()
    assert not np.array_equal(a_new, a_old)
    assert a_new.mean() >= a_old.mean()


# ---------------- Gydan CMT70 baseline carbon balance (optional data) ---------
GYDAN70 = os.environ.get('RTS_GYDAN70_DIR', '/workspace/russia-inputs/inputs/gydan_cmt70')

@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'vegetation.nc')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt70_baseline_carbon_balance(tmp_path):
    """RTS-off CMT70 closed with ~-11 gC/m2/month gaps in Aug/Sep: soil-integrator
    'blackhole' fallback (negative AVLN/ORGN) discarded the month's soil update."""
    o = r.run(BIN, str(tmp_path / 'c70'), False, tk_file='thermokarst-none.nc', data_dir=GYDAN70,
              extra={'baseline_start': 1940, 'baseline_end': 1970})
    C, F = _total_c(o), _net_flux(o)
    res = np.diff(C) - F[1:]
    assert np.max(np.abs(res[12:])) < 0.01


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'vegetation.nc')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt70_deepc_spinup_n_and_closure(tmp_path):
    """deep_carbon_v2 burial drains organic-layer N into frozen mineral layers;
    root N uptake beyond AVLN was then borrowed from an exhausted layer's ORGN
    (negative ORGN by eq year ~590) and an all-frozen month divided by
    totnetnmin=0 (NaN in every N pool by spin-up). Full pr/eq/sp spin-up."""
    import netCDF4 as nc4, subprocess
    w = str(tmp_path / 'c70eq')
    os.makedirs(w)
    cfgp = r.make_config(w, False, tk_file='thermokarst-none.nc', data_dir=GYDAN70,
                         extra={'baseline_start': 1940, 'baseline_end': 1970})
    p = subprocess.run([BIN, '-f', cfgp, '-l', 'fatal', '-p', '100', '-e', '1000', '-s', '200',
                        '-t', '76', '-n', '0'], cwd=r.REPO, capture_output=True, text=True, timeout=900)
    assert p.returncode == 0, p.stderr[-2000:]
    tot = {}
    for st in ('eq', 'sp'):
        d = nc4.Dataset(os.path.join(w, 'output', 'restart-%s.nc' % st))
        o = np.array(d['orgn'][:], dtype=float).reshape(-1)
        a = np.array(d['avln'][:], dtype=float).reshape(-1)
        d.close()
        o, a = o[o > -9000], a[a > -9000]
        assert o.size > 0 and np.all(np.isfinite(o)) and np.all(np.isfinite(a)), st
        assert np.all(o >= 0.0), (st, o.min())
        tot[st] = o.sum() + a.sum()
    # spin-up only redistributes/recycles soil N here: column N stays put
    assert abs(tot['sp'] - tot['eq']) < 0.02 * tot['eq'], tot
    # monthly C closure over the transient, incl. integrator-retry months (float
    # rounding accumulated over their many steps: up to 0.066 g/month before the
    # double-shadow state in Integrator::adapt)
    o = os.path.join(w, 'output')
    res = np.diff(_total_c(o)) - _net_flux(o)[1:]
    assert np.all(np.isfinite(res)) and np.max(np.abs(res)) < 0.01, np.max(np.abs(res))


def _gydan_cmt_two_stage(tmp_path, cmt, extra=None, return_dir=False, regrow=False):
    """Two-stage RTS run (slump 2021, scenario to 2045) on the Gydan CMT70 cell
    with the vegetation class set to `cmt`; asserts the monthly N balance
    (NBAL_DIAG) closes within 0.01 g and returns stage2b monthly C residuals."""
    import json, sys, shutil, netCDF4 as nc4
    data = str(tmp_path / ('in%02d' % cmt)); os.makedirs(data)
    for f in os.listdir(GYDAN70):
        if f != 'vegetation.nc':
            os.symlink(os.path.join(GYDAN70, f), os.path.join(data, f))
    shutil.copy(os.path.join(GYDAN70, 'vegetation.nc'), os.path.join(data, 'vegetation.nc'))
    v = nc4.Dataset(os.path.join(data, 'vegetation.nc'), 'a'); v['veg_class'][:] = cmt; v.close()
    w = str(tmp_path / ('c%02drts' % cmt)); os.makedirs(w)
    sf = os.path.join(w, 'output_spec.csv')
    r.write_outspec(sf, {**r.BASE_OUTPUTS, **r.RTS_OUTPUTS, **(r.REGROW_OUTPUTS if regrow else {})})
    ov = json.load(open(os.path.join(data, 'two_stage_overrides.json')))
    for k, val in ov['IO'].items():
        if isinstance(val, str) and val.endswith('.nc'):
            ov['IO'][k] = os.path.join(data, val)
    ov['IO']['output_spec_file'] = sf
    ov['IO']['parameter_dir'] = os.path.join(r.REPO, 'parameters/')
    if extra: ov['model_settings'].update(extra)
    json.dump(ov, open(os.path.join(w, 'ov.json'), 'w'))
    nbf = os.path.join(w, 'nbal.txt')
    wrap = os.path.join(w, 'nbwrap.sh')
    open(wrap, 'w').write('#!/bin/bash\nNBAL_DIAG=1 %s "$@" 2>>%s\n' % (os.path.abspath(BIN), nbf))
    os.chmod(wrap, 0o755)
    p = subprocess.run([sys.executable, os.path.join(r.REPO, 'scripts', 'run_rts_two_stage.py'),
                        '--binary', wrap, '--workdir', w, '--config-overrides', os.path.join(w, 'ov.json'),
                        '--hist', '76', '--sc', '30', '--mode', 'rts'],
                       capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]
    # monthly N balance of the scenario stage from the NBAL_DIAG check:
    # yr m N0 N1 avlninput orgninput avlnlost orgnlost tk_lat tk_resp vegin
    # (+ last column: RTS regrow N to veg/organic layer)
    nb = np.array([[float(x) for x in l.split()[1:]] for l in open(nbf) if l.startswith('NB ')])[-360:]
    _, _, n0, n1, ai, oi, al, ol, tl, tv, vi, rn = nb.T
    nres = (n1 - n0) - (ai + oi + vi + rn - al - tl - tv)
    assert np.max(np.abs(nres)) < 0.01, ('N', np.max(np.abs(nres)))
    def g(vn):
        d = nc4.Dataset(os.path.join(w, 'stage2b_rts', '%s_monthly_sc.nc' % vn))
        a = np.ma.filled(np.array(d[vn][:], dtype=float), np.nan); d.close()
        return np.nan_to_num(a.reshape(a.shape[0], -1)[:, 0])
    C = g('VEGC') + g('SOC') + g('DEADC') + g('DWDC')
    F = g('GPP') - g('RG') - g('RM') - g('RHSOM') - g('RHDWD') - g('TKLATC') - g('TKRESPC')
    if regrow: F = F + g('TKREGROWC')     # external reseed + organic-layer input
    assert g('TKRESPC').max() > 0          # the slump happened
    if return_dir: return np.diff(C) - F[1:], os.path.join(w, 'stage2b_rts')
    return np.diff(C) - F[1:]


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt21_rts_on_closure_double_blackhole(tmp_path):
    """CMT21 after the slump: some Augusts blackholed again in the retry; the
    pools stopped part-way through the month while the whole month's litter
    was reported (-11..-14 gC/month, -122 g over 30 yr)."""
    res = _gydan_cmt_two_stage(tmp_path, 21)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt70_deepc_rts_on_monthly_closure(tmp_path):
    """CMT70 two-stage RTS run (slump 2021, scenario to 2045). C: root litter
    with no receiving layer was dropped after the slump (up to 0.067 g/month).
    N: the slump restructuring moved layer C but not ORGN/AVLN (-91.5 g),
    rootless layers lost their net mineralization (-22 g/month) and uptake
    from layers without root litter was never removed (+0.3 g/month)."""
    res = _gydan_cmt_two_stage(tmp_path, 70)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt52_rts_on_no_fibric_layer(tmp_path):
    """CMT52 (heath-lichen): the slump leaves no fibric layer and, with the
    vegetation removed, no litter builds a new one. Soil_Env::updateDailySM
    dereferenced ground->fstshlwl (NULL) -> segfault the month after the slump;
    the month's litter was also dropped with no fibric layer to take it."""
    res = _gydan_cmt_two_stage(tmp_path, 52)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt55_rts_on_slump_month_root_removal(tmp_path):
    """CMT55 slump month: root fraction removed with the soil plus the retained
    living fraction exceeded 1, so TKLATC/TKRESPC counted 0.028 g more root C
    than left the root pool (negative dead-root C)."""
    res = _gydan_cmt_two_stage(tmp_path, 55)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))


FLOOR_ON = {'thermokarst_rh_cap': True, 'thermokarst_floor_export': True, 'thermokarst_floor_anox': True,
            'thermokarst_floor_rmin': True}

@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
def test_cmt70_floor_on_c_n_closure_and_rh_limit(tmp_path):
    """RTS floor (cap + export + anoxia, v0.2 defaults): monthly C and N closure < 0.01 g
    (lateral export in TKLATC/TKLATN); first post-slump year RH limited (cap0 = 4.5%/yr of floor
    mineral C, Lena chamber yr-1 160-184 gC/m2/yr) and floor export continues after the slump month."""
    import netCDF4 as nc4
    res, d = _gydan_cmt_two_stage(tmp_path, 70, extra=FLOOR_ON, return_dir=True)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))
    def g(vn):
        x = nc4.Dataset(os.path.join(d, '%s_monthly_sc.nc' % vn))
        a = np.ma.filled(np.array(x[vn][:], dtype=float), np.nan); x.close()
        return np.nan_to_num(a.reshape(a.shape[0], -1)[:, 0])
    m = int(np.argmax(g('TKRESPC') > 0))
    lat = g('TKLATC')
    assert np.all(lat[m + 1:m + 13] > 0)                        # floor export after the slump month
    assert g('RHSOM')[m + 1:m + 13].sum() < 300.0               # was > 2000 g/m2 in v0.1 (thick ice v2)


@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
@pytest.mark.parametrize('cmt', [5, 7])
def test_floor_rmin_dry_upland(tmp_path, cmt):
    """Floor minimum respiration (v0.2 default, rmin0 = 0.3): on dry CMT05/07 floors the first post-slump
    year RH is >= rmin0 x the pre-slump RH of the same months (C from existing pools), with monthly C and N
    closure < 0.01 g; with the switch off the floor RH is far lower."""
    import netCDF4 as nc4
    res, d = _gydan_cmt_two_stage(tmp_path / 'on', cmt, return_dir=True)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))
    res0, d0 = _gydan_cmt_two_stage(tmp_path / 'off', cmt, extra={'thermokarst_floor_rmin': False}, return_dir=True)
    def g(dd, vn):
        x = nc4.Dataset(os.path.join(dd, '%s_monthly_sc.nc' % vn))
        a = np.ma.filled(np.array(x[vn][:], dtype=float), np.nan); x.close()
        return np.nan_to_num(a.reshape(a.shape[0], -1)[:, 0])
    m = int(np.argmax(g(d, 'TKRESPC') > 0))
    rh = g(d, 'RHSOM')
    pre = np.mean([rh[m + 1 - 12 * k:m + 13 - 12 * k] for k in range(1, 6)], axis=0)  # 5-yr same-month mean
    yr1 = rh[m + 1:m + 13]
    assert np.all(yr1 >= 0.3 * pre - 1e-6)
    assert 0.25 <= yr1.sum() / pre.sum() <= 0.55
    assert g(d0, 'RHSOM')[m + 1:m + 13].sum() < yr1.sum()


# ---------------- post-slump vegetation regrowth (thermokarst_regrow) ----------------
@pytest.mark.parametrize('case,run', [('rts', 'rts_rgoff'), ('rts_ice', 'rts_ice_rgoff')])
def test_regrow_off_bit_identical_to_v02(runs, case, run):
    """thermokarst_regrow = false (the default): bitwise identical to the v0.2 references."""
    assert r.compare_dirs(os.path.join(REF02, case), runs[run], 0.0, 0.0) == []

def test_regrow_default_is_off(runs):
    assert r.compare_dirs(runs['rts'], runs['rts_rgoff'], 0.0, 0.0) == []

def test_regrow_synthetic_closure_and_inputs(runs):
    """Regrow on: identical before initiation; reseed + OL input in TKREGROWC only from stabilization
    (12 months after initiation); monthly C balance closes with TKREGROWC as an input."""
    o, b = runs['rts_rg'], runs['rts']
    for v in ['SOC', 'VEGC', 'GPP']:
        assert np.array_equal(r.read(o, v)[:M], r.read(b, v)[:M]), v
    rc = _g(o, 'TKREGROWC'); rn = _g(o, 'TKREGROWN')
    assert np.all(rc[:M + 12] == 0) and rc[M + 12] > 0 and np.all(rc >= 0) and np.all(rn >= 0)
    C, F = _total_c(o), _net_flux(o) + rc
    res = np.diff(C)[M + 1:] - F[M + 2:]
    assert np.max(np.abs(res)) < 0.01
    # floor N supply enters as an input (NINPUT) only after stabilization
    ni = _g(o, 'NINPUT')
    assert np.all(ni[:M + 12] == 0) and ni[M + 12:].sum() > 0

@pytest.mark.skipif(not os.path.exists(os.path.join(GYDAN70, 'two_stage_overrides.json')),
                    reason='Gydan CMT70 inputs not available')
@pytest.mark.parametrize('cmt', [70, 5, 7])
def test_regrow_gydan_closure_and_recovery(tmp_path, cmt):
    """Regrow v2 on (v0.2 floor defaults): monthly N balance (NBAL_DIAG incl. tapered floor N supply in
    avlninput and reseed/OL N column) and C balance (TKREGROWC input) close < 0.01 g; post-slump GPP
    higher than v0.2. v2 has no logistic growth limit (removed by design, flagged in commit), so
    CMT05/07 -- which recover in v0.2 already and were harmed by the v1 limit -- use the same defaults."""
    import netCDF4 as nc4
    res, d = _gydan_cmt_two_stage(tmp_path / 'rg', cmt, extra={**r.REGROW_ON},
                                  return_dir=True, regrow=True)
    assert np.max(np.abs(res)) < 0.01, np.max(np.abs(res))
    res0, d0 = _gydan_cmt_two_stage(tmp_path / 'v02', cmt, return_dir=True)
    def g(dd, vn):
        x = nc4.Dataset(os.path.join(dd, '%s_monthly_sc.nc' % vn))
        a = np.ma.filled(np.array(x[vn][:], dtype=float), np.nan); x.close()
        return np.nan_to_num(a.reshape(a.shape[0], -1)[:, 0])
    m = int(np.argmax(g(d, 'TKRESPC') > 0))
    assert np.array_equal(g(d, 'GPP')[:m + 1], g(d0, 'GPP')[:m + 1])
    assert g(d, 'TKREGROWC')[m + 12] > 0                       # reseed at stabilization
    assert g(d, 'GPP')[m + 13:].sum() > g(d0, 'GPP')[m + 13:].sum()
