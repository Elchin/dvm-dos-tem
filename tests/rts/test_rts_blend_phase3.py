"""Tests for scripts/rts_sus2case/blend_phase3.py (phase-3 offline cohort blend).

Synthetic closed runs for the unit tests; set RTS_PHASE2_ROOT (run_phase2.py tree) to also run the
checks on a real cell (c86). rts_postproc (RTS_POSTPROC_PATH, default /workspace/rts-postproc) is used
as an independent reference when importable.
"""
import os, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
from rts_sus2case import blend_phase3 as bp  # noqa: E402

NM = bp.NM


def closed(seed, trig=None, base=None):
    rng = np.random.default_rng(seed)
    t = rng.normal(0, 5, NM)
    if base is not None:
        t = base['TKNETC'].copy(); k = (trig - bp.FIRST) * 12 + bp.AUG
        t[k] += 3000.0; t[k + 1:] += rng.normal(2, 1, NM - k - 1)
    e = np.empty(NM); e[0] = 8e4; e[1:] = e[0] - np.cumsum(t[1:])
    lat = np.zeros(NM)
    if base is not None: lat[k] = 2500.0
    nee = t - lat
    return dict(ECOC=e, TKNETC=t, TKLATC=lat, NEE=nee)


@pytest.fixture(scope='module')
def synth():
    base = closed(1)
    runs = {y: closed(10 + y, y, base) for y in bp.RUN_YEARS}
    return base, runs, bp.kernels(base, runs)


def test_single_2021_cohort_A025_reproduces_direct_run(synth):
    base, runs, K = synth
    frac = np.zeros(len(bp.TAUS)); frac[0] = 0.25
    Bl = bp.blend(base, K, frac)
    for v in base:
        direct = 0.75 * base[v] + 0.25 * runs[2021][v]
        np.testing.assert_allclose(Bl[v], direct, rtol=1e-13, atol=1e-9)
    assert np.allclose(bp.area_series(frac)[5:], 0.25) and np.all(bp.area_series(frac)[:5] == 0)


def test_matches_rts_postproc_blend(synth):
    sys.path.insert(0, os.environ.get('RTS_POSTPROC_PATH', '/workspace/rts-postproc'))
    core = pytest.importorskip('rts_postproc.core')
    xr = pytest.importorskip('xarray')
    base, runs, K = synth
    frac = np.zeros(len(bp.TAUS)); frac[0] = 0.25
    Bl = bp.blend(base, K, frac)
    ds = lambda s: xr.Dataset({v: (('time',), s[v]) for v in s})
    ref = core.blend(ds(base), ds(runs[2021]), np.full(bp.NY, 0.25),
                     pools=('ECOC',), fluxes=('TKNETC', 'TKLATC', 'NEE'))
    for v in base:
        np.testing.assert_allclose(Bl[v], ref[v].values, rtol=1e-12, atol=1e-8)


def test_A0_is_baseline_bit_for_bit(synth):
    base, runs, K = synth
    Bl = bp.blend(base, K, np.zeros((3, len(bp.TAUS))))
    for v in base:
        for m in range(3): assert np.array_equal(Bl[v][m], base[v])


def test_shifted_kernel(synth):
    base, runs, K = synth
    i = list(bp.TAUS).index(2025); k = (2025 - bp.FIRST) * 12 + bp.AUG
    D = runs[2021]['ECOC'] - base['ECOC']
    assert np.all(K['ECOC'][i][:k] == 0) and K['ECOC'][i][k] != 0
    assert np.array_equal(K['ECOC'][i][48:], D[:-48])
    assert bp.run_year_for(2030) == 2021 and bp.run_year_for(2031) == 2031 and bp.run_year_for(2100) == 2091


def test_carbon_closure_over_cohorts(synth):
    base, runs, K = synth
    rng = np.random.default_rng(3)
    frac = rng.poisson(2.0, (5, len(bp.TAUS))) * 1e-5
    frac[:, 0] += 0.01
    Bl = bp.blend(base, K, frac)
    cm, cc = bp.closure(Bl['ECOC'], Bl['TKNETC'])
    assert cm < bp.TOL_M and cc < bp.TOL_CUM
    # identity also holds for the anomaly alone (closure is linear)
    am, ac = bp.closure(Bl['ECOC'] - base['ECOC'], Bl['TKNETC'] - base['TKNETC'])
    assert am < 1e-6 and ac < 1e-6


def test_rates_and_seeds():
    N = {2030: 40.0, 2050: 0.0, 2070: 10.0, 2090: 20.0}
    lam = bp.window_rates(N)
    assert lam.sum() == pytest.approx(70.0) and lam[bp.TAUS == 2025][0] == pytest.approx(2.0)
    ald = np.linspace(0.5, 1.5, bp.NY)
    l2, fb = bp.t2_rates(N, ald, 0.9)
    assert l2.sum() == pytest.approx(70.0) and fb == 1        # 2021-2040 below the reference -> uniform
    l3, fb3 = bp.t2_rates(N, ald, 0.0)
    w = (bp.TAUS >= 2081); assert fb3 == 0 and l3[w][-1] > l3[w][0]
    a = bp.poisson(lam, 4, 1, 7); b = bp.poisson(lam, 4, 1, 7)
    assert np.array_equal(a, b) and not np.array_equal(a[0], a[1])
    assert abs(bp.poisson(lam, 2000, 1, 7).sum(1).mean() - 70.0) < 1.0


def test_case_rules():
    row = dict(cell_area_km2=1000.0, nesterova_count=10.0,
               **{f'case_a_F_mean_{k}': 1e-4 for k in bp.WINDOWS}, **{f'case_a_F_gt050_{k}': 2e-4 for k in bp.WINDOWS},
               **{f'case_b_ratio_s3_{k}': r for k, r in zip(bp.WINDOWS, (1.0, 0.5, float('nan'), 2.0))})
    NA = bp.case_a_N(row)
    assert NA[2030] == pytest.approx(1e-4 * 1e9 / 1.44e4)
    assert bp.case_a_N(row, mf=bp.MF['high'])[2030] == pytest.approx(NA[2030] * 2e-3 / 7.5e-4)
    NB, n0 = bp.case_b_N(row)
    assert n0 == 10 and NB == {2030: 10.0, 2050: 5.0, 2070: 0.0, 2090: 20.0}
    sc, fb = bp.scenarios(row, 0, np.ones(bp.NY), 0.5, 3)
    assert sc['B_T1_f08'][1].shape == (3, len(bp.TAUS)) and sc['B_T1_f08'][2] == 10
    assert sc['B_legacy_f08'][1].sum() == 0 and sc['B_legacy_f08'][2] == 10 and sc['A_T1_f08'][2] == 0
    v1, _ = bp.scenarios(row, 0, np.ones(bp.NY), 0.5, 3, legacy_pulse=True)   # v1: legacy folded into 2021
    assert v1['B_legacy_f08'][1].sum() == 10 and v1['B_legacy_f08'][2] == 0
    assert np.array_equal(v1['B_T1_f08'][1][:, 0] - 10, sc['B_T1_f08'][1][:, 0])


def test_legacy_cohort_drops_initiation_pulse(synth):
    """Legacy (pre-2016) cohort: no initiation-month response (TKLATC/TKNETC pulse, ECOC step), the response
    from month+1 on unchanged, area kept, and dECOC = -TKNETC still closes."""
    base, runs, K = synth
    KL = bp.legacy_kernel(K); m0 = (2021 - bp.FIRST) * 12 + bp.AUG
    for v in ('TKNETC', 'TKLATC', 'NEE'):
        assert np.all(KL[v][:m0 + 1] == 0) and np.array_equal(KL[v][m0 + 1:], K[v][0][m0 + 1:])
    assert K['TKLATC'][0][m0] > 0 and K['TKNETC'][0][m0] > 0                 # the full kernel has the pulse
    np.testing.assert_allclose(KL['ECOC'][m0:], K['ECOC'][0][m0:] - K['ECOC'][0][m0], rtol=0, atol=1e-9)
    assert np.all(KL['ECOC'][:m0 + 1] == 0)
    A = 0.25
    Bl = bp.blend(base, K, np.zeros(len(bp.TAUS)), A, KL)
    assert Bl['TKLATC'][m0] == base['TKLATC'][m0] and Bl['TKNETC'][m0] == base['TKNETC'][m0]
    direct = 0.75 * base['TKNETC'] + 0.25 * runs[2021]['TKNETC']               # month+1 on = direct run
    np.testing.assert_allclose(Bl['TKNETC'][m0 + 1:], direct[m0 + 1:], rtol=1e-13, atol=1e-9)
    am, ac = bp.closure(Bl['ECOC'] - base['ECOC'], Bl['TKNETC'] - base['TKNETC'])
    assert am < 1e-9 and ac < 1e-9
    full = bp.blend(base, K, np.r_[A, np.zeros(len(bp.TAUS) - 1)])
    assert (full['TKNETC'] - Bl['TKNETC']).sum() == pytest.approx(A * K['TKNETC'][0][m0], rel=1e-9)
    Z = bp.blend(base, K, np.zeros((2, len(bp.TAUS))), np.zeros(2), KL)       # zero legacy area = baseline
    for v in base: assert np.array_equal(Z[v][1], base[v])


ROOT = os.environ.get('RTS_PHASE2_ROOT')


@pytest.mark.skipif(not ROOT, reason='RTS_PHASE2_ROOT not set')
def test_real_cell_c86():
    S, ald, ald_h = bp.load_cell(os.path.join(ROOT, 'c86'))
    base, runs = S['f08']; K = bp.kernels(base, runs)
    frac = np.zeros(len(bp.TAUS)); frac[0] = 0.25
    Bl = bp.blend(base, K, frac)
    for v in base:
        np.testing.assert_allclose(Bl[v], 0.75 * base[v] + 0.25 * runs[2021][v], rtol=1e-12, atol=1e-9)
    Z = bp.blend(base, K, np.zeros(len(bp.TAUS)))
    for v in base: assert np.array_equal(Z[v], base[v])
    cm, cc = bp.closure(Bl['ECOC'], Bl['TKNETC']); assert cm < bp.TOL_M and cc < bp.TOL_CUM
    assert len(ald) == bp.NY and len(ald_h) == 66
    KL = bp.legacy_kernel(K); m0 = (2021 - bp.FIRST) * 12 + bp.AUG
    L = bp.blend(base, K, np.zeros(len(bp.TAUS)), 0.25, KL)
    assert L['TKLATC'][m0] == base['TKLATC'][m0] and K['TKLATC'][0][m0] > 1000
    cm, cc = bp.closure(L['ECOC'], L['TKNETC']); assert cm < bp.TOL_M and cc < bp.TOL_CUM


def test_plot_smoke(tmp_path):
    """plot_phase3 on a tiny synthetic blend output (2 cells)."""
    import json, csv
    from rts_sus2case import plot_phase3 as pp
    NY = bp.NY; names = ['A_T1_f08', 'B_T1_f08', 'A_T2_f08', 'B_T2_f08', 'A_det_mflow_f08', 'A_det_mfhigh_f08',
                         'A_T1_f05', 'B_T1_f05', 'A_T1_gt050_f08', 'B_legacy_f08']
    z = {}
    for n in names:
        for k in ('area_frac', 'dTKLATC', 'dTKNETC', 'dNEE', 'dECOC', 'n_events'):
            z[f'{n}__{k}'] = np.abs(np.random.default_rng(0).normal(1, 0.1, (4, NY))).cumsum(-1)
    np.savez(tmp_path / 'domain_series.npz', **z)
    pc = {c: {'sc': {n: dict(area_frac_2100=1e-4, dTKNETC_cum_g=1e9) for n in names}} for c in ('c00', 'c01')}
    json.dump(pc, open(tmp_path / 'percell.json', 'w'))
    st = dict(p50=1.0, p05=0.5, p95=2.0, mean=1.0)
    S = dict(scenarios={n: {k: st for k in ('area_km2_2100', 'n_events', 'TKLATC_cum_Tg', 'dTKNETC_cum_Tg',
                                             'dNEE_cum_Tg', 'dECOC_2100_Tg')} for n in names},
             baseline={'f08': dict(ECOC_2016_Tg=1, ECOC_2100_Tg=0.9, TKNETC_cum_Tg=0.1)}, cells=2,
             positive_effect=dict(n=1, n_alt=1, cells_base_alt_ge3=1, cells=['c00']),
             closure=dict(c_month_max=0.0, c_cum_absmax=0.0), t2_fallback_windows=dict(A=0, B=0), seconds_total=1)
    for n in names: S['scenarios'][n]['dTKNETC_cum_peak'] = dict(Tg=1.0, year=2094, from_peak_to_2100_Tg=-0.1)
    json.dump(S, open(tmp_path / 'summary.json', 'w'))
    with open(tmp_path / 't.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['cell', 'j', 'i', 'lat', 'lon', 'cell_area_km2'])
        w.writerow(['c00', 0, 0, 67.75, 75.75, 1000]); w.writerow(['c01', 0, 1, 67.75, 76.25, 1000])
    sys.argv = ['x', str(tmp_path), str(tmp_path / 't.csv'), '--compare', str(tmp_path)]; pp.main()
    for f in ('fig_p3_maps.png', 'fig_p3_timeseries.png', 'fig_p3_cumulative.png', 'summary.md'):
        assert (tmp_path / f).stat().st_size > 0
    txt = (tmp_path / 'summary.md').read_text()
    assert 'peaks in 2094' in txt and 'Change from v1' in txt and 'initiation pulse is dropped' in txt
