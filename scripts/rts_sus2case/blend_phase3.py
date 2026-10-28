#!/usr/bin/env python3
"""Phase-3 offline cohort blend of the phase-2 single-start-year runs (Gydan 10x10, SSP5-8.5).

Linear cohort model (exact for the undisturbed part, no in-model coupling):
    X_cell(t) = X_base(t) + sum_tau f_tau * D_tau(t),   D_tau(t) = X_run(c(tau))(t - 12 s) - X_base(t - 12 s)
  * tau = initiation year 2021..2100 (Aug, as the runs); c(tau) = the phase-2 start year at or below tau
    (2021, 2031, ..., 2091) and s = tau - c(tau) whole years (response time-shifted; the anomaly is zero
    before the trigger, so a shifted cohort is still bit-for-bit baseline before its own Aug).
  * f_tau = n_tau * a_event / cell_area, a_event = 1.44 ha (ARTS v1.0 mean area, doi:10.18739/A2PK0738B);
    each patch carries its full area from initiation and keeps it (no repeat slumps, standing area only
    accumulates). Constant weights after initiation keep carbon closure exact: the blend of closed runs
    is closed (dECOC = -TKNETC month by month).
Area rules (per Jia et al. 2026 period k = 2030/2050/2070/2090, doi:10.1038/s41597-026-08284-y,
windows 2021-2040 / 2041-2060 / 2061-2080 / 2081-2100, 20 yr each):
  * Case A: F_k = mf (s_k - 0.5)/0.5 (cell-mean s; gt050 as sensitivity), mf 7.5e-4 [1.3e-4, 2e-3];
    expected initiations N_k = F_k cell_area / a_event over the window.
  * Case B: A0 = Nesterova et al. 2025 count (doi:10.1594/PANGAEA.974406) x 1.44 ha / cell_area as a
    legacy cohort initiating Aug 2021 (pre-2016 scars); new initiations RATE_k = count x r_k / 20 yr,
    r_k = SSP5 gt050_k / gt050_2030, 3x3 smoothed (cell-table case_b_ratio_s3_k).
Timing: T1 uniform hazard within each window (deterministic = expected counts; stochastic = Poisson
draws per year, fixed seeds); T2 hazard proportional to the exceedance of the baseline ALD over its
1950-1979 maximum (stage 2a), renormalised per window (uniform fallback if no exceedance).
"""
import argparse, json, os, sys, time
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from rts_sus2case.check_phase2 import read, outdir  # noqa: E402

FIRST, LAST, AUG = 2016, 2100, 7
NY = LAST - FIRST + 1; NM = 12 * NY
RUN_YEARS = (2021, 2031, 2041, 2051, 2061, 2071, 2081, 2091)
TAUS = np.arange(2021, LAST + 1)
WINDOWS = {2030: (2021, 2040), 2050: (2041, 2060), 2070: (2061, 2080), 2090: (2081, 2100)}
ANCHOR_M2 = 1.44e4
MF = dict(central=7.5e-4, low=1.3e-4, high=2.0e-3)
MF_TABLE = 7.5e-4
SEED = 20261028
TOL_M, TOL_CUM = 0.01, 1.0
VARS = ('ECOC', 'TKNETC', 'TKLATC', 'NEE')
POOLS = ('ECOC',)


# ---------------------------------------------------------------- series and kernels
def run_series(d):
    g = lambda v: (lambda a: None if a is None else a.reshape(a.shape[0], -1).sum(axis=1))(read(d, v))
    pool = g('TKCH4POOL')
    ecoc = g('VEGC') + g('SOC') + g('DEADC') + g('DWDC') + (0 if pool is None else pool)
    soil = g('TKRHCO2') + g('TKCH4EMIS') if g('TKRHCO2') is not None else g('RHSOM')
    nee = g('RG') + g('RM') + soil + g('RHDWD') - g('GPP')
    lat = g('TKLATC')
    return dict(ECOC=ecoc, TKNETC=g('TKNETC'), TKLATC=np.zeros(NM) if lat is None else lat, NEE=nee)


def shift(x, months):
    if months == 0: return x.copy()
    out = np.zeros_like(x); out[months:] = x[:-months]; return out


def run_year_for(tau):
    return max(y for y in RUN_YEARS if y <= tau)


def kernels(base, runs, taus=TAUS):
    """K[v] (len(taus), NM): anomaly of a unit-area cohort initiating in Aug tau."""
    K = {}
    for v in base:
        D = {y: runs[y][v] - base[v] for y in runs}
        K[v] = np.array([shift(D[run_year_for(t)], 12 * (t - run_year_for(t))) for t in taus])
    return K


def blend(base, K, frac):
    """frac (..., len(taus)) area fractions per initiation year -> blended series (..., NM)."""
    frac = np.asarray(frac, float)
    return {v: base[v] + frac @ K[v] for v in base}


def area_series(frac, taus=TAUS):
    """End-of-year disturbed fraction (..., NY)."""
    yrs = np.arange(FIRST, LAST + 1)
    M = (taus[:, None] <= yrs[None, :]).astype(float)
    return np.asarray(frac, float) @ M


def closure(ecoc, tknetc):
    m = np.diff(ecoc, axis=-1) + tknetc[..., 1:]
    return float(np.max(np.abs(m))), float(np.max(np.abs(tknetc[..., 1:].sum(-1) + ecoc[..., -1] - ecoc[..., 0])))


def yearly(x, pool):
    x = x.reshape(x.shape[:-1] + (NY, 12))
    return x[..., -1] if pool else x.sum(-1)


# ---------------------------------------------------------------- area rules and timing
def window_rates(N_by_k):
    """Expected initiations per year (len(TAUS)) from window totals {k: N_k}, uniform in window."""
    lam = np.zeros(len(TAUS))
    for k, (a, b) in WINDOWS.items():
        w = (TAUS >= a) & (TAUS <= b); lam[w] += N_by_k.get(k, 0.0) / w.sum()
    return lam


def t2_rates(N_by_k, ald_base, ald_ref):
    """T2: hazard ~ max(ALD_base(tau) - ald_ref, 0), renormalised per window; uniform fallback."""
    lam = np.zeros(len(TAUS)); fb = 0
    e = np.maximum(np.asarray(ald_base)[TAUS - FIRST] - ald_ref, 0.0)
    for k, (a, b) in WINDOWS.items():
        w = (TAUS >= a) & (TAUS <= b); N = N_by_k.get(k, 0.0)
        if e[w].sum() > 0: lam[w] = N * e[w] / e[w].sum()
        else:
            lam[w] = N / w.sum(); fb += int(N > 0)
    return lam, fb


def case_a_N(row, sus='mean', mf=MF['central']):
    ca = row['cell_area_km2'] * 1e6
    return {k: max(float(np.nan_to_num(row[f'case_a_F_{sus}_{k}'])), 0.0) * (mf / MF_TABLE) * ca / ANCHOR_M2
            for k in WINDOWS}


def case_b_N(row):
    n = float(np.nan_to_num(row['nesterova_count']))
    return {k: n * float(np.nan_to_num(row[f'case_b_ratio_s3_{k}'])) for k in WINDOWS}, n


def poisson(lam, members, case_code, cell_idx):
    return np.array([np.random.default_rng([SEED, case_code, cell_idx, m]).poisson(lam) for m in range(members)],
                    dtype=float)


def scenarios(row, cell_idx, ald_base, ald_ref, members):
    """name -> (folr tag, counts (M, len(TAUS)))."""
    out = {}; fbs = {}
    NA = case_a_N(row); NB, n0 = case_b_N(row)
    legacy = np.zeros(len(TAUS)); legacy[0] = n0
    lamA, lamB = window_rates(NA), window_rates(NB)
    for f in ('f08', 'f05'):
        out[f'A_T1_{f}'] = (f, poisson(lamA, members, 1, cell_idx))
        out[f'B_T1_{f}'] = (f, legacy + poisson(lamB, members, 2, cell_idx))
        out[f'A_det_{f}'] = (f, lamA[None])
        out[f'B_det_{f}'] = (f, (legacy + lamB)[None])
    for b in ('low', 'high'):
        out[f'A_det_mf{b}_f08'] = ('f08', window_rates(case_a_N(row, mf=MF[b]))[None])
    out['A_det_gt050_f08'] = ('f08', window_rates(case_a_N(row, 'gt050'))[None])
    out['A_T1_gt050_f08'] = ('f08', poisson(window_rates(case_a_N(row, 'gt050')), members, 3, cell_idx))
    la, fbs['A'] = t2_rates(NA, ald_base, ald_ref); lb, fbs['B'] = t2_rates(NB, ald_base, ald_ref)
    out['A_T2_f08'] = ('f08', la[None]); out['B_T2_f08'] = ('f08', (legacy + lb)[None])
    out['B_legacy_f08'] = ('f08', legacy[None])
    return out, fbs


# ---------------------------------------------------------------- per cell
def load_cell(cdir):
    S = {}
    for f in ('f08', 'f05'):
        base = run_series(outdir(os.path.join(cdir, f'base_{f}')))
        runs = {y: run_series(outdir(os.path.join(cdir, f'y{y}_{f}'))) for y in RUN_YEARS}
        S[f] = (base, runs)
    ald = read(outdir(os.path.join(cdir, 'base_f08')), 'ALD', 'yearly').reshape(-1)
    hist = os.path.join(cdir, 'stage2a', 'stage2a_hist', 'ALD_yearly_tr.nc')
    import netCDF4
    with netCDF4.Dataset(hist) as ds: ald_h = np.array(ds['ALD'][:], float).reshape(-1)
    return S, ald, ald_h


def do_cell(args):
    cdir, row, cell_idx, members = args
    S, ald, ald_h = load_cell(cdir)
    ald_ref = float(ald_h[:30].max())
    K = {f: kernels(*S[f]) for f in S}
    sc, fbs = scenarios(row, cell_idx, ald, ald_ref, members)
    a_ev = ANCHOR_M2 / (row['cell_area_km2'] * 1e6); ca = row['cell_area_km2'] * 1e6
    res = dict(cell=row['cell'], ald_ref=ald_ref, t2_fallback=fbs, sc={})
    worst = [0.0, 0.0]
    for name, (f, counts) in sc.items():
        base = S[f][0]; frac = counts * a_ev
        Bl = blend(base, K[f], frac)
        cm, cc = closure(Bl['ECOC'], Bl['TKNETC']); worst = [max(worst[0], cm), max(worst[1], cc)]
        r = dict(area_frac=area_series(frac), n_events=counts.sum(-1))
        for v in VARS:
            anom = Bl[v] - base[v]
            r[f'd{v}'] = yearly(anom, v in POOLS) * ca          # g per cell, per year (pool: Dec)
        r['TKNETC_abs'] = yearly(Bl['TKNETC'], False) * ca
        r['ECOC_abs'] = yearly(Bl['ECOC'], True) * ca
        res['sc'][name] = r
    res['closure'] = dict(c_month_max=worst[0], c_cum_absmax=worst[1])
    for f in S:
        base, runs = S[f]
        res[f'base_{f}'] = dict(ECOC_2016=float(base['ECOC'][0] * ca), ECOC_2100=float(base['ECOC'][-1] * ca),
                                TKNETC_cum=float(base['TKNETC'].sum() * ca))
        e = runs[2021]['ECOC'][-1] - base['ECOC'][-1]
        res[f'y2021_effect_{f}'] = dict(dECOC_2100_gm2=float(e),
                                        base_dECOC_gm2=float(base['ECOC'][-1] - base['ECOC'][0]))
    res['base_ALT_2091_2100'] = float(ald[-10:].mean())
    return res


# ---------------------------------------------------------------- domain
def summarize(results, table):
    names = results[0]['sc'].keys(); D = {}
    for n in names:
        acc = {}
        for r in results:
            ca_km2 = table[r['cell']]['cell_area_km2']
            for k, v in r['sc'][n].items():
                vv = v * ca_km2 if k == 'area_frac' else v          # area -> km2
                acc[k] = acc.get(k, 0) + vv
        D[n] = acc
    return D


def stats(x):
    x = np.atleast_2d(x)
    return dict(mean=float(x.mean()), p05=float(np.percentile(x, 5)), p50=float(np.median(x)),
                p95=float(np.percentile(x, 95)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root'); ap.add_argument('table'); ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--members', type=int, default=50); ap.add_argument('-j', '--jobs', type=int, default=16)
    ap.add_argument('--cells', default=None)
    a = ap.parse_args()
    import csv
    rows = list(csv.DictReader(open(a.table)))
    for r in rows:
        for k in r:
            if k != 'cell': r[k] = float(r[k]) if r[k] not in ('', 'nan') else float('nan')
    table = {r['cell']: r for r in rows}
    cells = a.cells.split(',') if a.cells else [r['cell'] for r in rows]
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    jobs = [(os.path.join(a.root, c), table[c], int(c[1:]), a.members) for c in cells]
    with Pool(a.jobs) as p: results = p.map(do_cell, jobs)
    t1 = time.time()
    D = summarize(results, table)
    np.savez_compressed(os.path.join(a.out, 'domain_series.npz'),
                        **{f'{n}__{k}': v for n, d in D.items() for k, v in d.items()})
    percell = {r['cell']: {**{k: r[k] for k in r if k != 'sc'},
                           'sc': {n: dict(area_frac_2100=float(np.mean(v['area_frac'][..., -1])),
                                          dTKNETC_cum_g=float(np.mean(v['dTKNETC'].sum(-1))),
                                          dECOC_2100_g=float(np.mean(v['dECOC'][..., -1])),
                                          TKLATC_cum_g=float(np.mean(v['dTKLATC'].sum(-1))),
                                          n_events=float(np.mean(v['n_events'])))
                                  for n, v in r['sc'].items()}} for r in results}
    json.dump(percell, open(os.path.join(a.out, 'percell.json'), 'w'), indent=1)
    S = {}
    for n, d in D.items():
        S[n] = dict(area_km2_2100=stats(d['area_frac'][..., -1]), n_events=stats(d['n_events']),
                    TKLATC_cum_Tg=stats(d['dTKLATC'].sum(-1) / 1e12),
                    dTKNETC_cum_Tg=stats(d['dTKNETC'].sum(-1) / 1e12),
                    dNEE_cum_Tg=stats(d['dNEE'].sum(-1) / 1e12),
                    dECOC_2100_Tg=stats(d['dECOC'][..., -1] / 1e12))
    base = {f: dict(ECOC_2016_Tg=sum(r[f'base_{f}']['ECOC_2016'] for r in results) / 1e12,
                    ECOC_2100_Tg=sum(r[f'base_{f}']['ECOC_2100'] for r in results) / 1e12,
                    TKNETC_cum_Tg=sum(r[f'base_{f}']['TKNETC_cum'] for r in results) / 1e12) for f in ('f08', 'f05')}
    pos = [r['cell'] for r in results if r['y2021_effect_f08']['dECOC_2100_gm2'] > 0]
    pos_alt = [c for c in pos if percell[c]['base_ALT_2091_2100'] >= 3.0]
    alt3 = [r['cell'] for r in results if r['base_ALT_2091_2100'] >= 3.0]
    summ = dict(cells=len(results), members=a.members, seconds_cells=t1 - t0, seconds_total=time.time() - t0,
                closure=dict(c_month_max=max(r['closure']['c_month_max'] for r in results),
                             c_cum_absmax=max(r['closure']['c_cum_absmax'] for r in results)),
                t2_fallback_windows=dict(A=sum(r['t2_fallback']['A'] for r in results),
                                         B=sum(r['t2_fallback']['B'] for r in results)),
                baseline=base, scenarios=S,
                positive_effect=dict(cells=pos, n=len(pos), with_base_alt_ge3=pos_alt, n_alt=len(pos_alt),
                                     cells_base_alt_ge3=len(alt3)))
    json.dump(summ, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1)
    print(json.dumps(dict(closure=summ['closure'], seconds=summ['seconds_total'],
                          positive=summ['positive_effect']['n'], pos_alt=summ['positive_effect']['n_alt']), indent=1))


if __name__ == '__main__':
    main()
