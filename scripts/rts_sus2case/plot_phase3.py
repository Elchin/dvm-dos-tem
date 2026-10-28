#!/usr/bin/env python3
"""Figures + summary markdown for the phase-3 blend (blend_phase3.py output dir).

fig_p3_maps.png        2100 cumulative RTS-disturbed area and RTS net C loss (cumulative TKNETC anomaly
                       2016-2100) per cell, Case A vs Case B (T1 ensemble mean, folr 0.8)
fig_p3_timeseries.png  domain disturbed area, TKLATC, TKNETC anomaly and NEE anomaly; T1 ensemble median
                       and 5-95 % band, Case A mf band (1.3e-4 / 2e-3), T2 dotted
fig_p3_cumulative.png  cumulative RTS net C loss with the sensitivities (gt050, folr 0.5, T2, mf band)
"""
import argparse, csv, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

YRS = np.arange(2016, 2101)
COL = dict(A='tab:blue', B='tab:orange')


def load(d):
    z = np.load(os.path.join(d, 'domain_series.npz'))
    D = {}
    for k in z.files:
        n, v = k.split('__'); D.setdefault(n, {})[v] = np.atleast_2d(z[k])
    return D, json.load(open(os.path.join(d, 'percell.json'))), json.load(open(os.path.join(d, 'summary.json')))


def grid(table, percell, fn):
    js = sorted({int(float(r['j'])) for r in table}); is_ = sorted({int(float(r['i'])) for r in table})
    G = np.full((10, 10), np.nan); lat = np.full(10, np.nan); lon = np.full(10, np.nan)
    for r in table:
        j, i = int(float(r['j'])), int(float(r['i'])); lat[j] = float(r['lat']); lon[i] = float(r['lon'])
        if r['cell'] in percell: G[j, i] = fn(percell[r['cell']], r)
    return G, lat, lon


def edges(c):
    c = np.asarray(c); d = np.nanmedian(np.diff(c)) if len(c) > 1 else 0.5
    return np.concatenate([c - d / 2, [c[-1] + d / 2]])


def maps(D, P, table, out):
    fig, ax = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for r_, case in enumerate(('A', 'B')):
        n = f'{case}_T1_f08'
        G1, lat, lon = grid(table, P, lambda p, r: p['sc'][n]['area_frac_2100'] * float(r['cell_area_km2']))
        G2, _, _ = grid(table, P, lambda p, r: p['sc'][n]['dTKNETC_cum_g'] / 1e9)
        lone, late = edges(np.nan_to_num(lon, nan=np.nanmean(lon))), edges(np.nan_to_num(lat, nan=np.nanmean(lat)))
        m1 = ax[r_, 0].pcolormesh(lone, late, G1, cmap='viridis', vmin=0, vmax=np.nanmax([G1.max(), 1e-9]))
        fig.colorbar(m1, ax=ax[r_, 0], label='disturbed area 2100 (km$^2$)')
        lim = np.nanmax(np.abs(G2)) or 1
        m2 = ax[r_, 1].pcolormesh(lone, late, G2, cmap='RdBu_r', vmin=-lim, vmax=lim)
        fig.colorbar(m2, ax=ax[r_, 1], label='RTS net C loss 2016-2100 (Gg C, + = loss)')
        tot1 = np.nansum(G1); tot2 = np.nansum(G2)
        ax[r_, 0].set_title(f'Case {case}: disturbed area, total {tot1:.1f} km$^2$', fontsize=10)
        ax[r_, 1].set_title(f'Case {case}: net C loss, total {tot2 / 1e3:.3f} Tg C', fontsize=10)
        for a in ax[r_]: a.set_xlabel('lon'); a.set_ylabel('lat')
    fig.suptitle('Gydan 10x10, SSP5-8.5, T1 ensemble mean (50 members), folr 0.8', fontsize=11)
    fig.savefig(out, dpi=300); plt.close(fig)


def band(ax, x, Y, c, label):
    lo, md, hi = np.percentile(Y, [5, 50, 95], axis=0)
    ax.fill_between(x, lo, hi, color=c, alpha=0.25, lw=0); ax.plot(x, md, color=c, lw=1.5, label=label)


def timeseries(D, out):
    keys = [('area_frac', 'disturbed area (km$^2$)', 1), ('dTKLATC', 'TKLATC (Gg C yr$^{-1}$)', 1e9),
            ('dTKNETC', 'TKNETC anomaly (Gg C yr$^{-1}$, + = loss)', 1e9), ('dNEE', 'NEE anomaly (Gg C yr$^{-1}$)', 1e9)]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True, sharex=True)
    for a, (k, lab, sc) in zip(ax.ravel(), keys):
        for case in ('A', 'B'):
            band(a, YRS, D[f'{case}_T1_f08'][k] / sc, COL[case], f'Case {case} T1 (median, 5-95 %)')
            a.plot(YRS, D[f'{case}_T2_f08'][k][0] / sc, ':', color=COL[case], lw=1.2, label=f'Case {case} T2')
        lo, hi = D['A_det_mflow_f08'][k][0] / sc, D['A_det_mfhigh_f08'][k][0] / sc
        a.plot(YRS, lo, '--', color=COL['A'], lw=0.8, label='Case A mf 1.3e-4 / 2e-3')
        a.plot(YRS, hi, '--', color=COL['A'], lw=0.8)
        if k == 'area_frac': a.set_yscale('log'); a.set_ylim(1, None)
        a.set_ylabel(lab, fontsize=9); a.axhline(0, color='k', lw=0.4)
    ax[0, 0].legend(fontsize=7)
    for a in ax[1]: a.set_xlabel('year')
    fig.suptitle('Gydan domain totals, SSP5-8.5, folr 0.8 (e2000 spin-up)', fontsize=11)
    fig.savefig(out, dpi=300); plt.close(fig)


def cumulative(D, out):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True, sharey=True)
    for a, case in zip(ax, ('A', 'B')):
        band(a, YRS, np.cumsum(D[f'{case}_T1_f08']['dTKNETC'], -1) / 1e12, COL[case], 'T1 folr 0.8 (main)')
        sens = [(f'{case}_T1_f05', 'T1 folr 0.5', '-.'), (f'{case}_T2_f08', 'T2 folr 0.8', ':')]
        if case == 'A':
            sens += [('A_T1_gt050_f08', 'T1 gt050 (s := gt050)', '-'), ('A_det_mflow_f08', 'mf 1.3e-4', '--'),
                     ('A_det_mfhigh_f08', 'mf 2e-3', '--')]
        else:
            sens += [('B_legacy_f08', 'legacy 2021 cohort only', '--')]
        for n, lab, ls in sens:
            y = np.cumsum(D[n]['dTKNETC'], -1) / 1e12
            a.plot(YRS, np.median(y, 0), ls, lw=1, color='k' if 'gt050' not in n else 'tab:green', label=lab)
        a.set_title(f'Case {case}: cumulative RTS net C loss', fontsize=10); a.set_xlabel('year'); a.legend(fontsize=7)
    ax[0].set_ylabel('Tg C (+ = loss)')
    fig.savefig(out, dpi=300); plt.close(fig)


def fmt(s, k, scale=1.0, f='{:.3g}'):
    x = s[k]
    return (f + ' [' + f + '–' + f + ']').format(x['p50'] * scale, x['p05'] * scale, x['p95'] * scale)


def md(S, out):
    sc = S['scenarios']; b = S['baseline']['f08']
    L = ['# Phase 3: offline cohort blend, Case A vs Case B (Gydan 10x10, SSP5-8.5 MRI-ESM2-0, e2000)', '',
         'Domain totals, T1 = 50-member Poisson ensemble, median [5–95 %]; Tg C, + = loss to the atmosphere/laterally.', '',
         '| scenario | area 2100 (km²) | events | TKLATC cum (Tg) | RTS net C loss (ΔTKNETC cum, Tg) | ΔNEE cum (Tg) | ΔECOC 2100 (Tg) |',
         '|---|---|---|---|---|---|---|']
    for n in sc:
        s = sc[n]
        L.append(f"| {n} | {fmt(s, 'area_km2_2100')} | {fmt(s, 'n_events', f='{:.0f}')} | {fmt(s, 'TKLATC_cum_Tg')} | "
                 f"{fmt(s, 'dTKNETC_cum_Tg')} | {fmt(s, 'dNEE_cum_Tg')} | {fmt(s, 'dECOC_2100_Tg')} |")
    L += ['', f"Baseline (no RTS, folr 0.8): ECOC 2016 {b['ECOC_2016_Tg']:.1f} Tg → 2100 {b['ECOC_2100_Tg']:.1f} Tg "
              f"(Δ {b['ECOC_2100_Tg'] - b['ECOC_2016_Tg']:+.1f} Tg; cumulative TKNETC {b['TKNETC_cum_Tg']:+.1f} Tg).", '',
          '**Caveat:** absolute losses depend on the spin-up length: e1000 vs e2000 shifts the absolute baseline and RTS '
          'losses by ~40 % (c86/c04 sensitivity: start ECOC +44 %, RTS ECOC effect ~+20 %, TKLATC cum +25 %, '
          'baseline loss at c04 +40 %); neither spin-up is equilibrated (deep/mineral drift ~35 gC/m²/yr).', '',
          f"Positive single-cohort ECOC effect (y2021 − base at Dec 2100 > 0, folr 0.8): {S['positive_effect']['n']} of "
          f"{S['cells']} cells; {S['positive_effect']['n_alt']} of them with baseline ALT(2091–2100) ≥ 3 m "
          f"({S['positive_effect']['cells_base_alt_ge3']} cells have baseline ALT ≥ 3 m). Cells: "
          f"{', '.join(S['positive_effect']['cells'])}.", '',
          f"Closure over cohorts: max monthly |dECOC+TKNETC| {S['closure']['c_month_max']:.2e} gC/m², max cumulative "
          f"{S['closure']['c_cum_absmax']:.3f} gC/m² (tolerances 0.01 / 1.0). T2 uniform-fallback windows: "
          f"A {S['t2_fallback_windows']['A']}, B {S['t2_fallback_windows']['B']}. Compute {S['seconds_total']:.0f} s.", '',
          'Refs: Jia et al. 2026 doi:10.1038/s41597-026-08284-y; Nesterova et al. 2025 doi:10.1594/PANGAEA.974406; '
          'ARTS v1.0 doi:10.18739/A2PK0738B (1.44 ha).', '']
    open(out, 'w').write('\n'.join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('blend_dir'); ap.add_argument('table'); a = ap.parse_args()
    D, P, S = load(a.blend_dir); table = list(csv.DictReader(open(a.table)))
    maps(D, P, table, os.path.join(a.blend_dir, 'fig_p3_maps.png'))
    timeseries(D, os.path.join(a.blend_dir, 'fig_p3_timeseries.png'))
    cumulative(D, os.path.join(a.blend_dir, 'fig_p3_cumulative.png'))
    md(S, os.path.join(a.blend_dir, 'summary.md'))


if __name__ == '__main__':
    main()
