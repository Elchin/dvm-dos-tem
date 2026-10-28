#!/usr/bin/env python3
"""Figures + summary markdown for the phase-3 blend (blend_phase3.py output dir).

fig_p3_maps.png        2100 cumulative RTS-disturbed area and RTS net C loss (cumulative TKNETC anomaly
                       2016-2100) per cell, Case A vs Case B (T1 ensemble mean, folr 0.8)
fig_p3_timeseries.png  domain disturbed area, TKLATC, TKNETC anomaly and NEE anomaly; T1 ensemble median
                       and 5-95 % band, Case A mf band (1.3e-4 / 2e-3), T2 dotted
fig_p3_cumulative.png  cumulative RTS net C loss with the sensitivities (gt050, folr 0.5, T2, mf band)
summary.md             domain tables, peak/reversal, method notes; --compare DIR adds a v1-vs-current table
Headline numbers are at 2100; the cumulative-loss peak (year of max of the median) is annotated.
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
        m1 = ax[r_, 0].pcolormesh(lone, late, G1, cmap='viridis', vmin=0, vmax=max(float(np.nanmax(G1)), 1e-9))
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


def peak(Y):
    """(year, value) of the maximum of the ensemble-median cumulative series."""
    m = np.median(np.cumsum(Y, -1), 0); i = int(np.argmax(m)); return int(YRS[i]), float(m[i]), float(m[-1])


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
        if k in ('dTKNETC', 'dNEE'):
            yp = peak(D['B_T1_f08']['dTKNETC'])[0]
            a.axvline(yp, color='grey', lw=0.6, ls='--')
            yt, va = (0.97, 'top') if k == 'dTKNETC' else (0.03, 'bottom')
            a.text(yp - 1, yt, f'{yp}: cum. loss peak;\nlater the unslumped baseline\nloses more (deep thaw, N cells)',
                   transform=a.get_xaxis_transform(), ha='right', va=va, fontsize=6.5, color='0.3')
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
        yp, vp, ve = peak(D[f'{case}_T1_f08']['dTKNETC'] / 1e12)
        a.plot([yp], [vp], 'v', color=COL[case], ms=6)
        a.annotate(f'peak {vp:.2f} Tg ({yp})\n2100: {ve:.2f} Tg (headline)\nreversal: baseline deep thaw\nin N cells exceeds slumped',
                   (yp, vp), xytext=(0.03, 0.97), textcoords='axes fraction', fontsize=7, va='top',
                   arrowprops=dict(arrowstyle='->', color='0.4', lw=0.6))
        a.set_title(f'Case {case}: cumulative RTS net C loss', fontsize=10); a.set_xlabel('year'); a.legend(fontsize=7, loc='lower right')
    ax[0].set_ylabel('Tg C (+ = loss)')
    fig.savefig(out, dpi=300); plt.close(fig)


def fmt(s, k, scale=1.0, f='{:.3g}'):
    x = s[k]
    return (f + ' [' + f + '–' + f + ']').format(x['p50'] * scale, x['p05'] * scale, x['p95'] * scale)


def md(S, out, S1=None):
    sc = S['scenarios']; b = S['baseline']['f08']
    pk = lambda n: sc[n].get('dTKNETC_cum_peak') or {}
    L = ['# Phase 3: offline cohort blend, Case A vs Case B (Gydan 10x10, SSP5-8.5 MRI-ESM2-0, e2000)', '',
         'Domain totals, T1 = 50-member Poisson ensemble, median [5–95 %]; Tg C, + = loss to the atmosphere/laterally.', '',
         '| scenario | area 2100 (km²) | events | TKLATC cum (Tg) | RTS net C loss (ΔTKNETC cum, Tg) | ΔNEE cum (Tg) | ΔECOC 2100 (Tg) |',
         '|---|---|---|---|---|---|---|']
    for n in sc:
        s = sc[n]
        L.append(f"| {n} | {fmt(s, 'area_km2_2100')} | {fmt(s, 'n_events', f='{:.0f}')} | {fmt(s, 'TKLATC_cum_Tg')} | "
                 f"{fmt(s, 'dTKNETC_cum_Tg')} | {fmt(s, 'dNEE_cum_Tg')} | {fmt(s, 'dECOC_2100_Tg')} |")
    pa, pb = pk('A_T1_f08'), pk('B_T1_f08')
    if pa and pb:
        L += ['', f"**Headline = 2100.** The cumulative RTS net C loss (T1 median) peaks in {pa['year']} at "
                  f"{pa['Tg']:.3f} Tg (Case A) and in {pb['year']} at {pb['Tg']:.3f} Tg (Case B), then falls by "
                  f"{-pa['from_peak_to_2100_Tg']:.3f} / {-pb['from_peak_to_2100_Tg']:.3f} Tg to 2100. Late-century "
                  "reversal: after ~2094 the unslumped baseline in the northern cells (e.g. c66, c76, c86) starts losing "
                  "0.7–2 kg C m⁻² yr⁻¹ to deep thaw (ALD deepening), while slumped patches (shallower ALD) do not, so the "
                  "RTS anomaly of TKNETC/NEE turns negative. This is in the runs themselves (also in the legacy-only "
                  "cohort), not a blend artefact."]
    if S1 is not None:
        L += ['', 'Change from v1 (v1: legacy cohort with its initiation pulse; T1 median, Tg C):', '',
              '| scenario | TKLATC v1 → now | net C loss v1 → now | ΔNEE v1 → now | ΔECOC 2100 v1 → now |', '|---|---|---|---|---|']
        for n in sc:
            if n not in S1['scenarios']: continue
            o = S1['scenarios'][n]; c = sc[n]
            L.append(f"| {n} | " + ' | '.join(f"{o[k]['p50']:.3g} → {c[k]['p50']:.3g}" for k in
                     ('TKLATC_cum_Tg', 'dTKNETC_cum_Tg', 'dNEE_cum_Tg', 'dECOC_2100_Tg')) + ' |')
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
          '**Method notes and assumptions** (confirmed by E. Jafarov, 2026-10-28):', '',
          '- *Legacy cohort (Case B):* pre-2016 scars (Nesterova count × 1.44 ha) keep their area from 2021 but '
          'their initiation pulse is dropped: the cohort uses the Aug-2021 response from month+1 on (fluxes zero in '
          'and before the initiation month; ECOC anomaly rebased to zero at the initiation month), i.e. only '
          'post-initiation exposure/erosion and regrowth. Closure dECOC = −TKNETC is exact; the legacy ECOC '
          'anomaly is the change since 2021. v1 (with pulse) is reproducible with --legacy-pulse.',
          '- *Main Case B assumption:* base initiation rate = one Nesterova inventory\'s worth of area per 20-yr '
          'window × Jia SSP5 gt050 ratio r_k (3x3 smoothed, relative to 2030).',
          '- Jia 2030/2050/2070/2090 maps apply to the windows 2021–2040 / 2041–2060 / 2061–2080 / 2081–2100.',
          '- Fixed 1.44 ha per event from initiation (no in-model area growth used); standing area only accumulates.',
          '- *Limitation:* initiation years between run start years reuse the lower start-year run (2021, 2031, …, '
          '2091), time-shifted, so the response carries the climate of up to 9 yr earlier.', '',
          'Refs: Jia et al. 2026 doi:10.1038/s41597-026-08284-y; Nesterova et al. 2025 doi:10.1594/PANGAEA.974406; '
          'ARTS v1.0 doi:10.18739/A2PK0738B (1.44 ha).', '']
    open(out, 'w').write('\n'.join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('blend_dir'); ap.add_argument('table')
    ap.add_argument('--compare', default=None, help='previous blend dir (adds a v1-vs-current table)')
    a = ap.parse_args()
    D, P, S = load(a.blend_dir); table = list(csv.DictReader(open(a.table)))
    maps(D, P, table, os.path.join(a.blend_dir, 'fig_p3_maps.png'))
    timeseries(D, os.path.join(a.blend_dir, 'fig_p3_timeseries.png'))
    cumulative(D, os.path.join(a.blend_dir, 'fig_p3_cumulative.png'))
    S1 = json.load(open(os.path.join(a.compare, 'summary.json'))) if a.compare else None
    md(S, os.path.join(a.blend_dir, 'summary.md'), S1)


if __name__ == '__main__':
    main()
