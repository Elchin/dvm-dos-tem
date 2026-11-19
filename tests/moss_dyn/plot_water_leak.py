#!/usr/bin/env python3
"""Task 7: whole-column water budget (water_daily.csv, energy_diag 1) for the Gydan T1 diagnostic runs.
usage: plot_water_leak.py RUNS_DIR OUT.png   (RUNS_DIR holds D_off, D_on, D_fix, each with stage2a_hist/water_daily.csv)
Daily residual r = dS - (rain + snow - eet - qover - qdrain), S = swe + moss + soil (liq+ice) + puddle + canopy water/snow.
r > 0 = water created (outflows over-reported), r < 0 = destroyed."""
import sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def rd(p):
    h = open(p).readline().strip().split(','); d = np.genfromtxt(p, delimiter=',', skip_header=1)
    return {k: d[:, i] for i, k in enumerate(h)}
def resid(D):
    S = D['swe'] + D['W_moss_liq'] + D['W_moss_ice'] + D['W_soil_liq'] + D['W_soil_ice'] + D['puddle'] + D['veg_rwater'] + D['veg_snow']
    F = D['rain'] + D['snow'] - D['eet'] - D['qover'] - D['qdrain']
    return np.diff(S) - F[1:], {k: v[1:] for k, v in D.items()}
R, out = sys.argv[1], sys.argv[2]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
lab = {'D_off': 'legacy (moss schemes off)', 'D_on': 'moss_water_scheme 1 + tc 2, before fix', 'D_fix': 'moss_water_scheme 1 + tc 2, fixed'}
col = {'D_off': 'k', 'D_on': 'tab:red', 'D_fix': 'tab:blue'}
for n in ['D_off', 'D_on', 'D_fix']:
    r, X = resid(rd('%s/%s/stage2a_hist/water_daily.csv' % (R, n))); ny = len(r) / 365.
    yr = 1950 + np.arange(len(r)) / 365.
    ax[0].plot(yr, np.cumsum(r), col[n], lw=1.3, label='%s: %+.2f mm/yr' % (lab[n], r.sum() / ny))
    if 'rich_err' in X and n == 'D_on':
        e = X['rich_err']; b2 = r - e
        ax[1].bar([0], [e.sum() / ny], color='tab:orange', label='R1: Richards rejected-substep dt (final - target column liquid)')
        ax[1].bar([1], [b2.sum() / ny], color='tab:purple', label='B2: early-return infil lost (partly frozen column, moss drain)')
        ax[1].bar([2], [r.sum() / ny], color='tab:red', label='total before fix')
    if n == 'D_fix': ax[1].bar([3], [r.sum() / ny], color='tab:blue', label='total after fix')
ax[0].axhline(0, color='0.6', lw=0.6); ax[0].set_xlabel('year (stage2a, 66 yr)'); ax[0].set_ylabel('cumulative column water residual (mm)')
ax[0].set_title('Gydan T1 (rmf 0.05, kdc x0.5): daily whole-column water budget'); ax[0].legend(fontsize=8)
ax[1].axhline(0, color='0.6', lw=0.6); ax[1].set_xticks(range(4)); ax[1].set_xticklabels(['R1', 'B2', 'before', 'after'])
ax[1].set_ylabel('mm/yr (+ = water created)'); ax[1].set_title('scheme-1 leak decomposition'); ax[1].legend(fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(out, dpi=110); print('wrote', out)
