#!/usr/bin/env python3
"""Work-notes figure 1: moss_dyn test runs at Toolik (T, CMT05) and Gydan (G, CMT70), RTS removal in Aug 2021.
  plot_tests.py OUT.png"""
import sys, os, json, subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = os.path.expanduser('~/tem-dev/moss-work'); H = os.path.join(W, 'harness')
def A(r):
    return json.loads(subprocess.check_output([sys.executable, os.path.join(H, 'analyze.py'), os.path.join(W, 'runs', r), '--moss-pfts', '5,6,7']))
runs = {r: A(r) for r in ('T0', 'T1', 'T2', 'G0', 'G1')}
def yrs(d): return np.arange(d['years'][0], d['years'][1] + 1)
lab = {'T0': 'T moss_dyn 0 (legacy 0.01 m)', 'T1': 'T moss_dyn 1 (dmax 0.30 = cmt maxdmossthick, t_half 5)',
       'T2': 'T moss_dyn 2 (moss PFT C / 11.25 kgC m$^{-3}$)', 'G0': 'G moss_dyn 0', 'G1': 'G moss_dyn 1 (dmax 0.035, t_half 5)'}
sty = {'T0': ('k', '-'), 'T1': ('tab:orange', '-'), 'T2': ('tab:green', '-'), 'G0': ('k', '--'), 'G1': ('tab:blue', '--')}
fig, ax = plt.subplots(2, 2, figsize=(12, 8)); ax = ax.ravel()
for r, d in runs.items():
    y = yrs(d); c, ls = sty[r]
    ax[0].plot(y, np.array(d['MOSSDZ']) * 100, color=c, ls=ls, label=lab[r])
    ax[2].plot(y, d['ALD'], color=c, ls=ls, label=r)
for r in ('T0', 'T1', 'T2'):
    d = runs[r]; ax[1].plot(yrs(d), d['moss_c_dec'], color=sty[r][0], label=r)
ax[0].set_ylabel('live moss thickness MOSSDZ (cm)'); ax[0].set_yscale('log')
ax[1].set_ylabel('moss PFT C, Dec (gC m$^{-2}$), CMT05 PFTs 5-7')
ax[2].set_ylabel('ALD (m, from the moss top)')
# C balance residuals (monthly)
for r in ('T0', 'T2', 'G1'):
    d = runs[r]
    ax[3].bar(r, d['cbal']['max_abs_resid'], color=sty[r][0])
ax[3].set_ylabel('max |monthly C residual| (gC m$^{-2}$)\nd(VEGC+SOC+DEADC+DWDC)+TKNETC, 1901/1950-2045')
for a in ax[:3]:
    a.set_xlim(1990, 2045); a.axvline(2021.6, color='r', lw=0.8, ls=':'); a.grid(alpha=0.3)
ax[0].legend(fontsize=7, loc='lower left'); ax[1].legend(fontsize=8); ax[2].legend(fontsize=8, ncol=2)
fig.suptitle('moss_dyn test runs (short spin-up p10/e300/s30, v4b settings; RTS slump Aug 2021, red dotted)', fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110); print('wrote', sys.argv[1])
