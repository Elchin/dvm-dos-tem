#!/usr/bin/env python3
"""Work-notes figure 2 + table: Gydan (G, CMT70) moss thickness sensitivity, legacy thermal (L) vs org_tc_scheme 1 S21H (O).
  plot_sens.py OUT.png OUT.json"""
import sys, os, json, subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = os.path.expanduser('~/tem-dev/moss-work'); H = os.path.join(W, 'harness')
def A(r): return json.loads(subprocess.check_output([sys.executable, os.path.join(H, 'analyze.py'), os.path.join(W, 'runs', r)]))
cases = [('G0', 'L', 'static 1 cm (legacy)'), ('S_L_d001', 'L', 'static 0.1 cm'), ('S_L_d03', 'L', 'static 3 cm'),
         ('S_L_d10', 'L', 'static 10 cm'), ('S_L_rg05', 'L', 'regrow dmax 5 cm, t_half 5'),
         ('S_O_ctl', 'O', 'static 1 cm (legacy)'), ('S_O_d001', 'O', 'static 0.1 cm'), ('S_O_d03', 'O', 'static 3 cm'),
         ('S_O_d10', 'O', 'static 10 cm'), ('S_O_rg05', 'O', 'regrow dmax 5 cm, t_half 5')]
col = {'static 0.1 cm': 'tab:red', 'static 1 cm (legacy)': 'k', 'static 3 cm': 'tab:orange', 'static 10 cm': 'tab:purple',
       'regrow dmax 5 cm, t_half 5': 'tab:blue'}
tab = []; fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for r, th, lab in cases:
    d = A(r); y = np.arange(d['years'][0], d['years'][1] + 1)
    ald = np.array(d['ALD']); mz = np.array(d['MOSSDZ']); below = ald - mz
    i0, i1 = list(y).index(1986), list(y).index(2015) + 1
    j0, j1 = list(y).index(2022), list(y).index(2026) + 1
    k0 = list(y).index(2036)
    tab.append(dict(case=r, thermal=th, moss=lab, moss_cm_2015=round(100 * mz[i1 - 1], 2),
                    ald_1986_2015=round(float(ald[i0:i1].mean()), 3), ald_below_moss_1986_2015=round(float(below[i0:i1].mean()), 3),
                    ald_2022_2026=round(float(ald[j0:j1].mean()), 3), ald_2036_2045=round(float(ald[k0:].mean()), 3),
                    cbal_max_abs=d['cbal']['max_abs_resid']))
    ls = '-' if th == 'L' else '--'
    ax[0].plot(y, 100 * mz, color=col[lab], ls=ls)
    ax[1].plot(y, ald, color=col[lab], ls=ls, label='%s %s' % (th, lab))
x = np.arange(5)
for k, th in enumerate('LO'):
    v = [t['ald_1986_2015'] for t in tab if t['thermal'] == th]
    ax[2].bar(x + 0.4 * k - 0.2, v, 0.4, label={'L': 'legacy ln-Kersten', 'O': 'org_tc_scheme 1 (S21H)'}[th])
ax[2].set_xticks(x); ax[2].set_xticklabels(['1 cm\n(legacy)', '0.1 cm', '3 cm', '10 cm', 'regrow\n5 cm'], fontsize=8)
ax[2].set_ylabel('mean ALD 1986-2015 (m, from moss top)'); ax[2].legend(fontsize=8)
ax[0].set_ylabel('MOSSDZ (cm)'); ax[1].set_ylabel('ALD (m, from moss top)'); ax[1].legend(fontsize=6, ncol=2)
for a in ax[:2]: a.set_xlim(1980, 2045); a.axvline(2021.6, color='r', lw=0.8, ls=':'); a.grid(alpha=0.3)
fig.suptitle('Gydan c55 (CMT70), moss_dyn 1, short spin-up test settings (not the calibrated A0/A2 set); solid = legacy thermal, dashed = org_tc_scheme 1 S21H; RTS Aug 2021', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
json.dump(tab, open(sys.argv[2], 'w'), indent=1); print(json.dumps(tab, indent=0))
