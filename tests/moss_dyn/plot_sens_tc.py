#!/usr/bin/env python3
"""Work-notes figure 3 + table: Gydan (G, CMT70) moss thickness sensitivity with moss_tc_scheme 1
(Chadburn et al. 2015 moss conductivity) vs moss_tc_scheme 0, under legacy thermal (L) and org_tc_scheme 1 S21H (O).
  plot_sens_tc.py OUT.png OUT.json"""
import sys, os, json, subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = os.path.expanduser('~/tem-dev/moss-work'); H = os.path.join(W, 'harness')
def A(r): return json.loads(subprocess.check_output([sys.executable, os.path.join(H, 'analyze.py'), os.path.join(W, 'runs', r)]))
geo = [('d01', '1 cm'), ('d03', '3 cm'), ('d10', '10 cm'), ('rg05', 'regrow 5 cm')]
tc0 = {('L', 'd01'): 'G0', ('L', 'd03'): 'S_L_d03', ('L', 'd10'): 'S_L_d10', ('L', 'rg05'): 'S_L_rg05',
       ('O', 'd01'): 'S_O_ctl', ('O', 'd03'): 'S_O_d03', ('O', 'd10'): 'S_O_d10', ('O', 'rg05'): 'S_O_rg05'}
col = {'d01': 'k', 'd03': 'tab:orange', 'd10': 'tab:purple', 'rg05': 'tab:blue'}
tab = []; fig, ax = plt.subplots(1, 4, figsize=(19, 4.6)); ts = {}
for th in 'LO':
    for g, lab in geo:
        for sch, r in ((0, tc0[(th, g)]), (1, 'C_%s_%s' % (th, g))):
            d = A(r); y = np.arange(d['years'][0], d['years'][1] + 1); yl = list(y)
            ald = np.array(d['ALD']); mz = np.array(d['MOSSDZ']); below = ald - mz
            i0, i1 = yl.index(1986), yl.index(2015) + 1; j0, j1 = yl.index(2022), yl.index(2026) + 1; k0 = yl.index(2036)
            tab.append(dict(case=r, thermal=th, moss=lab, moss_tc_scheme=sch, moss_cm_2015=round(100 * mz[i1 - 1], 2),
                            ald_1986_2015=round(float(ald[i0:i1].mean()), 3),
                            ald_below_moss_1986_2015=round(float(below[i0:i1].mean()), 3),
                            ald_2022_2026=round(float(ald[j0:j1].mean()), 3), ald_2036_2045=round(float(ald[k0:].mean()), 3),
                            cbal_max_abs=d['cbal']['max_abs_resid']))
            if sch == 1:
                ls = '-' if th == 'L' else '--'
                ax[0].plot(y, 100 * mz, color=col[g], ls=ls)
                ax[1].plot(y, ald, color=col[g], ls=ls, label='%s %s' % (th, lab))
x = np.arange(4); lbl = {('L', 0): 'legacy, moss_tc 0', ('L', 1): 'legacy, Chadburn',
                         ('O', 0): 'S21H, moss_tc 0', ('O', 1): 'S21H, Chadburn'}
for p, key in ((2, 'ald_1986_2015'), (3, 'ald_below_moss_1986_2015')):
    for k, (th, sch) in enumerate([('L', 0), ('L', 1), ('O', 0), ('O', 1)]):
        v = [t[key] for t in tab if t['thermal'] == th and t['moss_tc_scheme'] == sch]
        ax[p].bar(x + 0.2 * k - 0.3, v, 0.2, label=lbl[(th, sch)], hatch='//' if sch else None)
    ax[p].set_xticks(x); ax[p].set_xticklabels([l for _, l in geo], fontsize=8); ax[p].set_yscale('log')
    ax[p].legend(fontsize=7); ax[p].grid(alpha=0.3, axis='y')
ax[2].set_ylabel('mean ALD 1986-2015 (m, from moss top)'); ax[3].set_ylabel('mean ALD below moss 1986-2015 (m)')
ax[0].set_ylabel('MOSSDZ (cm)'); ax[1].set_ylabel('ALD (m, from moss top), moss_tc_scheme 1'); ax[1].legend(fontsize=6, ncol=2)
for a in ax[:2]: a.set_xlim(1980, 2045); a.axvline(2021.6, color='r', lw=0.8, ls=':'); a.grid(alpha=0.3)
fig.suptitle('Gydan c55 (CMT70), moss_dyn 1, short spin-up test settings; moss_tc_scheme 1 = Chadburn et al. 2015 (0.06-0.5 W/m/K linear); '
             'solid = legacy thermal, dashed = org_tc_scheme 1 S21H; RTS Aug 2021', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
json.dump(tab, open(sys.argv[2], 'w'), indent=1); print(json.dumps(tab, indent=0))
