#!/usr/bin/env python3
"""Work-notes figures 4/5 + table: Gydan (G, CMT70) A3b_snow300_otc thermal set, moss 1/3/10 cm,
legacy (L) vs moss_water_scheme 1 (W) vs moss_water_scheme 1 + moss_tc_scheme 2 (N).
Reads moss_daily.csv / energy_daily.csv (energy_diag 1) and the NetCDF outputs.
  plot_moss_water.py OUT_seasonal.png OUT_ald.png OUT.json"""
import sys, os, json, subprocess, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = os.path.expanduser('~/tem-dev/moss-work'); H = os.path.join(W, 'harness')
def A(r): return json.loads(subprocess.check_output([sys.executable, os.path.join(H, 'analyze.py'), os.path.join(W, 'runs', r)]))
def csv(p):
    with open(p) as f: hdr = f.readline().strip().split(',')
    d = np.genfromtxt(p, delimiter=',', skip_header=1)
    return {h: d[:, i] for i, h in enumerate(hdr)}
def keff(ku, kf, liq, ice):
    w = liq + ice; f = np.where(w > 0, ice / np.where(w > 0, w, 1), 0.0)
    return np.where(w > 0, ku ** (1 - f) * kf ** f, ku)
geo = [('d01', 1), ('d03', 3), ('d10', 10)]; sch = [('L', 'legacy (S21H on moss, fibric copy)'), ('W', 'moss water 1, S21H k'), ('N', 'moss water 1 + Johansen k (tc 2)')]
col = {'L': 'tab:gray', 'W': 'tab:blue', 'N': 'tab:red'}; lsd = {'d01': ':', 'd03': '--', 'd10': '-'}
tab = []; mon = {}
for g, cm in geo:
    for s, lab in sch:
        r = 'A_%s_%s' % (s, g); D = os.path.join(W, 'runs', r)
        m = csv(os.path.join(D, 'stage2a_hist', 'moss_daily.csv'))
        yr = 1950 + (m['day'] // 365).astype(int); doy = m['doy'].astype(int)
        sel = (yr >= 1986) & (yr <= 2015)
        dz = m['dz']; thl = m['liq'] / (1000 * dz); thi = m['ice'] / (917 * dz)
        ke = keff(m['k_unf'], m['k_frz'], m['liq'], m['ice'])
        mo = np.minimum((doy * 12) // 365, 11)
        summ = sel & (mo >= 6) & (mo <= 7); wint = sel & (mo <= 1)
        mon[r] = dict(thl=[float(thl[sel & (mo == k)].mean()) for k in range(12)],
                      thi=[float(thi[sel & (mo == k)].mean()) for k in range(12)],
                      k=[float(ke[sel & (mo == k)].mean()) for k in range(12)])
        res = [float(np.nanmax(np.abs(csv(os.path.join(D, st, 'moss_daily.csv'))['resid']))) for st in ('stage1', 'stage2a_hist', 'stage2b_rts')]
        en = []
        for st in ('stage1', 'stage2a_hist', 'stage2b_rts'):
            e = csv(os.path.join(D, st, 'energy_daily.csv'))
            en.append((float(np.nanmax(np.abs(e['resday']))), float(np.nanmax(np.abs(e['unattr'])))))
        a = A(r); y = np.arange(a['years'][0], a['years'][1] + 1); yl = list(y); ald = np.array(a['ALD'])
        i0, i1 = yl.index(1986), yl.index(2015) + 1
        tab.append(dict(case=r, scheme=s, moss_cm=cm,
                        thl_summer=round(float(thl[summ].mean()), 3), thi_summer=round(float(thi[summ].mean()), 3),
                        thl_winter=round(float(thl[wint].mean()), 3), thi_winter=round(float(thi[wint].mean()), 3),
                        k_summer=round(float(ke[summ].mean()), 3), k_winter=round(float(ke[wint].mean()), 3),
                        T_moss_summer=round(float(m['tem'][summ].mean()), 2), T_moss_winter=round(float(m['tem'][wint].mean()), 2),
                        ald_1986_2015=round(float(ald[i0:i1].mean()), 3),
                        ald_below_moss=round(float(ald[i0:i1].mean()) - cm / 100.0, 3),
                        ald_2022_2026=round(float(ald[yl.index(2022):yl.index(2026) + 1].mean()), 3),
                        moss_water_resid_max=max(res), energy_resday_max_Jm2=max(x[0] for x in en),
                        energy_unattr_max_Jm2=max(x[1] for x in en), cbal_max_abs=a['cbal']['max_abs_resid']))
        mon[r]['ald'] = (y.tolist(), ald.tolist())
# Fig 1: seasonal cycle
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4)); x = np.arange(1, 13)
for g, cm in geo:
    for s, lab in sch:
        r = 'A_%s_%s' % (s, g)
        ax[0].plot(x, mon[r]['thl'], color=col[s], ls=lsd[g], label='%s %d cm' % (s, cm))
        ax[1].plot(x, mon[r]['thi'], color=col[s], ls=lsd[g])
        ax[2].plot(x, mon[r]['k'], color=col[s], ls=lsd[g])
ax[0].set_ylabel('moss theta_liq (m3/m3)'); ax[1].set_ylabel('moss theta_ice (m3/m3)')
ax[2].set_ylabel('moss k as used by the solver (W/m/K)\nlam_u^(1-f) lam_f^f'); ax[2].set_yscale('log')
for a_ in ax: a_.set_xlabel('month'); a_.set_xticks(x); a_.grid(alpha=0.3)
ax[0].legend(fontsize=7, ncol=3)
fig.suptitle('Gydan c55 (CMT70) A3b_snow300_otc, moss layer monthly climatology 1986-2015; gray = legacy, blue = moss_water_scheme 1, '
             'red = + moss_tc_scheme 2 (ks 0.25); dotted/dashed/solid = 1/3/10 cm', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
# Fig 2: ALD
fig, ax = plt.subplots(1, 2, figsize=(14, 4.4))
for g, cm in geo:
    for s, lab in sch:
        r = 'A_%s_%s' % (s, g); yy, aa = mon[r]['ald']
        ax[0].plot(yy, aa, color=col[s], ls=lsd[g], label='%s %d cm' % (s, cm))
ax[0].set_xlim(1980, 2045); ax[0].axvline(2021.6, color='k', lw=0.6, ls=':'); ax[0].set_ylabel('ALD (m, from moss top)')
ax[0].legend(fontsize=7, ncol=3); ax[0].grid(alpha=0.3)
xb = np.arange(3)
for k, (s, lab) in enumerate(sch):
    v = [t['ald_1986_2015'] for g, cm in geo for t in tab if t['case'] == 'A_%s_%s' % (s, g)]
    ax[1].bar(xb + 0.27 * k - 0.27, v, 0.27, color=col[s], label=lab)
ax[1].set_xticks(xb); ax[1].set_xticklabels(['1 cm', '3 cm', '10 cm']); ax[1].set_ylabel('mean ALD 1986-2015 (m)')
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, axis='y')
fig.suptitle('Gydan c55 (CMT70) A3b_snow300_otc, moss_dyn 1 static thickness, short spin-up test settings; RTS Aug 2021', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[2], dpi=110)
json.dump(dict(table=tab, monthly={k: {kk: vv for kk, vv in v.items() if kk != 'ald'} for k, v in mon.items()}),
          open(sys.argv[3], 'w'), indent=1)
for t in tab: print(t)
