#!/usr/bin/env python3
"""RTS port figures + slump metrics from a port_check.py work dir (no model runs).

  python3 tests/rts/port_figs.py WORK OUTDIR

Writes OUTDIR/{rts_effect_vs_v03.png, rts_closure.png, rts_regrow.png, rts_port_metrics.json}.
Effect = run(thermokarst 1) - run(off_auto: thermokarst 0 + the three auto-enabled switches);
v0.3 effect = tests/rts/reference/v03/rts - reference/baseline (rts-regrow 7c9fd634 binary, m2d on).
"""
import json, os, sys
import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); REF = os.path.join(HERE, 'reference')
W, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
INIT = 127   # tr month index of the initiation (year 10, Aug)

def rd(d, v, ts='monthly'):
    f = os.path.join(d, '%s_%s_tr.nc' % (v, ts))
    if not os.path.exists(f): return None
    with nc.Dataset(f) as x: a = np.array(x[v][:], dtype=float)
    a = np.where(np.abs(a) > 1e20, np.nan, a); a = np.where(a <= -9999, np.nan, a)
    a = a.reshape(a.shape[0], -1)
    return np.nansum(a, axis=1) if a.shape[1] > 1 else a[:, 0]

O = lambda n: os.path.join(W, n, 'output')
on, off = O('on'), O('off_auto'); von, voff = os.path.join(REF, 'v03', 'rts'), os.path.join(REF, 'baseline')
M = {}
for v in ('SOC', 'VEGC', 'SHLWC', 'DEEPC', 'ORGN'):
    a, b, ra, rb = rd(on, v), rd(off, v), rd(von, v), rd(voff, v)
    e, re = a - b, ra - rb
    M['slump_drop_' + v] = dict(ours=float(e[INIT] - e[INIT - 1]), v03=float(re[INIT] - re[INIT - 1]))
for v in ('TKLATC', 'TKRESPC', 'TKLATN'):
    a, ra = rd(on, v), rd(von, v)
    M['slump_month_' + v] = dict(ours=float(a[INIT]), v03=float(ra[INIT]))
    M['total_' + v] = dict(ours=float(np.nansum(a)), v03=float(np.nansum(ra)))
def recov(a, b, frac=0.9):
    ya, yb = a.reshape(-1, 12).sum(1), b.reshape(-1, 12).sum(1)
    y0 = INIT // 12 + 1
    for y in range(y0, len(ya)):
        if ya[y] >= frac * yb[y]: return y - INIT // 12
    return None
g, gb, rg, rgb = rd(on, 'GPP'), rd(off, 'GPP'), rd(von, 'GPP'), rd(voff, 'GPP')
M['gpp_recovery_years_to_90pct_of_off'] = dict(ours=recov(g, gb), v03=recov(rg, rgb))
M['gpp_first_full_year_after_slump_frac_of_off'] = dict(
    ours=float(g.reshape(-1, 12).sum(1)[INIT // 12 + 1] / gb.reshape(-1, 12).sum(1)[INIT // 12 + 1]),
    v03=float(rg.reshape(-1, 12).sum(1)[INIT // 12 + 1] / rgb.reshape(-1, 12).sum(1)[INIT // 12 + 1]))
M['gpp_last_year_frac_of_off'] = dict(ours=float(g[-12:].sum() / gb[-12:].sum()), v03=float(rg[-12:].sum() / rgb[-12:].sum()))

t = np.arange(len(g)) / 12.0
fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for a_, v in zip(ax.flat, ('SOC', 'VEGC', 'GPP', 'TKLATC')):
    e = rd(on, v) - (rd(off, v) if v != 'TKLATC' else 0)
    re = rd(von, v) - (rd(voff, v) if v != 'TKLATC' else 0)
    a_.plot(t, re, color='0.5', lw=2.5, label='rts-regrow v0.3 (on - off)')
    a_.plot(t, e, color='C3', lw=1, label='this port (on - off_auto)')
    a_.axvline(INIT / 12, color='k', ls=':', lw=0.8); a_.set_title(v + (' (on - off)' if v != 'TKLATC' else ''))
    a_.set_ylabel('g/m2' + ('/month' if v in ('GPP', 'TKLATC') else ''))
ax[0, 0].legend(fontsize=8); ax[1, 0].set_xlabel('tr year'); ax[1, 1].set_xlabel('tr year')
fig.suptitle('RTS effect, 1-cell synthetic case (-p 10 -e 300 -s 30 -t 30): port vs v0.3')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'rts_effect_vs_v03.png'), dpi=110); plt.close(fig)

fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
for n, c, lab in (('on', 'C3', 'on (auto switches)'), ('on_fx', 'C2', 'on + soil_giveup_retry'), ('on_noauto', 'C0', 'on, auto switches forced 0'), ('off_auto', 'k', 'off_auto')):
    z = np.load(os.path.join(W, n, 'closure.npz')); tt = np.arange(1, len(z['C']) + 1) / 12.0
    ax[0].plot(tt, z['C'], c, lw=0.8, label=lab)
    if z['N'].size: ax[1].plot(np.arange(len(z['N'])) / 12.0, z['N'], c, lw=0.8, label=lab)
    ax[2].plot(tt, z['W'], c, lw=0.8, label=lab)
ax[0].set_ylabel('C residual g/m2/month'); ax[1].set_ylabel('N residual g/m2/month'); ax[2].set_ylabel('water residual mm/month')
ax[0].set_yscale('symlog', linthresh=1e-3); ax[1].set_yscale('symlog', linthresh=1e-3)
for a_ in ax: a_.axvline(INIT / 12, color='k', ls=':', lw=0.8)
ax[0].legend(fontsize=8, ncol=2); ax[2].set_xlabel('tr year')
fig.suptitle('Monthly closure residuals, synthetic case (tr)'); fig.tight_layout()
fig.savefig(os.path.join(OUT, 'rts_closure.png'), dpi=110); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
gy, gby, rgy, rgby = [x.reshape(-1, 12).sum(1) for x in (g, gb, rg, rgb)]
ax[0].plot(gy / gby, 'C3', label='port (config.js: Pizano n_rate1 2.2, t_active 1.0)'); 
_v3 = O('on_ta09')
_g3 = rd(_v3, 'GPP') if os.path.isdir(_v3) else None
if _g3 is not None:
    _g3y = _g3.reshape(-1, 12).sum(1); ax[0].plot(_g3y / gby, 'C0--', label='port, t_active 0.9 (not adopted)'); M['gpp_ratio_yearly_ta09'] = (_g3y / gby).tolist()
M['gpp_ratio_yearly'] = {'ours': (gy / gby).tolist(), 'v03': (rgy / rgby).tolist()}
ax[0].plot(rgy / rgby, color='0.5', lw=2.5, label='v0.3', zorder=0)
ax[0].axvline(INIT // 12, color='k', ls=':'); ax[0].set_title('annual GPP on / off'); ax[0].set_xlabel('tr year'); ax[0].legend()
for v, c in (('TKREGROWC', 'C1'), ('TKREGROWN', 'C4')):
    a = rd(on, v)
    if a is not None: ax[1].plot(t, np.cumsum(a), c, label=v + ' (cumulative)')
ax[1].set_title('RTS regrow external inputs (port)'); ax[1].set_xlabel('tr year'); ax[1].set_ylabel('g/m2'); ax[1].legend()
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'rts_regrow.png'), dpi=110); plt.close(fig)
json.dump(M, open(os.path.join(OUT, 'rts_port_metrics.json'), 'w'), indent=1)
print(json.dumps(M, indent=1))
