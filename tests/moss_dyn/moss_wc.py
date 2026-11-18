#!/usr/bin/env python3
"""Moss-layer liquid water and implied Chadburn k (thawed months), 1986-2015 monthly, stage2a_hist.
  moss_wc.py RUN [RUN...]"""
import sys, os, numpy as np, netCDF4
W = os.path.expanduser('~/tem-dev/moss-work/runs')
for r in sys.argv[1:]:
    D = lambda v: np.ma.filled(netCDF4.Dataset(os.path.join(W, r, 'stage2a_hist', '%s_monthly_tr.nc' % v))[v][:].astype(float), np.nan).reshape(-1, netCDF4.Dataset(os.path.join(W, r, 'stage2a_hist', '%s_monthly_tr.nc' % v))[v].shape[1])
    lwc, typ, dz, tl = D('LWCLAYER'), D('LAYERTYPE'), D('LAYERDZ'), D('TLAYER')
    u = netCDF4.Dataset(os.path.join(W, r, 'stage2a_hist', 'LWCLAYER_monthly_tr.nc'))['LWCLAYER'].units
    t0 = (1986 - 1950) * 12; t1 = (2016 - 1950) * 12
    m = typ[t0:t1, 0] == 0
    l, T = lwc[t0:t1, 0], tl[t0:t1, 0]
    th = T > 0
    jja = np.array([(i % 12) in (5, 6, 7) for i in range(t1 - t0)])
    print(r, 'units', u, 'top layer moss frac', m.mean().round(3), 'dz', np.nanmean(dz[t0:t1, 0]).round(4),
          'liq JJA mean', np.nanmean(l[jja]).round(4), 'liq thawed-months mean', np.nanmean(l[th]).round(4),
          'liq frozen-months mean', np.nanmean(l[~th]).round(4),
          'k_chad(thawed, S=liq/0.95)', (0.06 + 0.44 * np.clip(np.nanmean(l[th]) / 0.95, 0, 1)).round(3))
