import sys, netCDF4 as nc, numpy as np, csv, os
W=os.path.expanduser('~/tem-dev/ssf-work/runs/')
print('| run | ALD 1980-89 | ALD 1990-2009 | ALD 2010-24 | max ALD | worst full-day residual (W m-2) | 99.9th pct | whole-frozen dump days (bwater>1e5, nfront 0) | mean solves/day |')
print('|---|---|---|---|---|---|---|---|---|')
for v in sys.argv[1:]:
  o=W+v+'/main/output/'
  if not os.path.exists(o+'ALD_yearly_sc.nc'): print('|',v,'| missing |'); continue
  a=nc.Dataset(o+'ALD_yearly_sc.nc')['ALD'][:].squeeze()
  R=list(csv.DictReader(open(o+'energy_daily.csv')))
  rd=np.array([abs(float(r['resday'])) for r in R])/86400
  ev=sum(1 for r in R if int(r['nfront'])==0 and float(r['bwater'])>1e5)
  so=np.mean([float(r['solves']) for r in R])
  print(f"| {v} | {a[30:40].mean():.2f} | {a[40:60].mean():.2f} | {a[60:75].mean():.2f} | {a.max():.2f} | {rd.max():.1e} | {np.percentile(rd,99.9):.1e} | {ev} | {so:.0f} |")
