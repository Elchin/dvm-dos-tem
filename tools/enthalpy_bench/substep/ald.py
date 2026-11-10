import sys, netCDF4 as nc, numpy as np, csv, os
W=os.path.expanduser('~/tem-dev/ssf-work/runs/')
print('run      ALD80-89 ALD90-09 ALD10-24 maxALD  firstyr>0.6  maxres(full-day) maxres_main  bwater_sum_main(J/m2)')
for v in sys.argv[1:]:
  o=W+v+'/main/output/'
  try: a=nc.Dataset(o+'ALD_yearly_sc.nc')['ALD'][:].squeeze()
  except Exception as e: print(v,'ERR',e); continue
  R=list(csv.DictReader(open(o+'energy_daily.csv')))
  rd=np.array([abs(float(r['resday'])) for r in R]); bw=np.array([float(r['bwater']) for r in R])
  big=np.where(a>0.6)[0]
  print(f"{v:8s} {a[30:40].mean():.2f} {a[40:60].mean():.2f} {a[60:75].mean():.2f} {a.max():.2f} {1950+big[0] if len(big) else '-':>6} {rd.max():.2e} {np.percentile(rd,99.9):.2e} {np.abs(bw).max():.2e}")
