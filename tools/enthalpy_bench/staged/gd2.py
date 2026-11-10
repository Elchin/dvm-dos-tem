import pandas as pd, numpy as np, netCDF4 as nc, sys
for s in sys.argv[1:]:
  o=f"{s}/HPC_Warming_nfeed/main/output/"
  d=pd.read_csv(o+"integrator_giveup_diag.csv"); g=d[d.giveup==1]
  T=nc.Dataset(o+"TLAYER_monthly_sc.nc").variables["TLAYER"][:,:,0,0] if True else None
  Z=nc.Dataset(o+"LAYERDEPTH_monthly_sc.nc").variables["LAYERDEPTH"][:,:,0,0]
  DZ=nc.Dataset(o+"LAYERDZ_monthly_sc.nc").variables["LAYERDZ"][:,:,0,0]
  rows=[]
  for _,r in g.iterrows():
    t=r.year_index*12+r.month; l=int(r.layer)
    th=np.where(T[t]>0)[0]; deep=th.max() if len(th) else -1
    rows.append((l, float(Z[t,l]), float(DZ[t,l]), float(T[t,l]), deep, l==deep))
  x=pd.DataFrame(rows,columns=["lay","ztop","dz","T","deepest_thawed","is_front_layer"])
  print(s, "n=",len(x)); print(x.describe().round(3).loc[["mean","50%","min","max"]].to_string()); print(" frac event layer = deepest thawed layer: %.2f"%x.is_front_layer.mean(), " frac T<1C: %.2f"%(x["T"]<1).mean())
