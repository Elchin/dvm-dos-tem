"""Restart continuity: RUNS/<run>_split/B (restarted at year 35) vs the continuous RUNS/<run>/main for E1on_fb3, E1on_fb2, E1_fb3.
Writes FIGDIR/split_summary.json. Options: see common.py."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
_a = common.finish(common.parser(__doc__, "fire_balance").parse_args())
import sys, numpy as np, netCDF4 as nc, json
from pathlib import Path
R=_a.runs; S=35
KEY=["DEADC","DEADN","DWDC","DWDN","BURNAIR2SOILN","ORGN","AVLN","VEGC","SOC","NPP","RHSOM","GPP","TLAYER","SNOWTHICK","LAI"]
def get(o,v):
    for k,m in (("monthly",12),("yearly",1)):
        f=o/f"{v}_{k}_sc.nc"
        if f.exists():
            x=np.ma.filled(nc.Dataset(f)[v][:].astype(float),np.nan); return x.reshape(x.shape[0],-1), m
    return None,None
out={}
for c in ("E1on_fb3","E1on_fb2","E1_fb3"):
    C=R/c/"main/output"; B=R/(c+"_split")/"B/output"; d={}
    ident=0; tot=0
    for f in sorted(B.glob("*_sc.nc")):
        if f.name.startswith("restart"): continue
        k=365 if "_daily_" in f.name else 12 if "_monthly_" in f.name else 1
        v=f.name.split("_"+f.name.split("_")[-2])[0]
        with nc.Dataset(f) as a, nc.Dataset(C/f.name) as cc:
            x=np.ma.filled(a[v][:].astype(float),np.nan); y=np.ma.filled(cc[v][:].astype(float),np.nan)[S*k:S*k+x.shape[0]]
        tot+=1; ident+= np.array_equal(x,y,equal_nan=True)
    d["vars_identical"]=f"{ident}/{tot}"
    for v in KEY:
        x,m=get(B,v); y,_=get(C,v)
        if x is None: continue
        y=y[S*m:S*m+x.shape[0]]; xs=x.sum(1); ys=y.sum(1)
        dd=np.abs(xs-ys); d[v]=dict(maxabs=float(np.nanmax(dd)),maxrel=float(np.nanmax(dd)/max(np.nanmax(np.abs(ys)),1e-30)),first_month_abs=float(dd[0]),yr1_maxabs=float(np.nanmax(dd[:m])))
    # restart-sc.nc at the end (year 75)
    with nc.Dataset(B/"restart-sc.nc") as a, nc.Dataset(C/"restart-sc.nc") as cc:
        n=0; same=0
        for v in a.variables:
            if v in a.dimensions: continue
            n+=1; same+= np.array_equal(np.ma.getdata(a[v][:]), np.ma.getdata(cc[v][:]))
    d["restart_end_vars_identical"]=f"{same}/{n}"
    out[c]=d
print(json.dumps(out,indent=1))
json.dump(out,open(_a.figdir/"split_summary.json","w"),indent=1)
