"""Spin-up -> main re-deposition overshoot for RUNS/{F1c_on (19af750f), F1c_on_fb3, E1on_fb2, E1on_fb3}.
Writes FIGDIR/overshoot.json and WORK/overshoot_series.json (monthly series, for plot_restart_continuity.py)."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
_a = common.finish(common.parser(__doc__, "fire_balance").parse_args())
import netCDF4 as nc, numpy as np, json
from pathlib import Path
R=_a.runs
def rd(c,st,v):
    x=np.ma.filled(nc.Dataset(R/c/st/"output"/f"{v}_monthly_sc.nc")[v][:].astype(float),0.0); return x.reshape(x.shape[0],-1).sum(1)
out={}
for c in ("F1c_on","F1c_on_fb3","E1on_fb2","E1on_fb3"):
    a=np.concatenate([rd(c,"spin000","BURNAIR2SOILN"),rd(c,"main","BURNAIR2SOILN")])
    em=rd(c,"spin000","BURNVEG2AIRN")+rd(c,"spin000","BURNSOIL2AIRN")
    fm=int(np.nonzero(em)[0][0])            # spin fire month index
    nxt=20*12+15*12+6                       # main fire, year 15 July (Dsb after Bgc)
    dep=a[fm:nxt]/12.0                      # due in months fm+1..nxt (rate of the previous month)
    months=int(np.count_nonzero(dep))
    out[c]=dict(spin_fire_month=f"spin yr {fm//12} m{fm%12}", emitted_N=round(float(em[fm]),3), redeposited_N=round(float(dep.sum()),3),
                ratio=round(float(dep.sum()/em[fm]),4), months_depositing=months, main_redep_total=round(float((np.concatenate([[a[20*12-1]],rd(c,"main","BURNAIR2SOILN")[:-1]])/12).sum()),2),
                series=(a/12).tolist())
print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk!="series"} for k,v in out.items()},indent=1))
json.dump({k:{kk:vv for kk,vv in v.items() if kk!="series"} for k,v in out.items()},open(_a.figdir/"overshoot.json","w"),indent=1)
_a.work.mkdir(parents=True, exist_ok=True); json.dump(out,open(_a.work/"overshoot_series.json","w"))
