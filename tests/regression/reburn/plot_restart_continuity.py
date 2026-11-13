"""Figure FIGDIR/restart_continuity.png from the split runs and WORK/overshoot_series.json (run redep_overshoot.py first)."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
_a = common.finish(common.parser(__doc__, "fire_balance").parse_args())
import netCDF4 as nc, numpy as np, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
R=_a.runs; FIG=_a.figdir; S=35
def get(o,v):
    for k,m in (("monthly",12),("yearly",1)):
        f=o/f"{v}_{k}_sc.nc"
        if f.exists(): x=np.ma.filled(nc.Dataset(f)[v][:].astype(float),np.nan); return x.reshape(x.shape[0],-1).sum(1), m
f,ax=plt.subplots(1,2,figsize=(10.5,3.8),dpi=100)
a=ax[0]
for c,lab,col in (("E1on_fb2","old (19af750f)","tab:orange"),("E1on_fb3","new (fire state in restart)","k")):
    y,m=get(R/c/"main/output","DEADC"); t=1950+np.arange(len(y))/m
    yb,_=get(R/(c+"_split")/"B/output","DEADC"); tb=1950+S+np.arange(len(yb))/m
    a.plot(t,y,color=col,lw=2.6 if col!="k" else 1.0,label=lab+": continuous"); a.plot(tb,yb,color=col,ls="--",lw=1.4,label=lab+": restarted at 1985")
a.axvline(1950+S,color="gray",ls=":"); a.set_xlim(1975,2002); a.set_ylim(-1,30); a.set_ylabel("DEADC (g C m$^{-2}$)"); a.set_xlabel("year")
a.set_title("Standing dead after the 1980 fire, restart in 1985\n(TVC explicit fire, both fire fixes on)",fontsize=9); a.legend(fontsize=7); a.grid(alpha=.3)
o=json.load(open(_a.work/"overshoot_series.json")); a=ax[1]
for c,lab,col in (("F1c_on_fb3","new","k"),("F1c_on","old (19af750f)","tab:orange")):
    s=np.array(o[c]["series"]); fm=15*12+6; nxt=20*12+15*12+6
    dep=np.zeros_like(s); dep[1:]=s[:-1]; cum=np.cumsum(dep[fm+1:nxt+1]); t=(np.arange(fm+1,nxt+1))/12-20
    a.plot(t,cum,color=col,lw=1.4,label=f"{lab}: {o[c]['redeposited_N']:.0f} g N = {o[c]['ratio']:.2f} × emitted")
a.axhline(o["F1c_on_fb3"]["emitted_N"],color="gray",ls="--",lw=1,label=f"N emitted by the spin-up fire ({o['F1c_on_fb3']['emitted_N']:.0f} g)")
a.axvline(0,color="gray",ls=":"); a.text(0.3,10,"spin-up | main",fontsize=8)
a.set_xlabel("years from the start of the main stage"); a.set_ylabel("cumulative re-deposited N (g N m$^{-2}$)")
a.set_title("Re-deposition of the spin-up fire's N (classic FRI 15, fixes on)",fontsize=9); a.legend(fontsize=7,loc="lower right"); a.grid(alpha=.3)
f.tight_layout(); f.savefig(FIG/"restart_continuity.png"); print("ok")
