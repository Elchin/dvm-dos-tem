import netCDF4 as nc, numpy as np, os, json, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=os.path.expanduser("~/tem-dev/wet-work/rbsmoke/runs/"); FIG=os.path.expanduser("~/tem-dev/rb/docs/figures/reburn/"); os.makedirs(FIG,exist_ok=True)
CASES=[("F0g","fire 0 (stage asks for fire)","k"),("F1c","fire 1, classic","tab:orange"),("F1p_c","fire 1, process (Reburn)","tab:red")]
def rd(o,v,kind):
    f=o+f"{v}_{kind}_sc.nc"; x=np.ma.filled(nc.Dataset(f)[v][:].astype(float),0.0); x=x.reshape(x.shape[0],-1).sum(1)  # sum pft/layer dims
    return x
def yr_flux(o,v): m=rd(o,v,"monthly"); return m.reshape(-1,12).sum(1)
def yr_stock_m(o,v):
    if os.path.exists(o+f"{v}_yearly_sc.nc"): return rd(o,v,"yearly")
    m=rd(o,v,"monthly"); return m.reshape(-1,12)[:,11]
res={}
for c,_,_ in CASES:
  for st in ("spin000","main"):
    o=R+c+"/"+st+"/output/"
    burnC=yr_flux(o,"BURNVEG2AIRC")+yr_flux(o,"BURNSOIL2AIRC")
    if st=="spin000":
        res[(c,st)]=dict(burnC=burnC,fire_years=[int(i) for i in np.nonzero(burnC)[0]],burnN=yr_flux(o,"BURNVEG2AIRN")+yr_flux(o,"BURNSOIL2AIRN"),VEGC=rd(o,"VEGC","yearly"),SOC=yr_stock_m(o,"SOC")); continue
    C=yr_stock_m(o,"VEGC")+yr_stock_m(o,"DEADC")+yr_stock_m(o,"DWDC")+yr_stock_m(o,"SOC")  # DEEPC (amorphous SOM) is part of SOC
    burnC=yr_flux(o,"BURNVEG2AIRC")+yr_flux(o,"BURNSOIL2AIRC")
    fC=yr_flux(o,"NPP")-yr_flux(o,"RHSOM")-yr_flux(o,"RHDWD")-burnC
    rC=np.diff(C)-fC[1:]
    N=yr_stock_m(o,"VEGNTOT")+yr_stock_m(o,"ORGN")+yr_stock_m(o,"AVLN")+yr_stock_m(o,"DEADN")+yr_stock_m(o,"DWDN")
    burnN=yr_flux(o,"BURNVEG2AIRN")+yr_flux(o,"BURNSOIL2AIRN")-yr_flux(o,"BURNAIR2SOILN")/12.  # monthly output holds the annual re-deposition rate (deposited /12 per month)
    fN=yr_flux(o,"NINPUT")-yr_flux(o,"NLOST")-burnN
    rN=np.diff(N)-fN[1:]
    fy=[int(i) for i in np.nonzero(burnC)[0]]
    res[(c,st)]=dict(C=C,N=N,burnC=burnC,burnN=burnN,rC=rC,rN=rN,fire_years=fy,
        VEGC=rd(o,"VEGC","yearly"),SOC=yr_stock_m(o,"SOC"))
summ={}
for (c,st),d in res.items():
    if st=="spin000":
        summ[f"{c}/{st}"]=dict(fire_years=d["fire_years"],burnC=[round(float(d["burnC"][y]),2) for y in d["fire_years"]],burnN=[round(float(d["burnN"][y]),3) for y in d["fire_years"]],VEGC_end=round(float(d["VEGC"][-1]),1),SOC_end=round(float(d["SOC"][-1]),1)); continue
    fy=[y for y in d["fire_years"] if y>0]; nf=[y for y in range(1,len(d["C"])) if y not in d["fire_years"]]
    summ[f"{c}/{st}"]=dict(fire_years=d["fire_years"], burnC=[round(float(d["burnC"][y]),2) for y in d["fire_years"]],
        burnN=[round(float(d["burnN"][y]),3) for y in d["fire_years"]],
        C_resid_fire_yrs=[round(float(d["rC"][y-1]),3) for y in fy], C_resid_nonfire_maxabs=round(float(np.max(np.abs(d["rC"][[y-1 for y in nf]]))),3),
        N_resid_fire_yrs=[round(float(d["rN"][y-1]),4) for y in fy], N_resid_nonfire_maxabs=round(float(np.max(np.abs(d["rN"][[y-1 for y in nf]]))),4),
        C_end=round(float(d["C"][-1]),1), VEGC_end=round(float(d["VEGC"][-1]),1), SOC_end=round(float(d["SOC"][-1]),1), N_end=round(float(d["N"][-1]),2))
print(json.dumps(summ,indent=1))
json.dump(summ,open(FIG+"../../../tests/regression/reburn_smoke_summary.json" if False else "/tmp/rb_smoke_summary.json","w"),indent=1)
# Fig 1: burn emissions per year (spin + main concatenated)
f,a=plt.subplots(figsize=(6.4,3.4),dpi=100)
for i,(c,l,k) in enumerate(CASES):
    y=np.concatenate([res[(c,"spin000")]["burnC"],res[(c,"main")]["burnC"]]); x=np.arange(len(y))
    a.bar(x+(i-1)*0.3,y,width=0.3,color=k,label=l)
a.axvline(19.5,color="gray",ls=":"); a.text(20.5,a.get_ylim()[1]*0.85,"spin-up | main (1950–2024)",fontsize=8)
a.set_yscale("symlog",linthresh=10); a.set_xlabel("run year (20-yr spin-up, then 75-yr main)"); a.set_ylabel("fire C emission (g C m$^{-2}$ yr$^{-1}$)")
a.set_title("TVC (CMT51) smoke: fire events, FRI 15 yr, July, severity 3"); a.legend(fontsize=8); a.grid(alpha=.3,axis="y")
f.tight_layout(); f.savefig(FIG+"fire_events.png"); plt.close(f)
# Fig 2: C pools on vs off (main)
f,ax=plt.subplots(1,2,figsize=(7.2,3.2),dpi=100)
for c,l,k in CASES:
    d=res[(c,"main")]; yrs=1950+np.arange(len(d["VEGC"]))
    ax[0].plot(yrs,d["VEGC"],color=k,label=l); ax[1].plot(yrs,d["SOC"],color=k,label=l)
ax[0].set_title("VEGC (g C m$^{-2}$)",fontsize=9); ax[1].set_title("SOC, December (g C m$^{-2}$)",fontsize=9)
for a in ax: a.grid(alpha=.3); a.set_xlabel("year")
ax[0].legend(fontsize=7); f.suptitle("C pools, main stage: fire off vs classic vs process",fontsize=10)
f.tight_layout(); f.savefig(FIG+"c_pools_on_off.png"); plt.close(f)
# Fig 3: C balance residual
f,a=plt.subplots(figsize=(6.4,3.0),dpi=100)
for c,l,k in CASES:
    r=res[(c,"main")]["rC"]; a.plot(1951+np.arange(len(r)),r,color=k,lw=1,label=l)
for y in res[("F1c","main")]["fire_years"]: a.axvline(1950+y,color="gray",ls=":",lw=.8)
a.set_xlabel("year (main stage; dotted = fire years)"); a.set_ylabel("ΔC − (NPP − RH − fire C) (g C m$^{-2}$)"); a.set_title("Annual ecosystem C balance residual",fontsize=10); a.legend(fontsize=7); a.grid(alpha=.3)
f.tight_layout(); f.savefig(FIG+"c_balance_residual.png"); plt.close(f)
