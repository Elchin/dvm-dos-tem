import netCDF4 as nc, numpy as np, os, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=os.path.expanduser("~/tem-dev/wet-work/rbsmoke/runs/"); FIG=os.path.expanduser("~/tem-dev/rb/docs/figures/fire_balance/"); os.makedirs(FIG,exist_ok=True)
def rd(o,v,kind):
    x=np.ma.filled(nc.Dataset(o+f"{v}_{kind}_sc.nc")[v][:].astype(float),0.0); return x.reshape(x.shape[0],-1).sum(1)
def yr_flux(o,v): return rd(o,v,"monthly").reshape(-1,12).sum(1)
def stock(o,v):
    if os.path.exists(o+f"{v}_yearly_sc.nc"): return rd(o,v,"yearly")
    return rd(o,v,"monthly").reshape(-1,12)[:,11]
def ana(c):
    o=R+c+"/main/output/"; s=R+c+"/spin000/output/"
    C=stock(o,"VEGC")+stock(o,"DEADC")+stock(o,"DWDC")+stock(o,"SOC")
    burnC=yr_flux(o,"BURNVEG2AIRC")+yr_flux(o,"BURNSOIL2AIRC")
    rC=np.diff(C)-(yr_flux(o,"NPP")-yr_flux(o,"RHSOM")-yr_flux(o,"RHDWD")-burnC)[1:]
    N=stock(o,"VEGNTOT")+stock(o,"ORGN")+stock(o,"AVLN")+stock(o,"DEADN")+stock(o,"DWDN")
    burnN=yr_flux(o,"BURNVEG2AIRN")+yr_flux(o,"BURNSOIL2AIRN")
    a2s=rd(o,"BURNAIR2SOILN","monthly"); a2s_prev=np.concatenate([[rd(s,"BURNAIR2SOILN","monthly")[-1]],a2s[:-1]])
    dep=(a2s_prev/12.).reshape(-1,12).sum(1)   # re-deposition due in each month (rate set by the previous month's Dsb step)
    rN=np.diff(N)-(yr_flux(o,"NINPUT")-yr_flux(o,"NLOST")-burnN+dep)[1:]
    fy=[int(i) for i in np.nonzero(burnC)[0]]
    return dict(C=C,N=N,rC=rC,rN=rN,dep=dep,burnC=burnC,fy=fy,VEGC=stock(o,"VEGC"),DEADC=stock(o,"DEADC"),DWDC=stock(o,"DWDC"),SOC=stock(o,"SOC"),ORGN=stock(o,"ORGN"),AVLN=stock(o,"AVLN"))
CASES={"F1c":"classic, fixes off","F1c_N":"classic, N fix only","F1c_C":"classic, C fix only","F1c_on":"classic, both fixes","F1c_L":"classic, layer fixes only","F1c_onL":"classic, both fixes + layer fixes","F1p_c":"process, fixes off","F1p_on":"process, both fixes","F1p_L":"process, layer fixes only","F1p_onL":"process, both fixes + layer fixes","F0g":"fire off"}
D={c:ana(c) for c in CASES}
summ={}
for c,d in D.items():
    fy=[y for y in d["fy"] if y>0]; nf=[y for y in range(1,len(d["C"])) if y not in d["fy"]]
    summ[c]=dict(label=CASES[c],fire_years=d["fy"],
      C_resid_fire_yrs=[round(float(d["rC"][y-1]),3) for y in fy], C_resid_maxabs_all=round(float(np.max(np.abs(d["rC"]))),3),
      C_resid_nonfire_maxabs=round(float(np.max(np.abs(d["rC"][[y-1 for y in nf]]))),3), C_resid_sum=round(float(d["rC"].sum()),2),
      N_resid_fire_yrs=[round(float(d["rN"][y-1]),4) for y in fy], N_resid_maxabs_all=round(float(np.max(np.abs(d["rN"]))),4),
      N_resid_sum=round(float(d["rN"].sum()),2), N_redep_due_sum=round(float(d["dep"][1:].sum()),2),
      C_end=round(float(d["C"][-1]),1),VEGC_end=round(float(d["VEGC"][-1]),1),DEADC_end=round(float(d["DEADC"][-1]),2),DWDC_end=round(float(d["DWDC"][-1]),1),SOC_end=round(float(d["SOC"][-1]),1),N_end=round(float(d["N"][-1]),2),ORGN_end=round(float(d["ORGN"][-1]),2))
print(json.dumps(summ,indent=1)); json.dump(summ,open(FIG+"fire_balance_summary.json","w"),indent=1)
yrs=lambda n:1951+np.arange(n)
st={"F1c":("tab:orange","-"),"F1c_N":("tab:blue","--"),"F1c_C":("tab:green","--"),"F1c_on":("k","-"),"F1c_onL":("tab:cyan",":"),"F1p_c":("tab:red","-"),"F1p_on":("tab:purple","-"),"F1p_onL":("tab:cyan",":")}
for key,lab,fn,cs in (("rC","ΔC − (NPP − RH − fire C)  (g C m$^{-2}$ yr$^{-1}$)","c_balance_gap.png",["F1c","F1c_C","F1c_on","F1c_onL","F1p_c","F1p_on","F1p_onL"]),
                     ("rN","ΔN − (NINPUT − NLOST − fire N + re-dep.)  (g N m$^{-2}$ yr$^{-1}$)","n_balance_gap.png",["F1c","F1c_N","F1c_on","F1c_onL","F1p_c","F1p_on","F1p_onL"])):
    f,ax=plt.subplots(2,1,figsize=(7.0,5.0),dpi=100,sharex=True)
    for c in cs:
        a=ax[0] if c.startswith("F1c") else ax[1]; k,ls=st[c]; r=D[c][key]; a.plot(yrs(len(r)),r,color=k,ls=ls,lw=1.2,label=CASES[c])
    for a,c in zip(ax,("F1c","F1p_c")):
        for y in D[c]["fy"]: a.axvline(1950+y,color="gray",ls=":",lw=.8)
        a.axhline(0,color="gray",lw=.5); a.grid(alpha=.3); a.legend(fontsize=7)
    ax[0].set_title("TVC smoke, classic fire (FRI 15): annual balance gap, fixes off vs on",fontsize=9); ax[1].set_title("process (Reburn) fire",fontsize=9)
    ax[1].set_xlabel("year (main stage; dotted = fire years)"); f.text(0.01,0.5,lab,rotation=90,va="center",fontsize=8)
    f.tight_layout(rect=(0.03,0,1,1)); f.savefig(FIG+fn); plt.close(f)
f,ax=plt.subplots(2,3,figsize=(9.6,5.0),dpi=100,sharex=True)
for a,v in zip(ax.flat,("VEGC","DEADC","DWDC","SOC","ORGN","AVLN")):
    for c in ("F0g","F1c","F1c_on"):
        k={"F0g":"gray","F1c":"tab:orange","F1c_on":"k"}[c]; d=D[c][v]; a.plot(1950+np.arange(len(d)),d,color=k,lw=1.1,label=CASES[c])
    a.set_title(v+(" (g C m$^{-2}$)" if v.endswith("C") else " (g N m$^{-2}$)"),fontsize=9); a.grid(alpha=.3)
ax[0,0].legend(fontsize=7); f.suptitle("Pools, classic fire: fixes off vs on (fire off for reference)",fontsize=10)
f.tight_layout(); f.savefig(FIG+"pools_on_off.png"); plt.close(f)
