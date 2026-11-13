"""N leak fixes (dwd_n_leak_fix, deadroot_n_route_fix): N/C balance gaps before/after from RUNS/{F1c_onL,F1c_onLK,F1c_onLD,F1c_onLR,F1p_onL,F1p_onLK,E1on,E1onL}_fb4.
Writes FIGDIR/{n_leak_summary.json,n_leak_gaps.png}. Options: see common.py."""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
_a = common.finish(common.parser(__doc__, "fire_balance").parse_args())
# N leak fixes (dwd_n_leak_fix, deadroot_n_route_fix): N and C balance gaps before/after, from the
# rbsmoke runs (dvmdostem_fb4). Balance definitions as analyze_fire_balance.py.
import netCDF4 as nc, numpy as np, os, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R=str(_a.runs)+"/"; FIG=str(_a.figdir)+"/"
def rd(o,v,kind):
    x=np.ma.filled(nc.Dataset(o+f"{v}_{kind}_sc.nc")[v][:].astype(float),0.0); return x.reshape(x.shape[0],-1).sum(1)
def has(o,v,kind): return os.path.exists(o+f"{v}_{kind}_sc.nc")
def yr_flux(o,v): return rd(o,v,"monthly").reshape(-1,12).sum(1)
def stock(o,v):
    if has(o,v,"yearly"): return rd(o,v,"yearly")
    return rd(o,v,"monthly").reshape(-1,12)[:,11]
def ana(c):
    o=R+c+"/main/output/"; s=R+c+"/spin000/output/"
    C=stock(o,"VEGC")+stock(o,"DEADC")+stock(o,"DWDC")+stock(o,"SOC")
    burnC=yr_flux(o,"BURNVEG2AIRC")+yr_flux(o,"BURNSOIL2AIRC")
    rC=np.diff(C)-(yr_flux(o,"NPP")-yr_flux(o,"RHSOM")-yr_flux(o,"RHDWD")-burnC)[1:]
    N=stock(o,"VEGNTOT")+stock(o,"ORGN")+stock(o,"AVLN")+stock(o,"DEADN")+stock(o,"DWDN")
    burnN=yr_flux(o,"BURNVEG2AIRN")+yr_flux(o,"BURNSOIL2AIRN")
    a2s=rd(o,"BURNAIR2SOILN","monthly"); a2s_prev=np.concatenate([[rd(s,"BURNAIR2SOILN","monthly")[-1]],a2s[:-1]])
    dep=(a2s_prev/12.).reshape(-1,12).sum(1)
    rN=np.diff(N)-(yr_flux(o,"NINPUT")-yr_flux(o,"NLOST")-burnN+dep)[1:]
    # expected DWD N leak (legacy): sum_m RHDWD_m * DWDN/DWDC (monthly stocks, if output)
    dwdleak=None
    if has(o,"DWDN","monthly") and has(o,"DWDC","monthly"):
        rh=rd(o,"RHDWD","monthly"); wc=rd(o,"DWDC","monthly"); wn=rd(o,"DWDN","monthly")
        wcp=np.concatenate([[rd(s,"DWDC","monthly")[-1]],wc[:-1]]); wnp=np.concatenate([[rd(s,"DWDN","monthly")[-1]],wn[:-1]])
        r=np.where(wcp>0,wnp/np.where(wcp>0,wcp,1),0); dwdleak=(rh*r).reshape(-1,12).sum(1)
    fy=[int(i) for i in np.nonzero(burnC)[0]]
    return dict(rC=rC,rN=rN,fy=fy,dwdleak=dwdleak,N=N,C=C,ORGN=stock(o,"ORGN"),SOC=stock(o,"SOC"),VEGC=stock(o,"VEGC"),DWDN=stock(o,"DWDN"))
PAIRS=[("F1c_onL_fb4","F1c_onLK_fb4","classic FRI 15"),("F1p_onL_fb4","F1p_onLK_fb4","process (Reburn)"),("E1on_fb4","E1onL_fb4","explicit classic (E1)")]
out={}; D={}
for a,b,lab in PAIRS:
    for c in (a,b):
        if os.path.exists(R+c+"/main/output/NPP_monthly_sc.nc"): D[c]=ana(c)
    if a not in D or b not in D: continue
    da,db=D[a],D[b]
    out[lab]=dict(before=a,after=b,fire_years=da["fy"],
        N_gap_sum_before=round(float(da["rN"].sum()),3),N_gap_sum_after=round(float(db["rN"].sum()),3),
        N_gap_maxabs_before=round(float(np.abs(da["rN"]).max()),4),N_gap_maxabs_after=round(float(np.abs(db["rN"]).max()),4),
        N_gap_fireyrs_before=[round(float(da["rN"][y-1]),4) for y in da["fy"] if y>0],
        N_gap_fireyrs_after=[round(float(db["rN"][y-1]),4) for y in db["fy"] if y>0],
        C_gap_sum_before=round(float(da["rC"].sum()),3),C_gap_sum_after=round(float(db["rC"].sum()),3),
        C_gap_maxabs_before=round(float(np.abs(da["rC"]).max()),4),C_gap_maxabs_after=round(float(np.abs(db["rC"]).max()),4),
        C_gap_fireyrs_before=[round(float(da["rC"][y-1]),3) for y in da["fy"] if y>0],
        C_gap_fireyrs_after=[round(float(db["rC"][y-1]),3) for y in db["fy"] if y>0],
        expected_dwd_leak_sum_before=None if da["dwdleak"] is None else round(float(da["dwdleak"][1:].sum()),3),
        N_end_before=round(float(da["N"][-1]),2),N_end_after=round(float(db["N"][-1]),2),
        ORGN_end_before=round(float(da["ORGN"][-1]),2),ORGN_end_after=round(float(db["ORGN"][-1]),2),
        SOC_end_before=round(float(da["SOC"][-1]),1),SOC_end_after=round(float(db["SOC"][-1]),1),
        VEGC_end_before=round(float(da["VEGC"][-1]),1),VEGC_end_after=round(float(db["VEGC"][-1]),1),
        DWDN_end_before=round(float(da["DWDN"][-1]),3),DWDN_end_after=round(float(db["DWDN"][-1]),3))
single={}
for c,lab in (("F1c_onL_fb4","none"),("F1c_onLD_fb4","dwd_n_leak_fix only"),("F1c_onLR_fb4","deadroot_n_route_fix only"),("F1c_onLK_fb4","both")):
    d=ana(c); single[lab]=dict(run=c,N_gap_sum=round(float(d["rN"].sum()),3),N_gap_fireyrs=[round(float(d["rN"][y-1]),4) for y in d["fy"] if y>0],
        N_gap_nonfire_sum=round(float(sum(d["rN"][y-1] for y in range(1,len(d["N"])) if y not in d["fy"])),3),
        C_gap_fireyrs=[round(float(d["rC"][y-1]),4) for y in d["fy"] if y>0])
out["classic FRI 15, one switch at a time"]=single
print(json.dumps(out,indent=1)); json.dump(out,open(FIG+"n_leak_summary.json","w"),indent=1)
labs=[p for p in PAIRS if p[2] in out]
f,ax=plt.subplots(len(labs),2,figsize=(10,2.6*len(labs)),dpi=100,sharex=True,squeeze=False)
for i,(a,b,lab) in enumerate(labs):
    for j,(key,un) in enumerate((("rN","g N"),("rC","g C"))):
        x=1952+np.arange(len(D[a][key]))
        ax[i,j].plot(x,D[a][key],color="tab:orange",lw=1.2,label="leak fixes off (before)")
        ax[i,j].plot(x,D[b][key],color="k",lw=1.0,ls="--",label="dwd_n_leak_fix + deadroot_n_route_fix = 1 (after)")
        for y in D[a]["fy"]:
            if y>0: ax[i,j].axvline(1951+y,color="0.85",lw=0.6,zorder=0)
        ax[i,j].axhline(0,color="0.5",lw=0.5)
        ax[i,j].set_title(f"{lab}: {'N' if key=='rN' else 'C'} balance gap (sum {out[lab][('N' if key=='rN' else 'C')+'_gap_sum_before']:+.2f} → {out[lab][('N' if key=='rN' else 'C')+'_gap_sum_after']:+.2f} {un}/m²)",fontsize=8.5)
        ax[i,j].set_ylabel(f"{un} m$^{{-2}}$ yr$^{{-1}}$",fontsize=8); ax[i,j].tick_params(labelsize=7.5)
ax[0,0].legend(fontsize=7,loc="lower left"); ax[-1,0].set_xlabel("year"); ax[-1,1].set_xlabel("year")
f.suptitle("Annual balance gap ΔStock − (inputs − outputs), main stage 75 yr; grey lines = fire years",fontsize=9)
f.tight_layout(); f.savefig(FIG+"n_leak_gaps.png"); print("fig ok")
