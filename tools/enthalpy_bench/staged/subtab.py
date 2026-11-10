import json,os,sys,numpy as np,netCDF4 as nc,datetime as dt,glob
# Stage D sub-step study: TVC/HPC Warming_nfeed, thermal_substeps 1/2/4/8 vs 24 (reference) and legacy (-O2). T stats: main years 1-40 monthly; ALD: last 15 yr.
NY=15
def rd(d,v):
    with nc.Dataset(os.path.join(d,v)) as f:
        k=[x for x in f.variables if x not in('time','x','y','albers_conical_equal_area','lat','lon')][0]; a=np.array(f.variables[k][:],dtype=float)
    return a[...,0,0] if a.ndim>=3 else a
def tprof(d,depths):
    T=rd(d,'TLAYER_monthly_sc.nc');Z=rd(d,'LAYERDEPTH_monthly_sc.nc');DZ=rd(d,'LAYERDZ_monthly_sc.nc');o=np.full((T.shape[0],len(depths)),np.nan)
    for t in range(T.shape[0]):
        ok=(DZ[t]>0)&(DZ[t]<100)&(T[t]>-100);zc=(Z[t]+0.5*DZ[t])[ok];o[t]=np.where(np.array(depths)<=zc.max(),np.interp(depths,zc,T[t][ok]),np.nan)
    return o[:480]
D=os.path.expanduser('~/tem-dev/esc-work/stageD'); dep=[0.05,0.2,0.5,1.0]
for c in ['TVC_Warming_nfeed','HPC_Warming_nfeed']:
    R={}
    for s in ['dleg_o2','dsub1','dsub2','dsub4','dsub8','dsub24']:
        r=os.path.join(D,s,c); m=json.load(open(r+'/meta.json')); o=r+'/main/output'
        f=lambda k:dt.datetime.strptime(m[k],'%Y-%m-%d %H:%M:%S %z'); w=(f('finished')-f('started')).total_seconds()
        ed=glob.glob(o+'/energy_daily.csv'); sv=None
        if ed:
            import pandas as pd; e=pd.read_csv(ed[0]); sv=e['solves'].mean() if 'solves' in e else None
        R[s]=dict(w=w,ald=rd(o,'ALD_yearly_sc.nc')[-NY:],T=tprof(o,dep),sv=sv,sn=m['main']['soil_safetynet_events'])
    ref=R['dsub24']
    print('==',c,'(T: main yr1-40; ALD col: last 15 yr; ref = 24 substeps)'); print('run wall_s solves/day ALD dALD_vs24 meanT(5,20,50,100cm) rmseT_monthly_vs24(5,20,50,100cm) safetynet')
    for s,x in R.items():
        print(s,x['w'],None if x['sv'] is None else round(x['sv'],2),round(x['ald'].mean(),3),round((x['ald']-ref['ald']).mean(),3),
              np.round(np.nanmean(x['T'],0),2).tolist(),np.round(np.sqrt(np.nanmean((x['T']-ref['T'])**2,0)),3).tolist(),x['sn'])
