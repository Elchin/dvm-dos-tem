import pandas as pd, sys, os, glob, json
for s in sys.argv[1:]:
  for st in ["main","spin000"]:
    f=f"{s}/HPC_Warming_nfeed/{st}/output/integrator_giveup_diag.csv"
    if not os.path.exists(f): print(s,st,"no csv"); continue
    d=pd.read_csv(f); g=d[d.giveup==1]
    print(f"== {s} {st}: giveups={len(g)} rows={len(d)}")
    if len(g)==0: continue
    print(" reason/pool:", g.groupby(["reason","pool"]).size().to_dict())
    print(" layer:", g.layer.value_counts().sort_index().to_dict())
    print(" month:", g.month.value_counts().sort_index().to_dict())
    print(" year_index range:", g.year_index.min(), g.year_index.max(), " t0 median %.3f"%g.t0.median())
    print(" median avln %.3g avln0 %.3g rtn %.3g netnmin %.3g sum_rtn %.3g sum_avln0 %.3g nscale %.3g totnextract %.3g totnetnmin %.3g"%tuple(g[c].median() for c in ["avln","avln_monthstart","rtnextract","netnmin","sum_rtnextract","sum_avln_monthstart","nscale","totnextract","totnetnmin"]))
