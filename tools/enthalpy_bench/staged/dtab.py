import json,sys,pandas as pd
d=pd.DataFrame(json.load(open(sys.argv[1]))); pat=sys.argv[2] if len(sys.argv)>2 else "NoWarming"
pd.set_option("display.width",300); pd.set_option("display.max_columns",60); pd.set_option("display.max_rows",200)
order=["leg_dc0","ts1_dc0","s0_dc0","leg_dc1","ts1_dc1","s0_dc1","legprod_dc1","ts1prod_dc1"]
for site in ["TVC","HPC"]:
  x=d[(d.site==site)&(d.case.str.contains(pat))]
  cols=[c for c in x.columns if x[c].dtype!=object and (c.endswith("bias") or c.endswith("rmse") or c in("ALD","peakSWE","snowoff"))]
  g=x.groupby("run")[cols].mean().reindex([o for o in order if o in set(x.run)])
  print("==",site,pat,"(mean over",x.case.nunique(),"cases)"); print(g.round(2).T.to_string())
