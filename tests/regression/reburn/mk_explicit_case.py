import json, shutil, sys
from pathlib import Path
import netCDF4 as nc
SRC=Path.home()/"tem-dev/srd-work/cs/leg/TVC_Warming_nfeed"; d=Path.home()/"tem-dev/wet-work/rbsmoke/cs/E1"
assert not d.exists(); shutil.copytree(SRC,d)
cj=json.loads((d/"case.json").read_text())
BURN=["BURNAIR2SOILN","BURNSOIL2AIRC","BURNSOIL2AIRN","BURNTHICK","BURNVEG2AIRC","BURNVEG2AIRN","BURNVEG2DEADC","BURNVEG2DEADN","BURNVEG2SOILABVC","BURNVEG2SOILABVN","BURNVEG2SOILBLWC","BURNVEG2SOILBLWN"]
for st in ("spinup","main"):
    old=cj["stage_root"][st]; new=str(d/st); cj["stage_root"][st]=new
    c=d/st/"config/config.js"; j=json.loads(c.read_text().replace(old,new))
    j["model_settings"]["fire"]=1            # read by the new binary only
    j["stage_settings"]["sc"]["dsb"]=True    # legacy key; no "fire" key -> new binary falls back to dsb
    c.write_text(json.dumps(j,indent=2)+"\n")
    pf=Path(j["IO"]["proj_exp_fire_file"])
    with nc.Dataset(pf,"a") as x:
        n=x["exp_burn_mask"].shape[0]
        for i in (15,30,45,60):
            if i<n: x["exp_burn_mask"][i]=1; x["exp_jday_of_burn"][i]=196; x["exp_fire_severity"][i]=3; x["exp_area_of_burn"][i]=1
    s=d/st/"config/output_spec.csv"; L=s.read_text().splitlines()
    for i,l in enumerate(L):
        p=l.split(",")
        if p[0] in BURN: p[3]="y"; p[4]="m"; L[i]=",".join(p)
    s.write_text("\n".join(L)+"\n")
    print(st, pf, n)
(d/"case.json").write_text(json.dumps(cj,indent=1)+"\n")
