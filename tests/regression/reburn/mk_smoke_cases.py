# Build fire smoke-test case copies of the TVC_Warming_nfeed regression case (originals untouched).
# usage: mk_smoke_cases.py [--regcases DIR] [--cases OUTDIR] [--fireweed DIR] [--only F0g,F1c,F1p]
#   source case: REGCASES/leg/TVC_Warming_nfeed (calib-elchin only for now); writes CASES/{F0g,F1c,F1p}
#   (variants of these: mk_variants.py with variants.json)
import json, re, shutil, sys
from pathlib import Path
import netCDF4 as nc
sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
_p = common.parser("Build the F0g/F1c/F1p fire smoke cases.")
_p.add_argument("--fireweed", type=Path, default=common.REPO/"external/Fireweed/Input_Files")
_p.add_argument("--only", default="")
_a = common.finish(_p.parse_args())
SRC=_a.regcases/"leg/TVC_Warming_nfeed"
OUT=_a.cases; OUT.mkdir(parents=True, exist_ok=True)
FW=_a.fireweed.resolve()
BURN=["BURNAIR2SOILN","BURNSOIL2AIRC","BURNSOIL2AIRN","BURNTHICK","BURNVEG2AIRC","BURNVEG2AIRN","BURNVEG2DEADC","BURNVEG2DEADN","BURNVEG2SOILABVC","BURNVEG2SOILABVN","BURNVEG2SOILBLWC","BURNVEG2SOILBLWN"]
REBURN_CMT51=["123    //cmt2fm:// from Reburn ce80b344 cmt_firepar.txt, FM GS3 (123)",
              "-1.0   //cbd://",
              "-1.0   //cbh://"]
CASES={
 "F0g": dict(fire=0, proc=False, fp21=False),   # master switch off, stage asks for fire -> must equal fire-off
 "F1c": dict(fire=1, proc=False, fp21=False),   # classic fire, FRI 15 yr, 18-line firepar
 "F1p": dict(fire=1, proc=True,  fp21=True),    # process-based (Fireweed/Burnup), 21-line firepar
}
def strip_js(t): return re.sub(r"(?m)^\s*//.*$|(?<=[,{\[\s])//[^\n\"]*$","",t)
for name,o in CASES.items():
    if _a.only and name not in _a.only.split(","): continue
    d=OUT/name
    if d.exists(): sys.exit(f"{d} exists")
    shutil.copytree(SRC,d)
    cj=json.loads((d/"case.json").read_text())
    for st in ("spinup","main"):
        old=cj["stage_root"][st]; new=str(d/st); cj["stage_root"][st]=new
        c=d/st/"config/config.js"; t=c.read_text().replace(old,new)
        j=json.loads(strip_js(t))
        j["model_settings"]["fire"]=o["fire"]
        j["stage_settings"]["sc"]["fire"]=True
        j.setdefault("module_settings",{}).setdefault("dsb",{})["fire"]={"ignition_tr":1,"ignition_sc":1,"fire_process_wildfire":o["proc"]}
        for k,f in [("fire_fuel_model_file","Standard_Fuel_Models_Table.csv")]+[(f"fire_fosberg_{x}_file",f"NWGC_Fosberg_Table_{x.upper()}.txt") for x in "abcd"]:
            j["IO"][k]=str(FW/f)
        c.write_text(json.dumps(j,indent=2)+"\n")
        # FRI input copy: 15-yr return interval, burn on day 196 (mid-July), severity class 3, full area
        f=d/st/"input/fri-fire.nc"
        with nc.Dataset(f,"a") as x:
            x["fri"][:]=15; x["fri_jday_of_burn"][:]=196; x["fri_severity"][:]=3; x["fri_area_of_burn"][:]=1
        s=d/st/"config/output_spec.csv"; L=s.read_text().splitlines()
        for i,l in enumerate(L):
            p=l.split(",")
            if p[0] in BURN: p[3]="y"; p[4]="m"; L[i]=",".join(p)
        s.write_text("\n".join(L)+"\n")
    if o["fp21"]:
        fp=d/"parameters/cmt_firepar.txt"; L=fp.read_text().splitlines()
        i0=[i for i,l in enumerate(L) if l.startswith("// CMT51")][0]
        i=i0+1
        while i<len(L) and "r_retain_n" not in L[i]: i+=1
        L[i+1:i+1]=REBURN_CMT51; fp.write_text("\n".join(L)+"\n")
    (d/"case.json").write_text(json.dumps(cj,indent=1)+"\n")
    print("made",d)
