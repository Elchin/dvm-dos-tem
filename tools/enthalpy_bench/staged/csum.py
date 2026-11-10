import json,sys,glob,os
print("set n status_ok maxEnthDay maxUnattr nonconv Wpct(min..max) Cmaxann Nmean safetynet(main,list)")
for s in sys.argv[1:]:
    fs=sorted(glob.glob(s+"/*.closure.json")); R=[]
    for f in fs:
        try: R.append(json.load(open(f)))
        except Exception as e: R.append({"status":"BAD","case":os.path.basename(f)})
    ok=sum(r.get("status")=="OK" for r in R)
    g=lambda k:[r.get(k) for r in R if r.get(k) is not None]
    mx=lambda k:max(map(abs,g(k))) if g(k) else None
    w=g("W_resid_pctP")
    sn={r["case"]:r.get("main_safetynet") for r in R if r.get("main_safetynet")}
    print(os.path.basename(s.rstrip('/')),len(R),ok,mx("Enth_resday_max_abs_Wm2"),mx("Enth_unattr_max_abs_Wm2"),mx("Enth_fallback_failures"),
      "%.2f..%.2f"%(min(w),max(w)) if w else "-",mx("C_resid_max_abs_annual"),mx("N_resid_mean_gm2yr"),sn)
