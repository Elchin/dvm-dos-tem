"""split_run.py BIN CASEDIR CONTRUN OUTDIR S N : main stage split at year S: A = years 0..S-1 from the spin restart
of CONTRUN, B = years S..N-1 from A's restart-sc.nc with projected climate and explicit fire sliced by S years."""
import json, shutil, subprocess, sys, os
from pathlib import Path
import netCDF4 as nc
sys.path.insert(0, str(Path(__file__).resolve().parents[1])); import run_case as rc   # this repo's tests/regression/run_case.py
b, case, cont, out, S, N = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve(), Path(sys.argv[4]).resolve(), int(sys.argv[5]), int(sys.argv[6])
cm = json.loads((case/"case.json").read_text()); assert not out.exists(); out.mkdir(parents=True)
rc.ENV["LD_LIBRARY_PATH"] = os.environ.get("REG_LD_LIBRARY_PATH", rc.ENV["LD_LIBRARY_PATH"])
def slice_nc(src, dst, n0, dim="time"):
    with nc.Dataset(src) as a, nc.Dataset(dst, "w", format=a.data_model) as z:
        z.setncatts({k: a.getncattr(k) for k in a.ncattrs()})
        for k, d in a.dimensions.items(): z.createDimension(k, (len(d) - n0) if k == dim else (None if d.isunlimited() else len(d)))
        for k, v in a.variables.items():
            fv = v.getncattr("_FillValue") if "_FillValue" in v.ncattrs() else None
            w = z.createVariable(k, v.datatype, v.dimensions, fill_value=fv)
            w.setncatts({x: v.getncattr(x) for x in v.ncattrs() if x != "_FillValue"})
            w[:] = v[n0:] if (v.dimensions and v.dimensions[0] == dim) else v[:]
A = out/"A"; rc.stage(case, cm, "main", "main", A, cont/"spin000/output/restart-sc.nc")
ra = rc.run(A, S, b); print("A", ra, flush=True)
B = out/"B"; rc.stage(case, cm, "main", "main", B, A/"output/restart-sc.nc")
j = json.loads((B/"config/config.js").read_text())
for key, n0 in (("proj_climate_file", S*365), ("proj_exp_fire_file", S)):
    f = Path(j["IO"][key]); tmp = f.with_suffix(".full.nc"); f.rename(tmp); slice_nc(tmp, f, n0); tmp.unlink()
rb = rc.run(B, N - S, b); print("B", rb, flush=True)
