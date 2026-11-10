#!/usr/bin/env python3
"""Stage C variant case sets: mk_variant.py SRC DST CASES(comma|ALL) SPIN_MS_JSON MAIN_MS_JSON [--geomap VALUE]
Copies SRC (prepared case set) to DST, rewrites absolute paths, updates model_settings per stage.
--geomap writes main/input/geothermal_flux.nc (from run-mask.nc grid, constant VALUE W/m2) into both stages
and sets IO.geothermal_flux_file. SRC is not modified."""
import sys, json, shutil, argparse
from pathlib import Path
import numpy as np, netCDF4 as nc
ap = argparse.ArgumentParser(); ap.add_argument("src"); ap.add_argument("dst"); ap.add_argument("cases")
ap.add_argument("spin"); ap.add_argument("main"); ap.add_argument("--geomap", type=float)
a = ap.parse_args(); src = Path(a.src).resolve(); dst = Path(a.dst).resolve(); assert not dst.exists(), dst
cs = json.load(open(src / "caseset.json"))
cases = cs["cases"] if a.cases == "ALL" else a.cases.split(",")
dst.mkdir(parents=True)
ms = {"spinup": json.loads(a.spin), "main": json.loads(a.main)}
for c in cases:
    shutil.copytree(src / c, dst / c)
    cd = dst / c
    for st in ("spinup", "main"):
        cfg = cd / st / "config/config.js"; j = json.loads(cfg.read_text().replace(str(src), str(dst)))
        j["model_settings"].update(ms[st])
        if a.geomap is not None:
            rm = Path(j["IO"]["runmask_file"]); g = rm.parent / "geothermal_flux.nc"
            with nc.Dataset(rm) as r, nc.Dataset(g, "w") as w:
                for d, dim in r.dimensions.items(): w.createDimension(d, len(dim))
                v = w.createVariable("geothermal_flux", "f8", ("Y", "X") if "Y" in r.dimensions else tuple(r.dimensions))
                v[:] = np.full(v.shape, a.geomap); v.units = "W m-2"
            j["IO"]["geothermal_flux_file"] = str(g)
        t = json.dumps(j, indent=2) + "\n"; assert str(src) not in t; cfg.write_text(t)
    meta = json.loads((cd / "case.json").read_text().replace(str(src), str(dst))); meta["model_settings"].update(ms["main"])
    meta["derived_from"] = str(src); (cd / "case.json").write_text(json.dumps(meta, indent=1) + "\n")
cs.update({"cases": cases, "derived_from": str(src), "set": dst.name, "variant": {"spinup": ms["spinup"], "main": ms["main"], "geomap": a.geomap}})
(dst / "caseset.json").write_text(json.dumps(cs, indent=1) + "\n"); print("wrote", dst)
