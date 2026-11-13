#!/usr/bin/env python3
"""identdir.py DIR_A DIR_B : data-level identity of two dvmdostem output directories.

Every file in either directory is compared: NetCDF files variable by variable (raw bytes of the data,
so NaN/fill patterns count; `total_runtime` ignored), other files byte for byte; `config_record.js`
is skipped. Restart files (restart-*.nc) are included. Prints one line
  files=N diff=D identical=True|False [first differences]
and exits 0 if identical, 1 otherwise. For whole run directories (spin000/ and main/) use ident.py.
"""
import sys
from pathlib import Path
import numpy as np, netCDF4 as nc
IGN = {"config_record.js"}; IGV = {"total_runtime"}

def cmp_nc(a, b):
    with nc.Dataset(a) as x, nc.Dataset(b) as y:
        bad = [v for v in x.variables if v not in IGV and (v not in y.variables or np.ma.getdata(x.variables[v][:]).tobytes() != np.ma.getdata(y.variables[v][:]).tobytes())]
        bad += [v for v in y.variables if v not in x.variables and v not in IGV]
        return bad

def compare_dirs(A, B):
    A, B = Path(A), Path(B); n = 0; d = []
    for name in sorted(set(p.name for p in A.iterdir()) | set(p.name for p in B.iterdir())):
        if name in IGN: continue
        f, g = A / name, B / name
        if not (f.exists() and g.exists()): d.append("missing " + name); continue
        n += 1
        if f.suffix == ".nc":
            b = cmp_nc(f, g)
            if b: d.append(f"{name}:{b[:3]}")
        elif f.read_bytes() != g.read_bytes(): d.append(name)
    return n, d

if __name__ == "__main__":
    if len(sys.argv) != 3: sys.exit(__doc__)
    n, d = compare_dirs(sys.argv[1], sys.argv[2])
    print(f"files={n} diff={len(d)} identical={not d}", d[:8])
    sys.exit(0 if not d else 1)
