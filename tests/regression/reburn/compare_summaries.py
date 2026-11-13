#!/usr/bin/env python3
"""compare_summaries.py NEW.json REF.json [NEW REF ...] [--rtol 1e-9] [--atol 0]

Recursive comparison of summary JSONs (numbers within rtol/atol; the summaries are already rounded, so the
default is exact up to float formatting; NaN equals NaN). Keys only in REF are reported as missing; keys only in NEW are listed
as extra (not a failure). Prints PASS/FAIL per file; exits 1 if any file fails.
"""
import json, math, sys

def walk(a, b, path, rtol, atol, out):
    if isinstance(b, dict):
        if not isinstance(a, dict): out.append(f"{path}: type {type(a).__name__} != dict"); return
        for k in b:
            if k not in a: out.append(f"{path}/{k}: missing"); continue
            walk(a[k], b[k], f"{path}/{k}", rtol, atol, out)
    elif isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b): out.append(f"{path}: list length {len(a) if isinstance(a, list) else '-'} != {len(b)}"); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f"{path}[{i}]", rtol, atol, out)
    elif isinstance(b, (int, float)) and not isinstance(b, bool) and isinstance(a, (int, float)):
        if not (math.isclose(a, b, rel_tol=rtol, abs_tol=atol) or (math.isnan(a) and math.isnan(b))): out.append(f"{path}: {a} != {b}")
    elif a != b:
        out.append(f"{path}: {a!r} != {b!r}")

def main(argv):
    rtol, atol, files = 1e-9, 0.0, []
    it = iter(argv)
    for x in it:
        if x == "--rtol": rtol = float(next(it))
        elif x == "--atol": atol = float(next(it))
        else: files.append(x)
    if not files or len(files) % 2: sys.exit(__doc__)
    bad = 0
    for new, ref in zip(files[::2], files[1::2]):
        try: A, B = json.load(open(new)), json.load(open(ref))
        except FileNotFoundError as e: print("FAIL", new, "vs", ref, e); bad += 1; continue
        d = []; walk(A, B, "", rtol, atol, d)
        print("PASS" if not d else "FAIL", new, "vs", ref, f"({len(d)} differences)" if d else "")
        for x in d[:10]: print("   ", x)
        bad += bool(d)
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main(sys.argv[1:])
