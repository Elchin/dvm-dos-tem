#!/usr/bin/env python3
"""ident.py: bit identity of pairs of run directories (as written by run_case.py).

  ident.py [--json OUT] RUN_A RUN_B [RUN_A RUN_B ...]
  ident.py [--json OUT] --pairs FILE          (FILE: one "RUN_A RUN_B" per line, # comments)

For every pair, each stage folder (spin000/output, main/output) is compared with identdir.compare_dirs:
all NetCDF outputs and restart files, variable by variable. Prints one line per stage folder and a total
"identical stage folders: k/n"; exits 0 only if all are identical.
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from identdir import compare_dirs

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="*"); p.add_argument("--pairs", type=Path); p.add_argument("--json", type=Path)
    p.add_argument("--stages", default="spin000,main")
    a = p.parse_args()
    pairs = [tuple(a.runs[i:i + 2]) for i in range(0, len(a.runs), 2)]
    if a.pairs:
        pairs += [tuple(l.split()[:2]) for l in a.pairs.read_text().splitlines() if l.strip() and not l.startswith("#")]
    if not pairs or any(len(x) != 2 for x in pairs): sys.exit(__doc__)
    out = {}; ok = n = 0
    for A, B in pairs:
        for st in a.stages.split(","):
            da, db = Path(A) / st / "output", Path(B) / st / "output"; key = f"{B}/{st}"
            if not (da.is_dir() and db.is_dir()):
                out[key] = dict(a=str(da), identical=False, error="missing output dir"); n += 1
                print("MISSING", key); continue
            nf, d = compare_dirs(da, db); n += 1; ok += not d
            out[key] = dict(a=str(da), files=nf, diff=len(d), identical=not d, examples=d[:6])
            print("IDENTICAL" if not d else "DIFF     ", key, f"files={nf}", d[:3])
    print(f"identical stage folders: {ok}/{n}")
    if a.json: a.json.write_text(json.dumps(dict(identical=ok, total=n, folders=out), indent=1) + "\n")
    sys.exit(0 if ok == n else 1)

if __name__ == "__main__":
    main()
