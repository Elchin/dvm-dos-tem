#!/usr/bin/env python3
"""Build the regression and smoke case variants from variants.json (no hand-made copies).

  mk_variants.py [--regcases DIR] [--cases DIR] [--only NAME,...] [--table variants.json]
      casesets (e.g. s4_f1): REGCASES/<from>  -> REGCASES/<name>
      smoke variants:        CASES/<from>     -> CASES/<name>
  mk_variants.py --out DIR --verify-against DIR [...]
      build into DIR instead (scratch) and compare every built variant with the existing one of the same
      name under the --verify-against roots (parsed config.js with paths normalised, every other file by md5).

A variant is a hard-linked copy (`cp -a --link`) of its base; only the files that change are rewritten (unlinked
first, so the base is never modified): config.js (model_settings overlay, paths), case.json (stage_root,
model_settings), caseset.json, output_spec.csv (`output_ym`: set Yearly=y, Monthly=m), plus VARIANT.json.
Bases: leg (calibrated regression cases), F0g/F1c/F1p (mk_smoke_cases.py), E1 (mk_explicit_case.py).
Existing variants are left alone (same overlay -> "exists", different -> error).
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common

def unlink_write(p, text):
    p = Path(p)
    if p.exists(): p.unlink()
    p.write_text(text)

def stamp_ok(dst, spec):
    s = dst / "VARIANT.json"
    if s.exists() and json.loads(s.read_text()).get("spec") == spec:
        print("exists", dst); return True
    raise SystemExit(f"{dst} exists and is not this variant (VARIANT.json missing or different)")

def build_caseset(src, dst, spec):
    """As srd-work/scripts/mk_variants.py (which made s4_f1): config.js indent 4, case.json/caseset.json indent 1."""
    overlay = spec["model_settings"]
    if dst.exists(): return stamp_ok(dst, spec)
    subprocess.check_call(["cp", "-a", "--link", str(src), str(dst)])
    for dirpath, _, files in os.walk(dst):
        for fn in files:
            if fn not in ("config.js", "case.json", "caseset.json"): continue
            p = Path(dirpath) / fn; t = p.read_text().replace(str(src), str(dst))
            j = json.loads(t)
            if fn == "config.js":
                j["model_settings"].update(overlay); t = json.dumps(j, indent=4) + "\n"
            elif fn == "case.json":
                j.setdefault("model_settings", {}).update(overlay)
                j["stage_root"] = {"spinup": str(p.parent / "spinup"), "main": str(p.parent / "main")}
                t = json.dumps(j, indent=1) + "\n"
            else:
                j["derived_from"] = str(src); j["set"] = dst.name; j["ms_overlay"] = overlay
                t = json.dumps(j, indent=1) + "\n"
            assert fn == "caseset.json" or str(src) not in t, p
            unlink_write(p, t)
    unlink_write(dst / "VARIANT.json", json.dumps({"src": str(src), "dst": str(dst), "label": spec.get("label", ""),
                                                   "overlay": overlay, "spec": spec}, indent=1) + "\n")
    print("built", dst)

def build_smoke(src, dst, spec):
    """Smoke variant: config.js indent 2 (as mk_smoke_cases.py), case.json indent 1."""
    if dst.exists(): return stamp_ok(dst, spec)
    subprocess.check_call(["cp", "-a", "--link", str(src), str(dst)])
    cj = json.loads((dst / "case.json").read_text())
    for st in ("spinup", "main"):
        old = cj["stage_root"][st]; new = str(dst / st); cj["stage_root"][st] = new
        c = dst / st / "config/config.js"
        j = json.loads(c.read_text().replace(old, new).replace(str(src), str(dst)))
        j["model_settings"].update(spec["model_settings"])
        unlink_write(c, json.dumps(j, indent=2) + "\n")
        if spec.get("output_ym"):
            s = dst / st / "config/output_spec.csv"; L = s.read_text().splitlines()
            for i, l in enumerate(L):
                p = l.split(",")
                if p[0] in spec["output_ym"]: p[3] = "y"; p[4] = "m"; L[i] = ",".join(p)
            unlink_write(s, "\n".join(L) + "\n")
    cj.setdefault("model_settings", {}).update(spec["model_settings"])
    unlink_write(dst / "case.json", json.dumps(cj, indent=1) + "\n")
    unlink_write(dst / "VARIANT.json", json.dumps({"src": str(src), "dst": str(dst), "spec": spec}, indent=1) + "\n")
    print("built", dst)

# ---- verification against existing (hand-made) variants ----
PATHRE = re.compile(r"/[^\"]*/[^/\"]+/(spinup|main)/")   # <anything>/<case>/{spinup,main}/ -> @/{spinup,main}/
def nsub(t): return PATHRE.sub(r"@/\1/", t)
def norm(v):
    if isinstance(v, dict): return {k: norm(x) for k, x in v.items()}
    if isinstance(v, list): return [norm(x) for x in v]
    if isinstance(v, str): return nsub(v)
    return v

def md5(p): return hashlib.md5(p.read_bytes()).hexdigest()

def compare(new, old):
    """Parsed configs (paths normalised) and the md5 of every other file except case/VARIANT/caseset/MANIFEST json."""
    res = {"config_identical": True, "config_bytes_identical": True, "model_settings_diff": {}, "other_config_diff": [],
           "files_compared": 0, "files_diff": []}
    skip = {"config.js", "case.json", "VARIANT.json", "caseset.json", "MANIFEST.json"}
    for c in sorted(new.rglob("config.js")):
        rel = c.relative_to(new); o = old / rel
        if not o.exists(): res["config_identical"] = False; res["other_config_diff"].append(f"missing {rel}"); continue
        a, b = norm(json.loads(c.read_text())), norm(json.loads(o.read_text()))
        if nsub(c.read_text()) != nsub(o.read_text()): res["config_bytes_identical"] = False
        if a != b:
            res["config_identical"] = False
            ma, mb = a.get("model_settings", {}), b.get("model_settings", {})
            d = {k: [ma.get(k, "<absent>"), mb.get(k, "<absent>")] for k in set(ma) | set(mb) if ma.get(k, "<absent>") != mb.get(k, "<absent>")}
            if d: res["model_settings_diff"][str(rel)] = d
            for k in set(a) | set(b):
                if k != "model_settings" and a.get(k) != b.get(k): res["other_config_diff"].append(f"{rel}:{k}")
    names = {p.relative_to(new) for p in new.rglob("*") if p.is_file() and p.name not in skip}
    names |= {p.relative_to(old) for p in old.rglob("*") if p.is_file() and p.name not in skip}
    for rel in sorted(names):
        a, b = new / rel, old / rel
        res["files_compared"] += 1
        if not (a.exists() and b.exists()) or (a.stat().st_ino != b.stat().st_ino and md5(a) != md5(b)):
            res["files_diff"].append(str(rel))
    res["identical"] = res["config_identical"] and not res["files_diff"]
    return res

def main():
    p = common.parser(__doc__)
    p.add_argument("--table", type=Path, default=Path(__file__).resolve().parent / "variants.json")
    p.add_argument("--only", default="", help="comma-separated variant names (default: all)")
    p.add_argument("--out", type=Path, default=None, help="build into OUT/{reg,smoke}/ instead of the case dirs")
    p.add_argument("--verify-against", action="store_true", help="compare built variants with the existing ones in --regcases/--cases")
    p.add_argument("--json", type=Path, default=None, help="write the verification result here")
    a = common.finish(p.parse_args())
    T = json.loads(a.table.read_text()); only = set(filter(None, a.only.split(",")))
    reg_out = (a.out / "reg") if a.out else a.regcases
    smk_out = (a.out / "smoke") if a.out else a.cases
    reg_out.mkdir(parents=True, exist_ok=True); smk_out.mkdir(parents=True, exist_ok=True)
    built = []
    for name, spec in T["casesets"].items():
        if only and name not in only: continue
        build_caseset(a.regcases / spec["from"], reg_out / name, spec); built.append((reg_out / name, a.regcases / name))
    done = set()
    def smoke(name):                 # bases first (F1cdbg before F1cdbg_on)
        if name in done: return
        spec = T["smoke"][name]; base = spec["from"]
        if base in T["smoke"]: smoke(base); src = smk_out / base
        else: src = a.cases / base
        build_smoke(src, smk_out / name, spec); done.add(name); built.append((smk_out / name, a.cases / name))
    for name in T["smoke"]:
        if not only or name in only: smoke(name)
    if a.verify_against:
        assert a.out, "--verify-against needs --out (never compares in place)"
        out = {}
        for new, old in built:
            if not old.exists(): out[new.name] = {"missing_existing": str(old)}; continue
            if new.name in T["casesets"]:
                out[new.name] = {c.name: compare(c, old / c.name) for c in sorted(new.iterdir()) if (c / "spinup").is_dir()}
            else:
                out[new.name] = compare(new, old)
        ok = all((all(x["identical"] for x in r.values()) if "missing_existing" not in r and "identical" not in r else r.get("identical", False)) for r in out.values())
        for k, r in out.items():
            if "identical" in r or "missing_existing" in r:
                print(f"{k:12s}", "IDENTICAL" if r.get("identical") else "DIFF", "" if r.get("identical") else json.dumps(r)[:400])
            else:
                n = sum(x["identical"] for x in r.values())
                print(f"{k:12s}", f"{n}/{len(r)} cases IDENTICAL", "" if n == len(r) else json.dumps({c: x for c, x in r.items() if not x["identical"]})[:400])
        if a.json: a.json.write_text(json.dumps(out, indent=1) + "\n")
        print("ALL IDENTICAL" if ok else "DIFFERENCES FOUND"); sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
