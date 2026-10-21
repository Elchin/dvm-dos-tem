#!/usr/bin/env python3
"""Summarise a run_regression.sh output directory -> report.json / report.txt, exit 0 if PASS.
usage: summarize.py OUTDIR --expect-identical 0|1 --repeat 0|1 CASE [CASE ...]"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_runs, closure  # noqa: E402

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("cases", nargs="+")
    ap.add_argument("--expect-identical", type=int, default=1); ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--wtol", type=float, default=1.0); ap.add_argument("--ctol", type=float, default=0.5); ap.add_argument("--ntol", type=float, default=1.0)
    a = ap.parse_args(); out = Path(a.out); rows = []; fails = []
    for c in a.cases:
        r = dict(case=c); A, B, R = out / "A" / c, out / "B" / c, out / "R" / c
        try:
            ea = closure.evaluate(A, a.wtol, a.ctol, a.ntol); eb = closure.evaluate(B, a.wtol, a.ctol, a.ntol)
        except Exception as e:
            fails.append(f"{c}: cannot evaluate ({e})"); rows.append(r); continue
        r["A"] = ea; r["B"] = eb
        if ea["status"] != "OK" or eb["status"] != "OK": fails.append(f"{c}: run status A={ea['status']} B={eb['status']}")
        if not eb["closure_pass"]: fails.append(f"{c}: B closure: {'; '.join(eb['closure_fail'])}")
        if ea["status"] == "OK" and eb["status"] == "OK":   # closure must not get worse than the reference
            if abs(eb["W_resid_pctP"]) > abs(ea["W_resid_pctP"]) + 0.1: fails.append(f"{c}: water closure worse than A ({eb['W_resid_pctP']:.3f} vs {ea['W_resid_pctP']:.3f} % P)")
            if eb.get("N_resid_mean_gm2yr") is not None and ea.get("N_resid_mean_gm2yr") is not None and \
               abs(eb["N_resid_mean_gm2yr"]) > abs(ea["N_resid_mean_gm2yr"]) + 0.02: fails.append(f"{c}: N closure worse than A ({eb['N_resid_mean_gm2yr']:.3f} vs {ea['N_resid_mean_gm2yr']:.3f} gN/m2/yr)")
        for s in ("spin", "main"):
            if eb.get(f"{s}_safetynet", 0) > ea.get(f"{s}_safetynet", 0) or eb.get(f"{s}_veg_giveups", 0) > ea.get(f"{s}_veg_giveups", 0):
                fails.append(f"{c}: B has more integrator safety-net/give-up events in {s} than A")
        if ea["status"] == "OK" and eb["status"] == "OK":
            ab = compare_runs.cmp(A, B); r["A_vs_B"] = ab; r["A_vs_B_identical"] = all(s["ndiff"] == 0 for s in ab)
            if a.expect_identical and not r["A_vs_B_identical"]: fails.append(f"{c}: A vs B not bit-identical: " + json.dumps([[s["stage"], s["ndiff"], s["diffs"][:3]] for s in ab]))
        if a.repeat:
            if (R / "meta.json").exists() and json.loads((R / "meta.json").read_text())["status"] == eb["status"] and (B / "spin000").exists():
                br = compare_runs.cmp(B, R); r["B_vs_R"] = br; r["deterministic"] = all(s["ndiff"] == 0 for s in br)
            else:
                r["deterministic"] = False
            if not r["deterministic"]: fails.append(f"{c}: B is not deterministic (B vs repeat differ)")
        rows.append(r)
    rep = dict(binaries=(out / "binaries.txt").read_text().splitlines(), expect_identical=bool(a.expect_identical),
               repeat=bool(a.repeat), tolerances=dict(water_pctP=a.wtol, C_max_annual=a.ctol, N_mean=a.ntol),
               cases=rows, failures=fails, result="PASS" if not fails else "FAIL")
    (out / "report.json").write_text(json.dumps(rep, indent=1, default=float) + "\n")
    L = [*rep["binaries"], f"expect A==B bit-identical: {rep['expect_identical']}; determinism check: {rep['repeat']}", "",
         f"{'case':24s} {'A':>4s} {'B':>4s} {'A==B':>5s} {'B==R':>5s} {'W%P(B)':>8s} {'Cmax':>8s} {'Nmean':>8s} {'Nmax':>8s} {'NLOST':>7s} {'sn A/B':>7s}"]
    for r in rows:
        if "B" not in r: L.append(f"{r['case']:24s} ERROR"); continue
        b = r["B"]; f = lambda k: (f"{b[k]:8.4f}" if isinstance(b.get(k), float) else f"{'-':>8s}")
        sn = f"{r['A'].get('main_safetynet', '-')}/{b.get('main_safetynet', '-')}"
        L.append(f"{r['case']:24s} {r['A']['status'][:4]:>4s} {b['status'][:4]:>4s} {str(r.get('A_vs_B_identical', '-'))[:5]:>5s} "
                 f"{str(r.get('deterministic', '-'))[:5]:>5s} {f('W_resid_pctP')} {f('C_resid_max_abs_annual')} {f('N_resid_mean_gm2yr')} {f('N_resid_max_abs_annual')} "
                 f"{(b.get('NLOST_mean_gm2yr') or 0):7.4f} {sn:>7s}")
    L += ["", *("FAIL: " + x for x in fails), f"RESULT: {rep['result']}"]
    (out / "report.txt").write_text("\n".join(L) + "\n"); print("\n".join(L))
    return 0 if not fails else 1

if __name__ == "__main__":
    sys.exit(main())
