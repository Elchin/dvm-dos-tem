#!/usr/bin/env python3
"""Materialise a self-contained regression case set (once per machine).

usage: prepare_caseset.py OUTDIR [--set calib8] [--cases C1,C2,...] [--ms key=value ...]
  --ms overrides/adds model_settings (e.g. --ms nfactor_sameday=1 --ms vp_hpa_fix=1) to test switches.

calib8 = the 8 calibrated daily-mode cases {HPC,TVC}_{Warming,NoWarming}_{nfeed,nonfeed}
exactly as the v5i experiments (/workspace/v5i/scripts/exp5i.py) set them up:
  cases      /workspace/v5e_calib/cases_v5b/<case>/{spinup,main}/{config,input,parameters}  (read only)
  TVC params /workspace/tvc_hpc_stageA_kit/results_box/v5h/stage1_tvc/data/candidate_lhs3_r003.json (read only)
  HPC params $HPC_CAND (default /workspace/v5i/data/hpc_k5006.json)
  restarts   /workspace/v5h/mineralC_hpc_org/restarts/{obs40_m27 (TVC), hpc_org40 (HPC)}/<case>.nc
It reuses the existing helpers (vrun_seb.mk/edit, exp5i.BASE) so the parameter edits, model_settings
and config rewriting are identical to the calibration runs. Extra yearly outputs DEADC, DWDC, NINPUT,
NLOST, DEADN, DWDN are switched on in the main stage for C/N closure (output-only).

Layout written:  OUTDIR/<case>/{case.json, restart0.nc, parameters/, spinup/, main/}
(config paths inside spinup/ and main/ point at OUTDIR; run_case.py rewrites them per run).
Nothing under the source directories is modified."""
import argparse, json, os, shutil, sys, copy
from pathlib import Path
os.environ.setdefault("HPC_CAND", "/workspace/v5i/data/hpc_k5006.json")
os.environ.setdefault("VRUN_EXTRA_Y", "DEADN,DWDN")
SCRIPTS = [os.environ.get("VRUN_SCRIPTS", "/workspace/v5h/seb_explore/scripts"), os.environ.get("EXP5I_SCRIPTS", "/workspace/v5i/scripts")]
for s in SCRIPTS: sys.path.insert(0, s)
import vrun_seb as V      # noqa: E402  (read-only use)
import exp5i as X          # noqa: E402  (read-only use: BASE params / model_settings / restarts)

CASESETS = {"calib8": X.CASES["all8"]}

def prepare(out, case, ms_over=None):
    site = case[:3]; base = X.BASE[site]
    d = out / case
    if d.exists(): shutil.rmtree(d)
    d.mkdir(parents=True)
    params = copy.deepcopy(base["params"]); ms = dict(base["ms"]); ms.update(ms_over or {})
    pdir = d / "parameters"; shutil.copytree(V.T / case / "spinup/parameters", pdir)
    for f, kv in params.items():
        kv2 = {k.split(":", 1)[-1]: v for k, v in kv.items() if ":" not in k or k.startswith(site + ":")}
        if kv2: V.edit(pdir / f, V.CMT[site], kv2)
    rst = Path(base["rst"].format(case=case)); shutil.copy2(rst, d / "restart0.nc")
    V.FORCING.clear()
    for stage, extra in (("spinup", False), ("main", True)):
        V.mk(d / stage, case, stage, d / "restart0.nc", pdir, ms, extra)
        (d / stage / "output" / "restart-tr.nc").unlink()     # supplied per run
    meta = dict(case=case, site=site, model_settings=ms, params=params, restart0_src=str(rst),
                stage_root={s: str(d / s) for s in ("spinup", "main")}, n_spin=20, n_main=75,
                hpc_cand=os.environ["HPC_CAND"], extra_yearly=V.EXTRA_Y)
    (d / "case.json").write_text(json.dumps(meta, indent=1) + "\n")
    return case

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--set", default="calib8"); ap.add_argument("--cases"); ap.add_argument("--ms", action="append", default=[])
    a = ap.parse_args(); out = Path(a.out).resolve(); out.mkdir(parents=True, exist_ok=True)
    cases = a.cases.split(",") if a.cases else CASESETS[a.set]
    msov = {}
    for kv in a.ms:
        k, v = kv.split("=", 1); msov[k] = json.loads(v)
    for c in cases: print("prepared", prepare(out, c, msov), flush=True)
    (out / "caseset.json").write_text(json.dumps(dict(set=a.set, cases=cases, ms_override=msov), indent=1) + "\n")
