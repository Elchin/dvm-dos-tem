"""Shared path options for the Reburn / fire-balance helper scripts.

Every script takes the same options (all optional):
  --cases DIR   smoke case directory          (default $TEM_SMOKE_CASES or $TEM_WORK/wet-work/rbsmoke/cs)
  --runs DIR    smoke run directory           (default $TEM_SMOKE_RUNS  or $TEM_WORK/wet-work/rbsmoke/runs)
  --regcases DIR regression case-set root, holding leg/ and s4_f1/
                                              (default $TEM_REG_CASES   or $TEM_WORK/srd-work/cs)
  --figdir DIR  where figures and summary JSONs are written (default: this repo's docs/figures/<sub>)
  --work DIR    scratch for intermediate files (default <runs>/_analysis)
TEM_WORK defaults to ~/tem-dev (the calib-elchin layout, where the case data lives for now).
"""
import argparse, os
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]          # tests/regression/reburn/common.py -> repo root
TEM_WORK = Path(os.environ.get("TEM_WORK", str(Path.home() / "tem-dev"))).expanduser()


def parser(desc, figsub="fire_balance"):
    p = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cases", type=Path, default=Path(os.environ.get("TEM_SMOKE_CASES", TEM_WORK / "wet-work/rbsmoke/cs")))
    p.add_argument("--runs", type=Path, default=Path(os.environ.get("TEM_SMOKE_RUNS", TEM_WORK / "wet-work/rbsmoke/runs")))
    p.add_argument("--regcases", type=Path, default=Path(os.environ.get("TEM_REG_CASES", TEM_WORK / "srd-work/cs")))
    p.add_argument("--figdir", type=Path, default=REPO / "docs/figures" / figsub)
    p.add_argument("--work", type=Path, default=None)
    return p


def finish(a):
    """Resolve defaults that depend on other options; create output dirs."""
    for k in ("cases", "runs", "regcases", "figdir"):
        setattr(a, k, Path(getattr(a, k)).expanduser().resolve())
    a.work = (Path(a.work).expanduser().resolve() if a.work else a.runs / "_analysis")
    a.figdir.mkdir(parents=True, exist_ok=True)
    return a
