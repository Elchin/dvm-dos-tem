# Reburn / fire-balance regression helpers

Scripts behind `docs/WORK_NOTES_reburn.md` and `docs/WORK_NOTES_fire_balance.md`: case builders, run tools,
analysis, and `run_all.sh`, which reruns everything from the three reference commits.

## Case data (calib-elchin only for now)

The only external input is the set of 8 calibrated regression cases, `leg/`
(`{HPC,TVC}_{Warming,NoWarming}_{nfeed,nonfeed}`, each `case.json`, `restart0.nc`, `parameters/`, `spinup/`,
`main/` with inputs; 144 MB). It is **not in the repo** and lives on calib-elchin in `~/tem-dev/srd-work/cs/leg`.
Every other case is generated from it:

| cases | built by | from |
|---|---|---|
| `s4_f1/*` (8 regression cases, new thermal solver) | `mk_variants.py` (`variants.json`, `casesets`) | `leg/*` |
| `F0g`, `F1c`, `F1p` (FRI fire smoke cases) | `mk_smoke_cases.py` (+ `external/Fireweed/Input_Files`) | `leg/TVC_Warming_nfeed` |
| `E1` (explicit classic fire) | `mk_explicit_case.py` | `leg/TVC_Warming_nfeed` |
| `F1c_*`, `F1cdbg*`, `F1p_*`, `E1on*` | `mk_variants.py` (`variants.json`, `smoke`) | `F1c`, `F1p`, `E1` |

`variants.json` holds the `model_settings` overlay of every variant (and, for `F1cdbg`, the outputs switched to
yearly + monthly). `E1onL` is E1 with both fire fixes **and both N leak fixes** (the `L` there means leak).

## Paths

All scripts take `--cases DIR` (smoke cases), `--runs DIR` (smoke runs), `--regcases DIR` (root holding
`leg/`, `s4_f1/`), `--figdir DIR` (figures and summary JSONs; default this repo's `docs/figures/...`) and
`--work DIR` (intermediates). Defaults come from `TEM_SMOKE_CASES`, `TEM_SMOKE_RUNS`, `TEM_REG_CASES`, or
`TEM_WORK` (default `~/tem-dev`, the calib-elchin layout). See `common.py`. Nothing is written to `/tmp`.

## Files

| file | what |
|---|---|
| `run_all.sh` | build the three reference binaries, make the cases, run the identity and balance sets, check |
| `variants.json`, `mk_variants.py` | the variant table and its builder; `--out DIR --verify-against` rebuilds into DIR and compares with the existing variants |
| `mk_smoke_cases.py`, `mk_explicit_case.py` | the smoke bases F0g/F1c/F1p and E1 |
| `analyze_smoke.py` | fire events and C/N residuals, fire 0 / classic / process (`docs/figures/reburn/`) |
| `analyze_fire_balance.py` | fire C/N balance gaps, `fire_n_redep_fix`, `fire_c_bookkeep_fix` |
| `analyze_n_leaks.py` | N leak fixes `dwd_n_leak_fix`, `deadroot_n_route_fix` |
| `split_run.py`, `compare_split.py` | restart continuity: main stage split at year 35 vs continuous |
| `redep_overshoot.py`, `plot_restart_continuity.py` | spin-up → main re-deposition overshoot; restart figure |
| `compare_summaries.py` | compare fresh summary JSONs with the committed ones |
| `plot_test_summary.py` | `docs/figures/fire_balance/test_summary.png` from the committed summary JSONs |
| `common.py` | shared path options |

Shared tools one level up, in `tests/regression/`: `run_case.py` (one case: 20-yr spin-up + 75-yr main),
`launch4.sh` (a job list in parallel), `identdir.py` (two output directories, every NetCDF variable and restart
file) and `ident.py` (pairs of run directories, spin000 + main; JSON report).

## run_all.sh

```
tests/regression/reburn/run_all.sh --work /scratch/dir --conda ~/miniconda3/envs/dvmdostem-build   # everything
tests/regression/reburn/run_all.sh --work DIR --dry-run                                             # show the plan
tests/regression/reburn/run_all.sh --work DIR --bins DIR_WITH_dvmdostem_fb2_fb3_fb4 --sets cases,identity,check \
    --reg-sets leg --reg-cases TVC_Warming_nfeed --smoke F0g --fixes-on ""                          # small check
```

1. **build**: detached `git worktree`s of `19af750f` (fb2), `06c41688` (fb3) and `f94fceaf` (fb4) in
   `WORK/src`, `scripts/get_external.sh` (Fireweed `31051b0`, Burnup `6d56e29`; needs GitHub), `make`,
   binaries in `WORK/bin/dvmdostem_fb{2,3,4}`. `--bins DIR` skips this.
2. **cases**: `WORK/cases/` from `--regcases` (s4_f1 is hard-linked, so WORK must be on the same file system).
3. **identity**: 16 regression cases and E1/F1c/F1p/F0g/F1cdbg with fb2, fb3 and fb4, plus F1c_on and E1on with
   fb3 and fb4 (67 runs). **balance**: the remaining fire-balance, N-leak and overshoot runs (16 new runs; F1c,
   F0g, F1p_c reuse the fb2 identity runs) and 3 restart splits.
4. **check**: `ident.py` on fb2 vs fb3 and fb3 vs fb4; the analysis scripts write `WORK/summaries/`, and
   `compare_summaries.py` compares them with the committed JSONs.

**Expected results** (as committed): identity **88/88** stage folders identical (fb2 vs fb3: 42, fb3 vs fb4: 46;
`docs/figures/fire_balance/identity_summary.json`); `compare_summaries.py` **PASS** for
`fire_balance_summary.json`, `n_leak_summary.json`, `overshoot.json`, `split_summary.json` and
`reburn/smoke_summary.json` (the runs are deterministic, so the rounded numbers match exactly). The original
`F0g`/`F1c`/`F1p_c` smoke runs used the pre-fix Reburn binaries; with every new switch 0 the fb2 binary gives the
same outputs (the identity chain in the work notes), so the smoke summary is expected to match too.

**Cost**: about **1.5 core-h** on calib-elchin (3 builds ≈ 0.25, ≈ 86 run-equivalents of ≈ 1 min each ≈ 1.3) and
≈ 2 GB in WORK. `-P` sets the parallel runs and `make -j` (default 4).
