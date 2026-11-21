# Restart state carry (bug fix), 2026-11-21

Branch `fix/restart-state-carry` (local only, not pushed), off `prod/deepc-default-off` e570a3e9. Worktree `~/tem-dev/restartfix/wt`. Diagnosis: `~/tem-dev/restartdiag/NOTES_restart_sp_mismatch.md`.

## Problem
A run restarted from a restart file differed from a chained run. The cause was two pieces of in-memory state that the restart file did not carry:
1. **BNF previous-year memory.** `bnf::YearMemory` in `Soil_Bgc::bnfmem` (prev_npp[12], prev_et[12], prev_pft_npp[NUM_PFT], have_prev) was not in RestartData. Every restarted stage therefore had **zero BNF for one year** (bnf_scheme >= 1), and the run stayed perturbed afterwards (noise level, up to 0.6% AVLN).
2. **Ponded water.** `magic_puddle` / `puddle_ice` were written only when `runon_energy_mode >= 1`, and never in MPI builds. The legacy ponding uses the puddle in every mode, so daily/SEB runs dropped it at the restart (S22 TVC: 3.18 mm, giving TSHLW up to 0.10 °C in tr year 1).

## Change (bug fix, no switch)
Chained runs already reload their own restart at every stage. Once the state round-trips, chained outputs are bit-identical to before, and only split runs change, to match chained.
- `RestartData`:
  - New members `bnf_prev_npp[12]`, `bnf_prev_et[12]`, `bnf_prev_pft_npp[NUM_PFT]`, `bnf_have_prev` (0/1) and `has_bnf_state`.
  - New `read_px_bnf_vars` / `write_px_bnf_vars`. NetCDF vars: `bnf_prev_npp`, `bnf_prev_et` (Y,X,prevtwelve), `bnf_prev_pft_npp` (Y,X,pft), `bnf_have_prev` (Y,X) int, all with _FillValue MISSING.
  - `create_empty_file` now defines these four plus `puddle_ice` / `magic_puddle` for every run. Under WITHMPI they are set to `NC_INDEPENDENT` access, as the other vars are.
  - `write_pixel_to_ncfile` always writes the puddle and BNF vars. The static `RestartData::write_puddle` and its `runon_energy_mode` gate in TEM.cpp are removed.
  - The writers use `nc_open_par` (MPI_COMM_SELF) under WITHMPI, like `write_px_prev_pft_vars`. Serial builds can also add the vars to a file made by an older binary (redef).
  - `read_pixel` calls `read_px_bnf_vars`. `has_bnf_state` is set only if all four vars are present.
- `Soil_Bgc`:
  - `bnf_memory_to_restart()` is called from `Cohort::set_restartdata_from_state`.
  - `bnf_memory_from_restart()` is called at the end of `Soil_Bgc::set_state_from_restartdata`.
  - **Old restarts** (no vars), a never-written pixel (fill values) or an invalid flag leave the in-memory memory unchanged. That is the legacy behaviour: a fresh process gets first-year BNF 0, and a chained process keeps its memory.
  - The `cur_*` accumulators are not stored; they are 0 at a year-end restart.
- The puddle reader was already unconditional (`Soil_Env.cpp`, `rdata.magic_puddle >= 0`).
- `include/Bnf.h` header comment updated.
- **MPI:** implemented, but **untested**. The VM has no MPI toolchain (no mpicxx). Before this change, puddle vars were never written under MPI.

## Tooling
Restart files now have 6 extra variables, so restart comparisons use **shared variables only**. A variable present in one file only is reported as info and is not counted as a difference.
- `tests/regression/identdir.py` (used by `ident.py`, the regression harness and the sync-check recipe): `cmp_nc(..., shared_only=True)` for `restart-*.nc`, and an `info` list (`restart-only:`). Pre-edit copy: `~/tem-dev/restartfix/identdir.py.bak_20261121`.
- `~/tem-dev/idcheck/idcmp.py` (identity check, not in git): restart one-sided vars go to `restart_vars_one_side`. Backup `idcmp.py.bak_20261121`.
- `~/tem-dev/sync_check.sh` does no file comparisons itself, so it is unchanged.

## Tests
- Unit, new suite `restart_carry` (`tests/unit/test_restart_carry.cpp`, 5 cases):
  - `create_empty_file` defines the six vars, and `reinitValue` clears them.
  - Exact per-pixel round trip of the BNF memory and puddle, including the have_prev=0 pixel and a never-written pixel (fill values).
  - Restored memory gives first-year BNF, month_parts and PFT shares equal to the chained memory (sum > 0, legacy 0), and it rolls on identically.
  - **Old-restart compatibility:** no vars means nothing is read, a fresh process keeps first-year 0 and a chained process keeps its memory. A fill or invalid flag is rejected.
  - The serial writers add the vars to an old-format file.
- Full suite: **482/482 test cases pass** (477 + 5 new). **8 failed assertions**, all the known ones (6 KNOWN BUG xfail cases plus the 2 richards_water_energy assertions, the same 8 as before).

Harness (fixed binary 972b5705, `~/tem-dev/restartfix/out`, summary `harness_summary.json`). Split legs are compared with the matching chained files (every output of the stages the leg ran, plus the restarts it wrote; restarts on shared vars):

| case | split | files compared | differing after fix | differing before fix (e570a3e9) |
|---|---|---|---|---|
| tile1 cell, monthly, bnf_scheme 2 (pr10 eq40 sp20 tr40) | eq->sp | 72 | **0** | 11 (of 14, lean spec) |
| tile1 cell, monthly, bnf_scheme 2 | sp->tr | 71 | **0** | 6 (of 13, lean spec) |
| S22 TVC, daily SEB, bnf_scheme 2 (eq5 sp5 tr5 / sp5 tr20) | eq->sp | 64 | **0** | 15 (of 64) |
| S22 TVC, daily SEB, bnf_scheme 2 | sp->tr | 63 | **0** | 14 (of 63) |

Straight-through, fixed binary vs production 29f26b5d (binA e45cf6ef), all files:
- tile1 cell chained: 76 files, **0 differing**.
- S22 TVC chained sp5 tr20, started from the production restart-eq (an old-format restart without the new vars): 66 files, **0 differing**.
- The only extra is the six new variables in the fixed binary's restart files (info).

Tile1 first transient year NINPUT (= BNF): chained 0.05167 gN/m2. Split before the fix: 0. Split after the fix: 0.05167, bit-identical (figure `restart_state_carry.png`).

## Cost
Runs: 16, **0.026 core-h** (ledger `~/tem-dev/restartfix/LEDGER.txt`, cap 0.35).

Rebuilds are listed separately: dvm 55 s plus unittest builds 104 s (a compile error, fixed) and 91 s, all at -j4. That is 250 s wall, at most about 0.28 core-h.

## Notes
- Existing restarts stay valid. New restarts are larger by 2×12 + NUM_PFT + 1 doubles/ints per cell plus 2 doubles.
- Production runs split at stage boundaries will now get the BNF year (tile: +0.024 to 0.051 gN/m2 per restart boundary) and keep the puddle, which removes the split-vs-chained noise.
