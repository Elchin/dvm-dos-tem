# Stage 4: Richards soil hydrology - step-halving fix (`richards_dt_fix`) and water/energy tests

[Back to the stage index](README.md) · Design doc: [§15f](../../tools/enthalpy_bench/enthalpy_solver_design.md)

## Purpose and physics summary
Bug R1 in `Richards::update`: `dt_dz = dtsub/dz` was computed once per sub-step, before the trial loop. When a trial is rejected and
`dtsub` is halved, `computeLHS`/`computeRHS` and the error estimate kept the old (longer) `dt_dz`, while `dtdone` advanced only by the
halved step. Each accepted retry therefore applied the fluxes for longer than it booked, so the day's net flux was applied for more
than 86,400 s and water was created. The fix recomputes `dt_dz` before every retry (as CTSM does). The stage also adds the analytic
Richards / water-balance / combined freeze-up tests that exposed R1, and defines the `enthalpy/ssf-rdt` branch (default-4 sub-steps +
`richards_dt_fix`), which is the Stage E baseline code.

## Switches
| key | default | values | meaning |
|---|---|---|---|
| `richards_dt_fix` | **0** (code and production configs) | 0 / 1, else FATAL | 1 = recompute `dt_dz` before each Richards retry; applies to both thermal solvers |

Elchin, 2026-11-10: "Yes: richards_dt_fix on for new-solver work and Stage E baseline, production stays off until recalibration".
The code default stays 0; new-solver case sets set 1 explicitly (`tools/enthalpy_bench/stagec/mk_variant.py` writes 1 when
`thermal_solver = 1` unless the caller sets the key; `extra.sh` and `stageD_vm.sh` overlays include it). `config/config.js` only carries a
commented new-solver example. Diagnostic counters (Richards calls, rejections, days with a rejection) are printed per stage on stderr as
`Richards substep summary: ...` and never affect results.

## How to run
```jsonc
"model_settings": { "thermal_solver": 1, "richards_dt_fix": 1 }   // production legacy: leave richards_dt_fix unset (0)
```

## Tests
| file | suite(s) | covers |
|---|---|---|
| `tests/unit/test_richards_water_energy.cpp` | `richards_analytic`, `water_validity`, `water_energy_combined` (10 cases) | no-flux hydrostatic profile, steady downward flux vs ODE, flux-BC infiltration bounds, ET/uptake/drainage balance; `checkSoilLiquidWaterValidity` conservation and the sub-step-fix guard; combined Richards + enthalpy freeze-up (sfcc 1 and 0); R1 and frozen-lens cases as `may_fail` |
| `tests/unit/test_richards_dtfix.cpp` | `richards_dt_fix` (4 cases) | defaults; infiltration column with rejections; rejection on vs off; through Soil_Env for `thermal_solver` 0 and 1 |

Run: `tests/unit/dvmdostem_unittest -ts=richards_dt_fix` (or `-ts=richards_analytic`, ...). The Happy Valley / R&O 2000 benchmark and its
Python driver from the same test campaign are on the local branch `tests/richards-water-energy` only, not on this branch.

## Key results
- Unit tests (richards_tests / richards_dtfix reports): infiltration column at 20 mm/d with default tolerances (11 rejections): **3.331 mm/d
  created with `dt_fix 0`, 2.8e-14 mm with 1**. Through Soil_Env: worst 5.7 mm/d (ts 0) and 3.3 mm/d (ts 1) with 0, 2.8e-14 mm/d with 1.
- Combined freeze-up test: worst daily energy 1.1e-8 (sfcc 1) / 6.2e-11 (sfcc 0) W/m2; water drift < 2e-13 kg.
- Production impact (`prod_t125_s23`, 8 cases, 2010-2024): retries are common (2,000-31,000 rejected trial sub-steps per run, on 2-29 % of
  Richards days). The bug creates water: -0.004...-0.44 % P (legacy) and -0.10...-0.53 % P (new solver). With the fix every run is within
  -0.001...+0.06 % P. Legacy impact is small (obs_regression 21/21 unchanged); new-solver ALD changes -0.04...+0.07 m.
- Merge to `enthalpy/ssf-rdt` (§15f): 322/322 unit test cases; bit-identity on `prod_t125_s23` 8/8 for (a) legacy with the key absent,
  (b) new solver default-4 with the key absent vs `f8617974`, (c) new solver s4 + fix vs the fix-branch runs. Closure on (c): full-day energy
  <= 4.6e-8 W/m2, water -0.001...+0.059 % P.

## Known open issues
- A small positive water residual (up to +0.06 % P) remains with the fix on; it is outside Richards and not investigated.
- Frozen lens inside fstsoill..drainl: the tridiagonal system is truncated (documented `may_fail`; not reachable from Soil_Env today).
- Turning `richards_dt_fix` on in production waits for the post-merge recalibration.

## Commits in this stage
| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `582cfc3b` | 2026-11-10 | Tests only: Richards analytic tests, checkSoilLiquidWaterValidity water conservation + sub-step-fix guard, combined Richards + enthalpy freeze-up closure (sfcc 0/1); bug R1 and frozen-lens cases documented as may_fail. | `tests/unit/test_richards_water_energy.cpp` | §15f |
| 2 | `8ffc55ed` | 2026-11-10 | Adds model_settings.richards_dt_fix (default 0 = legacy, bit-identical): recompute dt_dz before each Richards retry after a rejected sub-step is halved (bug R1); Richards sub-step counters on stderr; unit tests test_richards_dtfix.cpp. | `include/ModelData.h`, `include/Richards.h`, `src/Cohort.cpp`, `src/ModelData.cpp`, `src/Richards.cpp`, `src/Runner.cpp`, `tests/unit/test_richards_dtfix.cpp` | §15f |
| 3 | `05ef374c` | 2026-11-10 | Tooling/config comments: new-solver case-set tools (mk_variant.py, extra.sh, stageD_vm.sh) write richards_dt_fix=1; config.js only documents a commented new-solver example. | `config/config.js`, `tools/enthalpy_bench/stagec/extra.sh`, `tools/enthalpy_bench/stagec/mk_variant.py`, `tools/enthalpy_bench/stagec/stageD_vm.sh` | §15f |
| 4 | `5c08f890` | 2026-11-10 | Design doc §15f: enthalpy/ssf-rdt combines default-4 sub-steps with richards_dt_fix (Stage E baseline definition); docs only. | `tools/enthalpy_bench/enthalpy_solver_design.md` | §15f |
