# Stage 2: energy-closure diagnostic for the legacy thermal step (phases 2b-2d)

[Back to the stage index](README.md)

## Purpose
Measure whether the legacy Stefan-front + Crank-Nicolson thermal step conserves energy, and gate it in the regression harness. Elchin
asked for this before the wetland merge (2026-11-09). The result motivated the enthalpy solver (Stage 3). Output only: no physics change.

## What is measured
For each daily thermal step (`Soil_Env::updateDailyGroundT`):
- column enthalpy H = sum over snow + soil + rock of [c_solid*T + liq*c_liq*T + ice*(c_ice*T - L_f)] (J/m2, reference liquid water at 0 C)
  before and after the step (`Soil_Env::layerEnthalpy`);
- heat through the solver's virtual top/bottom boundary nodes (Gtop, Gbot), with the solver's own conductances;
- residual res = dH - (Gtop + Gbot) - snowmk, written daily to `output/energy_daily.csv`. `snowmk` (from `a6c718d4`) is the snow-layer
  enthalpy change of the `meltingSnowLayer` melt marker, whose energy belongs to Snow_Env.

## Switches
| key | default | values | meaning |
|---|---|---|---|
| `energy_diag` | 0 | 0 / 1, else FATAL | 1 writes `energy_daily.csv`; changes no state (unit-tested). Harness test copies set 1 (`prepare_caseset.py` HARNESS_DEFAULTS, override `--ms energy_diag=0`); repo `config.js` stays 0 |

Harness option: `closure.py` / `summarize.py --etol` (run-mean FAIL threshold, **3 W/m2** since `1dca209e`, was 5).

## Tests
- Unit: `tests/unit/test_energy_diag.cpp`, suite `energy_diag` (6 cases on this branch): parse/default/FATAL; layer enthalpy including L_f;
  Gtop + Gbot equals the solver's discrete budget sum 0.5*(cap_i + cap_i-1)*dT_i to 1e-6; on/off identical; single-layer shortcut booking.
  Run: `tests/unit/dvmdostem_unittest -ts=energy_diag`.
- Harness logic: `python3 tests/regression/test_closure.py` (synthetic energy test).
- Regression: any `run_regression.sh` run on a prepared case set reports E(B), E(A), Ef(B), Erms(B).

## Key results (merge plan §8.1-8.3)
- With `energy_diag` 0, diag builds are bit-identical to the frozen binaries; with 1, outputs are bit-identical to the frozen v5q baseline
  (3952 / 3976 variables identical).
- Phase 2b (`517ab675`): the legacy scheme creates about **1.3-2.7 W/m2** on average, almost all on freeze/thaw-front days
  (prod: front days +4.5...+6.0 W/m2, no-front days -0.68...-0.77 W/m2).
- Phase 2d (`a6c718d4`, corrected diagnostic): run-mean residuals calib8 dc0 0.16-1.65, dc1 0.53-1.68; prod dc0 1.93-2.45,
  dc1 1.88-2.46 W/m2. All PASS at 3 W/m2. Unit tests 281/281 (6 legacy should_fail assertions unchanged).
- Legacy attribution (design doc §1): Stefan latent heat exceeds conducted heat by +2.0...+2.9 W/m2 (HPC) and +4.5...+5.3 W/m2 (TVC)
  on front days; the staggered half-capacity stencil gives about -0.65...-0.86 W/m2 on no-front days.

## Known open issues
- The budget covers the thermal step only; heat carried by water, snowfall/melt mass and layer restructuring are reported, not closed.
  The enthalpy solver (Stage 3) extends this to a full-day budget.
- The energy-conserving front prototype (`proto/front-energy`, adf4be54) was stopped and is **not** on this branch.

## Commits in this stage
| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `517ab675` | 2026-11-09 | Adds model_settings.energy_diag (default 0, output only): daily column enthalpy budget of the legacy thermal step to energy_daily.csv; harness energy check in closure.py/summarize.py; unit tests test_energy_diag.cpp. | 12 files: `tests/` (5), `src/` (4), `include/` (3) | merge plan §8.1 (Phase 2b); design doc §1 |
| 2 | `1dca209e` | 2026-11-09 | Regression harness: energy FAIL threshold 5 -> 3 W/m2 run-mean; prepare_caseset.py writes energy_diag 1 into harness test copies only (repo config.js stays 0). | `tests/README.md`, `tests/regression/closure.py`, `tests/regression/prepare_caseset.py`, `tests/regression/summarize.py` | merge plan §8.2 (Phase 2c) |
| 3 | `a6c718d4` | 2026-11-09 | energy_diag fixes (output only): books boundary heat of the single-layer solver shortcut; excludes the snow-melt marker (new snowmk column). | `include/Soil_Env.h`, `src/Cohort.cpp`, `src/Soil_Env.cpp`, `src/TemperatureUpdator.cpp`, `tests/README.md`, `tests/regression/closure.py`, `tests/unit/test_energy_diag.cpp` | merge plan §8.3 (Phase 2d) |
