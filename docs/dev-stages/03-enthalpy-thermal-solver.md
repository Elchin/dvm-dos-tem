# Stage 3: enthalpy thermal solver (soil and snow heat, freezing curve, sub-steps, melt-out and snow-layer fixes)

[Back to the stage index](README.md) · Full design and measurements: [`enthalpy_solver_design.md`](../../tools/enthalpy_bench/enthalpy_solver_design.md)

## Purpose and physics summary
The legacy thermal step (Stefan fronts + Crank-Nicolson) does not conserve energy (Stage 2: about +1.3...+2.7 W/m2). This stage
replaces it, behind a default-off switch, with a fully implicit, globally convergent enthalpy solver for the snow + soil + rock column.
Energy closes to round-off (full-day budget <= 1e-6 W/m2 is a gate), there are no tracked fronts in the physics (fronts and ALD are
diagnosed from the profile), unfrozen water follows the Schaefer & Jafarov (2016) curve already used by the BGC code, and melt and
layer remaps are booked conservatively. This README covers all soil- and snow-thermal commits.

### A. Design (`53252d3c`, design doc §0-§19)
Implicit backward-Euler, cell-centred finite-volume enthalpy scheme solved with the nested Newton method of Casulli & Zanolli (NCZ),
following Tubini et al. (2021); a Katzenelson solver (Hötten et al. 2026) is the free-water cross-check. Prognostic layer enthalpy
H_i (J/m2), written back to `Layer::tem/liq/ice`, so the rest of the model keeps reading T, liq and ice. Elchin's decisions of
2026-11-09 (§0a): D1 production freezing curve `thermal_sfcc = 1` (Schaefer & Jafarov 2016); D2 harmonic face conductance;
D3 radiation/SEB melt becomes a heat source in the top snow cell and the "surface fixed at 0 C while melting" rule is dropped
(new solver only); D4 geothermal flux zero/value/map; D5 full-day energy balance <= 1e-6 W/m2 is a pass/fail gate; D6 core-hour
budget deferred; D7 Katzenelson cross-check in Stage B.

### B. Standalone library (`106ce93f`, `3fb63332`; design doc §15a)
`include/EnthalpySolver.h`, `src/EnthalpySolver.cpp`: curve, NCZ, Thomas solve, Katzenelson. Compiled into the model, never called.
Replay of 14,595 TVC/HPC column-days at -O2: NCZ SFCC 1 **27.9 µs/day** (3.42 solves/day) vs legacy 217 µs/day (734 solves/day);
max daily closure 3.4e-8 W/m2 (SFCC 1), 2.0e-11 (SFCC 0). Unit tests 290/290 (9 new). Benchmark tools: `tools/enthalpy_bench/` ([README](../../tools/enthalpy_bench/README.md)).

### C. Integration behind `thermal_solver` (`accddbfc`, `5d783561` via `9bfb8941`/`6826f0f3`, `d5b68b3b`, `e44c5433`, `457917e2`; §15b)
Daily enthalpy step in `src/Soil_Env_enthalpy.cpp`, dispatched from `Soil_Env::updateDailyGroundT`; Robin top boundary (cn 100);
SEB/radiation melt as snow heat source; heat of fully melted snow cells passed down; fronts, `frozen`, `frozenfrac`, ALD/ALC
reconstructed daily from the profile; switch-on/restart projection booked as `init`; `Ground::enthalpy_remap` (split/merge conserve H);
harness `closure.py --etol-enth` (default 1e-6 W/m2). `stefan_front_resist_fix` is ignored with a WARNING under `thermal_solver = 1`.
`5d783561` fixes a pre-existing one-past-end read in EnvData (bit-identical) so -O2 -Werror builds work.

### C2/C2b. Energy-conserving remap and restart label (`f8f518b9`, `15fcca85`; §15c)
Moved water carries its enthalpy (rain at max(T_air, 0 C), melt and ponded water at 0 C); restart attribute `thermal_solver`
(written only when 1; absent = legacy; a mismatch only warns). C2b shares direct-melt latent heat over the top 0.10 m of soil when there
are no snow cells.

### D. Site evaluation (`557462f8`; §15c D1-D3)
TVC/HPC vs observations (2010-2024), sub-step study, root cause of the integrator safety-net increase. Tools in
`tools/enthalpy_bench/staged/` (`stageD_eval.py`, obs in `obs_stageD.json`).

### Sub-step fix and default 4 sub-steps (`7c5e4759`, `f8617974`; §15d, §15e)
Fix A: with `thermal_solver = 1` a `frozen == 1` layer is inert in `checkSoilLiquidWaterValidity` (its freezing-curve water was being
dumped, whole column on fully frozen days). Fix B: the enthalpy tolerance is a per-day budget (tol/nsub per sub-step).
`thermal_substeps` also accepts 12 and 48. Then (Elchin, 2026-11-10, "Yes, make 4 sub-steps the default") an absent
`thermal_substeps` means 4 when `thermal_solver = 1`.

### Melt-out closure fix (`282287bc`; §15g)
Calibration-expert's r1 case (TVC_Warming_nfeed, CMT51 `snwdennew` 75, `nfactor(w)` 0.3, `snow_scheme` 0, Stage E baseline settings)
missed the daily gate on 1 Jun 2024: res = -57.36 J/m2 (6.6e-4 W/m2). Root cause (instrumented): a top snow cell melted out at
+1.52 C and its heat was handed to an **empty** snow cell (W = 0, dz = 0, zero heat capacity), where `equilibrate()` returned T = 0
and the heat vanished. Fix (`Soil_Env_enthalpy.cpp`, step 5, new solver only): the receiver is the next cell below that can hold heat
(Cs > 0 or liq + ice > 0). When the next cell has heat capacity the arithmetic is unchanged.

### Empty snow-layer cleanup (`60691f7a`, `d4f86733`; §15h, text reproduced here)
Root cause behind the melt-out miss. In `Snow_Env::updateDailyM`, when the thermal step melts a lower snow layer completely while the
layer above keeps ice, `meltSnowLayersAfterT()` leaves ice = liq = 0 with the old dz. `Ground::constructSnowLayers()` removes emptied
layers only from the top, `Ground::combineSnowLayers()` tests the pre-compaction dz (still > mindz), and `SnowLayer::updateThick()` sets
dz = 0 only afterwards, so the layer (W = 0, dz = 0) enters the next day's thermal step. In r1 this happened once in 75 years
(31 May 2024: the 29.66 mm pack melted to a 0.009 mm remnant over the emptied bottom layer).

Change (switch on): `Ground::removeEmptySnowLayers()` runs right after `updateSnowLayerPropertiesDaily()` and before
`updateSnowHorizon()`. On this branch the call also comes before `stratifySnowLayers()` (snow-class scheme). A snow layer is removed when |ice| < 1e-9 and |liq| < 1e-9 kg/m2 and |dz| <= 1e-9 m; any residual mass
goes to the snow layer above (else below, else `snow.extramass`), so water is conserved exactly. Removal relinks the layer list and
then runs `updateLayerIndex()` and `resortGroundLayers()`. Under `thermal_solver = 1` the call sits inside the snow-window enthalpy
booking. Counters `snow_empty_removed/_days/_mass`; with `energy_diag = 1` days with removals go to `output/snow_empty_daily.csv`.

Note: commit `d4f86733` says "design doc 15h", but only the §15g cross-reference and a ModelData.h comment landed in
`enthalpy_solver_design.md`; the §15h text is this subsection.

## Switches (all default-off / bit-identical; `thermal_solver = 0` takes the untouched legacy path)
| key | default | values | meaning |
|---|---|---|---|
| `thermal_solver` | **0** | 0 legacy Stefan + CN, 1 enthalpy NCZ | master switch |
| `thermal_sfcc` | 1 | 1 SJ2016 power law + 0 C ramp (production choice D1), 0 free water + ramp | non-default only with `thermal_solver = 1` |
| `thermal_substeps` | **4 when `thermal_solver = 1` and the key is absent**; must be absent or 1 with `thermal_solver = 0` | 1, 2, 4, 8, 12, 24, 48 | backward-Euler steps per day (approved default, Elchin 2026-11-10) |
| `thermal_face_mode` | 0 | 0 harmonic (D2), 1 Tubini max rule (fallback) | face conductance |
| `thermal_nl_solver` | 0 | 0 NCZ, 1 Katzenelson (SFCC 0 only; verification) | nonlinear solver |
| `geothermal_flux_mode` | 0 | 0 zero flux, 1 constant `geothermal_flux`, 2 map `IO.geothermal_flux_file` (NetCDF var `geothermal_flux`, W/m2) | 1/2 need `thermal_solver = 1` |
| `geothermal_flux` | 0.0 | W/m2 | used by mode 1 |
| `snow_empty_cleanup` | 0 | 0 / 1, else FATAL | 1 = remove emptied snow layers (ice, liq < 1e-9 kg/m2 and dz <= 1e-9 m) at the end of the daily snow update, both solvers (from `60691f7a`) |

Out-of-range values are FATAL, and non-default dependent keys with `thermal_solver = 0` are FATAL. Restart files written with
`thermal_solver = 1` carry a global attribute `thermal_solver = 1`; either solver reads either file (warning on mismatch).

## How to enable
```jsonc
"model_settings": {
  "thermal_solver": 1,      // enthalpy solver (default 0 = legacy)
  "thermal_sfcc": 1,        // SJ2016 freezing curve (default; 0 = free water)
  // "thermal_substeps": 4, // absent = 4 under thermal_solver 1
  "richards_dt_fix": 1,     // Stage 4; set explicitly for new-solver work
  "energy_diag": 1,
  "snow_empty_cleanup": 1  // Stage E baseline choice; default 0
}
```

## Tests
Unit tests (doctest): `make unittest` (or `./build-conda.sh -j8 unittest`), then e.g. `tests/unit/dvmdostem_unittest -ts=enthalpy_integration`.

| file | suite(s) | covers |
|---|---|---|
| `tests/unit/test_enthalpy_solver.cpp` | (no suite name; 9 cases) | Stage B library: SJ2016 curve = `sjLayerPhi`, Neumann analytic, conservation, large dt, Katzenelson vs NCZ |
| `tests/unit/test_enthalpy_integration.cpp` | `enthalpy_integration`, `enthalpy_remap_c2` (15 cases) | Stage C/C2: switch parsing, closure, remap, rock heat capacity, NaN, default-4 sub-steps |
| `tests/unit/test_enthalpy_substep_fix.cpp` | `enthalpy_substep_fix` (3 cases) | fix A whole-column freeze, legacy unchanged, fix B tolerance |
| `tests/unit/test_enthalpy_meltout.cpp` | `enthalpy_meltout` (4 cases) | r1 1 Jun 2024 column, cascade over two empty cells, empty cell between snow cells, unchanged path |
| `tests/unit/test_snow_empty_cleanup.cpp` | `snow_empty_cleanup` (8 cases) | switch parse/FATAL; middle/top/bottom/cascaded removal with list integrity and exact mass; residual mass; thin layers kept; `updateDailyM` path for both solvers; next-day closure |

Regression:
- `tests/regression/closure.py RUNDIR --etol-enth 1e-6`: daily solve and full-day residual gates, zero non-converged steps, non-finite = FAIL.
- `tests/regression/run_regression.sh -j 4 BIN_A BIN_B CASESET OUT` with `thermal_solver` 0 for the bit-identity gate (see [tests/README.md](../../tests/README.md)).
- Stage C2 verification driver: `tools/enthalpy_bench/stagec/c2_vm.sh <commit>`; case-set variants: `tools/enthalpy_bench/stagec/mk_variant.py`.
- `tests/regression/check_meltout_r1.sh BIN [CASESET] [RUNDIR]`: re-runs r1 and applies the 1e-6 W/m2 gate (default case set path is on calib-elchin).

## Key results
- **Bit-identity, `thermal_solver = 0`:** md5 / data identical to the legacy reference at every step (Stage C dc0/dc1 x calib8/prod; C2b; sub-step fix 24/24; default-4 on `dc0_calib8`).
- **Closure, `thermal_solver = 1`:** Stage C2b 32/32 PASS, max full-day residual 3.7e-8 W/m2 (SFCC 0: 1.1e-10), unattributed 0, remap fallbacks 0; water -0.42...+0.59 % P.
- **Runtime (Stage C, whole model):** 0.59x legacy at -O0 and 0.73x at -O2 (47.6 vs 64.8 core-h per century per 10,000 cells at -O2).
- **Physics vs legacy (Stage C, last 30 yr):** soils colder and ALD shallower (e.g. calib8 dc0 HPC ALD -0.26...-0.37 m). SFCC 0 gives nearly the same result as SFCC 1 (ALD within 0.04 m), so the change matches the removed legacy energy leak, not the freezing curve.
- **Stage D vs observations:** legacy is already too cold in winter (DJF -2.5...-6 C); the new solver roughly doubles the cold-season bias (DJF -7.5...-10 C) and worsens HPC ALD bias (-0.02...+0.10 -> -0.16...-0.31 m). Snow-off improves (HPC -17 -> -5 d). The cold-season error points at the snowpack (HPC peak SWE 55-59 % low).
- **Sub-steps (Stage D):** ALD within 0.05 m of the 24-step reference for every count; 4 sub-steps cost +17 % wall time and the -O2 model is still faster than legacy (21 vs 23 s).
- **Sub-step fix (§15d):** HPC_Warming_nfeed at 24 sub-steps went from ALD 2010-24 1.34 m (max 1.91) to 0.61 m; after the fix ALD converges across 1-48 sub-steps (0.61-0.66 m) and the worst residual falls with sub-steps (2.0e-8 W/m2 at 24).
- **Default 4 (§15e):** key absent = explicit 4 (bit-identical); explicit 1 = old 1-step result; legacy bit-identical.
- **Melt-out fix:** r1 worst day 6.6e-4 -> **5.6e-8 W/m2**, 1 Jun 2024 -57.36 -> 9.9e-8 J/m2; 3 of the 4 new unit tests fail without the fix (res -57.9, -419, -27 J/m2). Legacy 8/8 and new-solver s4_f1 8/8 TVC/HPC cases bit-identical (the zero-mass receiver never occurs there).
- **Empty-layer cleanup:** switch off bit-identical on all 8 TVC/HPC cases for both solvers and on r1. Switch on: no empty layer occurs in the 8 TVC/HPC cases or legacy r1 (bit-identical); r1 new solver: one removal (31 May 2024), outputs differ from 1 Jun 2024 only, max |dT| 2.6e-5 C, closure 5.6e-8 W/m2 on and off. Cleanup **without** the melt-out heat fix also closes r1 (6.6e-4 -> 5.6e-8 W/m2) and is bit-identical to cleanup + heat fix, so the empty layer was the only path to the lost heat; the §15g fix stays as a guard.
- **Full unit suite on this branch tip:** 346/346 test cases pass (snow_empty_fix report); the 8 failed assertions are the known xfail / may_fail cases.

## Known open issues
- **Calibration:** T125/S23 is not valid with `thermal_solver = 1` (soil T, ALD and thaw timing change). Recalibration is Calibration-expert's Stage E.
- **Cold winter bias** under the new solver (Stage D: DJF -7.5...-10 C at TVC/HPC; new solver passes 7/21 obs_regression gated checks in the snow-class report, TVC failures). Points at snowpack/insulation parameters.
- **Restart dependence (new solver only):** results depend on the restart split point, up to 1.1 C DJF T0.10 at Franklin Bluffs (snow-checks report); legacy is identical. Not investigated.
- **Integrator safety net:** more N-uptake give-ups in HPC nfeed cases (e.g. calib8 dc1 HPC 6/4 -> 64/60 after C2b) because the thaw is shallower and colder; `soil_giveup_retry = 1` removes them. A calibration-time decision, not a solver change.
- **Stage D scope not done:** `thermal_k_picard` (designed in §7/§12, not implemented as a key), multi-front-year analysis, HPC MPI tile timing.
- `snow_empty_cleanup` stays default 0. For the Stage E baseline Elchin agreed (2026-11-10) to run it on together with the melt-out fix; Calibration-expert switches only after their verification checks pass.

## Commits in this stage
| # | commit | date | summary | files touched | design doc / ref |
|---|---|---|---|---|---|
| 1 | `106ce93f` | 2026-11-09 | Stage B: standalone enthalpy solver library (NCZ nested Newton + Katzenelson cross-check), 9 unit tests and the compute-cost benchmark tools; compiled into the model but never called. | 10 files: `tools/enthalpy_bench/` (5), `Makefile`, `SConstruct`, `include/EnthalpySolver.h`, `src/EnthalpySolver.cpp`, `tests/unit/test_enthalpy_solver.cpp` | §15a |
| 2 | `3fb63332` | 2026-11-09 | Benchmark build only: -O2 bench build demotes aggressive-loop-optimizations to a warning (pre-existing EnvData.cpp out-of-bounds read). | `tools/enthalpy_bench/run_bench.sh` | §15a note (ii) |
| 3 | `53252d3c` | 2026-11-09 | Adds the enthalpy solver design document (Stage A design, Elchin's decisions D1-D7 of 2026-11-09, Stage B measured cost). | `tools/enthalpy_bench/enthalpy_solver_design.md` | §0-§19, §0a, §15a |
| 4 | `accddbfc` | 2026-11-09 | Stage C: integrates the enthalpy solver behind model_settings.thermal_solver (default 0 = legacy, bit-identical) with thermal_sfcc, thermal_face_mode, thermal_nl_solver, thermal_substeps, geothermal_flux_mode; shared SJ2016 curve (SJCurve); fronts/ALD rebuilt from the profile; closure.py --etol-enth gate; unit tests test_enthalpy_integration.cpp. | 21 files: `src/` (9), `include/` (6), `tests/` (4), `Makefile`, `SConstruct` | §15b |
| 5 | `5d783561` | 2026-11-09 | EnvData: fixes one-past-end reads in getTempAtDepthFromArray/getVWCAtDepthFromArray (bit-identical); the model now builds at -O2 -Werror. Committed on wetland_deep_C, merged into the solver branches. | `src/EnvData.cpp` | §15a note (ii), §15b |
| 6 | `9bfb8941` | 2026-11-09 | Merge of wetland_deep_C (brings 5d783561) into enthalpy/stageB. | merge, vs first parent: `src/EnvData.cpp` | §15b |
| 7 | `6826f0f3` | 2026-11-09 | Merge of enthalpy/stageB (with 5d783561) into enthalpy/stageC. | merge, vs first parent: `src/EnvData.cpp` | §15b |
| 8 | `d5b68b3b` | 2026-11-09 | Stage C: FATAL layer dump when the column enthalpy is non-finite; closure gate fails on NaN. | `include/Soil_Env.h`, `src/Soil_Env_enthalpy.cpp`, `tests/regression/closure.py` | §15b |
| 9 | `e44c5433` | 2026-11-09 | Stage C: rock heat capacity taken from ParentLayer (Layer::vhcsolid is shadowed and MISSING for rock); rock unit test; NaN closure test. | `src/Ground.cpp`, `src/Soil_Env_enthalpy.cpp`, `tests/regression/test_closure.py`, `tests/unit/test_enthalpy_integration.cpp` | §15b |
| 10 | `457917e2` | 2026-11-10 | Stage C docs/tools: design doc §15b (results, deviations), verification drivers and comparison tools under tools/enthalpy_bench/stagec (model sources unchanged). | `tools/enthalpy_bench/enthalpy_solver_design.md`, `tools/enthalpy_bench/stagec/cmp_ts.py`, `tools/enthalpy_bench/stagec/compare_ts.py`, `tools/enthalpy_bench/stagec/esc_vm.sh`, `tools/enthalpy_bench/stagec/extra.sh`, `tools/enthalpy_bench/stagec/mk_variant.py` | §15b |
| 11 | `f8f518b9` | 2026-11-10 | Stage C2: energy-conserving remap of the snow and soil-water windows (moved water carries its enthalpy), layer-window water check, thermal_solver attribute in restart files; unit tests. | 9 files: `src/` (4), `include/` (2), `tests/` (2), `tools/enthalpy_bench/stagec/c2_vm.sh` | §15c (C2) |
| 12 | `15fcca85` | 2026-11-10 | Stage C2b: with no snow cells, direct-melt latent heat is shared by the top 0.10 m of soil by heat capacity instead of the first cell (a thin dry moss cell had reached -2850 C); unit test. | `src/Soil_Env_enthalpy.cpp`, `tests/unit/test_enthalpy_integration.cpp` | §15c (C2b) |
| 13 | `557462f8` | 2026-11-10 | Stage D: TVC/HPC evaluation vs observations, thermal_substeps study, integrator safety-net root cause; design doc §15c (tools and docs only). | 10 files: `tools/enthalpy_bench/` (10) | §15c (D1-D3) |
| 14 | `7c5e4759` | 2026-11-10 | Sub-step fix: thermal_solver 1 keeps the freezing-curve water of frozen cells in checkSoilLiquidWaterValidity (fix A) and uses a per-day enthalpy tolerance (fix B); thermal_substeps also accepts 12/48; unit tests test_enthalpy_substep_fix.cpp. | 14 files: `tools/enthalpy_bench/` (9), `src/` (3), `include/ModelData.h`, `tests/unit/test_enthalpy_substep_fix.cpp` | §15d |
| 15 | `f8617974` | 2026-11-10 | Default thermal_substeps = 4 when thermal_solver = 1 and the key is absent (Elchin, 2026-11-10); legacy unchanged; unit test; config.js comment. | `config/config.js`, `include/ModelData.h`, `src/ModelData.cpp`, `tests/unit/test_enthalpy_integration.cpp`, `tools/enthalpy_bench/enthalpy_solver_design.md` | §15e |
| 16 | `282287bc` | 2026-11-10 | Melt-out closure fix (thermal_solver 1 only): heat of a melted-out snow cell skips zero-mass snow cells (was lost: r1, 1 Jun 2024, -57.36 J/m2); unit tests test_enthalpy_meltout.cpp; regression check_meltout_r1.sh; design doc §15g. | `src/Soil_Env_enthalpy.cpp`, `tests/regression/check_meltout_r1.sh`, `tests/unit/test_enthalpy_meltout.cpp`, `tools/enthalpy_bench/enthalpy_solver_design.md` | §15g |
| 17 | `60691f7a` | 2026-11-10 | Adds model_settings.snow_empty_cleanup (default 0 = legacy, bit-identical; both solvers): removes emptied snow layers at the end of Snow_Env::updateDailyM (root cause of the r1 melt-out miss); diag snow_empty_daily.csv; unit tests test_snow_empty_cleanup.cpp. | `config/config.js`, `include/Ground.h`, `include/ModelData.h`, `src/Cohort.cpp`, `src/Ground.cpp`, `src/ModelData.cpp`, `src/Snow_Env.cpp`, `tests/unit/test_snow_empty_cleanup.cpp` | §15h (text in this README) |
| 18 | `d4f86733` | 2026-11-10 | Docs: design-doc cross-reference to §15h in §15g and a ModelData.h comment matching the 1e-9 thresholds. Note: the §15h section text itself did not land in the design doc; it is reproduced in this README. | `include/ModelData.h`, `tools/enthalpy_bench/enthalpy_solver_design.md` | §15g/§15h |
