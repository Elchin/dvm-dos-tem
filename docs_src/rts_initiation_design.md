# RTS initiation in dvm-dos-tem: design (adapted to the `RTS` branch)
Local branch `rts-initiation` is based on origin/RTS (391b7a6). It has not been pushed.

## What already existed on `RTS` (reused, not duplicated)
- `Thermokarst` / `ThermokarstData`, a copy of WildFire/FirData. Methods: `should_initiate(yr, m, stage)`,
  `initiate(yr)`, `getThermokarstSoilRemoval()`, `getThermokarstAbgVegetation()`.
- Forcing: `hist_/proj_exp_thermokarst_file` holding `exp_thermokarst_mask`,
  `exp_jday_of_thermokarst`, `_severity` and `_area` as (time, Y, X) arrays. Initiation happens in tr/sc stages when
  mask[yr]==1 and the month of jday matches.
  **This serves as the "RTS initiation year per cell" input.** No initiation (the rts_init_year=-1 case) is an all-zero mask.
- `Cohort::updateMonthly_Thermokarst` is called from `updateMonthly_Dsb` when
  `model_settings.thermokarst` is true. It follows the fire sequence:
  initiate, then Bd2Layer, then `Ground::updateMineralProperties`, then `adjustSoilAfterburn`, then Layer2Bd, then reset dims.
- Flux and diagnostic structs (`*_thermokarst`), plus the restart field `thermokarsta2sorgn`.
- Previous behavior: it removed soil down to the thaw depth / previous-year ALD. Removed C was sent to the atmosphere
  through `firpar.r_retain_c`, which was **never initialized**. Vegetation removal was hardcoded at 0.99. The branch had no outputs.

## Added on `rts-initiation`
| Concern | Implementation |
|---|---|
| Partial organic removal | `Thermokarst::organicRemovalThickness(folr, moss, shlw, deep, thaw)`. When folr>=0, removal = folr x total organic thickness. When folr<0, the legacy thaw-depth behavior is kept. |
| Lateral export | `partitionLostCN()` splits removed soil C/N + woody debris and removed veg + standing dead. A `c_resp_frac` share goes to the atmosphere (`soi2a`, `v2a`). The rest goes to the new `soi2lat_thermokarst` struct (orgc, orgn, vegc, vegn) in ThermokarstData. Nothing is retained. |
| Thaw | `Ground::thawToDepth(depth, drained)` uses the pure helper `Ground::thawColumn`. Whole layers inside the depth become ice to liquid, liquid is capped at poro*dz*1000 and the excess is drained, then tem>=0.01 and frozen=-1. Fronts inside the zone are dropped, and a thaw front is added if frozen ground remains below. |
| Drainage | After initiation `cd.drainage_type = thermokarst_post_drainage`, which Soil_Env reads every step. It is carried in restart files as `drainage_type`; older restarts without it still load. |
| Params | config.js `model_settings` (PLACEHOLDER, tunable): `thermokarst_folr`, `thermokarst_c_resp_frac`, `thermokarst_veg_remove_frac`, `thermokarst_thaw_depth`, `thermokarst_post_drainage`. Defaults keep the legacy behavior. |
| Outputs | `TKLATC`, `TKLATN`, `TKREMTHICK`, `TKICEMELT`, `DRAINTYPE`. tkdata is now recorded every month, not only in the initiation month. |

## Tests
- `make rts-test` runs the C++ unit tests (`tests/rts/test_rts_units.cpp`) and the pytest regression suite (`tests/rts/`).
- They run on the 1-cell synthetic case `testing-data/rts-synthetic` (regenerate it with `make_rts_synthetic.py`).
- References are in `tests/rts/reference`. Regenerate them with `make_reference.py`; compare with `compare_outputs.py`.
- CI is `.github/workflows/rts-tests.yml`.

## Blended-run plan
For each cell, run once with thermokarst off (baseline) and once with it on (RTS), then blend afterward by RTS area fraction.
`exp_area_of_thermokarst` could carry that fraction.

## Next steps / open questions
Parameterize from ground-ice maps (excess ice and subsidence); the `feature/thermokarst-prototype` branch has subsidence code.
Mineral-soil removal for deep slumps. Per-CMT parameter files instead of global config.

## Excess ice and subsidence (ported from feature/thermokarst-prototype)
- **Ported:** the excess-ice profile `excess = DENICE * overlap * f/(1-f)`, the melt rule "excess melts only once the layer thaws, subsidence = melted / DENICE", and the water-table coupling (WTD reduced by each new subsidence increment).
- **Not ported:** the prototype's finite-volume thermal column (`thermokarst::Column`) and its topology remapping.
- **Placement:** `Layer::excess_ice` sits only in frozen **mineral** layers and is not counted in dz, so it does not affect heat or water calculations while it stays frozen.
- **Calls:** `Ground::initExcessIce` runs once at the first thermokarst step (January of transient year 0). `Ground::meltExcessIce` runs monthly and again right after the RTS thaw.
- **Config:** `thermokarst_excess_fraction`, `thermokarst_excess_top`, `thermokarst_excess_bottom`. These are placeholders; f = 0 turns it off.
- **Restart:** `tkexinit`, `tksubsidence`, `EXICEsoil`.
- **Outputs:** `TKSUBSID`, `TKEXICE`, `TKEXMELT`.

## Two-stage execution workflow (scripts/run_rts_two_stage.py)

Mirrors the original thermokarst branch. Ground ice is initialized only after the historical period, and the RTS initiates during the transient run.

* **Stage 1:** pr/eq/sp (default `-p 100 -e 1000 -s 200`), thermokarst OFF → `stage1/restart-sp.nc`.
  One cell takes about 71 s.
* **Stage 2a:** historical tr run restarted from stage 1, thermokarst OFF → `stage2a_hist/restart-tr.nc`.
* **Stage 2b:** sc run restarted from 2a with the thermokarst module ON.
  * Excess ground ice (`thermokarst_excess_*`) is initialized in the first month.
  * The sc years before the flagged RTS year act as the short spin-up with ice present.
  * The RTS initiates in the year and day flagged in `projected-explicit-thermokarst*.nc`.

Config switches:
* `model_settings.thermokarst`
* `thermokarst_initiation`: false gives a baseline run, with ice and subsidence but no slump.
* `thermokarst_ch4`: rts-ch4 only.

Driver options:
* `--mode baseline` runs without initiation.
* `--stage1-restart` reuses a stage-1 restart.
* `--ch4` turns on the methane path.
* The driver writes `timings.json`.

Tests: `tests/rts/test_rts_two_stage.py` covers the short CI workflow. Set `RTS_FULL_TIMING=1` to also run the full stage 1 and assert it takes ≤ 300 s.
