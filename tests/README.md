# dvm-dos-tem tests (daily-climate-mode)

Two layers:

* `tests/unit/` has unit tests in C++. They use [doctest](https://github.com/doctest/doctest) 2.4.11, which is vendored as the header-only `doctest.h`. The test binary links every model object except `obj/TEM.o`, which holds `main`.
* `tests/regression/` has a scripted regression harness. It compares two `dvmdostem` binaries on a case set and checks bit-identity, determinism, mass closure and run health.

None of these files change model behaviour. The only model-source change made for testing is the behaviour-neutral seam `include/NFactorSelect.h` (see "Seams" below). Since v5k the suite also covers the switchable fixes B1-B4 and the draft N fixes N1-N3 (all default 0 = legacy; see "Switchable fixes" below).

## Unit tests

```bash
./build-conda.sh -j8 unittest          # builds obj/*.o, then tests/unit/dvmdostem_unittest, then runs it
# or, with any toolchain:
make unittest CC=g++ SITE_SPECIFIC_INCLUDES=... SITE_SPECIFIC_LIBS=...
# re-run / filter:
tests/unit/dvmdostem_unittest -ts=seb              # one suite
tests/unit/dvmdostem_unittest -tc="*Kersten*"      # by name
tests/unit/dvmdostem_unittest -lts                 # list suites
tests/unit/dvmdostem_unittest --order-by=rand      # check that tests are independent
```

Test objects depend on `include/*.h`, so they rebuild when a model header changes. They have to: the tests reach private members, so a layout change needs a recompile.

How the tests are written:

* `test_access.h` first includes every std, boost, netcdf and jsoncpp header. It then defines `private`/`protected` as `public` and includes the model headers. This gives direct access to members and private helpers without editing model code, and it does not change any class layout.
* Code paths that call `exit()` are tested in a forked child: `ut::exit_status_of`, a "death test". A possible infinite loop is run in a forked child under `alarm()`.
* Climate and CO2 tests write tiny synthetic NetCDF files to `/tmp` (`nc_fixtures.h`). No other test uses input files.
* Expected values come from physics, the literature or hand calculation, and the source is given in each test. Examples: WMO/Tetens SVP, Côté & Konrad 2005, Sturm 1997, Goodrich 1982, Yi et al. 2009, FAO-56 Ra.
* Known bugs are not fixed here. Each one is a test marked `* doctest::should_fail(true)`, with `[KNOWN BUG: xfail]` in its name and a comment explaining it. Such a test "passes" while the bug exists, and doctest reports it as failing once the bug is fixed. When that happens, remove the decorator.

| suite | file | area |
|---|---|---|
| climate | test_climate.cpp | SVP, VPD (vp_hpa_fix 0/1), daily precip/drizzle/snow, DOY helpers, daily girr, daily-file length/calendar fail-fast |
| co2_nlost | test_climate.cpp | SC CO2 indexed by calendar year (offset), coverage/gap/calendar guards |
| kersten | test_soil_thermal.cpp | Kersten number ln vs log10 (kersten_form), frozen Ke = Sr, bounds |
| thermal_options | test_soil_thermal.cpp | mineral texture formulas, mineral_tcsolid/vhcsolid/tcdry_form, org_tc_scheme (Côté & Konrad) |
| snow | test_snow_organic.cpp, test_seb.cpp | Goodrich vs Sturm conductivity, the v5j snow_tc_form switch, densification, snow albedo |
| organic_thickness | test_snow_organic.cpp | Yi et al. 2009 thickness ↔ C (CMT01/CMT51/HPC), expelled water |
| stefan | test_stefan_water.cpp | front rebuild when the deque is empty, all-fronts-swept, combineExtraFronts |
| water | test_stefan_water.cpp | freeze-expansion clip → expelled water, partial clip, puddle kept |
| richards | test_stefan_water.cpp | MAX_RICH_IND sizing, clearRichardsArrays, index/substep fail-fast guards, early-return branch |
| integrator | test_integrator_bgc.cpp | soil integrator safety net (delivered C/N, counters) |
| root_n | test_integrator_bgc.cpp | N-conserving per-layer root N extraction (hand values plus a 200-case property test) |
| nfactor | test_nfactor.cpp | n-factor selection, nfactor_sameday 0/1, first day of month |
| seb | test_seb.cpp | v5j SEB scheme 2 (Soil_Env::solveSEB): energy closure, melt, brackets, vp units |
| fixes_v5k | test_fixes_v5k.cpp | v5k switches: B1 svp over ice (+SEB esat), B2 Richards early-return infiltration, B3 combineExtraFronts termination, B4 snow albedo cap; N1-N3 soil N conservation (deltastate), non-vascular N extraction; C1 soil integrator NMIN tolerance |

Tests in `test_seb.cpp` need the v5j commits (9dae595b and d41c7f79), which are cherry-picked onto this branch.

## Regression harness

```bash
# 1. once per machine: materialise a self-contained case set (reads the calibration sources read-only)
python3 tests/regression/prepare_caseset.py /path/calib8
python3 tests/regression/prepare_caseset.py /path/calib8_nfday --ms nfactor_sameday=1   # switch variants
# 2. compare two binaries (A = reference, B = candidate)
tests/regression/run_regression.sh -j 8 BIN_A BIN_B /path/calib8 /path/out_dir
#    --no-expect-identical  intended behaviour change: skip the A==B requirement, keep the other checks
#    --no-repeat            skip the determinism re-run of B
#    --strict               turn WARNs (relative closure, extra safety-net events) into FAILs
# diagnostic env vars for a diagnostic binary (recorded in meta.json), e.g.
REG_EXTRA_ENV="DVM_ET_DIAG=1" python3 tests/regression/run_case.py BIN CASEDIR RUNDIR
```

Protocol for each case is in `run_case.py`. It is the same launcher and argv as the calibration runs (`vrun_seb.py`), with one segment:

1. `spin000`: a 20-year `sc` spin-up on the 1950-69 forcing, started from the case's `restart-tr.nc`.
2. `main`: `sc -n 75` (1950-2024), started from the spin-up's `restart-sc.nc`.

Both use `-p 0 -e 0 -s 0 -t 0 -l disabled --no-output-cleanup`. The environment contains only `PATH=/usr/bin:/bin HOME=/tmp LANG=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 LD_LIBRARY_PATH=$REG_LD_LIBRARY_PATH`, and stdin is `/dev/null`. Runs A, B and R (a repeat of B) go in sibling directories whose paths have the same length.

Checks, in `summarize.py` (results in `out_dir/report.txt` and `report.json`):

| check | rule |
|---|---|
| A vs B bit-identity | Every variable of every NetCDF file in `spin000/output` and `main/output` must match (`run_status:total_runtime` ignored), and so must `integrator_safety_net.csv` (`compare_runs.py`, the logic of v5i `bitid.py`). Required unless `--no-expect-identical`. |
| determinism | B vs R must be bit-identical (always required). |
| run health | exit 0, `run_status == 100`, no `fail_log*`, for both A and B. |
| water closure (B) | \|P − EET − QRUNOFF − QDRAINAGE − ΔS\| < 1 % of P over 1951-2024 (ΔS = layer LWC+IWC plus SWE). |
| C closure (B) | max over years of \|Δ(SOC+VEGC+DEADC+DWDC) + NEE\| < 0.5 gC/m²/yr. |
| N closure (B) | \|mean Δ(VEGNTOT+ORGN+AVLN+DEADN+DWDN) − (NINPUT − NLOST)\| < 1 gN/m²/yr. A coarse bound: it catches NLOST double counting (10-1000s gN/m²/yr), but the nfeed cases carry an unexplained ~0.2 gN/m²/yr residual (see REPORT). |
| WARN | B has more integrator safety-net or veg give-up events than A; B's water or N closure is worse than A's. These become FAILs with `--strict`. |

**Bit-identity needs a fixed environment.** On the GCP VM, outputs were seen to depend on the size of the inherited environment. Run the harness on the box with `REG_LD_LIBRARY_PATH` left at its default `/workspace/wf4-runtime/lib`, and compare only runs made on the same machine with the same runtime libraries. The VM-built binaries carry an rpath to the conda env, which the box does not have, so `LD_LIBRARY_PATH` supplies the libraries there.

## Seams (model-source changes made for testability)

* `include/NFactorSelect.h`: `nfactor_sel::selector_tair`, `nfactor_sel::select_nfactor` and `TDRV_MONTH_START`. This is the n-factor selection code from `Cohort::updateMonthly_Env`, moved unchanged. Proven bit-identical on all 8 calibration cases (spin-up and main, every variable) with `nfactor_sameday` = 0 and = 1.

## Switchable fixes (v5k: merged in daily-climate-mode, tags v5k (B1-B4, N1-N3) and v5k1 (+C1); all default 0 = legacy, bit-identical when off)

| key (model_settings) | values | fixes | tests (switch on) | legacy xfail kept |
|---|---|---|---|---|
| `svp_ice_fix` | 0/1 | B1: svp at t <= 0 used the Tetens water coefficient with the WMO ice denominator (+51 % at -20 C); 1 = WMO ice formula 6.112 exp(22.46 t/(272.62+t)) hPa in `calculate_saturated_vapor_pressure` and in the SEB esat | values vs WMO/Murphy-Koop, continuity at 0 C, daily svp/vpd, SEB sublimation | test_climate.cpp (2) |
| `richards_infil_fix` | 0/1 | B2: `Richards::update` early return put the infiltration into Richards' private EnvData (lost); 1 = `Soil_Env::rebookUnappliedRichards` sends it to the puddle (<= ponding_max_mm) and runoff, books evaporation as unapplied, uptake as not applied | flag, closure soil+puddle+runoff, over-full puddle | test_stefan_water.cpp |
| `stefan_fronts_fix` | 0/1 | B3: `Stefan::combineExtraFronts` set `mindz` once (hang with > 1 surplus pair); 1 = reset per pass | termination (fork+alarm), nearest pair removed, rebuild with > MAX_NUM_FNT+2 boundaries | test_stefan_water.cpp |
| `snow_albedo_fix` | 0/1 | B4: `Snow_Env::getAlbedoVis` returned 0.9 below -10 C; 1 = snwalbmax | cap, continuity at -10 C | test_seb.cpp |
| `n_rootless_fix` | 0/1/2 | N1 (draft): rootless layers lost their net mineralisation N (and ignored extraction); 1 = no mineralisation there, 2 = mineralised N to AVLN (extraction applied) | deltastate conservation | (documented in the N1 test, fx = 0) |
| `n_som_move_fix` | 0/1 | N2 (draft): SOM N moved shallow->deep was overwritten ('=' not '+=') in the top humic layers (main HPC N sink), receivers got the full instead of the thickness share, donors gave N for the cumulative C | deltastate conservation, equal donors | N2 legacy xfail |
| `soil_nmin_tol_fix` | 0/1 | C1: the soil step control tested the signed net-N-mineralisation accumulator NMIN against ptol*\|NMIN\|; when NMIN ~ 0 every doubled step is rejected, adapt() crawls at dt ~ 2^-19 (264k steps, TVC_B W nfeed Nov 2012 with B1+B2+B4) and the float soil state loses the tiny C decrements to rounding (+2.74 gC unbooked, runtime 29 s vs 16 s); 1 = tolerance ptol*max(\|NMIN\|, \|NIMMOB\|, \|AVLN\|) | boundcon accept/reject cases | legacy-behaviour test (reject on NMIN ~ 0) |
| `n_nonvasc_extract` | 0/1 | N3 (draft): non-vascular N uptake had no soil extraction (implicit N source); 1 = extracted from moss/fibric layers by avln*liq and capped at 95 % of that avln | afterIntegration extraction | (documented, fx = 0) |

Proof of neutrality: `run_regression.sh` A = bin_v5j vs B = bin_v5k / bin_v5k1, all new switches 0, on the 8 calib8 cases (with determinism) and on TVC_B + HPC_H1 with the v5j combined-round settings: bit-identical (see /workspace/v5k/V5K_REPORT.md).
