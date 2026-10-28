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
| fixes_v5l | test_fixes_v5l.cpp | v5l ET-split switches: key defaults/validation/parsing (+ death tests), S1 per-cover stomatal light, S2 buried non-vascular snow interception, S3 LAI-based rain interception with storage, S4 moss rain as surface water (relabel / to ground), evap_moss accumulation |
| fixes_v5m | test_fixes_v5m.cpp | v5m L1 `n_safetynet_conserve`: default/parsing (+ death tests); legacy (off) state unchanged; on: soil N change = remaining inputs − (1 − t0)·root extraction (veg + soil N conserved), donor order ORGN (by extraction share) → ORGN → AVLN, rootless layers as in deltastate, zero pools → debt carried and paid, negative pools untouched, `debitProportional` bookkeeping |
| fixes_v5n | test_fixes_v5n.cpp | v5n G1 `n_extract_limiter` (+ `n_extract_tau`), G2 `soil_giveup_retry`, G3 `n_netnmin_rescale_guard`, `soil_giveup_diag`: defaults/parsing (+ death tests); G1 limiter shape, zero-order sink drives AVLN < 0 / limiter keeps it ≥ 0, off and ample-AVLN bit-identical, AVLN/ORGN split with soil N change unchanged, rooted/rootless application; G2 per-month flag; G3 guard 0 = legacy loop bit for bit, guard 1 skips a negative scale and applies a positive one, zero denominator |
| daily_c_split | test_daily_c_split.cpp | v5p `daily_c_split` (output only): key default/parsing (+ death tests for 2, −1, 0.5, true, "1"), `*_DS` request rules (daily only, FATAL with the switch off), proportional split and exact monthly sum, NaN/negative weights = 0, zero flux, all-zero-weight uniform fallback and its flags/month counter, RG weight max(GPP−RM, 0) → GPP weight → uniform, derived NPP/INNPP, RECO/NEE identities, 500-month random property test (|Σ−M| ≤ 1e-13·max(1,|M|)), model weights from Vegetation_Bgc/Soil_Bgc, layer remap helpers |
| daily_water_out | test_daily_water_out.cpp | v5p-daily-water-fix (output only): `dwout::record_tran` + `veg_endOfDay` daily TRANSPIRATION sums exactly to `m_v2a.tran` (12 months), out-of-range day indices ignored, `append_days` order/clamping, `atm_endOfDay` `daily_eet`/`daily_pet` sum exactly to `m_l2a.eet`/`pet`, `Cohort::getEd4land_daily` with 3 PFTs: cohort EET/PET = veg + ground while the old per-PFT sum = veg + 3 × ground |
| bnf_v5q | test_bnf_v5q.cpp | v5q biological N fixation (`bnf_scheme`, `bnf_*`): defaults (Meyerholt et al. 2016 Table A1), validation (+ death tests: scheme 4/−1/1.5/true/"1", non-numeric and wrong-sign parameters, optional keys with scheme 0, scheme > 0 without n_rootless_fix 2), PRO/AET values vs the paper and the Cleveland et al. 1999 kgN/ha form, monthly weights, previous-year memory, Soil_Bgc begin/record, deltastate N conservation with a rootless top layer, safety-net delivery (scheme 0 = legacy), output row; Phase 2 `bnf_to_plant`: parse/death, month_parts, PFT shares, Soil_Bgc split, Vegetation_Bgc LABN conservation, Integrator::monthFractionDone, BNF_VEG row |

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
| N closure (B) | \|mean Δ(VEGNTOT+ORGN+AVLN+DEADN+DWDN) − (NINPUT [+ BNF_VEG, v5q] − NLOST)\| < 1 gN/m²/yr, **enforced only when N1-N3 are all on** in B (`n_rootless_fix` 1/2, `n_som_move_fix` 1, `n_nonvasc_extract` 1; read from B's main config.js, missing keys = binary default: v5o 2/1/1, or 0/0/0 with `closure.py --n-defaults legacy` / `REG_N_DEFAULTS=legacy` for pre-v5o binaries). Otherwise a WARN reports the residual and names the three bugs (legacy nfeed residual ~0.2 gN/m²/yr mean, single years up to ~0.9). With the fixes on the nfeed cases close to ≤ 2.7e-4 gN/m²/yr in every year (float32 output round-off). Logic tested by `python3 tests/regression/test_closure.py`. |
| WARN | B has more integrator safety-net or veg give-up events than A; B's water or N closure is worse than A's; B's N closure is not enforced (N1-N3 not all on). These become FAILs with `--strict`. |

**Bit-identity needs a fixed environment.** On the GCP VM, outputs were seen to depend on the size of the inherited environment. Run the harness on the box with `REG_LD_LIBRARY_PATH` left at its default `/workspace/wf4-runtime/lib`, and compare only runs made on the same machine with the same runtime libraries. The VM-built binaries carry an rpath to the conda env, which the box does not have, so `LD_LIBRARY_PATH` supplies the libraries there.

## Seams (model-source changes made for testability)

* `include/NFactorSelect.h`: `nfactor_sel::selector_tair`, `nfactor_sel::select_nfactor` and `TDRV_MONTH_START`. This is the n-factor selection code from `Cohort::updateMonthly_Env`, moved unchanged. Proven bit-identical on all 8 calibration cases (spin-up and main, every variable) with `nfactor_sameday` = 0 and = 1.

## Switchable fixes (v5k: merged in daily-climate-mode, tags v5k (B1-B4, N1-N3) and v5k1 (+C1); all default 0 = legacy, bit-identical when off. **Since v5o N1-N3 default ON (2/1/1)**, see the v5o section)

| key (model_settings) | values | fixes | tests (switch on) | legacy xfail kept | references |
|---|---|---|---|---|---|
| `svp_ice_fix` | 0/1 | B1: svp at t <= 0 used the Tetens water coefficient with the WMO ice denominator (+51 % at -20 C); 1 = WMO ice formula 6.112 exp(22.46 t/(272.62+t)) hPa in `calculate_saturated_vapor_pressure` and in the SEB esat | values vs WMO/Murphy-Koop, continuity at 0 C, daily svp/vpd, SEB sublimation | test_climate.cpp (2) | WMO-No. 8 (2008) CIMO Guide Annex 4.B ([WMO e-Library](https://library.wmo.int/records/item/41650-guide-to-meteorological-instruments-and-methods-of-observation)); cf. Buck 1981 [doi:10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2](https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2); test values Murphy & Koop 2005 [doi:10.1256/qj.04.94](https://doi.org/10.1256/qj.04.94) |
| `richards_infil_fix` | 0/1 | B2: `Richards::update` early return put the infiltration into Richards' private EnvData (lost); 1 = `Soil_Env::rebookUnappliedRichards` sends it to the puddle (<= ponding_max_mm) and runoff, books evaporation as unapplied, uptake as not applied | flag, closure soil+puddle+runoff, over-full puddle | test_stefan_water.cpp | no published source; code bug (upstream 419a4eb3, fix 07f6c8d0); restores water mass conservation (Richards 1931 [doi:10.1063/1.1745010](https://doi.org/10.1063/1.1745010)) |
| `stefan_fronts_fix` | 0/1 | B3: `Stefan::combineExtraFronts` set `mindz` once (hang with > 1 surplus pair); 1 = reset per pass | termination (fork+alarm), nearest pair removed, rebuild with > MAX_NUM_FNT+2 boundaries | test_stefan_water.cpp | no published source; loop termination (upstream 419a4eb3, fix 07f6c8d0); algorithm Woo et al. 2004 [doi:10.1029/2004GL019475](https://doi.org/10.1029/2004GL019475) |
| `snow_albedo_fix` | 0/1 | B4: `Snow_Env::getAlbedoVis` returned 0.9 below -10 C; 1 = snwalbmax | cap, continuity at -10 C | test_seb.cpp | no published source; consistency with snwalbmax / continuity (upstream 419a4eb3, fix 07f6c8d0); cf. Warren 1982 [doi:10.1029/RG020i001p00067](https://doi.org/10.1029/RG020i001p00067) |
| `n_rootless_fix` | 0/1/2 (**default 2 since v5o**) | N1: rootless layers lost their net mineralisation N (and ignored extraction); 1 = no mineralisation there, 2 = mineralised N to AVLN (extraction applied) | deltastate conservation | (documented in the N1 test, fx = 0) | no published source; N mass conservation (Soil_Bgc::deltastate, fix 29d72af8) |
| `n_som_move_fix` | 0/1 (**default 1 since v5o**) | N2: SOM N moved shallow->deep was overwritten ('=' not '+=') in the top humic layers (main HPC N sink), receivers got the full instead of the thickness share, donors gave N for the cumulative C | deltastate conservation, equal donors | legacy loss documented with explicit 0 (former xfail now a passing default test, v5o) | no published source; N mass conservation, N moves with C (Soil_Bgc::deltastate, fix 20b2723e) |
| `soil_nmin_tol_fix` | 0/1 | C1: the soil step control tested the signed net-N-mineralisation accumulator NMIN against ptol*\|NMIN\|; when NMIN ~ 0 every doubled step is rejected, adapt() crawls at dt ~ 2^-19 (264k steps, TVC_B W nfeed Nov 2012 with B1+B2+B4) and the float soil state loses the tiny C decrements to rounding (+2.74 gC unbooked, runtime 29 s vs 16 s); 1 = tolerance ptol*max(\|NMIN\|, \|NIMMOB\|, \|AVLN\|) | boundcon accept/reject cases | legacy-behaviour test (reject on NMIN ~ 0) | no published source for the fix; error scaling: Press et al. 2007 ([Numerical Recipes 3rd ed.](https://numerical.recipes/book.html)) §17.2; Hairer et al. 1993 [doi:10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1) |
| `n_nonvasc_extract` | 0/1 (**default 1 since v5o**) | N3: non-vascular N uptake had no soil extraction (implicit N source); 1 = extracted from moss/fibric layers by avln*liq and capped at 95 % of that avln | afterIntegration extraction | (documented, fx = 0) | no published source for the scheme; N mass conservation; bryophyte N uptake Turetsky 2003 [doi:10.1639/05](https://doi.org/10.1639/05) |

Proof of neutrality: `run_regression.sh` A = bin_v5j vs B = bin_v5k / bin_v5k1, all new switches 0, on the 8 calib8 cases (with determinism) and on TVC_B + HPC_H1 with the v5j combined-round settings: bit-identical (see /workspace/v5k/V5K_REPORT.md).

## Switchable fixes (v5l: merged in daily-climate-mode, local tag v5l; all default 0 = legacy, bit-identical to v5k1 when off)

ET-partition fixes from the v5k ET-split diagnostic. Unit tests: suite `fixes_v5l` (15 cases).

| key (model_settings) | values | fixes | tests | references |
|---|---|---|---|---|
| `et_fpc_once_fix` | 0/1 | S1: transpiration was cover-scaled twice: the radiation/PAR given to the stomatal and Penman-Monteith terms is already FPC-scaled (`updateRadiation`), and Cohort multiplies `d_v2a.tran` by FPC again ("ALT fix 1", 28167031); 1 = stomatal PAR and transpiration radiation use the per-cover values (pardown/fpc, sw/fpc), so Cohort's single FPC multiply is the only one. Wet-canopy evaporation and tran_pet are unchanged | legacy formula, per-cover formula independent of FPC, identity at FPC = 1 | double scaling from 28167031; cover once: Gerten et al. 2004 [doi:10.1016/j.jhydrol.2003.09.029](https://doi.org/10.1016/j.jhydrol.2003.09.029), Schaphoff et al. 2018 [doi:10.5194/gmd-11-1343-2018](https://doi.org/10.5194/gmd-11-1343-2018), [CLM5 Tech Note](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Hydrology/CLM50_Tech_Note_Hydrology.html) eq. 2.7.68; Monsi & Saeki 1953 [doi:10.1093/aob/mci052](https://doi.org/10.1093/aob/mci052) |
| `nv_snow_bury_fix` | 0/1 | S2: moss/lichen intercepted snow while under the snowpack (`getSnowInterception` = min(snfl*fpc, iscoef*LAI)) and the snow then sublimated as "canopy" sublimation; 1 = no snow interception and zero snow-storage capacity for non-vascular PFTs when the start-of-day snow depth > `nv_snow_bury_depth`; stored snow drips | buried/not buried, threshold, vascular unaffected, mass | Hedstrom & Pomeroy 1998 [doi:10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4); Sturm et al. 2001 [doi:10.1175/1520-0442(2001)014<0336:SSIIAT>2.0.CO;2](https://doi.org/10.1175/1520-0442(2001)014<0336:SSIIAT>2.0.CO;2) |
| `nv_snow_bury_depth` | m, >= 0, default 0.0 | S2 threshold. 0: the moss/lichen layer is the top of the soil column, so any snowpack lies on top of it | as above | moss layer is the top of the soil column: Yi et al. 2009 [doi:10.1029/2008JG000841](https://doi.org/10.1029/2008JG000841) |
| `rain_int_lai_fix` | 0/1 | S3: rain interception was max_int = 0.21 P + 0.5 mm, not related to LAI, and storage was not carried over in a consistent way; 1 = capacity = ircoef*LAI (mm, per ground area; ircoef in mm per unit LAI, from the CMT parameter file), intercepted = min(rain*fpc, capacity - stored), stored water evaporates first next day, excess drips. Literature 0.1-0.25 mm/LAI (CLM5 0.1 per LAI+SAI after Dickinson et al. 1993; BATS 0.2); parameter files have 0.041 | legacy formula, capacity, storage, evaporate first, shrinking capacity drips, mass | Dickinson et al. 1993 [doi:10.5065/D67W6959](https://doi.org/10.5065/D67W6959) / [CLM5 Tech Note](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Hydrology/CLM50_Tech_Note_Hydrology.html) eq. 2.7.12; von Hoyningen-Huene 1981 (DVWK, no DOI; [LISFLOOD refs](https://ec-jrc.github.io/lisflood-model/5_1_references/)), van der Knijff et al. 2010 [doi:10.1080/13658810802549154](https://doi.org/10.1080/13658810802549154); TVC obs Zwieback et al. 2019 [doi:10.1088/1748-9326/ab1049](https://doi.org/10.1088/1748-9326/ab1049), Chang et al. 2022 [doi:10.1016/j.rse.2021.112747](https://doi.org/10.1016/j.rse.2021.112747) |
| `nv_rain_surface` | 0/1/2 | S4: rain intercepted by moss/lichen was evaporated and reported as canopy evaporation; 1 = same water flux, booked as the new moss-surface evaporation `EVAP_MOSS` (in EET; closure unchanged); 2 = no moss/lichen rain interception; the rain reaches the ground (soil/moss-layer evaporation and infiltration) | relabel, to ground with and without S3, accumulation | Chadburn et al. 2015 [doi:10.5194/gmd-8-1493-2015](https://doi.org/10.5194/gmd-8-1493-2015); Yi et al. 2009 [doi:10.1029/2008JG000841](https://doi.org/10.1029/2008JG000841); Heijmans et al. 2004 [doi:10.1016/j.agrformet.2003.12.006](https://doi.org/10.1016/j.agrformet.2003.12.006) |

New outputs (output_spec.csv, off by default; monthly/yearly cohort totals, mm/time): `EVAP_CANOPY`, `SUBL_CANOPY`, `EVAP_MOSS`, `EVAP_SOIL`, `SUBL_SNOW`, `RAIN_INTERCEPT`, `SNOW_INTERCEPT`. EET = EVAP_CANOPY + SUBL_CANOPY + EVAP_MOSS + TRANSPIRATION + EVAP_SOIL + SUBL_SNOW. EVAP_MOSS counts in EET only when `nv_rain_surface` = 1.

Proof of neutrality: `run_regression.sh` A = bin_v5k1 vs B = bin_v5l (rc and merged), all new switches off, on calib8 and on the 8 V-Nall baseline cases (with determinism), plus V-Nall with all new keys written explicitly as 0: bit-identical, closure OK (see /workspace/v5l/V5L_REPORT.md).

## Switchable fix (v5m, merged bf05797e from v5m-nleak 23e8077b; default 0 = legacy, bit-identical to v5l1 when off)

Soil N leak in the integrator safety net (be8a62fa). Unit tests: suite `fixes_v5m` (8 cases). Report: /workspace/v5m/NLEAK_REPORT.md (archived on the VM in ~/box_archive/v5m/).

| key (model_settings) | values | fixes | tests | references |
|---|---|---|---|---|
| `n_safetynet_conserve` | 0/1 | L1: when the soil integrator gives up (adapt() resets the soil to month fraction t0) `applySoilSafetyNet` delivers the remaining inputs but not the remaining root N extraction, although the vegetation (integrated first) already received the whole month's N uptake: each event creates (1 − t0)·Σ rtnextract g N/m2. 1 = debit that N from soil ORGN of the rooted layers (by extraction share, capped at the pool), then ORGN of any layer, then AVLN; any rest (all pools empty) is carried and debited at the start of the next soil month; the remaining N deposition is delivered where deltastate puts it (also into a rootless top layer with `n_rootless_fix` = 2). C: no leak, unchanged | off = legacy; on: veg + soil N conserved, donors, rootless layers, zero pools (debt carried, paid), negative pools untouched, helper bookkeeping, parse/death tests | no published source for the fix (code: Integrator::applySoilSafetyNet, be8a62fa); closed TEM N cycle, plant N uptake debited from soil N: McGuire et al. 1992 [doi:10.1029/92GB00219](https://doi.org/10.1029/92GB00219), Raich et al. 1991 [doi:10.2307/1941899](https://doi.org/10.2307/1941899), Euskirchen et al. 2009 [doi:10.1890/08-0806.1](https://doi.org/10.1890/08-0806.1), Yuan et al. 2012 [doi:10.1890/11-1957.1](https://doi.org/10.1890/11-1957.1); ORGN as first donor (plant organic N uptake): Näsholm et al. 1998 [doi:10.1038/31921](https://doi.org/10.1038/31921), Schimel & Bennett 2004 [doi:10.1890/03-8002](https://doi.org/10.1890/03-8002); state resets/clipping must book the removed mass: Shampine et al. 2005 [doi:10.1016/j.amc.2004.12.011](https://doi.org/10.1016/j.amc.2004.12.011); N closure in C-N spin-up: Thornton & Rosenbloom 2005 [doi:10.1016/j.ecolmodel.2005.04.008](https://doi.org/10.1016/j.ecolmodel.2005.04.008) |

Diagnostics with the switch on (stderr only; `integrator_safety_net.csv` format unchanged): a line `n_safetynet_conserve: root N extraction to debit=… debited=… carried_debt=…` after each safety-net warning.

## Switchable fixes (v5n: commits 904a4a7c, 5c597aec, merge 4bb7cd99, tag v5n; all default 0 = legacy, bit-identical to v5m when off)

Soil-integrator give-ups (adapt() reaching dt_min 100 times in a month, then the safety net). Unit tests: suite `fixes_v5n` (11 cases). Report: /workspace/v5n/GIVEUP_REPORT.md (archived on the VM in ~/box_archive/v5n/). Meant to be used with `n_safetynet_conserve` = 1. Defaults: `n_extract_limiter` = 0, `soil_giveup_retry` = 0, `n_netnmin_rescale_guard` = 0, `soil_giveup_diag` = 0; `n_extract_tau` = 0.2 is read only when G1 (or a G2 retry month) is on.

| key (model_settings) | values | fixes | tests | references |
|---|---|---|---|---|
| `n_extract_limiter` (+ `n_extract_tau`, default 0.2, range (0, 1]) | 0/1 | G1: root N extraction rtnextract is a zero-order (constant) sink on AVLN during the soil month, so AVLN in rooted layers 3-5 reaches 0 mid-month in N-stressed Jul/Aug and no step size can avoid it. 1 = the AVLN part is r·min(1, AVLN/(τ·r)), the rest of r comes from ORGN of the same layer (total extraction, plant uptake and N mass unchanged) | defaults/parse, limiter shape, mechanism, off/ample bit-identical, split, rooted/rootless | Saturating plant N uptake in TEM: Raich et al. 1991 [doi:10.2307/1941899](https://doi.org/10.2307/1941899), McGuire et al. 1992 [doi:10.1029/92GB00219](https://doi.org/10.1029/92GB00219); uptake down-regulated when mineral N cannot supply demand: Thornton et al. 2007 [doi:10.1029/2006GB002868](https://doi.org/10.1029/2006GB002868); organic N uptake: Näsholm et al. 1998 [doi:10.1038/31921](https://doi.org/10.1038/31921), Schimel & Bennett 2004 [doi:10.1890/03-8002](https://doi.org/10.1890/03-8002); positivity needs a state-dependent sink: Shampine et al. 2005 [doi:10.1016/j.amc.2004.12.011](https://doi.org/10.1016/j.amc.2004.12.011), Burchard et al. 2003 [doi:10.1016/S0168-9274(03)00101-6](https://doi.org/10.1016/S0168-9274(03)00101-6) |
| `soil_giveup_retry` | 0/1 | G2: the month runs as legacy; only if adapt() gives up is it re-integrated from the month start with the G1 limiter on for that month (all other months bit-identical to legacy) | per-month flag | as G1 |
| `n_netnmin_rescale_guard` | 0/1 | G3: when avlnlost < 0, deltan rescales every layer's net N mineralisation by nscale = nminadj/Σnetnmin. If Σnetnmin < 0 (net immobilisation) and nminadj > 0, nscale < 0 flips every sign: immobilising layers with AVLN ≈ 0 get a large negative net flux and AVLN goes negative. 1 = skip a negative nscale (positive scale as legacy) | guard 0 = legacy loop bit for bit, guard 1 skip/apply, zero denominator | **no published source for the fix** (code: Soil_Bgc::deltan rescale); net mineralisation as gross mineralisation minus microbial uptake: Schimel & Bennett 2004 [doi:10.1890/03-8002](https://doi.org/10.1890/03-8002); immobilisation limited by mineral N: Thornton et al. 2007 [doi:10.1029/2006GB002868](https://doi.org/10.1029/2006GB002868) |
| `soil_giveup_diag` | 0/1 | diagnostics only (no state change): per give-up / retry the pool, layer, cause (negative pool, error test; for a step accepted at dt_min the previous rejection), t0, AVLN, ORGN, rtnextract, netnmin, nimmob, nscale and the deltan budget terms; stderr `GIVEUP_DIAG` lines and `integrator_giveup_diag.csv`, plus a `v5n summary` line (retries, retry failures, forced months, RKF45 calls) | parse | none needed (diagnostic output) |

Result (n_safetynet_conserve = 1, 1950-2099, 5 cases): G1 (τ 0.2) + G3 = 0 give-ups in every case (legacy 45 in total), RKF45 calls on S2 HPC 9398 (legacy) → 3250, run time unchanged; see the report.

All-8 test (2026-10-27, 28 cases: calib8, V-Nall, adopted v5l, v5l+H9; 1950-2099): L1+G1+G3 = 0 give-ups (main and spin) vs 4 main + 2 spin for legacy and L1 alone; nonfeed cases bit-identical in all configs; differences from L1 alone at noise level (soil-T bias ≤ 0.003 °C, GPP ≤ 0.04 %). Switches off: bit-identical to v5m on calib8 and V-Nall (run_regression.sh). Report: /workspace/v5n/ALL8_REPORT.md (VM ~/box_archive/v5n/).

## v5o: N conservation fixes N1-N3 on by default (commit 5812ff19, tag v5o)

New defaults: `n_rootless_fix` = 2, `n_som_move_fix` = 1, `n_nonvasc_extract` = 1 (ModelData Json parse and default constructor, Soil_Bgc / Vegetation_Bgc members). The give-up switches L1 `n_safetynet_conserve`, G1 `n_extract_limiter`, G2 `soil_giveup_retry`, G3 `n_netnmin_rescale_guard` and `soil_giveup_diag` stay default 0. Legacy stays selectable: `"n_rootless_fix": 0, "n_som_move_fix": 0, "n_nonvasc_extract": 0` reproduces v5n bit-for-bit. Breaks bit-identity with v5m/v5n for configs that did not set the keys (only the nfeed cases change; nonfeed cases are bit-identical).

Why: with the legacy N code N does not close in the nfeed cases (closed N system, NINPUT = NLOST = 0): mean −0.24 gN/m²/yr at HPC NoWarming nfeed, single years up to 0.9 gN/m² (n75) and up to 6 gN/m² in 2025-2100 (v5n NFEED_RESIDUAL_REPORT.md: N2 sink, N3 source, N1).

Configs: `config/config.js` and `prepare_caseset.py` (constant `N_FIXES`) write the three keys explicitly; for the legacy code pass `--ms n_rootless_fix=0 --ms n_som_move_fix=0 --ms n_nonvasc_extract=0`. Unit tests: suite `fixes_v5o` (4 cases: defaults in members / ModelData() / ModelData(Json), legacy selectable and invalid values exit, default N1 and N3 conserve); the former N2 legacy xfail now runs with the default and passes (the legacy loss is still documented with explicit 0). Harness: `test_closure.py` (6 cases).

Proof (box, calib8, spin 20 + main 75, `run_regression.sh -j 4 --no-repeat`):
* explicit legacy 0/0/0: bin_v5n vs bin_v5o **bit-identical on all 8 cases** (RESULT PASS; N closure reported as WARN, not enforced).
* new defaults (no keys) vs bin_v5n: nonfeed 4 bit-identical, nfeed 4 differ (intended); RESULT PASS, N enforced. N residual mean ≤ 5.3e-5, max |annual| ≤ 1.7e-4 gN/m²/yr (legacy: mean −0.24…+0.05, max 0.86).
* defaults == explicit 2/1/1 (calib8 copy with the keys): data-identical on all 8 cases (spin-up and main).
* V-Nall nfeed (keys already explicit): bin_v5o data-identical to bin_v5n; N mean ≤ 3.2e-5, max |annual| ≤ 2.7e-4 gN/m²/yr (one TVC NoWarming year at 2.68e-4 ≈ 1 float32 ulp of the 2932 gN/m² total N).

## v5p: output-only daily split of the monthly C fluxes (`daily_c_split`, branch v5p-daily-split; default 0 = off, bit-identical to v5o)

Design: /workspace/v5o/daily_split/DESIGN_DAILY_SPLIT.md (approved 2026-10-28 with all six recommended options). Code: `include/DailyCSplit.h` (pure helpers), `Cohort::ds_*` (Cohort.cpp), daily driver copies in EnvData, weight wrappers `Vegetation_Bgc::ds_gpp_weight / ds_rm_weight`, `Soil_Bgc::ds_rh_weight`, output in Runner.cpp. Unit tests: suite `daily_c_split` (11 cases). Report: /workspace/v5p/V5P_REPORT.md (VM ~/box_archive/v5p/).

**This is a statistical split of a monthly flux, not daily physics.** The model still computes every C flux once per month. With the switch on, each monthly total M (per PFT, per soil layer) is divided over the days of the calendar month as X_d = M·w_d/Σw, with w_d the model's own environmental scalar evaluated with that day's drivers, so the daily values sum exactly (float round-off) to the monthly value. Nothing is fed back: pools, restarts and all monthly/yearly outputs are unchanged.

| key (model_settings) | values | what it does | tests | references |
|---|---|---|---|---|
| `daily_c_split` | 0/1 (anything else: FATAL) | 1 = write the daily-only outputs `GPP_DS, NPP_DS, RM_DS, RG_DS, RHSOM_DS, RHDWD_DS, RECO_DS, NEE_DS, INGPP_DS, INNPP_DS` (ecosystem totals, g/m2/day; files `<VAR>_DS_daily_<stage>.nc`) and the monthly flag output `DS_FALLBACK`. Weights: GPP/INGPP: ci/(kc+ci)·par/(ki+par)·ftemp·thaw of `getGPP` with daily Ta, PAR, EET/PET (gv), thaw; RM: `getRaq10`(Ta_d); RG (and RG* = INGPP−RM−INNPP): max(GPP_d − RM_d, 0), then the GPP weight, then uniform; NPP_d = GPP_d − RM_d − RG_d, INNPP_d = INGPP_d − RM_d − RG*_d; RHSOM per soil layer and RHDWD (first non-moss layer): rhmoist(SWS_d)·rhq10(Ts_d) (rhmoistfrozen below 0 °C); RECO_d = RHSOM_d + RHDWD_d + RM_d + RG_d; NEE_d = RECO_d − GPP_d. All-zero weights with M ≠ 0: uniform M/n, flagged in `DS_FALLBACK` (bits 1 GPP, 2 RM, 4 RG→GPP weight, 8 RG uniform, 16 RHSOM, 32 RHDWD; 64 = Env-loop and Bgc soil layers differ, Bgc layers mapped by midpoint depth, not a fallback). A `*_DS` row ticked while the switch is 0, or a `*_DS` row ticked Monthly/Yearly, is FATAL. The existing GPP/RM/… rows are unchanged (ticking Daily on them would drop their monthly file). | switch parsing (+ death tests), request rules, weights, zero flux, fallback + counter, sum exactness (hand cases and 500 random months), RG/NPP/RECO/NEE identities, model weight wrappers | **no published source for this exact method.** Driver-weighted downscaling of coarse model C fluxes: Olsen & Randerson 2004 [doi:10.1029/2003JD003968](https://doi.org/10.1029/2003JD003968); Fisher et al. 2016 [doi:10.5194/bg-13-4271-2016](https://doi.org/10.5194/bg-13-4271-2016). Exact preservation of the coarse total (temporal disaggregation; principle only): Denton 1971 [doi:10.1080/01621459.1971.10482227](https://doi.org/10.1080/01621459.1971.10482227); Chow & Lin 1971 [doi:10.2307/1928739](https://doi.org/10.2307/1928739) |

Proof (box, wf4 runtime, calib8 = 8 cases, spin-up 20 + main 75):
* switch 0: `run_regression.sh -j 4` bin_v5o vs bin_v5p **bit-identical on all 8 cases** (spin-up and main, with determinism; RESULT PASS).
* switch 1 vs switch 0 (same binary): every pre-existing variable of every spin-up/main output file, restarts included, **data-identical on all 8 cases**; the only extra files are the 11 DS files. Switch 1 without DS rows: data-identical to switch 0 (checked on TVC Warming nonfeed).
* Σ_days X_d − M over 900 months × 8 cases: |err| ≤ 5.7e-14 g/m2 (rel ≤ 2.5e-15) for GPP, RM, RG, RHSOM, RHDWD, RECO, NEE (vs RECO−GPP), INGPP, INNPP; NPP_DS sums to GPP−RM−RG to 5.7e-14, which differs from the model's NPP output by the model's own float32 identity residual (≤ 1.4e-5 g/m2/month).
* Cost: compute of the split itself not measurable (switch 1 without DS rows: same run time as off); writing the 10 daily files: main run ~21 s → ~26 s (serial, TVC Warming nonfeed); +2.2 MB per 75-year main run (+0.67 MB spin-up).

## Output metadata fix (branch v5p-transp-units, from v5p; metadata only, no switch, no numeric change)

`config/output_spec.csv`: the units of `TRANSPIRATION` were `mm/day` at every resolution, but the stored values are sums over the output step: `EnvData::veg_endOfDay` adds the daily `d_v2a.tran` (mm/day) into `m_v2a.tran`, and `EnvData::veg_endOfMonth` adds that into `y_v2a.tran`. Runner writes d_/m_/y_ for daily/monthly/yearly files. The label is now `mm/time`, so the NetCDF `units` attribute reads mm/day, mm/month or mm/year (ModelData replaces `time`), as for EVAP_CANOPY/EVAP_SOIL/EVAP_MOSS/SUBL_* (v5l). Evidence, v5p binary, TVC Warming nfeed, 3-year diagnostic runs:
* monthly July value 0.706 mm vs 0.0228 mm/day;
* Σ monthly = yearly to 2.2e-16;
* the monthly components add up to EET (v5o: T 29.1 + canopy evap 7.9 + soil evap 25.0 + moss evap 8.5 + canopy subl 5.7 + snow subl 13.2 = 89.3 mm/yr = EET 89.27).

Checked and left unchanged:
* `EET`, `PET` (`mm/m2/time`): the time part is correct at all three resolutions (Σ monthly = yearly). The `mm/m2` form is historic upstream wording and was not changed.
* `EVAP_*`, `SUBL_*`, `*_INTERCEPT` (`mm/time`): correct.
* `RAINFALL`/`SNOWFALL` (`mm`, step totals): no time unit, not wrong.

Known data issues found here (see /workspace/bnf/units/check.txt), fixed on branch v5p-daily-water-fix (next section):
* daily `TRANSPIRATION` output holds one value per year (Runner TODO "single daily value");
* daily `EET` is the per-PFT sum of `daily_eet`, which counts ground evaporation once per PFT (about 2.2× monthly at TVC);
* daily `PET` sums differ from monthly in the shoulder months.

Unit tests: suite `output_units` (2 cases: spec labels; day → month → year accumulation). Proof: run_regression bin_v5p vs bin_v5p_units on calib8: all data bit-identical; only the TRANSPIRATION `units` attribute differs (with the new spec). References: **no published source**; rationale: the units attribute must describe the stored value, CF Conventions v1.8 §3.1 "Units" [CF18](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html).

## Daily water output fix (branch v5p-daily-water-fix, from v5p-transp-units; output only, no switch)

Daily NetCDF files for `TRANSPIRATION`, `EET` and `PET` (cohort total, no PFT dimension) did not hold the values the model sums into the monthly outputs. Monthly values come from the cohort `edall`: `EnvData::atm_endOfDay` adds `d_l2a.eet`/`pet` into `m_l2a` and stores them in `daily_eet`/`daily_pet[day]`; `EnvData::veg_endOfDay` adds `d_v2a.tran` into `m_v2a.tran`.
* `EET`, `PET`: Runner summed `cohort.ed[ip].daily_eet/pet` over all PFTs. Each PFT's `d_l2a` holds the full ground terms (soil evaporation, snow sublimation; `assignGroundEd2pfts_daily`), so ground was counted once per active PFT, and inactive PFTs added stale arrays. Now Runner writes `edall->daily_eet/pet` (`dwout::append_days`).
* `TRANSPIRATION`: Runner pushed `edall->d_v2a.tran` once per month and wrote one value per year. Now `Cohort::updateMonthly_Env` records the final daily `d_v2a.tran` after `veg_endOfDay` into the new output-only array `EnvData::daily_tran[31]` (`dwout::record_tran`), and Runner writes all `DINY` days.
* Unchanged: model state, fluxes, monthly/yearly outputs, restarts, per-PFT `ed[ip]` arrays (used by v5p `daily_c_split` weights). The per-PFT daily branches of these variables are unreachable (PFT is "invalid" in output_spec) and were left as they were.

Unit tests: suite `daily_water_out` (5 cases). Proof (numbers in /workspace/v5p_dwater/DAILY_WATER_FIX_REPORT.md): run_regression bin_v5p_units vs bin_v5p_dwater on calib8: all monthly/yearly outputs and restarts bit-identical; daily diagnostic runs: Σ daily = monthly for the three variables. References: **no published source**; rationale: the daily output of an accumulated flux must be the term that is summed into the monthly output (output bookkeeping only).

## v5q: biological N fixation (`bnf_scheme`, branch v5q-bnf from v5p-transp-units; default 0 = off, bit-identical to v5p when off)

DVM-DOS-TEM had no BNF: with nfeed = 1 the only external soil N input `m_a2soi.avlninput` (output `NINPUT`) was always 0. Code: `include/Bnf.h` (pure helpers and parsing), `Soil_Bgc::bnfBeginMonth / bnfRecordMonth` and the `prepareIntegration` hook, the soil safety-net clause in `Integrator.cpp`, the calls in `Cohort::updateMonthly_Bgc`, output `BNF` in Runner.cpp. Unit tests: suite `bnf_v5q`. Report: /workspace/bnf/V5Q_BNF_REPORT.md.

| key (model_settings) | values | what it does | tests | references |
|---|---|---|---|---|
| `bnf_scheme` | 0/1/2/3 (anything else: FATAL) | 0 = off. With nfeed = 1 an annual rate B (gN/m2/yr) from the **previous model year** goes to AVLN of the top soil layer (`m_a2soi.avlninput`, so it is part of `NINPUT`), month m gets B·w_m with w_m = previous-year max(NPP_m, 0)/Σ (uniform 1/12 if Σ = 0). 1 PRO: B = c·(1 − exp(d·NPP)), NPP = previous-year ecosystem NPP (gC/m2/yr, < 0 → 0); 2 AET: B = max(0, a·ET + b), ET = previous-year transpiration + canopy evaporation (mm/yr; the paper's ET excludes soil evaporation); 3 CONST: B = `bnf_const`. The first model year of each process gets 0 (no memory; not in the restart file). **Needs `n_rootless_fix` = 2** (default since v5o; FATAL otherwise): then the top-layer input reaches AVLN whether the layer has roots or not, and the soil safety net also delivers the rest-of-month BNF there. nfeed = 0: no effect. Output `BNF` (g/m2/time, monthly/yearly; 0 when off) | values vs paper, validation/death tests, weights/memory, conservation (rootless top layer), safety net | Meyerholt et al. 2016 [doi:10.5194/bg-13-1491-2016](https://doi.org/10.5194/bg-13-1491-2016) Eqs. A1/A2, Table A1; Cleveland et al. 1999 [doi:10.1029/1999GB900014](https://doi.org/10.1029/1999GB900014); BNF into soil mineral N as in CLM-CN: Thornton et al. 2007 [doi:10.1029/2006GB002868](https://doi.org/10.1029/2006GB002868); caveat (no significant BNF–ET/NPP relation in field data): Davies-Barnard & Friedlingstein 2020 [doi:10.1029/2019GB006387](https://doi.org/10.1029/2019GB006387) |
| `bnf_pro_c`, `bnf_pro_d` | default 1.8 gN/m2/yr, −0.003 m2 yr/gC (c ≥ 0, d ≤ 0) | PRO parameters | values | Meyerholt et al. 2016 Table A1 |
| `bnf_aet_a`, `bnf_aet_b` | default 0.00234 gN/mm/m2, −0.0172 gN/m2/yr (a ≥ 0) | AET parameters (= Cleveland 0.0234·ET − 0.172 kgN/ha/yr) | values, unit equivalence | Meyerholt et al. 2016 Table A1; Cleveland et al. 1999 |
| `bnf_const` | default 0 gN/m2/yr (≥ 0; nonzero needs scheme > 0) | CONST rate (per CMT via the case config) | parse, Soil_Bgc | Cleveland et al. 1999 (biome rates) |
| `bnf_asym` | default 0 gN/m2/yr (≥ 0; nonzero needs scheme > 0) | constant free-living (e.g. moss-associated cyanobacteria) term added to B, same monthly weights; the temperature/moisture/light modulation of Meyerholt Eqs. A20-A23 is not implemented | memory test | Meyerholt et al. 2016; DeLuca et al. 2002 [doi:10.1038/nature01051](https://doi.org/10.1038/nature01051); Stewart et al. 2011 [doi:10.1657/1938-4246-43.2.267](https://doi.org/10.1657/1938-4246-43.2.267); Rousk et al. 2013 [doi:10.3389/fmicb.2013.00150](https://doi.org/10.3389/fmicb.2013.00150) |
| `bnf_avln_cap` | default 0 = no taper (≥ 0; nonzero needs scheme > 0) | factor max(0, 1 − AVLN/cap), AVLN = column total at the start of the month. **The linear form is ours (no published source)**; Meyerholt's FOR scheme only states BNF → 0 above 2 gN/m2 soil inorganic N | cap_factor, Soil_Bgc | Meyerholt et al. 2016 (FOR, Sect. 2.2.1) |
| `bnf_to_plant` (Phase 2) | 0/1 (else FATAL; 1 needs scheme > 0) | 0 = all BNF to soil AVLN (Phase 1). 1 = the scheme part (PRO/AET/CONST, symbiotic) goes to the labile N (LABN) of the PFTs with vegcov > 0, split by previous-year PFT NPP (positive part; equal shares if all ≤ 0; all PFTs incl. moss/lichen); `bnf_asym` stays in the soil; with no active PFT everything goes to the soil. Added in `Vegetation_Bgc::deltastate` (nfeed only); after a vegetation integrator give-up only the integrated fraction `bh_t0` is booked. Output `BNF_VEG` (g/m2/time; part of `BNF`, not of `NINPUT`); `closure.py` adds `BNF_VEG` to the N inputs when present. **The NPP-proportional PFT split is ours (no published source)** | parse/death, month_parts, pft_shares, per-PFT memory, Soil_Bgc split, deltastate LABN gain = BNF, monthFractionDone, output row, test_closure.py n_budget | symbiotic BNF to the plant, asymbiotic to the soil: Meyerholt et al. 2016 [doi:10.5194/bg-13-1491-2016](https://doi.org/10.5194/bg-13-1491-2016); split: no published source |

Design notes: the rate uses the previous year so that BNF is known before the month's vegetation step (no extra iteration); no calibrated parameter (incl. `fnloss`) was changed.

## References for every daily-mode switch (v5-daily … v5n)

Standing rule since 2026-10-22: every code change carries its scientific references in the code comment at the switch site, in the commit message and here. The full bibliography is in `docs/daily-mode-REFERENCES.md`, including the non-switch fixes, the SEB formulation details and the switches without a published source. All DOIs and URLs were checked on 2026-10-22.

| key | references |
|---|---|
| `climate_timestep` = daily | no published source for the mode; CF 365_day calendar ([CF 1.8](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html)); daily Ra: [FAO-56](https://www.fao.org/4/x0490e/x0490e00.htm) eq. 21 |
| `vp_hpa_fix` | input in hPa: Harris et al. 2020 [doi:10.1038/s41597-020-0453-3](https://doi.org/10.1038/s41597-020-0453-3); VPD units: [FAO-56](https://www.fao.org/4/x0490e/x0490e00.htm); Monteith 1965 ([PMID 5321565](https://pubmed.ncbi.nlm.nih.gov/5321565/)) |
| `nfactor_sameday` | Klene et al. 2001 [doi:10.1080/15230430.2001.12003416](https://doi.org/10.1080/15230430.2001.12003416); uninitialised-selector code bug (upstream 419a4eb3, fix 9b3a4840) |
| `nonvasc_trans_withdraw` | no published source; water mass conservation; context Heijmans et al. 2004 [doi:10.1016/j.agrformet.2003.12.006](https://doi.org/10.1016/j.agrformet.2003.12.006) |
| `seb_scheme` (+ `seb_*`) | bulk transfer Brutsaert 1982 [doi:10.1007/978-94-017-1497-6](https://doi.org/10.1007/978-94-017-1497-6); simplified stability after Louis 1979 [doi:10.1007/BF00117978](https://doi.org/10.1007/BF00117978) and Dyer 1974 [doi:10.1007/BF00240838](https://doi.org/10.1007/BF00240838); snow emissivity Warren 1982 [doi:10.1029/RG020i001p00067](https://doi.org/10.1029/RG020i001p00067); esat WMO-No. 8 (2008); [CLM5 Tech Note, Fluxes](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Fluxes/CLM50_Tech_Note_Fluxes.html) |
| `seb_lw_source` | Prata 1996 [doi:10.1002/qj.49712253306](https://doi.org/10.1002/qj.49712253306); Crawford & Duchon 1999 [doi:10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2](https://doi.org/10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2) |
| `seb_snow_tc` | Sturm et al. 1997 [doi:10.3189/S0022143000002781](https://doi.org/10.3189/S0022143000002781); legacy Goodrich 1982 [doi:10.1139/t82-047](https://doi.org/10.1139/t82-047) |
| `kersten_form` | Johansen 1975 [doi:10.21236/ADA044002](https://doi.org/10.21236/ADA044002); Lawrence & Slater 2008 [doi:10.1007/s00382-007-0278-1](https://doi.org/10.1007/s00382-007-0278-1) |
| `org_tc_scheme` (+ `org_tc_*`) | Côté & Konrad 2005 [doi:10.1139/t04-106](https://doi.org/10.1139/t04-106); Letts et al. 2000 [doi:10.1080/07055900.2000.9649643](https://doi.org/10.1080/07055900.2000.9649643); Chadburn et al. 2015 [doi:10.5194/gmd-8-1493-2015](https://doi.org/10.5194/gmd-8-1493-2015) |
| `org_drain_rate` (+ `org_drain_*`) | Quinton & Marsh 1999 [doi:10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D](https://doi.org/10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D); Letts et al. 2000 |
| `mineral_tcsolid` / `mineral_vhcsolid` / `mineral_tcdry_form` | Johansen 1975; Farouki 1981 [doi:10.21236/ADA111734](https://doi.org/10.21236/ADA111734); [CLM5 Tech Note, soil temperature](https://escomp.github.io/CTSM/release-clm5.0/tech_note/Soil_Snow_Temperatures/CLM50_Tech_Note_Soil_Snow_Temperatures.html) |
| `frozen_infil_eice` | Swenson et al. 2012 [doi:10.1029/2012MS000165](https://doi.org/10.1029/2012MS000165); CLM5 eq. 2.7.48 |
| `soil_evap_beta` | Lee & Pielke 1992 [doi:10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2](https://doi.org/10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2); Clapp & Hornberger 1978 [doi:10.1029/WR014i004p00601](https://doi.org/10.1029/WR014i004p00601) |
| TVC Sedges klai (parameter test, not adopted) | Wallace & Baltzer 2019 [doi:10.1007/s10021-019-00435-0](https://doi.org/10.1007/s10021-019-00435-0), data [doi:10.5683/SP2/F4PQRY](https://doi.org/10.5683/SP2/F4PQRY); Grünberg & Boike 2019 [doi:10.1594/PANGAEA.904270](https://doi.org/10.1594/PANGAEA.904270); Yazbeck et al. 2025 [doi:10.1029/2024JG008721](https://doi.org/10.1029/2024JG008721); Chang et al. 2022 [doi:10.1016/j.rse.2021.112747](https://doi.org/10.1016/j.rse.2021.112747); Zwieback et al. 2019 [doi:10.1088/1748-9326/ab1049](https://doi.org/10.1088/1748-9326/ab1049) |
| `daily_c_split` (v5p, output only) | **no published source for this exact method**; driver-weighted downscaling: Olsen & Randerson 2004 [doi:10.1029/2003JD003968](https://doi.org/10.1029/2003JD003968), Fisher et al. 2016 [doi:10.5194/bg-13-4271-2016](https://doi.org/10.5194/bg-13-4271-2016); exact coarse totals (principle only): Denton 1971 [doi:10.1080/01621459.1971.10482227](https://doi.org/10.1080/01621459.1971.10482227), Chow & Lin 1971 [doi:10.2307/1928739](https://doi.org/10.2307/1928739) |
| `bnf_scheme` (+ `bnf_*`, v5q) | Meyerholt et al. 2016 [doi:10.5194/bg-13-1491-2016](https://doi.org/10.5194/bg-13-1491-2016); Cleveland et al. 1999 [doi:10.1029/1999GB900014](https://doi.org/10.1029/1999GB900014); Thornton et al. 2007 [doi:10.1029/2006GB002868](https://doi.org/10.1029/2006GB002868); Davies-Barnard & Friedlingstein 2020 [doi:10.1029/2019GB006387](https://doi.org/10.1029/2019GB006387); free-living term: DeLuca et al. 2002 [doi:10.1038/nature01051](https://doi.org/10.1038/nature01051), Stewart et al. 2011 [doi:10.1657/1938-4246-43.2.267](https://doi.org/10.1657/1938-4246-43.2.267), Rousk et al. 2013 [doi:10.3389/fmicb.2013.00150](https://doi.org/10.3389/fmicb.2013.00150); taper form: no published source; `bnf_to_plant` PFT split: no published source |

The v5k, v5l and v5m switches are referenced in the "references" column of their tables above.
