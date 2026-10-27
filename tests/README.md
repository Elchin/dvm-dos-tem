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
| fixes_v5m | test_fixes_v5m.cpp | v5m draft L1 `n_safetynet_conserve`: default/parsing (+ death tests); legacy (off) state unchanged; on: soil N change = remaining inputs − (1 − t0)·root extraction (veg + soil N conserved), donor order ORGN (by extraction share) → ORGN → AVLN, rootless layers as in deltastate, zero pools → debt carried and paid, negative pools untouched, `debitProportional` bookkeeping |

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

| key (model_settings) | values | fixes | tests (switch on) | legacy xfail kept | references |
|---|---|---|---|---|---|
| `svp_ice_fix` | 0/1 | B1: svp at t <= 0 used the Tetens water coefficient with the WMO ice denominator (+51 % at -20 C); 1 = WMO ice formula 6.112 exp(22.46 t/(272.62+t)) hPa in `calculate_saturated_vapor_pressure` and in the SEB esat | values vs WMO/Murphy-Koop, continuity at 0 C, daily svp/vpd, SEB sublimation | test_climate.cpp (2) | WMO-No. 8 (2008) CIMO Guide Annex 4.B ([WMO e-Library](https://library.wmo.int/records/item/41650-guide-to-meteorological-instruments-and-methods-of-observation)); cf. Buck 1981 [doi:10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2](https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2); test values Murphy & Koop 2005 [doi:10.1256/qj.04.94](https://doi.org/10.1256/qj.04.94) |
| `richards_infil_fix` | 0/1 | B2: `Richards::update` early return put the infiltration into Richards' private EnvData (lost); 1 = `Soil_Env::rebookUnappliedRichards` sends it to the puddle (<= ponding_max_mm) and runoff, books evaporation as unapplied, uptake as not applied | flag, closure soil+puddle+runoff, over-full puddle | test_stefan_water.cpp | no published source; code bug (upstream 419a4eb3, fix 07f6c8d0); restores water mass conservation (Richards 1931 [doi:10.1063/1.1745010](https://doi.org/10.1063/1.1745010)) |
| `stefan_fronts_fix` | 0/1 | B3: `Stefan::combineExtraFronts` set `mindz` once (hang with > 1 surplus pair); 1 = reset per pass | termination (fork+alarm), nearest pair removed, rebuild with > MAX_NUM_FNT+2 boundaries | test_stefan_water.cpp | no published source; loop termination (upstream 419a4eb3, fix 07f6c8d0); algorithm Woo et al. 2004 [doi:10.1029/2004GL019475](https://doi.org/10.1029/2004GL019475) |
| `snow_albedo_fix` | 0/1 | B4: `Snow_Env::getAlbedoVis` returned 0.9 below -10 C; 1 = snwalbmax | cap, continuity at -10 C | test_seb.cpp | no published source; consistency with snwalbmax / continuity (upstream 419a4eb3, fix 07f6c8d0); cf. Warren 1982 [doi:10.1029/RG020i001p00067](https://doi.org/10.1029/RG020i001p00067) |
| `n_rootless_fix` | 0/1/2 | N1 (draft): rootless layers lost their net mineralisation N (and ignored extraction); 1 = no mineralisation there, 2 = mineralised N to AVLN (extraction applied) | deltastate conservation | (documented in the N1 test, fx = 0) | no published source; N mass conservation (Soil_Bgc::deltastate, fix 29d72af8) |
| `n_som_move_fix` | 0/1 | N2 (draft): SOM N moved shallow->deep was overwritten ('=' not '+=') in the top humic layers (main HPC N sink), receivers got the full instead of the thickness share, donors gave N for the cumulative C | deltastate conservation, equal donors | N2 legacy xfail | no published source; N mass conservation, N moves with C (Soil_Bgc::deltastate, fix 20b2723e) |
| `soil_nmin_tol_fix` | 0/1 | C1: the soil step control tested the signed net-N-mineralisation accumulator NMIN against ptol*\|NMIN\|; when NMIN ~ 0 every doubled step is rejected, adapt() crawls at dt ~ 2^-19 (264k steps, TVC_B W nfeed Nov 2012 with B1+B2+B4) and the float soil state loses the tiny C decrements to rounding (+2.74 gC unbooked, runtime 29 s vs 16 s); 1 = tolerance ptol*max(\|NMIN\|, \|NIMMOB\|, \|AVLN\|) | boundcon accept/reject cases | legacy-behaviour test (reject on NMIN ~ 0) | no published source for the fix; error scaling: Press et al. 2007 ([Numerical Recipes 3rd ed.](https://numerical.recipes/book.html)) §17.2; Hairer et al. 1993 [doi:10.1007/978-3-540-78862-1](https://doi.org/10.1007/978-3-540-78862-1) |
| `n_nonvasc_extract` | 0/1 | N3 (draft): non-vascular N uptake had no soil extraction (implicit N source); 1 = extracted from moss/fibric layers by avln*liq and capped at 95 % of that avln | afterIntegration extraction | (documented, fx = 0) | no published source for the scheme; N mass conservation; bryophyte N uptake Turetsky 2003 [doi:10.1639/05](https://doi.org/10.1639/05) |

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

## Switchable fix (v5m draft, branch v5m-nleak, NOT merged; default 0 = legacy, bit-identical to v5l1 when off)

Soil N leak in the integrator safety net (be8a62fa). Unit tests: suite `fixes_v5m` (8 cases). Report: /workspace/v5m/NLEAK_REPORT.md.

| key (model_settings) | values | fixes | tests | references |
|---|---|---|---|---|
| `n_safetynet_conserve` | 0/1 | L1: when the soil integrator gives up (adapt() resets the soil to month fraction t0) `applySoilSafetyNet` delivers the remaining inputs but not the remaining root N extraction, although the vegetation (integrated first) already received the whole month's N uptake: each event creates (1 − t0)·Σ rtnextract g N/m2. 1 = debit that N from soil ORGN of the rooted layers (by extraction share, capped at the pool), then ORGN of any layer, then AVLN; any rest (all pools empty) is carried and debited at the start of the next soil month; the remaining N deposition is delivered where deltastate puts it (also into a rootless top layer with `n_rootless_fix` = 2). C: no leak, unchanged | off = legacy; on: veg + soil N conserved, donors, rootless layers, zero pools (debt carried, paid), negative pools untouched, helper bookkeeping, parse/death tests | no published source for the fix (code: Integrator::applySoilSafetyNet, be8a62fa); closed TEM N cycle, plant N uptake debited from soil N: McGuire et al. 1992 [doi:10.1029/92GB00219](https://doi.org/10.1029/92GB00219), Raich et al. 1991 [doi:10.2307/1941899](https://doi.org/10.2307/1941899), Euskirchen et al. 2009 [doi:10.1890/08-0806.1](https://doi.org/10.1890/08-0806.1), Yuan et al. 2012 [doi:10.1890/11-1957.1](https://doi.org/10.1890/11-1957.1); ORGN as first donor (plant organic N uptake): Näsholm et al. 1998 [doi:10.1038/31921](https://doi.org/10.1038/31921), Schimel & Bennett 2004 [doi:10.1890/03-8002](https://doi.org/10.1890/03-8002); state resets/clipping must book the removed mass: Shampine et al. 2005 [doi:10.1016/j.amc.2004.12.011](https://doi.org/10.1016/j.amc.2004.12.011); N closure in C-N spin-up: Thornton & Rosenbloom 2005 [doi:10.1016/j.ecolmodel.2005.04.008](https://doi.org/10.1016/j.ecolmodel.2005.04.008) |

Diagnostics with the switch on (stderr only; `integrator_safety_net.csv` format unchanged): a line `n_safetynet_conserve: root N extraction to debit=… debited=… carried_debt=…` after each safety-net warning.

## References for every daily-mode switch (v5-daily … v5l)

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

The v5k, v5l and v5m (draft) switches are referenced in the "references" column of their tables above.
