# Version v5q-gpp-fleaf2: TVC S22 / HPC S21H MIP candidate (tag `tvc-hpc-mip-S22-S21H`)

Written 2026-11-02. This file is documentation only. No source, parameter or test file was changed when it was added.

## 1. What this version is

- The model code is branch `v5q-gpp-fleaf2` at `d776b0e0` (gpp_fleaf_mode), plus this document. It is 90 commits on top of upstream `master` at `419a4eb3`.
- The tag `tvc-hpc-mip-S22-S21H` marks the commit that adds this document. Its source tree is identical to `d776b0e0`.
- The two runs documented here are the current MIP candidate pair:
  - **TVC S22**, Trail Valley Creek, CMT51;
  - **HPC S21H**, CMT01.
- This is a **candidate, not adopted as production**. Adoption is Elchin's decision. Production parameter files in `parameters/` are unchanged. All site edits live in copied parameter folders outside the repo (section 4).
- Policy for every switch added on the branch:
  - With the key absent or set to its default, the model is **bit-identical** to the previous step.
  - Only N1–N3 (`n_rootless_fix` 2, `n_som_move_fix` 1, `n_nonvasc_extract` 1) are on by default, since v5o. The legacy values stay selectable.
  - Every switch has unit tests and carries its references at the switch site, in `tests/README.md` and in `docs/daily-mode-REFERENCES.md`.

## 2. Lineage

```
upstream master 419a4eb3
 └─ daily-climate-mode (daily forcing, water/snow/Stefan/Richards fixes, thermal and SEB options)
     ├─ merge unit-tests        3d62b90f  (doctest suite tests/unit, regression harness, NFactorSelect seam)
     ├─ merge v5k-fixes         95ada552  (B1–B4, N1–N3 drafts)          tag v5k
     ├─ merge v5k-ctol          e1fc88d5  (C1 soil_nmin_tol_fix)         tag v5k1
     ├─ merge v5l-et            e8389ca3  (S1–S4 ET split, ET outputs)
     ├─ merge v5m-nleak         bf05797e  (L1 n_safetynet_conserve)
     ├─ merge v5n-giveup        4bb7cd99  (G1–G3, give-up diagnostics)   tag v5n
     ├─ merge v5o-ndefaults     ba204048  (N1–N3 default on)
     ├─ merge v5p-daily-split   e41a727e  (daily_c_split, output only)
     ├─ v5p water/units fixes   24d3163e, 5116575f
     └─ v5q series (linear)     c86304c7 BNF … 67bf69a6 rh_frozen_mode → d776b0e0 gpp_fleaf_mode  (this branch)
```

Commits `419a4eb3..d776b0e0`, oldest first (`git log --reverse --date=short --format='%h %ad %s' 419a4eb3..d776b0e0`):

| commit | date | subject |
|---|---|---|
| `76bc0858` | 2026-10-20 | Add opt-in daily climate forcing mode for SC stage (model_settings.climate_timestep) |
| `48053527` | 2026-10-20 | Fix snow loss: keep Snow::extramass across reset() in Ground::updateSnowHorizon |
| `3a7cc57f` | 2026-10-20 | Daily climate mode: correct monthly output time coordinate in SC stage |
| `6044df4c` | 2026-10-20 | Stefan: rebuild missing freeze/thaw fronts when the front deque is empty |
| `d2c7d35e` | 2026-10-20 | ALD: compute active layer depth from the 0 degC isotherm; keep old value as ALD_FLAG |
| `503bb21b` | 2026-10-20 | Soil_Bgc: guard N-mineralization rescale against divide-by-zero; cap root N extraction |
| `0e22205a` | 2026-10-20 | Scenario CO2: index projected CO2 by calendar year, fail loudly if not covered |
| `bb93da53` | 2026-10-20 | Daily climate mode: fail with a clear message if the daily file is shorter than the scenario |
| `de49c6d5` | 2026-10-20 | Soil_Bgc: make the root N extraction cap N-conserving |
| `ed28e049` | 2026-10-20 | Soil_Bgc: remove root N extraction caps; keep only the netnmin rescale guard |
| `40529d67` | 2026-10-20 | Soil evaporation: reset days-since-rain only on days with >= 2 mm rain |
| `333d7197` | 2026-10-20 | Water fix 1/4: fix out-of-bounds write of root_water_up in Soil_Env::updateDailySM |
| `a671ad77` | 2026-10-20 | Water fix 2/4: snow bookkeeping - extramass melts on bare warm ground; melt not removed twice from snowpack |
| `337020ac` | 2026-10-20 | Water fix 3/4: supply-limited ET - floor top-up and unapplied soil evaporation removed from reported ET |
| `28167031` | 2026-10-20 | Water fix 4/4: plant-cover weighting - transpiration scaled by fpc once; reported EET and soil uptake both use sum(tran*fpc) |
| `d6eee2c6` | 2026-10-20 | Daily output fix 1/2: daily horizon soil-T outputs (TSHLW, TDEEP, TMINEA/B/C) written at day-of-year offset (were overwriting days 1..dinm each month -> only January slots, holding December) |
| `689d8b97` | 2026-10-20 | Daily output fix 2/2: other per-month daily outputs (VWCSHLW/DEEP/MINEA-C, TC*, HK*, FRONTS*, PERCOLATION, LAYER_DRAIN, ROOT_WATER_UPTAKE) written at day-of-year offset |
| `17d53b56` | 2026-10-20 | Memory fixes (startup UB): calculate_daily_prec 1-based arrays sized DINM+1; NUL-safe time-units read in get_timeseries_start_year; cohort.md points to Runner-owned ModelData (was dangling ptr to ctor by-value arg) |
| `7a7fe61a` | 2026-10-20 | Hang guard: Richards fail-fast - throw if a layer index exceeds MAX_SOI_LAY (topind<1 or drainind>MAX_SOI_LAY) or the substep loop becomes non-finite/endless (dtsub<=0, non-finite, >1e6 substeps); the cell goes to fail_log instead of hanging (A3 r1 s008) |
| `e3fcc91c` | 2026-10-20 | Hang fix: remove snow layers emptied by melt in Ground::constructSnowLayers (after water fix 2/4 only sublimation reached the ablation branch, so on full-melt days with zero sublimation empty snow layers stayed and shifted soil indices past MAX_SOI_LAY in Richards; fixes A3 r1 s008 TVC_W_nfeed hang) |
| `e86633fe` | 2026-10-20 | Carbon option A: N-conserving per-layer root N extraction restored (from de49c6d5) without the total cap - extraction above a layer available N is moved to layers with remaining available N, total = vegetation N uptake; if column demand exceeds column avln the vegetation per-layer extraction is used unchanged. In A3 TVC nfeed tests soil integrator give-ups drop from ~190 to 0-4 months per 75 yr and the C-closure residual from -30..-39 to >-0.5 gC/m2/yr. |
| `be8a62fa` | 2026-10-20 | Carbon option B: soil integrator safety net - when Integrator::adapt gives up (state reset to month fraction t0), the soil inputs for the remaining 1-t0 of the month (litterfall C/N to RAWC/ORGN, woody debris C/N, N deposition) are delivered so C is conserved; each event logs a warning, prints a WARNING line on stderr and appends a row to <output_dir>/integrator_safety_net.csv; per-stage summary (soil safety-net events, veg integrator give-ups) on stderr |
| `228b0f08` | 2026-10-21 | DRAFT structural 1/4: snow radiation-driven melt at Tair >= -1 C routed to melt water (was booked as sublimation -> atmosphere/EET); true sublimation (Tair < -1) unchanged |
| `939790f5` | 2026-10-21 | DRAFT structural 2/4: optional infiltration into frozen soil (model_settings.frozen_infil_eice, default off): ice impedance 10^(-E*Fice)*hksat, pore space, refreezing with latent heat release |
| `08a525ba` | 2026-10-21 | DRAFT structural 3/4: restore snow-surface n-factor (remove TESTING line in Snow_Env::updateDailyM) |
| `8053cabc` | 2026-10-21 | DRAFT structural 4/4: optional soil-moisture limitation of soil evaporation (model_settings.soil_evap_beta, default off): beta=0.25(1-cos(pi*theta/theta_fc))^2 on first shallow layer |
| `4b1041fd` | 2026-10-21 | Kersten-form switch for unfrozen soil thermal conductivity: model_settings.kersten_form (0 = legacy ke = ln(Sr)+1, default, bit-identical; 1 = ke = log10(Sr)+1, Johansen fine-soil / CLM form). With ln, ke=0 (k = tcdry) for Sr <= 1/e, which makes moist organic layers near-insulating (TVC fibric k ~0.03-0.04 W/m/K in July). |
| `d9c30fc1` | 2026-10-21 | Water-leak fix: book freeze-expansion clip and organic-shrink water as subsurface drainage; keep the surface puddle across years |
| `b7e43750` | 2026-10-21 | Stefan: when today's front sweeps out every existing front, set all soil layers to the new front's state |
| `032e8f13` | 2026-10-21 | Richards: size work arrays MAX_SNW_LAY+MAX_SOI_LAY+2 and fix clearRichardsArrays off-by-one |
| `0f76919c` | 2026-10-21 | dvmdostem: return exit code 2 when any cell fails or times out |
| `fe93d93c` | 2026-10-21 | NLOST output: report the soil available-N (AVLN) loss only |
| `79668759` | 2026-10-21 | Optional moss/organic thermal conductivity: model_settings.org_tc_scheme (default 0 = off) |
| `014dda0d` | 2026-10-21 | Optional lateral drainage of thawed near-surface organic layers: model_settings.org_drain_rate (default 0 = off) |
| `3e225cec` | 2026-10-21 | Optional settable mineral thermal properties: model_settings.mineral_tcsolid / mineral_vhcsolid / mineral_tcdry_form (default off) |
| `a83e6d7c` | 2026-10-21 | Water-leak fix: comment-only cleanup of leftover DRAFT / NOT MERGED wording (v5f) |
| `6c78bd4c` | 2026-10-21 | Fix (a) (v5i): optional vapor-pressure unit fix, model_settings.vp_hpa_fix (default 0 = off) |
| `9b3a4840` | 2026-10-21 | Fix (b) (v5i): defined tdrv and optional same-day n-factor selection, model_settings.nfactor_sameday (default 0 = off) |
| `9dae595b` | 2026-10-21 | Fix (c) (v5j): optional withdrawal of non-vascular transpiration from soil water, model_settings.nonvasc_trans_withdraw (default 0 = off) |
| `d41c7f79` | 2026-10-21 | Option 4 (v5j): SEB top boundary under snow (seb_scheme 2) and Sturm 1997 snow conductivity (seb_snow_tc), defaults off |
| `ecdce425` | 2026-10-21 | tests: doctest unit-test suite (tests/unit) and make unittest target |
| `d718934d` | 2026-10-21 | seam (behaviour-neutral): n-factor selection as pure inline helpers (include/NFactorSelect.h) |
| `b990efe1` | 2026-10-21 | tests: n-factor selection tests; fix Soil_Bgc test fixture (FirData); regression harness |
| `073cf9a8` | 2026-10-21 | Fix (c) (v5j): optional withdrawal of non-vascular transpiration from soil water, model_settings.nonvasc_trans_withdraw (default 0 = off) |
| `d0c73aec` | 2026-10-21 | Option 4 (v5j): SEB top boundary under snow (seb_scheme 2) and Sturm 1997 snow conductivity (seb_snow_tc), defaults off |
| `839e4c8b` | 2026-10-21 | tests: v5j SEB (solveSEB) and Sturm snow-conductivity switch tests; test objects depend on model headers |
| `93e4b0fa` | 2026-10-21 | tests: README; snow albedo tests (xfail: hard-coded 0.9 below -10 C); harness: WARN vs FAIL (--strict), A closure columns |
| `3d62b90f` | 2026-10-22 | Merge branch unit-tests: doctest unit tests (tests/unit, make unittest), regression harness (tests/regression), behaviour-neutral seam include/NFactorSelect.h |
| `07f6c8d0` | 2026-10-22 | v5k: switchable fixes B1-B4, all default off (svp_ice_fix, richards_infil_fix, stefan_fronts_fix, snow_albedo_fix) |
| `b9e9028b` | 2026-10-22 | tests: unit tests for v5k fixes B1-B4 (switch on); legacy xfail tests kept |
| `29d72af8` | 2026-10-22 | v5k: draft fix N1 (n_rootless_fix, default 0): N mineralised in rootless soil layers left ORGN without reaching AVLN |
| `20b2723e` | 2026-10-22 | v5k: draft fixes N2 (n_som_move_fix) and N3 (n_nonvasc_extract), default 0: unreported N sink/source in SOM N movement and non-vascular N uptake |
| `1da189b4` | 2026-10-22 | tests: unit tests for draft N fixes N1-N3 (switch on) + legacy N2 xfail |
| `d3d3fbba` | 2026-10-22 | v5k: N3 draft also caps non-vascular uptake at 95% of the available N it is extracted from |
| `3357cb0b` | 2026-10-22 | tests: README for v5k switches; run_case.py REG_EXTRA_ENV for diagnostic env vars |
| `95ada552` | 2026-10-22 | Merge branch v5k-fixes: switchable fixes B1-B4 (svp_ice_fix, richards_infil_fix, stefan_fronts_fix, snow_albedo_fix) and draft N fixes (n_rootless_fix, n_som_move_fix, n_nonvasc_extract), all default 0; unit tests |
| `d427d649` | 2026-10-22 | v5k: fix C1 soil_nmin_tol_fix (default 0): signed NMIN accumulator stalls soil step control -> float C rounding gain; unit tests |
| `e1fc88d5` | 2026-10-22 | Merge branch v5k-ctol: fix C1 soil_nmin_tol_fix (default 0), soil integrator NMIN error tolerance; unit tests |
| `7446dd40` | 2026-10-22 | tests: README for v5k1 (C1 soil_nmin_tol_fix); doc only |
| `798388f5` | 2026-10-22 | v5l: switchable ET-split fixes S1-S4 (et_fpc_once_fix, nv_snow_bury_fix(+nv_snow_bury_depth), rain_int_lai_fix, nv_rain_surface), all default 0; ET-component outputs EVAP_CANOPY/SUBL_CANOPY/EVAP_MOSS/EVAP_SOIL/SUBL_SNOW/RAIN_INTERCEPT/SNOW_INTERCEPT; unit tests |
| `e8389ca3` | 2026-10-22 | Merge branch v5l-et: switchable ET-split fixes S1-S4 (all default 0, bit-identical to v5k1 when off); ET-component outputs; unit tests |
| `64d55ffc` | 2026-10-22 | tests: README for v5l (ET-split switches S1-S4, suite fixes_v5l, ET-component outputs); doc only |
| `1b36084e` | 2026-10-22 | docs+comments only (v5l1): scientific references for every daily-mode switch and fix (v5-daily ... v5l) |
| `23e8077b` | 2026-10-27 | v5m draft (switch, default off): n_safetynet_conserve - conserve N in the soil integrator safety net (L1) |
| `bf05797e` | 2026-10-27 | Merge branch v5m-nleak: switchable fix L1 n_safetynet_conserve (default 0, bit-identical to v5l1 when off): conserve N in the soil integrator safety net; unit tests |
| `7b600a01` | 2026-10-27 | docs: tests/README and daily-mode-REFERENCES for v5m (L1 n_safetynet_conserve merged, default off) |
| `904a4a7c` | 2026-10-27 | v5n draft (switches, all default off): G1 n_extract_limiter (+n_extract_tau), G2 soil_giveup_retry, soil_giveup_diag; unit tests |
| `5c597aec` | 2026-10-27 | v5n draft: G3 n_netnmin_rescale_guard, give-up diagnostics v2, n_extract_tau default 0.2; docs (all switches default off) |
| `4bb7cd99` | 2026-10-28 | Merge branch v5n-giveup: switchable soil give-up fixes, all default off (bit-identical to v5m when off) |
| `e5afbe37` | 2026-10-28 | docs: v5n merged (G1 n_extract_limiter + n_extract_tau, G2 soil_giveup_retry, G3 n_netnmin_rescale_guard, soil_giveup_diag; all default off) |
| `5812ff19` | 2026-10-28 | v5o: N conservation fixes N1-N3 on by default (n_rootless_fix=2, n_som_move_fix=1, n_nonvasc_extract=1) |
| `273f9440` | 2026-10-28 | v5o docs: N1-N3 default-on section, closure rule, regression proof; REFERENCES v5o entry (McGuire 1992 doi:10.1029/92GB00219, Raich 1991 doi:10.2307/1941899, Turetsky 2003 doi:10.1639/05; Crossref 2026-10-27) |
| `ba204048` | 2026-10-28 | Merge v5o-ndefaults: N conservation fixes N1-N3 on by default (2/1/1; legacy 0 selectable, bit-identical to v5n); closure.py enforces N only with all three on |
| `a0190562` | 2026-10-28 | v5p: output-only daily split of the monthly C fluxes (daily_c_split, default 0) |
| `e41a727e` | 2026-10-28 | Merge v5p-daily-split: output-only daily split of the monthly C fluxes (daily_c_split, default 0; off bit-identical to v5o) |
| `24d3163e` | 2026-10-28 | output_spec: TRANSPIRATION units mm/day -> mm/time (metadata only) |
| `5116575f` | 2026-10-28 | v5p-daily-water-fix: daily TRANSPIRATION/EET/PET outputs equal the terms summed into monthly (output only) |
| `c86304c7` | 2026-10-28 | v5q Phase 1: biological N fixation into soil AVLN (model_settings.bnf_scheme, default 0 = off) |
| `128241d6` | 2026-10-28 | v5q Phase 2: symbiotic BNF to plant labile N (model_settings.bnf_to_plant, default 0 = off) |
| `46e425ab` | 2026-10-28 | v5q-ndrain: n_leach_fix switch - soil AVLN pools lose exactly the booked NLOST |
| `13cf9f04` | 2026-10-28 | v5q-negpool: soil_negpool_fix switch - no negative soil N pools after the RKF45 step |
| `d50f755b` | 2026-10-28 | v5q-d2m-gate: model_settings.d2m_transfer_gate (0 legacy / 1 thaw-gated / 2 off) |
| `982d79a3` | 2026-10-29 | v5q-layern-fix: layer N follows C in the dynamic soil-layer remap (layer_n_remap_fix, default 0 = legacy) |
| `0c11be6c` | 2026-10-30 | v5q-layer-bookkeep: dynamic soil-layer C/N bookkeeping fixes (layer_bookkeep_fix, default 0 = legacy) |
| `9ece944b` | 2026-10-30 | v5q-thawfrac-decomp: decomposition modifier weighted by layer thawed fraction (thawfrac_decomp, default 0 = legacy) |
| `99ae7e0c` | 2026-10-30 | v5q-leach-curn: N drainage leaching from the current drainage-zone AVLN (n_leach_curn_fix, default 0) |
| `6630518f` | 2026-10-30 | v5q-thawfrac-lmoist: optional layer-mean moisture for the thawed part (thawfrac_moist_mode, default 0 = unchanged) |
| `df47179d` | 2026-10-30 | v5q-humic-stefan: humic sublayer cap, Stefan front-layer resistance fix, ALD_FRONT output (all default off) |
| `67bf69a6` | 2026-11-01 | v5q-rh-frozen: rh_frozen_mode, frozen-part RH from unfrozen water above the permafrost table (default 0 = bit-identical) |
| `d776b0e0` | 2026-11-02 | v5q-gpp-fleaf: gpp_fleaf_mode, TEM leaf-phenology factor in GPP (default 0 = bit-identical) |

## 3. model_settings switches added on this branch, and the values used

- **Default** is the binary default from `src/ModelData.cpp` at `d776b0e0`.
- **TVC S22** and **HPC S21H** are read from the run configs:
  - `step22/runs/main/S22_n150/<case>/main/config/config.js`
  - `step21/runs/main/S21H_n150/<case>/main/config/config.js`
- The model_settings are identical across the 4 cases of each site. The cases differ only in forcing (Warming/NoWarming) and in `stage_settings.*.nfeed` (nfeed/nonfeed).
- "absent (=default)" means the key is not in the config, so the binary default applies.
- References are copied from `tests/README.md` (section "References for every daily-mode switch" and the per-version tables). Full citations and URLs are in `docs/daily-mode-REFERENCES.md`. No citation was added for this document.

Other model_settings keys in both configs (not added on this branch): `cell_timelimit` 0, `dynamic_lai` 1, `baseline_start` 1901, `baseline_end` 1931. In the nfeed cases `stage_settings.*.nfeed` is true in 4 of the 5 stage blocks; in the nonfeed cases it is false in all 5.

The following keys are absent in both configs, so their defaults apply:
- the BNF sub-parameters `bnf_pro_c` 1.8, `bnf_pro_d` −0.003, `bnf_aet_a` 0.00234, `bnf_aet_b` −0.0172, `bnf_const` 0, `bnf_asym` 0 and `bnf_avln_cap` 0 (no cap);
- `n_extract_tau` 0.2;
- the SEB extras `seb_lw_source` 0, `seb_diag` 0, `seb_z0_snow` 0.001, `seb_z0_soil` 0.01, `seb_zwind` 10, `seb_ztemp` 2, `seb_emis_snow` 0.99, `seb_emis_soil` 0.97, `seb_wind` −1, `seb_minwind` 0.5, `seb_stab_b` 10, `seb_tran_frac` 1 and `seb_pres` 100000;
- `org_tc_depth` 10.

| switch | default | TVC S22 | HPC S21H | what it does | references (from tests/README.md) | added in |
|---|---|---|---|---|---|---|
| `climate_timestep` | "monthly" | "daily" | "daily" | "daily" = drive the SC (and spin/main) stages with daily forcing instead of monthly | no published source for the mode; CF 1.8 365_day calendar; daily Ra: FAO-56 eq. 21 | `76bc0858` |
| `frozen_infil_eice` | -1 (off) | 6 | 6 | ≥ 0 = infiltration into a frozen top layer with ice impedance 10^(−E·F_ice), E = value (off: no infiltration while frozen) | Swenson et al. 2012 doi:10.1029/2012MS000165; CLM5 eq. 2.7.48 | `939790f5` |
| `soil_evap_beta` | 0 | 1 | 1 | 1 = beta-function limit of soil evaporation by top-layer moisture | Lee & Pielke 1992 doi:10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2; Clapp & Hornberger 1978 doi:10.1029/WR014i004p00601 | `8053cabc` |
| `kersten_form` | 0 | 1 | 1 | unfrozen Kersten number: 0 = ln(Sr)+1 (legacy), 1 = log10(Sr)+1 | Johansen 1975 doi:10.21236/ADA044002; Lawrence & Slater 2008 doi:10.1007/s00382-007-0278-1 | `4b1041fd` |
| `org_tc_scheme` | 0 | 1 | 1 | 1 = Côté & Konrad organic-soil thermal conductivity (uses org_tc_kunf/kfrz/dry/thre_f/thre_h/depth) | Côté & Konrad 2005 doi:10.1139/t04-106; Letts et al. 2000 doi:10.1080/07055900.2000.9649643; Chadburn et al. 2015 doi:10.5194/gmd-8-1493-2015 | `79668759` |
| `org_tc_kunf` | 0.6 | 0.3547 | 0.27 | org_tc_scheme 1: unfrozen conductivity parameter | as org_tc_scheme | `79668759` |
| `org_tc_kfrz` | 0.25 | 0.3764 | 0.5134 | org_tc_scheme 1: frozen conductivity term | as org_tc_scheme | `79668759` |
| `org_tc_dry` | -1 | 0.0524 | 0.08571 | org_tc_scheme 1: dry conductivity (−1 = built-in) | as org_tc_scheme | `79668759` |
| `org_tc_thre_f` | 0.275 | 0.1931 | 0.1795 | org_tc_scheme 1: fibric threshold | as org_tc_scheme | `79668759` |
| `org_tc_thre_h` | 0.62 | 0.6977 | 0.5638 | org_tc_scheme 1: humic threshold | as org_tc_scheme | `79668759` |
| `org_drain_rate` | 0 | 0.2337 | 0.2661 | > 0 = lateral fill-and-spill drainage (1/day) of thawed fibric/humic layers with top < org_drain_depth towards theta_re, booked as qdrain | Quinton & Marsh 1999 doi:10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D; Letts et al. 2000 | `014dda0d` |
| `org_drain_thre_f` | 0.275 | 0.3133 | 0.2622 | org_drain_rate: fibric threshold | as org_drain_rate | `014dda0d` |
| `org_drain_thre_h` | 0.62 | 0.6584 | 0.4919 | org_drain_rate: humic threshold | as org_drain_rate | `014dda0d` |
| `org_drain_depth` | 0.3 | absent (=default) | 0.3103 | org_drain_rate: depth (m) | as org_drain_rate | `014dda0d` |
| `mineral_tcsolid` | -1 (off) | 3.109 | 4.044 | mineral solid conductivity override (W/m/K) | Johansen 1975; Farouki 1981 doi:10.21236/ADA111734; CLM5 Tech Note, soil temperature | `3e225cec` |
| `mineral_vhcsolid` | -1 (off) | absent (=default) | 1834000.0 | mineral solid volumetric heat capacity override (J/m3/K) | as mineral_tcsolid | `3e225cec` |
| `mineral_tcdry_form` | 0 | 1 | 1 | mineral dry conductivity: 0 = legacy geometric mean, 1 = alternative form | as mineral_tcsolid | `3e225cec` |
| `vp_hpa_fix` | 0 | 1 | 1 | 1 = vapour pressure input treated in hPa (VPD units fix) | Harris et al. 2020 doi:10.1038/s41597-020-0453-3; FAO-56; Monteith 1965 (PMID 5321565) | `6c78bd4c` |
| `nfactor_sameday` | 0 | 1 | 1 | 1 = n-factor selector from the same day (fixes an uninitialised selector) | Klene et al. 2001 doi:10.1080/15230430.2001.12003416 | `9b3a4840` |
| `nonvasc_trans_withdraw` | 0 | 1 | 1 | 1 = non-vascular transpiration withdrawn from soil water (water mass conservation) | no published source; context Heijmans et al. 2004 doi:10.1016/j.agrformet.2003.12.006 | `9dae595b` |
| `seb_scheme` | 0 | 2 | 2 | 2 = daily skin energy balance on snow days (n-factor boundary on snow-free days); daily mode only; seb_* extras at defaults | Brutsaert 1982 doi:10.1007/978-94-017-1497-6; Louis 1979 doi:10.1007/BF00117978; Dyer 1974 doi:10.1007/BF00240838; Warren 1982 doi:10.1029/RG020i001p00067; WMO-No. 8 (2008); CLM5 Tech Note, Fluxes | `d41c7f79` |
| `seb_snow_tc` | 0 | 1 | 1 | 1 = Sturm et al. 1997 snow thermal conductivity | Sturm et al. 1997 doi:10.3189/S0022143000002781; legacy Goodrich 1982 doi:10.1139/t82-047 | `d41c7f79` |
| `svp_ice_fix` | 0 | 1 | 1 | B1: WMO ice saturation vapour pressure at t ≤ 0 | WMO-No. 8 (2008) Annex 4.B; Buck 1981 doi:10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2; Murphy & Koop 2005 doi:10.1256/qj.04.94 | `07f6c8d0` |
| `richards_infil_fix` | 0 | 1 | 1 | B2: Richards early-return infiltration sent to puddle/runoff (not lost) | no published source; code bug; Richards 1931 doi:10.1063/1.1745010 | `07f6c8d0` |
| `stefan_fronts_fix` | 0 | 0 | 0 | B3: Stefan::combineExtraFronts mindz reset per pass (hang fix) | no published source; Woo et al. 2004 doi:10.1029/2004GL019475 | `07f6c8d0` |
| `snow_albedo_fix` | 0 | 1 | 1 | B4: snow albedo below −10 C = snwalbmax (continuity) | no published source; cf. Warren 1982 doi:10.1029/RG020i001p00067 | `07f6c8d0` |
| `n_rootless_fix` | 2 (since v5o) | 2 | 2 | N1: mineralised N of rootless layers to AVLN (2) | no published source; N mass conservation | `29d72af8` |
| `n_som_move_fix` | 1 (since v5o) | 1 | 1 | N2: SOM N moved shallow→deep added, not overwritten | no published source; N mass conservation | `20b2723e` |
| `n_nonvasc_extract` | 1 (since v5o) | 1 | 1 | N3: non-vascular N uptake extracted from moss/fibric AVLN (95 % cap) | no published source for the scheme; Turetsky 2003 doi:10.1639/05 | `20b2723e, d3d3fbba` |
| `soil_nmin_tol_fix` | 0 | 1 | 1 | C1: soil step-control tolerance ptol·max(|NMIN|,|NIMMOB|,|AVLN|) | no published source for the fix; Press et al. 2007 §17.2; Hairer et al. 1993 doi:10.1007/978-3-540-78862-1 | `d427d649` |
| `et_fpc_once_fix` | 0 | 1 | 1 | S1: transpiration cover-scaled once (not twice) | Gerten et al. 2004 doi:10.1016/j.jhydrol.2003.09.029; Schaphoff et al. 2018 doi:10.5194/gmd-11-1343-2018; CLM5 Tech Note eq. 2.7.68; Monsi & Saeki 1953 doi:10.1093/aob/mci052 | `798388f5` |
| `nv_snow_bury_fix` | 0 | 1 | 1 | S2: no snow interception by moss/lichen under the snowpack | Hedstrom & Pomeroy 1998 doi:10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4; Sturm et al. 2001 doi:10.1175/1520-0442(2001)014<0336:SSIIAT>2.0.CO;2 | `798388f5` |
| `nv_snow_bury_depth` | 0.0 | 0.0 | 0.0 | S2 threshold (m) | Yi et al. 2009 doi:10.1029/2008JG000841 | `798388f5` |
| `rain_int_lai_fix` | 0 | 1 | 1 | S3: rain interception capacity = ircoef·LAI with storage | Dickinson et al. 1993 doi:10.5065/D67W6959 / CLM5 eq. 2.7.12; van der Knijff et al. 2010 doi:10.1080/13658810802549154; Zwieback et al. 2019 doi:10.1088/1748-9326/ab1049; Chang et al. 2022 doi:10.1016/j.rse.2021.112747 | `798388f5` |
| `nv_rain_surface` | 0 | 1 | 1 | S4: moss/lichen-intercepted rain booked as EVAP_MOSS (1) | Chadburn et al. 2015 doi:10.5194/gmd-8-1493-2015; Yi et al. 2009; Heijmans et al. 2004 | `798388f5` |
| `n_safetynet_conserve` | 0 | 1 | 1 | L1: soil integrator safety net debits the remaining root N extraction | no published source for the fix; McGuire et al. 1992 doi:10.1029/92GB00219; Raich et al. 1991 doi:10.2307/1941899; Shampine et al. 2005 doi:10.1016/j.amc.2004.12.011; Thornton & Rosenbloom 2005 doi:10.1016/j.ecolmodel.2005.04.008 | `23e8077b` |
| `n_extract_limiter` | 0 | absent (=default) | absent (=default) | G1: AVLN-limited root N extraction (rest from ORGN), with n_extract_tau (0.2) | Raich et al. 1991; McGuire et al. 1992; Thornton et al. 2007 doi:10.1029/2006GB002868; Näsholm et al. 1998 doi:10.1038/31921; Shampine et al. 2005 | `904a4a7c` |
| `soil_giveup_retry` | 0 | absent (=default) | absent (=default) | G2: re-integrate a given-up soil month with the G1 limiter | as n_extract_limiter | `904a4a7c` |
| `soil_giveup_diag` | 0 | absent (=default) | absent (=default) | diagnostics only (GIVEUP_DIAG, integrator_giveup_diag.csv) | none (diagnostic output) | `904a4a7c` |
| `n_netnmin_rescale_guard` | 0 | absent (=default) | absent (=default) | G3: skip a negative netnmin rescale factor | no published source for the fix; Schimel & Bennett 2004 doi:10.1890/03-8002; Thornton et al. 2007 | `5c597aec` |
| `daily_c_split` | 0 | absent (=default) | absent (=default) | output only: daily split *_DS of the monthly C fluxes + DS_FALLBACK | no published source for this exact method; Olsen & Randerson 2004 doi:10.1029/2003JD003968; Fisher et al. 2016 doi:10.5194/bg-13-4271-2016; Denton 1971 doi:10.1080/01621459.1971.10482227; Chow & Lin 1971 doi:10.2307/1928739 | `a0190562` |
| `bnf_scheme` | 0 | 2 | 2 | BNF to top-layer AVLN from previous-year rate; 2 = AET: B = max(0, a·ET + b) (bnf_aet_a 0.00234, bnf_aet_b −0.0172; no cap) | Meyerholt et al. 2016 doi:10.5194/bg-13-1491-2016; Cleveland et al. 1999 doi:10.1029/1999GB900014; Thornton et al. 2007; Davies-Barnard & Friedlingstein 2020 doi:10.1029/2019GB006387; taper form: no published source | `c86304c7` |
| `bnf_to_plant` | 0 | 0 | 0 | 1 = scheme BNF to PFT labile N (NPP split) | as bnf_scheme; PFT split: no published source | `128241d6` |
| `n_leach_fix` | 0 | 1 | 1 | AVLN pools lose exactly the booked avlnlost (drainage-zone layers) | no published source (code bug; N mass conservation); McNamara et al. 2008 doi:10.1002/hyp.6920; Khosh et al. 2017 doi:10.1002/2016jg003377 | `46e425ab` |
| `soil_negpool_fix` | 0 | 1 | 1 | negative soil N pools clamped with the deficit booked to n_debt | no published source for the fix; Shampine et al. 2005 doi:10.1016/j.amc.2004.12.011 | `13cf9f04` |
| `d2m_transfer_gate` | 0 | 1 (spin phase A: 2) | 1 (spin phase A: 2) | SOM transfer to mineral soil gated by ALD (1); 2 = transfer off (spin phase A) | Koven et al. 2009 doi:10.1029/2009GL040150; Koven et al. 2013 doi:10.5194/bg-10-7109-2013; Lawrence et al. 2019 doi:10.1029/2018MS001583; Melton et al. 2019 doi:10.5194/gmd-12-4443-2019; f-scaling: no published source | `d50f755b` |
| `layer_n_remap_fix` | 0 | 1 | 1 | layer orgn/avln remapped with the C pools on re-division | Yi et al. 2009 doi:10.1029/2008JG000841; Yi et al. 2010 doi:10.1029/2010JG001302; remap rule: no published source | `982d79a3` |
| `layer_bookkeep_fix` | 0 | 0 | 0 | layer bookkeeping conservation fixes | Yi et al. 2009; Yi et al. 2010; the fixes: no published source | `0c11be6c` |
| `thawfrac_decomp` | 0 | 1 | 1 | decomposition per layer blended by thawed fraction f from the Stefan fronts | Koven et al. 2013; Lawrence et al. 2019; Yi et al. 2009; linear blend: no published source | `9ece944b` |
| `n_leach_curn_fix` | 0 | 1 | 1 | drainage leaching uses current drainage-zone AVLN (needs n_leach_fix 1) | no published source for the fix; Shampine et al. 2005; Burchard et al. 2003 doi:10.1016/S0168-9274(03)00101-6; McNamara et al. 2008; Khosh et al. 2017 | `99ae7e0c` |
| `thawfrac_moist_mode` | 0 | 0 | 0 | 1 = layer-mean moisture in the thawed part (needs thawfrac_decomp 1) | Koven et al. 2013; Lawrence et al. 2019; Yi et al. 2009; mode choice: no published source | `6630518f` |
| `humic_max_thick` | 0 (off) | 0.15 | 0.15 | cap (m) on humic sublayer thickness (≤ 3 sublayers) | Yi et al. 2009; Yi et al. 2010; the cap: no published source | `df47179d` |
| `stefan_front_resist_fix` | 0 | 1 | 1 | thawing front resistance counts only the thawed part of the front layer | Woo et al. 2004 doi:10.1029/2004GL019475; Yi et al. 2009; the correction: no published source | `df47179d` |
| `rh_frozen_mode` | 0 | 0 | 0 | 1 = SJ2016 frozen-soil RH above the permafrost table (needs thawfrac_decomp 1) | Schaefer & Jafarov 2016 doi:10.5194/bg-13-1991-2016; Romanovsky & Osterkamp 2000; Nicolsky et al. 2009 doi:10.1016/j.coldregions.2008.03.003; Mikan et al. 2002 doi:10.1016/S0038-0717(02)00168-2; Natali et al. 2019 doi:10.1038/s41558-019-0592-8; table restriction: no published source | `67bf69a6` |
| `gpp_fleaf_mode` | 0 | 1 | absent (=default) | 1 = INGPP of every PFT × clamp(fleaf,0,1); 2 = non-woody PFTs only | Raich et al. 1991 doi:10.2307/1941899; McGuire et al. 1992 doi:10.1029/92GB00219; use for moss/lichen/sedge: no published source | `d776b0e0` |

Notes:
- `gpp_fleaf_mode` is absent at HPC because S21H was run with the `67bf69a6` binary, which predates the switch. On `d776b0e0` the absent key means 0. With 0 the model is bit-identical to `67bf69a6`: S22 bit-identity test, 86 files / 287 variables, 0 differing.
- `d2m_transfer_gate` is 1 in the mains. The spin phase A uses 2 (transfer off) until ALD is stable; phase B uses 1.
- `daily_c_split` is 0 in the documented mains. The *_DS outputs come from separate reruns with `daily_c_split` = 1 (section 7).

## 4. Parameter changes vs the repo parameter files

The run parameter folders are copies outside the repo. The repo `parameters/` files are unchanged; this branch changes no parameter file (`git diff 419a4eb3 d776b0e0 -- parameters` is empty). The tables compare the run folders with `parameters/` at `d776b0e0`, value by value, for the CMT each site uses.

**TVC S22.** Folder `step22/params/K075RN_C105/<case>`, identical for the 4 cases (folder md5 2fe29de4…). CMT51, Tussock Tundra TVC.

| file | parameter | repo `parameters/` | TVC S22 |
|---|---|---|---|
| cmt_calparbgc.txt | cmax (Sedges/Sphagnum/Lichen) | 436.780059 / 284.232584 / 84.321964 | **226.57131 / 171.68256 / 41.96052** (S22 = S20 K075RN × 1.05; S20: 215.7822 / 163.5072 / 39.9624) |
| cmt_calparbgc.txt | krb(0) | −0.756569 / −0.577634 / −0.311025 | −1.506569 / −1.327634 / −1.061025 |
| cmt_calparbgc.txt | cfall(2) (Sedges) | 0.232827 | 0.03 |
| cmt_calparbgc.txt | nfall(2) (Sedges) | 0.006953 | 0.000896 |
| cmt_bgcsoil.txt | fnloss | 0.0 | 0.10 |
| cmt_bgcsoil.txt | initminec | 38120 | 44800 |
| cmt_bgcsoil.txt | initsoln | 2788 | 3121.45 |
| cmt_bgcvegetation.txt | knuptake (Sedges) | 0.0042 | 0.008 |
| cmt_bgcvegetation.txt | cnevenbl (Sphagnum/Lichen) | 22.58 / 29.79 | 37 / 50 |
| cmt_dimvegetation.txt | klai (Sedges) | 0.15 | 1.0 |
| cmt_envcanopy.txt | ircoef (Sedges) | 0.041 | 0.2 |
| cmt_dimground.txt | snwdennew / snwdenmax | 230 / 300 | 150 / 268 |
| cmt_envground.txt | tcsolid(f) / tcsolid(h) | 0.33 / 0.4 | 0.2229 / 0.2229 |
| cmt_envground.txt | porosity(f) / porosity(h) | 0.85 / 0.8 | 0.881 / 0.805 |
| cmt_envground.txt | hksat(f) | 0.001 | 0.0003908 |
| cmt_envground.txt | nfactor(s) / nfactor(w) | 1.5 / 1.0 | 0.78 / 0.65 |

S22 vs S20 (K075RN) differs only in the CMT51 cmax row. cmt_calparbgc md5 went from 17e240d0… to 226395821a8f1ea6adc786c40771e991, and the folder md5 from 999eb938… to 2fe29de4….

**HPC S21H.** Folder `step21/params/S21H/<case>`, identical for the 4 cases (folder md5 16238b8115e890c8a438231b5f414b28). It is byte-identical to `step9/runs/spinH/<case>/par` (H9) and carries no TVC edits. CMT01.

| file | parameter | repo `parameters/` | HPC S21H |
|---|---|---|---|
| cmt_bgcsoil.txt | fsompr | 0.6112 | 0.5683 |
| cmt_bgcsoil.txt | fsomcr | 0.0171 | 0.060 |
| cmt_bgcsoil.txt | fnloss | 0.0 | 0.030 |
| cmt_bgcvegetation.txt | knuptake (PFT0–3) | 0.0042 ×4 | 0.012 ×4 |
| cmt_bgcvegetation.txt | cnevenbl (PFT0) | 46.23 | 72 |
| cmt_calparbgc.txt | micbnup | 1.168 | 0.1 |
| cmt_dimground.txt | snwdennew / snwdenmax | 120 / 250 | 70 / 210 |
| cmt_dimground.txt | coefshlwa / coefdeepa | 0.012 / 0.008 | 0.01675803 / 0.01237385 |
| cmt_envcanopy.txt | ircoef (PFT0–2) | 0.041 / 0.041 / 0.041 | 0.25 / 0.1 / 0.1 |
| cmt_envground.txt | tcsolid(f) / tcsolid(h) | 0.25 / 0.25 | 0.5185 / 0.5185 |
| cmt_envground.txt | porosity(f) / porosity(h) | 0.95 / 0.8 | 0.9348 / 0.7159 |
| cmt_envground.txt | nfactor(s) | 1.5 | 0.395 |

The run folders hold a subset of CMTs, CMT00–CMT74. Both also differ from the repo in CMT50 (initminec, initsoln) and CMT65 (micbnup, kdc*), which neither site uses. In the HPC folder, CMT51 initminec/initsoln also differ.

## 5. Which binary and commit ran each site

| site | run | binary | md5 | built from |
|---|---|---|---|---|
| TVC | S22 spins + n150 mains | `step22/bins/bin_v5q_gppfleaf2` | 564b48f4cda2125b4c565c5ad9f53f58 | `d776b0e0` (this branch) |
| HPC | S21H spins + n150 mains | `step21/bins/bin_v5q_rhfrozen` | 60570d4bf8a4d365f08a89de754a53ea | `67bf69a6` (parent of `d776b0e0`) |

- The S22 bit-identity check ran `d776b0e0` with `gpp_fleaf_mode` off against the S20 binary `60570d4b`: 86 files / 287 variables, 0 differing. So the HPC settings are reproducible on this tag with `gpp_fleaf_mode` absent or 0.
- This was inferred from the TVC test. HPC itself was not rerun on `d776b0e0`.
- Both binaries contain `daily_c_split`.

## 6. Unit tests and how to run them

- The suite is doctest, in `tests/unit`. The per-switch test files are listed in `tests/README.md`.
- Build and run: `make unittest`, or with the conda environment used here, `./build-conda.sh -j16 unittest`. The binary is `tests/unit/dvmdostem_unittest`.
- The regression harness is `tests/regression/run_regression.sh`.
- Result at `d776b0e0`, from a clean build in a fresh clone on the VM, 2026-11-02:
  - `[doctest] test cases: 269 | 269 passed | 0 failed`;
  - `assertions: 315608 | 315602 passed | 6 failed`;
  - `Status: SUCCESS!`.
- The 6 failed assertions are the documented legacy-behaviour xfail cases (`[KNOWN BUG: xfail]`, "Failed as expected"). They are in test_climate svp (2), test_seb snow albedo (2) and test_stefan_water (2: combineExtraFronts, Richards early return).
- Baseline `67bf69a6`: 265 cases / 315571 assertions. The gpp_fleaf tests (`test_gppfleaf_v5q.cpp`, 4 cases) are the difference.

## 7. Run layout and outputs (box `/workspace/bnf`, step21 and step22 also exist on the VM under `~/box_archive/bnf`)

**Cases.** `<SITE>_{Warming,NoWarming}_{nfeed,nonfeed}`. The forcing is daily and SEB, from `step4/inputs/forcing_seb/<case>_projected-climate_seb.nc`. Mains use `-p 0 -e 0 -s 0 -t 0 -n 150` (1950–2099) from the spin restart.

**TVC S22.**
- Mains: `step22/runs/main/S22_n150/<case>/{meta.json, main/{config,input,parameters,output,run.log}}`.
- Spins: `step22/runs/spinP2/<case>`, driver `phase22.py`.
- Scripts: `launch22.sh → lanefull22.sh → lanep22.sh → mkjobs22.py → run_main22.sh`.
- Analysis: `step22/data`, plots `mipplots_s22`. Report: `V5Q_TVC_S22.md`.

**HPC S21H.**
- Mains: `step21/runs/main/S21H_n150/<case>`.
- Spins: `step21/runs/spinP2/<case>`.
- Scripts: `lanep21h.sh`, `phase21.py`, `mkjobs21h.py`, `run_main21.sh`.
- Analysis: `step21/data`, plots `mipplots_s21h`. Report: `V5Q_HPC_S21H.md`.

**Spin stop rule (relcyc).**
- Thresholds: relative drift < 1e-5/yr (SOC) and < 1e-4/yr (VEGC), |dN| < 0.01, on 20-yr cycle means. Cap 6000 yr. DEEPC is report-only.
- Phase A: gate 2 until ALD is stable (cap 2000). Phase B: gate 1, 600-yr launches.

**Main outputs.** These are `main/output/*_sc.nc`, 82–83 files.
- Monthly C, N and water fluxes and pools, including `BNF`, `NLOST`, `EVAP_MOSS`, `TRANSPIRATION`.
- Yearly `ALD`, `ALD_FLAG`, `ALD_FRONT`, `VEGC`, `SOM*`.
- Daily `SWE`, `TSHLW`, `TDEEP`, `TMINE*`, `VWC*`.
- `integrator_safety_net.csv`.

**Daily C split (*_DS).** These come from reruns of the same mains with `daily_c_split` = 1:
- Location: `step22/runs/main/S22_n150_DS/<case>/main` and `step21/runs/main/S21H_n150_DS/<case>/main`.
- Method: README `_ds_rerun/README.md`. Same binary, same restart, the 11 DS rows appended to output_spec.csv, spins not rerun.
- Files are `<V>_DS_daily_sc.nc`, for V = GPP, NPP, RM, RG, RHSOM, RHDWD, RECO, NEE, INGPP, INNPP.
  - Units: g/m²/day, 54750 days (1950–2099, 365-day calendar).
  - Relations: NEE_DS = RECO_DS − GPP_DS; RH = RHSOM_DS + RHDWD_DS; RA = RM_DS + RG_DS.
- `DS_FALLBACK_monthly_sc.nc` holds the per-month fallback flag bits.
- `verify.py` gives ALL PASS for all 8 cases:
  - every shared file is data-identical to the non-DS main (only `run_status.total_runtime` differs);
  - max |Σ daily − monthly| ≤ 8.53e-14 g/m² (relative ≤ ~1e-15) over all variables and cases;
  - no NaN.

## 8. Key results (copied verbatim from the box reports `/workspace/bnf/V5Q_TVC_S22.md` §4–§7 and `/workspace/bnf/V5Q_HPC_S21H.md` §3–§7)

Units: fluxes gC m⁻² yr⁻¹, stocks gC m⁻², ALD m, BNF/NLOST/dNtot/dt gN m⁻² yr⁻¹.

### TVC S22: 1981–2010 and 2080–99 (cells are S22 / S20)

| case | window | GPP | RECO | RH | NEE | VEGC (non-woody = total at TVC) | SOC | ALD / ALD_FRONT | BNF | NLOST | dNtot/dt |
|---|---|---|---|---|---|---|---|---|---|---|---|
| W_nfeed | 1981–2010 | 297.5 / 346.5 | 295.6 / 342.7 | 113.0 / 129.7 | −1.89 / −3.83 | 207 / 241 | 55635 / 58164 | 0.557 / 0.557 (S20 0.553) | 0.069 / 0.080 | 0.090 / 0.095 | −0.024 / −0.017 |
| W_nfeed | 2080–99 | 397.9 / 456.7 | 390.0 / 452.5 | 133.0 / 155.8 | −7.91 / −4.24 | 250 / 288 | 55536 / 58011 | 0.666 / 0.662 (S20 0.660) | 0.086 / 0.096 | 0.137 / 0.142 | −0.051 / −0.044 |
| NW_nfeed | 1981–2010 | 239.2 / 283.4 | 241.7 / 282.5 | 98.3 / 112.5 | +2.49 / −0.91 | 171 / 201 | 54497 / 57438 | 0.539 / 0.532 (0.527) | 0.051 / 0.064 | 0.073 / 0.085 | −0.021 / −0.020 |
| NW_nfeed | 2080–99 | 235.8 / 301.1 | 233.6 / 298.8 | 91.9 / 117.7 | −2.21 / −2.33 | 161 / 204 | 54485 / 57614 | 0.538 / 0.533 (0.533) | 0.050 / 0.066 | 0.071 / 0.088 | −0.024 / −0.024 |
| W_nonfeed | 1981–2010 | 289.1 / 337.8 | 288.9 / 336.2 | 111.3 / 128.8 | −0.13 / −1.57 | 202 / 236 | 55563 / 57728 | 0.558 / 0.559 (0.548) | 0 / 0 | 0 / 0 | 0 / 0 |
| W_nonfeed | 2080–99 | 391.1 / 448.6 | 387.4 / 443.2 | 135.3 / 152.4 | −3.67 / −5.47 | 246 / 284 | 55697 / 57943 | 0.716 / 0.689 (0.681) | 0 / 0 | 0 / 0 | 0 / 0 |
| NW_nonfeed | 1981–2010 | 234.0 / 275.8 | 232.2 / 274.5 | 92.2 / 109.2 | −1.83 / −1.33 | 167 / 196 | 54737 / 56872 | 0.532 / 0.531 (0.527) | 0 / 0 | 0 / 0 | 0 / 0 |
| NW_nonfeed | 2080–99 | 237.3 / 293.7 | 236.1 / 292.2 | 93.6 / 115.7 | −1.21 / −1.59 | 163 / 200 | 54869 / 57042 | 0.539 / 0.532 (0.527) | 0 / 0 | 0 / 0 | 0 / 0 |

### TVC S22: AmeriFlux CA-TVC 2017 (Warming mains; "cov" = weighted by the tower's valid-day coverage)

| quantity | **S22 W_nfeed** | S20 W_nfeed | S22 W_nonfeed | S20 W_nonfeed | tower NT / DT |
|---|---|---|---|---|---|
| annual GPP | **312.3** | 344.7 | 299.9 | 328.0 | 275.2 / 295.8 |
| annual RECO (cov) | **328.7 (321.9)** | 368.0 (360.4) | 335.8 (329.3) | 367.8 (360.6) | 289.1 / 252.6 |
| annual NEE (cov) | **+16.3 (+9.6)** | +23.3 (+15.7) | +35.9 (+29.3) | +39.8 (+32.6) | +6.1 (USTAR 5–95 %: −17.0…+55.8) |
| GPP Jun / Jul / Aug / Sep | **75.0 / 123.6 / 94.4 / 12.9** | 95.9 / 117.2 / 93.6 / 19.1 | 62.3 / 122.4 / 92.8 / 15.3 | 80.2 / 116.1 / 91.8 / 21.2 | NT 58.4 / 112.6 / 86.9 / 24.4; DT 49.2 / 123.8 / 77.0 / 26.5 |
| Jun/Jul GPP ratio | **0.61** | 0.82 | 0.51 | 0.69 | 0.52 / 0.40 |
| RECO Jun / Jul / Aug / Sep | 53.6 / 72.0 / 63.3 / 41.9 | 64.0 / 77.7 / 68.8 / 47.7 | 61.5 / 69.0 / 58.2 / 38.6 | 69.6 / 73.6 / 62.7 / 42.2 | NT 35.2 / 55.9 / 84.4 / 48.5 |
| cold-season NEE 2017 | **66.8** | 75.1 | 66.9 | 73.8 | 65.5 |
| cold-season NEE 2005–24 mean [min–max] | **49.6 [37.8–66.8]** (below band) | 58.9 [42.8–81.7] (in) | 46.9 [34.7–66.9] (below) | 54.4 [39.9–76.4] (in) | band 51.9–106.1 |

### TVC S22: spins

| case | A_stable Y (S20) | B launches / yrs | relcyc met | cap verdict (proj150_cyc_ok) | proj150_cyc SOC / VEGC / Ntot / MINEC % | DEEPC %/150 yr (report-only; S20) | B_dDEEPC_rising |
|---|---|---|---|---|---|---|---|
| W_nfeed | 160 (240) | 10 / 6000 (cap) | no | ok | 0.155 / 0.004 / 0.101 / 0.063 | 0.79 (0.84–0.92 range) | yes |
| NW_nfeed | 140 (140) | 10 / 6000 (cap) | no | ok | 0.162 / 0.043 / 0.104 / 0.066 | 0.82 | yes |
| W_nonfeed | 140 (140) | 10 / 6000 (cap) | no | ok | 0.398 / 0.060 / 0.000 / 0.371 | 0.79 | no |
| NW_nonfeed | 140 (140) | 10 / 6000 (cap) | no | ok | 0.378 / 0.013 / 0.000 / 0.353 | 0.74 | no |

- Safety-net events, spin A+B: 0 / 1 / 0 / 0. Main: 1 / 0 / 0 / 0.
- MIP in-range: **18/32** (S20 21/32) for 1981–2010 and **16/32** (S20 17/32) for 2080–99.

### HPC S21H: 1981–2010 and 2080–2099 (cells are S21H / S10H)

| case | window | GPP | NPP | RECO | RH | NEE | VEGC | SOC | ALD | ALD_FRONT | BNF | NLOST | NUP | dNtot/dt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| W_nfeed | 1981-2010 | 304 / 261 | 175 / 163 | 299 / 257 | 171 / 159 | −4.3 / −4.8 | 2304 / 1674 | 26497 / 28407 | 0.88 / 0.80 | 0.87 / — | 0.110 / 0.107 | 0.115 / 0.102 | 1.60 / 1.40 | −0.0046 / +0.0049 |
| W_nfeed | 2080-2099 | 379 / 334 | 191 / 180 | 374 / 330 | 186 / 176 | −4.6 / −4.5 | 2921 / 2301 | 26344 / 28246 | 1.09 / 1.00 | 1.07 / — | 0.131 / 0.128 | 0.123 / 0.116 | 1.77 / 1.59 | +0.0116 / +0.0171 |
| NW_nfeed | 1981-2010 | 264 / 233 | 162 / 154 | 264 / 233 | 162 / 154 | −0.0 / −0.7 | 1947 / 1444 | 25900 / 27919 | 0.82 / 0.76 | 0.82 / — | 0.100 / 0.097 | 0.106 / 0.099 | 1.46 / 1.31 | −0.0077 / −0.0028 |
| NW_nfeed | 2080-2099 | 265 / 233 | 164 / 155 | 265 / 233 | 164 / 154 | +0.0 / −0.2 | 1941 / 1436 | 25898 / 27941 | 0.83 / 0.76 | 0.82 / — | 0.100 / 0.097 | 0.109 / 0.100 | 1.49 / 1.31 | −0.0108 / −0.0058 |
| W_nonfeed | 1981-2010 | 368 / 367 | 210 / 215 | 357 / 352 | 199 / 199 | −11.3 / −15.4 | 2988 / 2890 | 22592 / 23889 | 0.95 / 0.95 | 0.94 / — | 0 / 0 | 0 / 0 | — | 0 / 0 |
| W_nonfeed | 2080-2099 | 504 / 497 | 242 / 240 | 490 / 484 | 227 / 226 | −14.3 / −13.4 | 4448 / 4376 | 22471 / 23785 | 1.17 / 1.13 | 1.13 / — | 0 / 0 | 0 / 0 | — | 0 / 0 |
| NW_nonfeed | 1981-2010 | 293 / 288 | 175 / 176 | 296 / 290 | 179 / 178 | +3.1 / +2.1 | 2351 / 2254 | 21250 / 23605 | 0.91 / 0.90 | 0.89 / — | 0 / 0 | 0 / 0 | — | 0 / 0 |
| NW_nonfeed | 2080-2099 | 298 / 292 | 180 / 180 | 298 / 292 | 180 / 180 | +0.2 / +0.1 | 2349 / 2257 | 21224 / 23597 | 0.91 / 0.91 | 0.89 / — | 0 / 0 | 0 / 0 | — | 0 / 0 |

### HPC S21H: spins

All S21H spins and mains have safety-net counts of 0. The relcyc stop rule was met for W_nfeed, NW_nfeed (at 6000 yr) and W_nonfeed. NW_nonfeed reached the 6000-yr cap and missed relcyc narrowly (crSOC −1.13e-5 vs 1e-5); its B_cap_report proj150_cyc_ok is True.

### HPC S21H: MIP ranges and site observations

| window | S21H in | S10H in | S21H out-of-range |
|---|---|---|---|
| 1981-2010 | **28/32** | 27/32 | W_nonfeed VEGC above (2988 vs max 2180); NW_nonfeed GPP below (293 vs 302), ER below (296 vs 305), VEGC above (2351 vs 1000) |
| 2080-2099 | **18/32** | 16/32 | W_nfeed NPP, RH below, NEE above (−4.6 vs max −9.9), ALD below (1.09 vs 1.39); W_nonfeed GPP, NPP, ER, RH, ALD below, VEGC above; NW_nonfeed GPP, ER below, VEGC above, NEE above (+0.22 vs 0.04) |

| obs | value | S21H | S10H |
|---|---|---|---|
| Tower 2018 GPP (NT/DT), partial year (254/365 valid days) | 452 / 487 | W_nfeed 318, W_nonfeed 345 (NW 310 / 384) | 279, 348 (NW 266 / 364) |
| Tower 2018 ER (NT/DT) | 348 / 381 | W_nfeed 296, W_nonfeed 353 | 257, 377 |
| Tower 2018 NEE | −107.8 | W_nfeed −21.5, W_nonfeed +8.0 | −22.0, +29.0 |
| GSC 01TC02 3-5 cm, Aug 2016-Jul 2017: JJA / JFM / annual (°C) | 3.87 / −12.0 / −4.0 | W_nfeed 3.69 / −8.02 / −2.17; W_nonfeed 3.96 / −8.18 / −2.27; NW_nfeed 3.77 / −12.51 / −5.23 | 3.53 / −8.30 / −2.35; 3.86 / −8.29 / −2.38; 3.61 / −12.50 / −5.29 |
| model T5 2013-24 JJA (W_nfeed) | — | 4.23 | 4.07 |
| GSC thaw tube 0.75 m (2012/2016) | 0.75 | ALD W_nfeed 0.92 (front 0.94), NW_nfeed 0.81, W_nonfeed 1.06 (1.00), NW_nonfeed 0.88 | 0.90, 0.80, 1.03, 0.93 |

### Open issues (from the reports' assessments)

**TVC S22 vs the tower:**
- GPP is still 16 above the upper tower value, and RECO is 40 above.
- June is still too high, and September GPP is too low (13 vs 24–27).
- Cold-season NEE for 2005–24 is 49.6, just below the 51.9–106.1 band.
- VEGC is about 14 % lower than S20.
- The nonfeed NPP/RH are below the MIP minimum.

**HPC:** the model is below the partial-year 2018 tower GPP and NEE.

**Attribution:** the attribution of the HPC S21H vs S10H changes to individual switches was not tested.
