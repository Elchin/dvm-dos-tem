# RTS-light v0.2 CH4 prototypes (branches rts-ch4-a, rts-ch4-b)

Both prototypes start from rts-light-v0.1 (23322af0). Switch: `model_settings.thermokarst_ch4`
(default `false`; with `false` every output is bit-identical to v0.1, which the tests check). CH4 is computed
monthly in the thermokarst (dsb) step, so it runs in stage 2b of `scripts/run_rts_two_stage.py`, in both
`--mode rts` and `--mode baseline`. Vegetation physics, soil hydrology/thermal and interception are unchanged.

| | (a) reduced form (rts-ch4-a) | (b) monthly upstream port (rts-ch4-b) |
|---|---|---|
| production | `min(1, ch4_frac*rel_q10^((T-10)/10)) * sat * RH_layer`, a partition of RH | `kfrac * sat * prodq10^((T-10)/10) * sum_k kd_k C_k`, limited to `avail_frac*C_layer`; removed from the SOM pools (extra anaerobic pathway, as upstream `Soil_Bgc::CH4Flux`) |
| state | none | CH4 pool per soil layer (`tkch4pool` in restart) |
| oxidation | `(prod-plant)*f_ox(WTD_eff)` | Michaelis-Menten in the oxic part of each layer (upstream Vmax/km/Q10, converted to monthly units) + rhizosphere (50% of plant transport) + transit oxidation of CH4 escaping from below the WT (`f_ox`) |
| transport | plant fraction `plant_frac*min(1,LAI/lai_ref)` bypasses oxidation | plant: first-order `plant_k*fLAI*froot`; diffusion+ebullition: first-order escape `tau_sat`/`tau_unsat` (no solver) |
| frozen layers | no CH4 | stored, released after thaw |
| RH_CO2 (`TKRHCO2`) | `RHSOM - prod + oxid` | `RHSOM + oxid` |

`f_ox(z) = oxid_min + (oxid_max - oxid_min)(1 - exp(-z/oxid_zscale))`, `sat` = saturated fraction of the layer below
the model water table (`Soil_Env::getWaterTable`), raised to `sat + (1-sat)*s` with `s = 1-(1-pond)(1-bg_sat)`,
and `WTD_eff = WTD*(1-s)`.

**C balance (both schemes):** with C = VEGC+SOC+DEADC+DWDC+TKCH4POOL,
`dC = GPP - RG - RM - TKRHCO2 - RHDWD - TKCH4EMIS - TKLATC - TKRESPC` (test: < 0.01 g/m2/month, CMT70, RTS on).
In (a) `TKRHCO2 + TKCH4EMIS = RHSOM`, so RHSOM must not be added to the CH4 flux (that was the double count in the
old rts-ch4 diagnostic).

## RTS-specific terms
- **Post-slump ponding:** `pond(m) = pond_sat0 * exp(-m/(12*pond_tau))`, m = months since initiation (0 before).
  v0.2 defaults `pond_sat0 = 0.2`, `pond_tau = 5 yr`. The emission is insensitive to them: on CMT70/56
  (Gydan, thick ice v2), sat0 between 0 and 0.2 changes the 2021-25 CH4 by <8%, because the thawed slump column
  already sits below the model water table. Moderate values were chosen because young exposed slump floors are
  mostly well drained, with CH4 uptake on the Tibetan Plateau (Nat. Commun. 2025, doi:10.1038/s41467-025-66869-4) and
  CH4 near undisturbed 3-12 yr after thaw (Mu et al. 2018, doi:10.1021/acs.est.7b04979).
- **Post-slump methanogen recovery:** `mg(m) = 1 - (1-mg0) exp(-m/(12*mg_tau))` (1 before initiation) scales
  ch4_frac. Freshly exposed Pleistocene/mineral sediment lacks an active methanogen community (Knoblauch et al.
  2021, doi:10.1029/2021JG006543: CH4 only 3.7-5.1 mg C/m2/d against CO2 0.24-2.6 g C/m2/d). CH4 rises again
  over decades as the floor wets and revegetates (Mu et al. 2018: 3.5x undisturbed after 20 yr).
  `mg0 = 0.05`, `mg_tau = 10 yr`, CALIBRATED so that the first-5-yr post-slump emission of CMT70/56 is ~4-5
  gC/m2/yr. That is near the log-centre of the chamber range for slump floors: 0.4 (Lena, July 3.7-5.1 mg/d
  x ~100 d) to ~16 gC/m2/yr (Baldwin Peninsula RTS 159 mg/d x ~100 d, EGU26-19294). The mg factor leaves CO2
  physics unchanged: suppressed CH4 C stays in TKRHCO2.
  Full data table: russia-inputs/ch4_compare/chamber_data.md/.csv.
- **Post-slump floor RH (v0.2, `include/RTSFloor.h`; three switches, all default ON since v0.2; all false (+ `thermokarst_ch4_mg0 = 0.05`) = bitwise v0.1, tested).**
  Without them, v0.1 RHSOM at CMT70 (Gydan, thick ice v2) reaches 908 gC/m2/month and ~3500 gC/m2 in the
  first post-slump year. Chamber data give 160-184 gC/m2/yr (Lena, yr 1) and floor Reco at or below undisturbed
  (ch4_compare/chamber_data.md).
  1. `thermokarst_rh_cap`: labile-limited mineralization of the exposed mineral-layer C.
     `RH_l <= k(m)/12 * C_l`, with `k(m) = cap_inf + (cap0-cap_inf) exp(-m/(12 cap_tau))`.
     `cap0 = 0.045/yr`, `cap_inf = 0.002/yr`, `cap_tau = 1.5 yr`.
     Knoblauch et al. 2021 (doi:10.1029/2021JG006543) mineralized ~0.8% of thawed C in yr 1, with rates falling
     as labile C is depleted. Calibrated: first-year floor RH 178 (CMT70) and 126 (CMT56) gC/m2.
  2. `thermokarst_floor_export`: floor mass wasting (mudflows, debris tongue) moves exposed mineral C/N laterally.
     The monthly fraction is `lat0/12 exp(-m/(12 lat_tau))` of each mineral SOM pool and orgN, added to
     TKLATC/TKLATN. `lat0 = 0.02/yr`, `lat_tau = 5 yr`.
     Sources: slumps export most of the ecosystem C off-site (Abbott & Jones 2015, doi:10.1111/gcb.13069;
     Abbott et al. 2015, doi:10.5194/bg-12-3725-2015), and sediment export decays as the floor stabilises
     (Kokelj et al. 2013, JGR Earth Surf., doi:10.1002/jgrf.20063).
     The data constrain this only weakly: 0.05-0.1/yr exported 13-23 kgC/m2 more over 25 yr.
  3. `thermokarst_floor_anox`: oxygen limitation on the saturated floor, `RH_l *= 1 - (1-anox_f) sat_l`.
     sat_l is the saturated fraction below the water table, raised by the ponding factor. `anox_f = 0.3`:
     anaerobic C release is ~3.4x smaller than aerobic (Schädel et al. 2016, Nat. Clim. Change,
     doi:10.1038/nclimate3054).
  4. `thermokarst_floor_rmin` (default true): minimum floor respiration. Column RHSOM in each month after
     initiation is at least `r(m) * RH_pre`, with `r(m) = 1 - (1-rmin0) exp(-m/(12 rmin_tau))`,
     `rts_floor_rmin0 = 0.3` and `rts_floor_rmin_tau = 10 yr`. RH_pre is a running mean of column RHSOM for the
     same calendar month over the last 5 years before the trigger (`rts_floor_rmin_years`), so the seasonal cycle
     is kept. It is column-level, not per layer, because the slump restructures the layers. The deficit is
     respired from the existing SOM pools. Layers are weighted by SOM C x moisture x Q10 scaler, pools in
     proportion to their C, and a layer gives at most 5% of its C per month. No C is created, so C/N closure
     holds (tested). If there is no pre-trigger history in the current stage, the term does nothing.
     Data: dry/bare floors have floor/undisturbed ratios of 0.28-0.55, and stabilised slumps return to
     undisturbed (Jensen et al. 2014, doi:10.1088/1748-9326/9/2/025001; Cassidy et al. 2016; Tibet 2025;
     russia-inputs/ch4_compare/chamber_data_dryupland.md).
     Result with 0.3 / 10 yr (no recalibration needed): yr-1 floor/pre-slump ratio is 0.34 (CMT05) and 0.34
     (CMT07), up from 0.12 and 0.04. CMT70 yr-1 RH is unchanged at 178 gC/m2 (the minimum is not binding there).
     A sweep of rmin0 0.3-0.4 and tau 5-10 gives ratios 0.34-0.47.
  v0.2 defaults (include/RTSFloor.h, include/RTSMethane.h, config/config.js): all three switches true plus
  `thermokarst_ch4_mg0 = 1.0`. With RH constrained,
  the methanogen suppression is no longer needed: first-5-yr CH4 is 2.5 (CMT70) and 2.0 (CMT56) gC/m2/yr, inside
  the chamber range of 0.4-16. With the floor switches off, use `mg0 = 0.05` (the v0.1-era calibration).
  All terms apply only after initiation, so baseline runs are unchanged (tested). The C/N balance closes
  (< 0.01 g/month, tested).
- **Exposed thawed carbon:** no new term. The v0.1 slump thaws the deep column and RH from the newly thawed
  layers feeds CH4 production (a) or the SOM-based production (b).
- **Lateral export split:** `TKLATCH4 = lat_ch4_frac * TKLATC` with `lat_ch4_frac = 0.02`, the part of the exported
  C assumed to be mineralized downstream as CH4. Thaw-slump DOC/POC is mostly mineralized aerobically to CO2
  (e.g. Abbott et al. 2014; Littlefair et al. 2017), so this share is small. The split is a diagnostic and leaves
  the cell C balance unchanged.
- **Area blend:** `scripts/rts_area_blend.py`: `X = (1-A) X_baseline + A X_rts` for pools and fluxes
  (A = 0 returns the baseline exactly).

## Parameters (same for all tundra CMTs; `thermokarst_ch4_<name>` model_settings keys override)
| name | (a) | (b) | basis |
|---|---|---|---|
| ch4_frac | **0.135** | - | calibrated: 5.0 gC/m2/yr for wet sedge tundra (CMT06, Gydan cell, bg_sat=1) |
| kfrac | - | **0.29** | calibrated to the same target |
| rel_q10 / prodq10 | 1.5 | 3.0 | methanogenesis is more temperature sensitive than RH. Q10 ~2.8 in the BAWLD-CH4 wetland synthesis (Kuhn et al. 2021); Yvon-Durocher et al. 2014 |
| tref | 10 °C | 10 °C | upstream prodTref |
| oxid_min / oxid_max / oxid_zscale | 0.2 / 0.95 / 0.1 m | same (transit) | 20-90%+ of produced CH4 is oxidized, increasing with WT depth (Le Mer & Roger 2001; Preuss et al. 2013) |
| plant_frac, lai_ref | 0.3, 1 | - | sedge-mediated transport is a large share of wet-tundra emissions (Kwon et al. 2017; Knoblauch et al. 2015) |
| oxid_vmax, oxid_km, oxidq10 | - | 5.76 gC/m3/d, 0.06 gC/m3, 1.9 | upstream (Fan et al. 2013 values: 20 umol/L/h, 5 umol/L) converted |
| plant_k, froot_ref, rhizo_ox | - | 0.1 /d, 0.1, 0.5 | rhizo_ox from upstream (50%) |
| tau_sat / tau_unsat | - | 20 d / 1 d | escape time scales replacing the diffusion solver |
| avail_frac | - | 0.05 /month | cap: production <= 5% of layer SOC per month |
| bg_sat | 0 | 0 | background anoxic/inundated fraction of the cell (e.g. a BAWLD tundra-wetland fraction); 1 = WT at the surface (calibration case) |
| pond_sat0, pond_tau | 0.2, 5 yr | 0.5, 10 yr (b unchanged) | chamber-constrained, weakly sensitive (see above) |
| mg0, mg_tau | 1.0 (v0.2 default; 0.05 with floor off), 10 yr | - | calibrated to slump-floor chamber CH4 (see above) |
| lat_ch4_frac | 0.02 | same | see above |
| rh_cap: cap0, cap_inf, cap_tau | 0.045 /yr, 0.002 /yr, 1.5 yr | same | calibrated (floor section) |
| floor_export: lat0, lat_tau | 0.02 /yr, 5 yr | same | Abbott & Jones 2015; Kokelj et al. 2013 |
| floor_anox: anox_f | 0.3 | same | Schädel et al. 2016 |
| floor_rmin: rmin0, rmin_tau, years | 0.3, 10 yr, 5 yr | same | Jensen et al. 2014; dry-upland compilation |

**Calibration targets:**
- FLUXNET-CH4 (Knox et al. 2019, BAMS; Delwiche et al. 2021, ESSD): most freshwater wetlands emit
  0-20 gC/m2/yr, and wet tundra is among the lowest-emitting wetland classes.
- Arctic wet tundra sites: typically ~1-10 gC/m2/yr (e.g. Samoylov, Wille et al. 2008; Utqiagvik/Atqasuk,
  Zona et al. 2016, where the cold season contributes >= 50%).
- BAWLD-CH4 (Kuhn et al. 2021, ESSD): tundra wetlands have the highest medians, while dry tundra
  (WT > 5 cm below the surface) is low and includes net uptake.

We use 5 gC/m2/yr for fully saturated wet sedge tundra. The dry case (CMT07 heath, model WT) then gives
~0.1-0.2 gC/m2/yr. Neither scheme has CH4 uptake by upland soils (no atmospheric CH4 oxidation sink).

## Post-slump vegetation regrowth (`thermokarst_regrow`, default off)

This is a reduced-form recovery rule v2 (`include/RTSRegrow.h`). Spec: `russia-inputs/rts_recovery_equations.md`. Vegetation equations are unchanged. The rule acts only through parameters (cmax/nmax per PFT), tracked external inputs and the soil N input hook. With the switch off, outputs are bitwise identical to v0.2 (and to v0.1 when the floor switches are also off). **v2 change (by design): the v1 per-PFT logistic growth ceiling is removed entirely** — it held GPP below natural recovery for 12–36 yr and reduced CMT05/07 to 0.49/0.27 × RTS-off. No factor below 1 is ever applied to cmax/nmax. Legacy keys `rts_regrow_limit`, `rts_regrow_t50_*`, `rts_regrow_k_*`, `rts_regrow_release` are ignored.

t = years since stabilization = (months since initiation − 12·`rts_regrow_t_active`) / 12, with `rts_regrow_t_active` = 1 yr.

| process | rule | keys (default) | tracked as |
|---|---|---|---|
| active slump | the 99% kill is kept; no growth ceiling (v1's hold removed in v2) | `rts_regrow_t_active` (1) | – |
| reseed at t = 0 | each PFT raised to f·C_pre,p (pre-slump tissue split and C:N) | `rts_regrow_f_reseed` (0.07) | TKREGROWC/N (external) |
| site timescale (diagnostic) | tau = 1.35·GPP_pre^−1.68 (Xia et al. 2026: yr, mean annual GPP in kgC m−2 yr−1; GPP_pre = model GPP over ≤ 5 yr before initiation), clamped 2–150 yr; logged at initiation, no longer used to cap growth in v2 | `rts_regrow_tau_ref` (7) | – |
| shrub boost | deciduous-shrub cmax and nmax × B(t) = 1 + (B0 − 1)·exp(−t/Tb) ≥ 1 (raises only, never caps) | `rts_regrow_B0` (1.3), `rts_regrow_Tb` (80) | – |
| organic layer | up to 32 (t ≤ 60) then 6.5 gC m−2 yr−1 added to the top fibric/humic layer as slow SOM (Pizano's rate is net accumulation), at C:N 30; only while organic C is below its target; cumulative input capped at the C the slump removed | `rts_regrow_ol_rate1/2`, `_ol_T`, `_ol_cn`, `_ol_frac_raw` (0), `_ol_maxmult` (1) | TKREGROWC/N (external) |
| floor N supply (tapered in v2) | base 2.2 (t ≤ 60) then 0.4 gN m−2 yr−1, minus the nominal OL N, into AVLN of the first rooted layer through `m_a2soi.avlninput` (an input like deposition). Months with no rooted layer carry it over. The AVLN part is scaled by `1 − r^pow` clamped to [0,1], with r = trailing-12-month GPP / GPP_ref(t) and GPP_ref(t) = GPP_pre · INGPP(t)/INGPP_pre (N-unlimited GPP, same trailing window), so the reference tracks the undisturbed climate/CO₂ trend: essentially full while GPP is well below the contemporary pre-slump level, then falling steeply to zero once it reaches it, preventing the v1 no-limit over-fertilisation (end GPP 1.07–1.25 × RTS-off). Static-reference tapers (against GPP_pre alone) were tested and rejected: they starved slow cells, which plateaued at the absolute pre-slump GPP while undisturbed GPP kept rising. | `rts_regrow_n_rate1/2`, `_n_T`, `rts_regrow_n_taper` (true), `rts_regrow_n_taper_pow` (4.0) | NINPUT; NBAL_DIAG avlninput |
| RH-cap exemption | net N mineralization uses the RH before the floor cap (cap applies to CO2 only) | `rts_regrow_navln_exempt` (true) | internal transfer |

NBAL_DIAG lines have one more column: TKREGROWN. A double integrator fallback now also adds the unintegrated part of `avlninput` (always 0 without regrow). Limitation: the regrow state (pre-slump C, GPP_pre) is not in the restart file. A run restarted after initiation does no regrowth. Single-cell results and recommended defaults: `russia-inputs/regrowth_diag/regrow/README.md`.
