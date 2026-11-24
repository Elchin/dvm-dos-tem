# Hg tracer, tier T2: design (explicit total-Hg budget)

Branch `feat/hg-t2` (local) off `feat/hg-tracer` 9e0f657c (T1 final). Author: Hg tracer development, 2026-11-24.
Sources: `HG_REVIEW_AND_EFFORT.md` §3-§5 (T2 definition), ORCHIDEE-Hg v1.0 `src_stomate/stomate_Hg.f90`
(line numbers below refer to that file), T1 work notes (`docs_src/WORK_NOTES_hg_tracer.md` §1-§10).

## 1. Scope

T2 replaces the parts of T1 that are fixed fractions or "follow carbon" assumptions with an explicit
total-Hg budget at the atmosphere-surface interface, in the soil solution and at the export boundary:

* where Hg enters (foliar Hg(0)/Hg(II) uptake tied to LAI, snowpack storage and re-emission, wet deposition
  timed with precipitation, optional photoreduction of deposited Hg(II));
* where Hg leaves (process-based Hg(0) re-emission from soil instead of a fixed `hg_f_evade`, an explicit
  dissolved pool and daily export with overland flow and drainage instead of a fixed `hg_f_export`, fire
  retention, particulate RTS export);
* how it is validated (Toolik flux budget, vegetation and soil Hg, river export and its spring-pulse timing).

Not in T2: a real DOC pool (the "faithful" option, +14-22 pw), MeHg process chemistry (T3; the T1 MeHg
index is kept), Hg feedback on C/N/water/energy. T2 stays diagnostic: it reads C, water, snow and
radiation states and never writes them, so C/N/water/energy stay bit-identical with Hg on. The only planned
element that would feed back is the optional DOC pool (doclost is 0 today, Soil_Bgc line 1828); it is out
of scope and would need its own default-off switch and a decision (D6).

What T1 already delivered from the review's T2 list (so T2 does not redo it): three-component deposition
series from GEOS-Chem with pre-industrial start (T1b/T1c), an effective Kd-partitioned dissolved fraction
with downward transport and lateral drainage export (T1c), fire emission of veg + soil Hg (T1), RTS removal
booked as export (T1), spin-up with an exact affine jump and horizon-aware closure (T1b-T1d).

## 2. Processes added beyond T1

| # | Process | T1 | T2 | Phase |
|---|---|---|---|---|
| a | Foliar uptake | NPP x tissue Hg:C (hg_r_leaf/stem/root) x veg_scale | leaf: f_fix(LAI) x canopy fraction x (dry Hg(0) + dry Hg(II)) (ORCHIDEE); stem/root keep the NPP x ratio pathway (translocation / root uptake) | 1 |
| b | Wet deposition timing | annual / 12 | in proportion to the month's precipitation | 1 |
| c | Snowpack Hg | none (deposition goes straight to soil) | snow Hg pool: receives Hg(II) on the snow-covered fraction, photoreduction in the photic layer, release with melt (to soil) and sublimation (to air) | 1 |
| d | Photoreduction of deposited Hg(II) on canopy/ground | none | optional, ORCHIDEE f = 0.67 exp(-0.016 SW) | 1 (default off, D2) |
| e | Moss Hg | moss C in the leaf compartment, Hg follows C | moss gets foliar uptake through its LAI like other PFTs | 1 |
| f | Explicit dissolved pool hg_diss[] | instantaneous Kd partition inside the monthly transport sub-step | prognostic hg_diss[l] with first-order sorption toward Kd equilibrium | 2 |
| g | Daily export | monthly qinfl/qdrain + layer_drain | daily sub-steps driven by Soil_Env daily qover/qdrain/layer_drain, leachLayers zone; overland flow export from the top layers (Toolik) | 2 |
| h | Vertical mixing | none | cryoturbation D = 1e-3 m2/yr in the active layer to 0 at 3 m, bioturbation 1e-4 (ORCHIDEE) | 2 |
| i | Soil Hg(0) re-emission | fixed hg_f_evade of decay-released Hg | fd_Hg (1 - CUE) form (ORCHIDEE, fd_Hg 0.053) with temperature / snow-cover modifiers; soil-air Hg(0) winter sink replaces the hg0_soil column | 3 |
| j | Fire | all burned Hg emitted | veg/soil split with ash retention (Zolkos 2024; Kumar & Wu 2019) | 4 |
| k | RTS particulate export | RTS lateral C x layer Hg:C | same, booked separately as particulate (St Pierre 2018) with mineral Hg of slumped layers | 4 |
| l | Canopy wash-off | none | leaf Hg(II) fraction washed to the ground with throughfall (rain months) | 4 |

## 3. Equations (phase 1 in full; later phases as design)

Units: ug Hg m-2 per month unless noted. Monthly step (Cohort BGC sub-step, label `bgc`).

**Deposition components** (hg_dep_file; T2 adds the optional column `hg0_dry` = gross GEOS-Chem dry
Hg(0); if absent, hg0_dry = hg0_soil / 0.37, the factor the T1 files were built with):

* wet: W_m = W_yr x P_m / P_yr (P from the year's daily forcing, Climate::prec_d; W_yr / 12 if P_yr = 0);
* dry Hg(II): D2_m = D2_yr / 12; dry Hg(0) gross: D0_m = D0_yr / 12; soil Hg(0) sink S0_m = S0_yr / 12.

**Snow cover** f_s = min(1, SWE_m / SWE_full), SWE_m = monthly-mean snow water (mm).

**Foliar uptake** (ORCHIDEE l.1046-1056): per PFT p with vegcov_p > 0 and LAI_p > 0
f_fix,p = F_max x min(1, LAI_p / LAI_sat);  U_p = a_fol x f_fix,p x vegcov_p x (1 - f_s) x (D0_m + D2_m),
with Sum_p f_fix,p vegcov_p (1 - f_s) capped at 1 (scaled). U_p goes to the leaf Hg of PFT p.
The intercepted dry Hg(II) share F_c = Sum_p U_p(D2) / D2_m is removed from the ground deposition; the
non-fixed dry Hg(0) is re-emitted (bidirectional exchange; not booked, GEOS-Chem dry Hg(0) is gross).
Stem and root keep T1's NPP_part x r_part x veg_scale (Zhou & Obrist 2021: most stem/root Hg is
translocated or root-derived at low concentration). With a_fol = 1 this is ORCHIDEE's scheme.

**Ground and snow deposition** of Hg(II): G_m = W_m + D2_m (1 - F_c).
Photoreduction (optional, ORCHIDEE l.852): f_pr = 0.67 exp(-0.016 SW_m) on the snow-free part; re-emitted
(booked as gross deposition and as evasion).
To the snowpack: G_snow = f_s G_m; to the soil: G_soil = (1 - f_s) G_m (1 - f_pr) + S0_m.

**Snowpack Hg** H_s (ug m-2, new state; pool in the conservation total, restart and spin-up vector):
1. H_s += G_snow;
2. photoreduction in the photic layer: L_red = H_s f_ph (1 - exp(-k_red 24 <cos Z>_m n_m)), with
   f_ph = min(1, z_ph / z_snow), z_snow = SWE_m / rho_snow, k_red = 0.090 h-1 (ORCHIDEE l.845, from the
   Elementa 11:00129 halogen-Hg model, 2023), <cos Z>_m the 24-h mean cosine of the solar zenith angle at
   the month's mid-day declination and the cohort latitude (computed, no new forcing), n_m days;
   booked as evasion;
3. release: r = (M_m + U_m) / (M_m + U_m + SWE_m) with melt M_m and sublimation U_m (mm/month);
   all released if SWE_m < 0.5 mm. Melt share to the top C soil layer (internal transfer), sublimation
   share to the atmosphere (evasion). ORCHIDEE uses the same concentration x melt / sublimation form (l.1018-1020).

Every phase-1 operation is linear in Hg with coefficients from C, water, snow and radiation, so the step
stays affine and the T1d spin-up jump stays exact (snow pool added to the affine vector).

**Phase 2 (dissolved pool, daily export).** hg_diss[l] (ug m-2); sorbed = soil[l] - hg_diss[l].
Equilibrium dissolved fraction phi_l = liq_l / (liq_l + Kd_l M_l) (T1c). Sorption relaxation
d hg_diss/dt = (phi_l H_l - hg_diss_l) / tau_s, tau_s ~ days (ORCHIDEE equilibrates each step, Ota et al. 2013 form;
ORCHIDEE Kd 10^-5.3 m3 g-1 x bulk density, which must be converted and tested against T1's 10^3.3 / 10^3.8 L kg-1).
Daily: export_l = hg_diss_l q_lat,l / liq_l (Richards layer_drain), overland: hg_diss of the top
layers (z < 4.5 cm, ORCHIDEE runoff depth) x qover / (liq + qover); downward advection at a fraction
v_Hg = 0.2 of the water velocity (ORCHIDEE). Vertical mixing: d/dz (D dC/dz) with D = 1e-3 m2 yr-1 in the
active layer, linearly to 0 at 3 m (cryoturbation), 1e-4 m2 yr-1 bioturbation (ORCHIDEE cryoturbate_Hg).
Operator splitting: the daily Soil_Env water fluxes are accumulated in a per-day buffer (Hg-independent),
the Hg sub-steps run at the monthly BGC point in daily order (keeps the C/Hg snapshot discipline).

**Phase 3 (soil re-emission).** Hg released by decay R_l = rh_l x Hg:C_l (T1; rh is the
respired, i.e. (1 - CUE) part of decomposed C). Evasion E_l = fd_Hg x R_l x g(T_l, f_s), with
fd_Hg 0.053 (ORCHIDEE l.2873-3014, after Smith-Downey et al. 2010) and g = 0 below 0 C / under snow
(Obrist 2017: winter soil is a net Hg(0) sink); the rest re-binds or enters hg_diss. Soil-air Hg(0) uptake
S0 (T1 hg0_soil column) becomes D0 x v_soil x f_s ratio, calibrated to Obrist's 37 % winter share.

**Phase 4.** Fire: emitted = (1 - r_ash) x burned layer Hg + burned veg Hg, r_ash from Zolkos 2024 / Kumar
& Wu 2019; RTS: particulate export = slumped C x Hg:C plus mineral Hg of removed layers (St Pierre 2018);
wash-off: leaf Hg(II) x f_wash x P_rain / (P_rain + P_ref).

## 4. Parameters (model_settings; no parameters/ edits)

| Name | Default | Source / note | Phase |
|---|---|---|---|
| hg_t2 | 0 | master switch (needs hg_tracer = 1) | 1 |
| hg_t2_foliar | 1 | LAI foliar uptake replaces leaf NPP uptake (with hg_t2) | 1 |
| hg_fol_fmax | 0.25 | ORCHIDEE-Hg l.1051 | 1 |
| hg_fol_lai | 1.25 m2/m2 | ORCHIDEE-Hg l.1049 | 1 |
| hg_fol_scale (a_fol) | 1.0 | calibration factor (GEOS-Chem dry Hg(0) at Toolik 3.0 vs Obrist 2017 net 6.5) | 1 |
| hg_t2_snow | 1 | snowpack Hg pool (with hg_t2) | 1 |
| hg_snow_kred | 0.090 h-1 | ORCHIDEE-Hg l.845 (Elementa 11:00129, 2023) | 1 |
| hg_snow_photic | 0.10 m | photic depth of snow (King & Simpson 2001, JGR 106, 12499: UV e-folding 5-20 cm) | 1 |
| hg_snow_rho | 250 kg m-3 | tundra snow density, only for the photic fraction | 1 |
| hg_swe_full | 10 mm | SWE for full snow cover | 1 |
| hg_t2_wetprec | 1 | wet deposition by precipitation | 1 |
| hg_t2_photored | 0 | ORCHIDEE photoreduction of deposited Hg(II) (D2) | 1 |
| hg_pr_a, hg_pr_b | 0.67, 0.016 m2/W | ORCHIDEE-Hg l.852 | 1 |
| hg_t2_diss, hg_tau_sorb | 0, 5 d | explicit dissolved pool (Ota et al. 2013 sorption form in ORCHIDEE) | 2 |
| hg_t2_daily | 0 | daily export sub-steps | 2 |
| hg_v_adv | 0.2 | ORCHIDEE Hg advection fraction of water velocity | 2 |
| hg_d_cryo, hg_d_bio | 1e-3, 1e-4 m2/yr | ORCHIDEE cryoturbate_Hg | 2 |
| hg_t2_evasion, hg_fd | 0, 0.053 | ORCHIDEE fd_Hg (Smith-Downey et al. 2010, doi:10.1029/2009JG001124) | 3 |
| hg_fire_ash | 0 | Zolkos et al. 2024; Kumar & Wu 2019 | 4 |
| hg_wash | 0 | canopy wash-off fraction | 4 |

Validation targets for tissue ratios: Zhou & Obrist 2021 (ES&T 55, 14245; global vegetation Hg
assimilation, leaf ~10-60 ng/g), Olson et al. 2018 (tundra veg Hg 20-60 ng/g), Obrist et al. 2017.

## 5. Data needs (open, no logins)

* GEOS-Chem v13.4.1 monthly wet / dry Hg(0) / dry Hg(II)+HgP, pre-industrial run: Sereni & Angot 2026, Zenodo
  doi:10.5281/zenodo.22810103 (CC-BY-4.0; already on the box). Phase 1 adds the gross `hg0_dry` column to the
  site files (from the same NetCDF); phase 2 adds a monthly climatology (spring AMDE Hg(II) dry deposition).
* Obrist et al. 2017 Toolik flux budget (Nature, open supplement), Olson et al. 2018 (GBC), Schuster et al. 2018
  (GRL; data on ScienceBase), Zhou & Obrist 2021, Zolkos et al. 2020 (GBC) and 2024, St Pierre et al. 2018,
  Fabre et al. 2024, Schaefer et al. 2020 (Yukon), Thompson et al. 2025. All open access or with open supplements.
* Model forcing already in the test sites: precipitation (daily split), SWE, melt, sublimation, nirr, LAI, latitude.

## 6. Coupling points (dvmdostem)

| Need | Where | How |
|---|---|---|
| LAI, vegcov per PFT | `cd.m_veg.lai[ip]`, `cd.m_veg.vegcov[ip]` (Cohort BGC block) | read into hgtr::Flux |
| SWE, melt, sublimation | `edall->m_snws.swesum` (monthly mean), `m_snw2soi.melt`, `m_snw2a.sublim` (monthly sums; EnvData 452-484) | read |
| Precipitation | `edall->m_a2l.prec`, year total from `climate.prec_d` | read |
| Shortwave | `edall->m_a2l.nirr` (monthly mean W m-2) | read |
| Latitude | `Cohort::lat` | <cos Z> computed in hgtr |
| Tracer step | `hgtr::step` (`include/HgTracer.h`); Cohort BGC sub-step at the `bgc` label | phase-1 block before the veg loop |
| Spin-up | `hgtr::NVEC`, `to_vec_/from_vec_` | snow pool appended |
| Restart | `hgtr::pack/unpack`, `RestartData::read_px_hg_vars` | layout v3 (+1 value); v2 and v1 still read |
| Outputs | `Runner.cpp` Hg block, `config/output_spec.csv` | HGSNOW, HG_FOLIAR, HG_SNOWRED (rows default off) |
| Phase 2 daily water | `Soil_Env` daily qover/qdrain/layer_drain, `Soil_Bgc::leachLayers` (line 95) zone | per-day buffer filled in the env loop |
| Phase 3 | decomposition flux in the tracer step (T1 `f.decomp`) | evasion function |
| Phase 4 | `WildFire::burn` (line 295) `r_retain_c`; Thermokarst/RTS label `rts` | fractions |

## 7. Switches

All default 0 (off) at model level; hg_t2 needs hg_tracer = 1. Off: none of the T2 code runs; the T1 step is
unchanged, so Hg-off outputs are bit-identical to binM and T1 Hg-on results unchanged (the restart gains one
value and the spin-up vector one zero component, both exact). On: C/N/water/energy bit-identical (diagnostic).
No T2 process feeds back. The optional DOC pool would (D6) and stays out of T2.

## 8. Test plan

Unit tests (doctest, `tests/unit/test_hg_tracer.cpp`):
1. parse: hg_t2 defaults, ranges, ignored without the tracer;
2. f_fix(LAI) piecewise and the canopy cap; uptake with snow cover;
3. wet deposition by precipitation sums to the annual rate;
4. snow pool: accumulation, photic-layer photoreduction, melt and sublimation release, exact balance;
5. 24-h mean cos zenith vs numerical integration (polar day / night, equator);
6. hg_t2 = 0 step identical (bit-exact) to the T1 step on random states;
7. random multi-step runs with hg_t2 = 1 (geometry changes, fire, snow): balance to roundoff;
8. spin-up jump == explicit repetition with the snow pool;
9. restart v3 round trip, v2 / v1 still read.

Model checks (regression harness): Hg off vs binM (8-case prod or the 2-case subset), bit-identical;
hg_t2 on vs Hg off: C/N/water/energy outputs identical; HG_BALERR at roundoff; site runs TVC / HPC /
Toolik with the T1d protocol (75-yr segments, 4 cores).

## 9. Validation sites and data

| Site | Data | Phase-1 metric | Later |
|---|---|---|---|
| Toolik (Obrist 2017; Olson 2018) | net Hg(0) 6.5 +- 0.7 /yr, veg uptake 63 %, winter soil sink 37 %, snowmelt small; veg Hg 20-60 ng/g | foliar uptake vs 4.1, snow re-emission share, leaf Hg:C | export (Zolkos 2020 Kuparuk/Toolik inlet), profiles |
| TVC, HPC (Trail Valley / Havikpak, T0 inventory, Azaroff 2026 synthesis) | organic RHgC 0.32-0.35, mineral 0.97-1.69 | uptake drop vs T1 (21.7 / yr at TVC) | organic Hg:C (open item 1), export (Mackenzie, Zolkos 2020) |
| Schuster 2018 cores | RHgC 1.6 median | - | profiles (phase 5) |
| Yukon / Mackenzie rivers | Schaefer 2020, Zolkos 2020, Fabre 2024 yields, spring pulse | - | daily export timing (phase 2/5) |

## 10. Phase plan

| Phase | Content | Person-weeks | Core-h (runs) |
|---|---|---|---|
| 1 | Surface interface: foliar uptake (LAI), snowpack Hg, wet deposition timing, optional photoreduction, restart v3, outputs, unit tests | 2-3 | 0.15-0.2 (identity subset 0.07-0.10, sites 0.03-0.06) |
| 2 | Explicit hg_diss[], daily export (qover / qdrain / layer_drain, leachLayers zone), overland flow, cryo/bioturbation | 3-4 | 0.4-0.6 (full identity 0.27, sites + sensitivities 0.2) |
| 3 | Process-based soil re-emission (fd_Hg, temperature / snow) and soil-air Hg(0) sink | 2-3 | 0.3-0.4 |
| 4 | Fire retention, RTS particulate, canopy wash-off | 1-2 | 0.2-0.3 (fire and RTS cases) |
| 5 | Calibration / validation (Toolik budget, veg Hg, profiles, river export and timing), C spin-up protocol for TVC, sensitivity / identifiability (~200 member ensembles x 3 sites), docs | 3-5 | 2.5-4 |
| total T2 | | 11-17 (T1 delivered ~5-8; 16-25 cumulative, review range 15-25) | 3.6-5.5 |

The remaining T1 cap leaves 0.34 core-h, enough for phase 1 only. Phases 2-5 need a new budget (D1).

## 11. T1 open items and where T2 addresses them

| T1 open item | Cause (T1d) | T2 |
|---|---|---|
| TVC / HPC organic Hg:C 1.1-2.6 (T0 0.32-0.35) | vegetation pump NPP x r_leaf (TVC uptake 21.7 /yr, 6x Toolik) with only 2.5 % of decay-released Hg leaving | phase 1 ties leaf uptake to LAI and atmospheric dry deposition (bounded by GEOS-Chem dry flux); phase 3 replaces the fixed 2 % evasion; phase 2 lets dissolved Hg leave daily |
| HPC export 5.8 /yr, drainage-dominated | monthly Kd-equilibrium drainage of the dissolved fraction | phase 2: sorption kinetics, advection at 0.2 of water velocity, daily drainage vs overland split; phase 5 calibration to Zolkos 2020 yields |
| Toolik export 0.24, no lateral drainage | Toolik column has no Richards layer drainage; overland flow is not used | phase 2: overland flow (qover) export from the top 4.5 cm and the snowmelt pulse; phase 1 snow melt release feeds it |
| TVC carbon not at equilibrium (organic C -30 % in 300 yr in the 1950-69 stage) | C spin-up protocol, not Hg | phase 5: a longer C eq / sp stage or C-equilibrium check before the Hg jump; flagged D4, not a code change in Hg |
| TVC %MeHg with no Thompson class | wetland class missing | T3, not T2 |

## 12. Decisions for Elchin

* D1: a run budget for phases 2-5 (about 3.4-5.3 core-h beyond phase 1).
* D2: ORCHIDEE's canopy / ground photoreduction f = 0.67 exp(-0.016 SW) decreases with light; kept default off.
* D3: foliar scale a_fol: ORCHIDEE's scheme with GEOS-Chem dry Hg(0) gives less uptake than Obrist 2017 measured.
* D4: C spin-up protocol for TVC (longer eq stage) before the Hg validation.
* D5: whether T2 defaults (hg_t2 = 1) become default with the tracer once phase 5 passes.
* D6: DOC pool (faithful option) is the only element that would feed back on C; out of T2 unless approved.
* D7 (from phase 1): snow sublimation. ORCHIDEE loses the Hg of sublimated snow to the air, and that term gives most of
  the 80-86 % snow re-emission. Alternative: Hg stays in the pack (concentrates) and leaves with the melt.
* D8 (from phase 1): the stem/root NPP x r pathway still carries 4.1 of TVC's 5.3 uptake. Options: keep it, lower
  hg_r_stem/root to Zhou & Obrist 2021 wood/root values, or translocate it from leaf Hg.

## 13. Phase 1 status (2026-11-24)

Implemented and run (see WORK_NOTES_hg_t2.md section 1). Organic Hg:C went from 2.63 to 0.78 at TVC and from 1.06 to 0.35 at HPC.
Export went from 4.26 to 1.40 at TVC and from 5.80 to 1.91 at HPC. Toolik vegetation uptake fell from 3.56 to 0.53, which is far below
Obrist 2017 (D3). Phase-1 runs cost 0.179 core-h, leaving a ledger total of 2.237 of 2.4.

## 14. Phase 2 status (2026-11-24)
Done on local feat/hg-t2: c32cd335, ecc13285, ae1a200b (+ docs). hg_t2_diss (default 0) adds the dissolved pool, the
daily sorption / overland / lateral / vertical export, and the monthly cryo/bioturbation. Beyond section 3 it adds a moss
sorbent mass (hg_t2_rho_moss) and lets the removal fractions act together. hg_t2_log_kd_org was refit to 3.5 and stays
provisional until phase 5; D7 is the default. All identity, conservation and unit checks pass
(WORK_NOTES_hg_t2.md section 2). Next is phase 3 (D3, GPP/stomatal foliar uptake).

## 15. Phase 3 status (2026-11-24)
Done on local feat/hg-t2: e3f1c3c2, 4c0d6914, 98aa09ee (+ docs). Three switches, all default off:
- hg_t2_evasion: per-layer re-emission fd x (1 - f_snow) x [T > 0] x exp(-z / ev_z), plus a process ground Hg(0) sink.
- hg_t2_foliar 2 (D3): stomatal uptake from (N-unlimited) GPP, CLM5-Hg / Medlyn, calibrated with hg_gs_eff 0.75. The LAI scheme remains hg_t2_foliar 1.
- hg_t2_meltroute: snowmelt Hg delivered daily with the melt.
The D8 sensitivity is done; the stem/root ratios stay as they were. Toolik uptake (3.99 stomatal) and net Hg(0) (6.0) match
Obrist 2017. Organic Hg:C and export are now above target at all sites, which is the main phase-5 item.
Checks: WORK_NOTES_hg_t2.md section 3.

## 16. Phase 4 status (2026-11-24)
Done on local feat/hg-t2: a94f93bb (+ docs). Three items, all default off, hg_t2 = 1 only:
- hg_fire_ash: the burned soil Hg share kept in ash, on the new top layer.
- hg_t2_rts: RTS particulate export, including the Hg of the respired C share.
- hg_wash / hg_wash_pref: wash-off of the month's intercepted dry Hg(II).
Tested in a Toolik case with an explicit 1926 fire and a 1956 RTS. Ash retention adds a fire-year export pulse; the RTS change
is +9 % in the slump year. Wash-off is 0.02-0.16 ug/m2/yr at the sites. Checks: WORK_NOTES_hg_t2.md section 4.

## 17. Phase 5a status (2026-11-24)
Phase 5a is done on local feat/hg-t2, docs and preset only. There are no code changes, so the binary is still binT2p4.
- D4: chained TVC carbon spin-up of 680 yr in new restarts, ending at -0.46 %/100 yr drift.
- Refit: a 64-member identifiability ensemble with an emulator fit.
- Fitted set (config/hg_t2_calibrated.json): t2_log_kd_org 3.46, v_ground 1.1e-4, fd 0.10, ev_z 0.20, r_stem/r_root 0.0044/0.0067 (D8a), snow_kred 0.090 (not identifiable) and k_meth_proxy 71.
- Targets met: Toolik uptake 4.06, net Hg(0) 6.70, the %MeHg classes (0.70 / 0.100 / 2.06) and mineral Hg:C.
- Targets missed (structural): export at TVC (6.5) and Toolik (4.9); organic Hg:C at TVC and Toolik (+10-15 %) and HPC (-5 %); and the winter share (28 % vs 37 %, tied to v_ground).
- Defaults are proposed, not flipped.
- All identity, conservation and unit checks pass.
- Details are in WORK_NOTES_hg_t2.md section 5. The optional 5b round is a wider ev_z range or structural options.

## 18. Phase 5b status (2026-11-24)
Done on local feat/hg-t2 (not pushed): d787f0fd (code + tests), 4339a0fc (preset values), + docs. Binary binT2p5c.
- New switches inside hg_t2, default off at component level: hg_t2_kdsplit + hg_t2_dlog_kd_hum (humic-layer Kd offset,
  layer type 2), hg_t2_xpath + hg_f_conn (organic-layer overland/lateral share retained to the first mineral layer),
  hg_t2_snowsink + hg_v_snow (ground Hg(0) sink under snow, weighted by snow fraction).
- hg_t2 = 1 now takes the calibrated preset (hg_t2_preset 1, explicit keys win). hg_t2 itself stays default 0 (D5).
- T2-only hg_t2_r_stem / hg_t2_r_root / hg_t2_k_meth_proxy (T1 names with hg_t2 = 1 are a parse error).
- D4: restarts/TVC_Warming_nfeed_D4_680yr_restart-sc.nc is the standard TVC start (cs/sites_d4, cs/STANDARD.md).
- Refit (E6 26 members, ev_z to 0.5 m, exact v_ground / v_snow solve): kd_org 3.453, dlog 0.242, f_conn 1 (xpath off),
  fd 0.083, ev_z 0.50, v_ground 4.88e-5, v_snow 1.283e-4; r_stem / r_root 0.0044 / 0.0067, k_meth_proxy 71, snow_kred 0.090.
- Met: Toolik uptake 4.06, net Hg(0) 6.40, winter share 0.367, Rm, %MeHg 0.709 / 0.100 / 2.05, TVC Ro 0.409.
  Missed: export TVC 6.15 / Toolik 4.77 (> 4), Toolik Ro 0.448, HPC Ro 0.268. Flow-path retention cuts export but raises Ro.
- All identity, conservation and unit checks pass. 5b cost 1.136 core-h; ledger 7.184 / 7.537.
- Details and the T2 summary: WORK_NOTES_hg_t2.md sections 6-7.

## 19. T2 final status
T2 phases 1-5b are complete on local feat/hg-t2. Decisions left for Elchin: D5 (default-on of hg_t2 with the tracer), D6
(DOC pool / deeper retention path for the export misfit), ev_z 0.5 vs 0.3, keeping or dropping kdsplit, and a T2 copy of k_meth_ch4.

## 20. Phase 5c status (2026-11-24)
Done on local feat/hg-t2 (not pushed): c19c2b9c (code + tests), 2972fc49 (preset values), + docs. Binary binT2p5e.
- Option A (approved): hg_t2_deeproute (needs hg_t2_xpath) with hg_f_deep / hg_z_deep sends the flow-path retention (1 - f_conn)
  to the layer at z_deep below the cryo/bioturbation zone instead of the first mineral layer. It is bit-identical to 5b when
  off, and conservative.
- The mechanism behaves as intended: with retention at 1.5 m, Ro is unchanged and export falls by 40 %. Retention inside the
  active layer (0.5 m, or the first mineral layer) mixes back and raises Ro.
- Refit by local Gauss-Newton / Broyden on site runs (12 runs): kd_org 3.377, f_conn 0.536, z_deep 1.5, fd 0.0688, ev_z 0.5,
  v_ground 4.37e-5, v_snow 1.304e-4; dlog_kd_hum 0.242 kept (no effect). All 12 bands met (cost 0.008, was 6.86).
- Export 3.83 / 1.85 / 2.64 at TVC / HPC / Toolik. Ro 0.407 / 0.291 / 0.391. Rm 1.545 / 1.689 / 1.477. Toolik net Hg(0) 6.50,
  winter share 0.370. %MeHg 0.709 / 0.100 / fen 2.03.
- ev_z 0.3 fits equally (V8b).
- New caveat: the routed Hg accumulates as a permafrost sink of 1-2.3 ug/m2/yr (column storage up 2-3 ug/m2/yr at TVC / HPC).
- All identity, conservation and unit checks pass. 5c cost 0.683 of 1.2 core-h; ledger 7.867.
- Details: WORK_NOTES_hg_t2.md section 8.

## 21. Phase 5d status (2026-11-24)
Done on local feat/hg-t2 (not pushed): 91b456de (hg_t2_k_meth_ch4 copy + test), f60c7170 (preset), + docs. Binary binT2p5f.
- Preset: ev_z 0.3, hg_t2_dlog_kd_hum 0 (switch kept). New T2-only hg_t2_k_meth_ch4 (3250; the T1 name with hg_t2 = 1 is an
  error).
- Refit from the 5c Jacobian (no exploration runs): kd_org 3.356, f_conn 0.534, f_deep 1 (lower f_deep raises Ro), fd 0.0719,
  v_ground 4.14e-5, v_snow 1.305e-4, z_deep 1.5.
- 10 of 12 bands met. HPC Ro 0.288 (< 0.29) and HPC Rm 1.706 (> 1.69) are about 1 % out; with ev_z 0.3 the local model has
  no point meeting all 12. Export 4.0 / 2.1 / 2.6, Toolik net Hg(0) 6.50 and winter share 0.370, %MeHg 0.708 / 0.100 / fen 2.02.
- Deep store below 1 m: +3.5 (TVC) and +1.8 (HPC) ug/m2/yr.
- All checks pass. The preset-vs-explicit, routing-off and hg_t2 = 0 checks ran on the Toolik column only, because of the
  budget. 5d cost 0.278 of 0.3; ledger 8.145.
- Details: WORK_NOTES_hg_t2.md section 9.
