# Hg tracer T2: work notes

Branch `feat/hg-t2` (local, not pushed) off 9e0f657c. Worktree `~/tem-dev/hg_t2/wt` on calib-elchin.
Design: `docs_src/HG_T2_design.md`. Ledger: `~/tem-dev/hg_t2/LEDGER.txt`, continuing the T1 total of 2.0586 core-h; cap 2.4.

## 1. Phase 1: surface interface (2026-11-24)

**Code** (`include/HgTracer.h`, `src/Cohort.cpp`, `src/Runner.cpp`, `src/RestartData.cpp`, `config/output_spec.csv`, `tests/unit/test_hg_tracer.cpp`):
* `hg_t2` (default 0) turns on the phase-1 block, `hgtr::surface_t2_`. It runs in the BGC sub-step only (`Flux::t2 = 1`).
  * Foliar uptake (`hg_t2_foliar`), following ORCHIDEE-Hg l.1046-1056: leaf input of PFT p is
    hg_fol_scale × f_fix(LAI_p) × vegcov_p × (1 − f_snow) × (dry Hg(0) + dry Hg(II)).
    * f_fix = 0.25 × min(1, LAI/1.25), and the canopy total is capped at 1.
    * This replaces the leaf NPP × r_leaf uptake. Stem and root keep NPP × r.
  * Wet deposition (`hg_t2_wetprec`) follows the month's precipitation: rain_d + snow_d of the year.
  * Snowpack Hg pool (`hg_t2_snow`).
    * Inputs: Hg(II) falling on the snow-covered fraction min(1, SWE/10 mm).
    * Photoreduction of the photic part min(1, 0.10 m / z_snow) at 0.090 h⁻¹ × 24-h mean cos Z (ORCHIDEE l.845). This loss is booked as evasion.
    * Release in proportion to melt, sublimation and SWE (ORCHIDEE l.1018-1020): the melt share goes to the top C soil layer, the sublimation share is evasion.
  * Optional ORCHIDEE photoreduction of ground Hg(II), 0.67·exp(−0.016·SW) (`hg_t2_photored`, default 0).
* Deposition file: new optional column `hg0_dry` (gross dry Hg(0)). When it is absent, hg0_dry = hg0_soil / 0.37.
  * For the GEOS-Chem files this is exactly the GEOS-Chem gross dry Hg(0), rounded to 4 decimals.
  * For the Obrist Toolik file it is the measured net Hg(0) deposition, 6.5.
  * Phase 1 needs no new data files.
* The snow pool is part of:
  * `State::total()`, so conservation is exact;
  * the spin-up affine vector (NVEC + 1);
  * the restart: layout v3 has 95 values; v2 (94) and v1 (92) are still read.
* New outputs: HGSNOW, HG_FOLIAR, HG_SNOWDEP, HG_SNOWRED, HG_SNOWMELT, HG_PHOTORED.
* No feedback: C, N, water and energy are only read.

**Checks:**
* Unit tests.
  * New tests: parse; f_fix and the 24-h cos Z against numerical integration; hg_t2 = 0 bit-exact against the T1 step; foliar balance and canopy cap; snow accumulation, reduction and release; photoreduction; 3000 random steps balanced to roundoff; spin-up jump with snow equal to explicit repetition; restart v3/v2.
  * Standalone Hg suite: 27/27. Full `make unittest`: 509/509 cases. The 8 known may-fail assertions are the same as in T1.
* T1 path in the model: binT2p1 with hg_t2 = 0 (P1T1) against binT1d_final (VF4), TVC/HPC spin + main and Toolik sp + tr.
  * Every common NetCDF variable is bit-identical.
  * The only exception is the restart `hg_state` length, 95 vs 94 (layout v3, expected).
* Hg off: binM against binT2p1 on hg_t1/cs/sub_off, run_regression --no-repeat: PASS (bit-identical).
* hg_t2 = 1: binM against binT2p1 on cs/sub_t2 (hg_tracer 1, hg_mehg 1, hg_t2 1): C/N/water/energy identical = True (cmp_on.py).
* HG_BALERR is at most 2e-10 ug/m2 in all site runs.

**Site runs.**
* Protocol: T1d VF4, with GEOS-Chem deposition at TVC/HPC and Obrist at Toolik. Spin-up is a 10-yr cycle jumped to 5000 yr; segments are ≤ 75 yr; 3 processes run in parallel.
* All values are means of the last 20 yr. Export targets are pan-Arctic river yields of 0.8-4.0 (Zolkos 2020). Organic Hg:C target: T0 0.32-0.35.

| Site | Variant | R org | R min | Uptake (foliar) | Export | Snow dep → re-emit / melt | Evasion |
|---|---|---|---|---|---|---|---|
| TVC | T1d (P1T1) | 2.63 | 1.46 | 21.73 | 4.26 | - | 6.93 |
| TVC | T2 p1 (P1T2) | 0.78 | 1.38 | 5.30 (1.22) | 1.40 | 2.51 → 2.16 / 0.35 | 4.56 |
| HPC | T1d | 1.06 | 2.00 | 12.67 | 5.80 | - | 2.99 |
| HPC | T2 p1 | 0.35 | 1.69 | 2.79 (1.75) | 1.91 | 1.70 → 1.59 / 0.11 | 2.61 |
| Toolik | T1d | 0.42 | 1.58 | 3.56 | 0.24 | - | 0.95 |
| Toolik | T2 p1 | 0.22 | 1.53 | 0.53 (0.22) | 0.21 | 1.80 → 1.49 / 0.31 | 2.32 |

**Sensitivities:**
* P1F4, hg_fol_scale 4 (the scheme's upper bound: all dry deposition fixed under a full canopy):
  * foliar uptake: Toolik 0.87, TVC 2.64, HPC 3.54;
  * organic Hg:C: TVC 0.86, HPC 0.43, Toolik 0.24.
* P1K10, hg_snow_kred / 10: snow re-emission falls by only 11-15 %, because the sublimation term dominates. Melt release to soil: TVC 0.64, HPC 0.34, Toolik 0.46.

**Findings:**
1. Tying leaf uptake to LAI × dry deposition removes most of the T1 vegetation pump.
   * TVC: 21.7 → 5.3. Of that, 4.1 is still the stem/root NPP × r pathway.
   * HPC organic Hg:C 0.35 is in the T0 range, and HPC export 1.9 falls within the Zolkos range. This closes T1 open item 2 for phase 1.
   * TVC organic Hg:C 0.78 is still about 2× T0. The remainder comes from the stem/root NPP pathway (shrub tundra) and the fixed 2 % evasion; it is addressed in phases 3 and 5.
2. The ORCHIDEE LAI scheme cannot reproduce the Toolik vegetation pump.
   * Obrist 2017 found ~4.1 µg/m²/yr growing-season uptake and a net 6.5. T2 gives 0.22 foliar; even the upper bound gives 0.87.
   * Toolik's LAI and canopy cover are low. The measured pump would need a GPP/stomatal-conductance scheme (CLM5-Hg, Yuan 2023), or a much larger effective dry Hg(0) flux. Decision D3.
3. Snowpack: 80-86 % of the Hg(II) deposited to snow is re-emitted, mostly through the ORCHIDEE sublimation-loss form.
   * If Hg stays in the pack when snow sublimates, which is likely and is how most snow-chemistry studies treat it, melt release rises.
   * This is a structural choice. Decision D7: proposed default "sublimation concentrates Hg".
4. Toolik export is unchanged at 0.21, because it has no lateral drainage. That is phase 2: overland flow and the snowmelt pulse.
5. Toolik organic Hg:C 0.22 is below the Toolik organic range of 0.29-0.41 that hg_log_kd_org was calibrated to. Kd_org should be recalibrated with the T2 inputs (phase 2/5).
6. The T1 "net Hg(0) uptake" diagnostic (uptake + 0.37·dep − evasion) is not comparable under T2: evasion now includes Hg(II) re-emitted from snow. Phase 3 should add an explicit Hg(0)/Hg(II) split of the fluxes.

**Core-h:** phase-1 runs 0.1786 (P1T1, P1T2, P1F4, P1K10 ≈ 0.027 each; I1 0.033; I2 0.034). Total 2.2372 of 2.4; 0.163 remains.

**Files on the VM, under ~/tem-dev/hg_t2/:**
* bins/binT2p1: md5 3b04ea7d38427e044e3f7e233bf07cb7, GIT_SHA g639e3b.
* src/: runP1.sh, runP1s.sh, runT2.sh, mk_t2.py, mk_cs_t2.py, cmp_t1.py, an_t2.py, an2.sh, plots_t2.py, ledger.py, build.sh.
* an2/: per-variant json/txt.
* Figures: docs_src/figures_hg_t2/hg_t2_fig1_surface_monthly.png, hg_t2_fig2_sites.png.

## 2. Phase 2: dissolved pool, explicit export, mixing; Kd_org refit (2026-11-24)

Commits (local feat/hg-t2, not pushed): c32cd335 code (hg_t2_diss, default 0; D7; restart v4; outputs; tests),
ecc13285 refit fix (moss sorbent, joint daily removal, tau_sorb 0.25 d), ae1a200b Kd_org default 3.5, then the docs commit.
Binary bins/binT2p2c (md5 fcc2bf38..., GIT_SHA ae1a200b). Runs: hg_t2/src/runP2.sh stages a-e (logs out/runP2_*.log).

### What it does (hg_t2 = 1, hg_t2_diss = 1)
- Separate dissolved pool diss[l] (in total, restart v4 = 117 values, spin-up vector NVEC + 22). Decay-released Hg
  that does not evade goes to diss (hg_f_export not used). Geometry changes re-sorb diss first.
- Monthly transport sub-step loops over the month's days with the model's daily water (qinfl, qover, Richards layer
  drain, daily layer liq / ice, buffered in Cohort after grnd_endOfDay): sorption relaxation toward
  phi = liq / (liq + Kd M) (tau 0.25 d); overland export from the layers starting above 4.5 cm, qover / (qover + liq);
  lateral drainage diss x drain / liq; vertical advection 0.2 q / liq; the three fractions act together on the same
  diss, sum capped at 0.9; stopped at the frost table. Then monthly cryo/bioturbation of sorbed Hg (ORCHIDEE).
- Moss layers (type 0) carry no soil C, so the T1c sorbent M = C / 0.5 was 0 there: phi = 1 and all deposition reaching
  the moss ran off (first pass P2A: Toolik organic Hg:C 0.02, Kd-insensitive). Fix: M_moss = dz x hg_t2_rho_moss.
- D7: hg_snow_sublim_loss 0 (default): sublimated snow Hg stays in the pack and leaves with the melt.

### Parameters (model_settings; no parameters/ edits)
| name | value | source |
|---|---|---|
| hg_t2_log_kd_org | 3.5 (log10 L/kg) | refit here (scan 2.8-4.3) to Toolik organic Hg:C 0.29-0.41 (Obrist et al. 2017 Nature 547:201) and the T0 0.32-0.35; provisional until the phase-5 joint refit |
| hg_log_kd_min | 3.8 (unchanged) | Allison & Allison 2005 EPA/600/R-05/074 |
| hg_tau_sorb | 0.25 d | Hg(II)-soil sorption near equilibrium within hours (Yin et al. 1996 J Environ Qual 25:837); scan: tau 0.02-0.5 d alike, 5 d changes Hg:C |
| hg_t2_rho_moss | 30 kg/m3 | live + dead moss dry bulk density 10-50 kg/m3 (O'Donnell et al. 2009 Soil Sci 174:646) |
| hg_v_adv | 0.2 | ORCHIDEE-Hg |
| hg_z_runoff | 0.045 m | ORCHIDEE-Hg runoff mixing depth |
| hg_d_cryo, hg_z_cryo, hg_d_bio | 1e-3 m2/yr, 3 m, 1e-4 m2/yr | ORCHIDEE-Hg cryoturbate_Hg |
| hg_snow_sublim_loss | 0 (D7) | approved design choice D7 |

### Sites (last 20 model years; ug/m2/yr; Hg:C ug/g). P2FT1 = hg_t2 0 (= T1d), P2N0 = phase 1 + D7, P2F = phase 2 defaults
| site | run | R org | R min | uptake (foliar) | export (overland / drain) | evasion | net storage | %MeHg JJAS |
|---|---|---|---|---|---|---|---|---|
| TVC | P2FT1 | 2.63 | 1.46 | 21.7 | 4.26 | 6.93 | +16.4 | 1.24 |
| TVC | P2N0 | 0.78 | 1.38 | 5.30 (1.22) | 1.40 | 4.55 | +4.67 | 1.11 |
| TVC | P2F | 0.43 | 1.28 | 5.30 (1.22) | 5.45 (1.53 / 3.92) | 4.51 | +0.65 | 0.98 |
| HPC | P2FT1 | 1.06 | 2.00 | 12.7 | 5.80 | 2.99 | +9.24 | 0.10 |
| HPC | P2N0 | 0.35 | 1.69 | 2.79 (1.75) | 1.92 | 2.60 | +3.19 | 0.07 |
| HPC | P2F | 0.33 | 1.53 | 2.79 (1.75) | 3.13 (0.64 / 2.49) | 3.76 | +0.81 | 0.12 |
| Toolik | P2FT1 | 0.42 | 1.58 | 3.56 | 0.24 | 0.95 | +7.48 | 0.12 |
| Toolik | P2N0 | 0.22 | 1.53 | 0.53 (0.22) | 0.21 | 2.26 | +3.11 | 0.11 |
| Toolik | P2F | 0.35 | 1.29 | 0.53 (0.22) | 3.07 (3.07 / 0) | 2.30 | +0.21 | 0.20 |

Targets: Toolik organic Hg:C 0.29-0.41 (met at Toolik and HPC; TVC 0.43 is high); Toolik plant uptake ~4.1
(Obrist 2017; model 0.53, not addressed until phase 3 D3); export vs Zolkos et al. 2020 river yields 0.8-4.0
(HPC, Toolik inside; TVC 5.45 above).

Kd scan (binT2p2b, fig3): log Kd 3.3 / 3.5 / 3.6 / 3.8 / 4.3 gives Toolik R org 0.23 / 0.35 / 0.43 / 0.57 / 0.86,
TVC 0.33 / 0.43 / 0.49 / 0.59 / 0.76, HPC 0.28 / 0.33 / 0.35 / 0.39 / 0.43; TVC export 6.4 / 5.5 / 4.9 / 3.8 / 1.6.
No single Kd fits TVC Hg:C and TVC export together with Toolik; 3.5 favours the Toolik / HPC Hg:C targets.
Sensitivity hg_t2_rho_moss 0 (fix off): Toolik R org 0.08, TVC export 8.8.

### Checks
- Hg-off: sub_off (2 cases) PASS and full 8-case prod vs binM PASS (bit-identical).
- hg_t2 = 0 (P2FT1, binT2p2c) vs T1d VF4: every output bit-identical; only the restart hg_state length differs
  (117 vs 94). (HG_DRAIN is now written monthly in these variant runs; its monthly sums equal VF4's yearly to 2e-16.)
- hg_t2 = 1 with hg_t2_diss 0 and hg_snow_sublim_loss 1 (P2S1, c32cd335) = phase-1 P1T2 bit-exact.
- hg_t2 + diss on (sub_t2d, 2 cases): C/N/water/energy identical to binM, ndiff 0; balance error max 9e-10 ug.
- Conservation: site balance errors <= 2e-9 ug on 1.4e5 ug (1e-14 relative).
- Unit tests: standalone Hg suite 35/35 (609 assertions); full unittest 517/517 cases (8 known may-fail assertions as T1).

### Findings and open items
- Export timing: peaks Jul-Aug (fig4), May-Jun only 25-40 %; observed Arctic river Hg export is freshet-dominated
  (Zolkos 2020). Melt delivers little Hg (snow re-emission takes 80-86 %); candidate for phase 4/5.
- D7 effect is small: Toolik snow re-emission 1.49 -> 1.42, melt release 0.31 -> 0.37 ug/m2/yr.
- At Toolik, cryo/bioturbation moves 1.6 ug/m2/yr of mineral Hg up into the organic horizon (mineral Hg per m is
  higher); this props up organic Hg:C with uptake still at 0.53.
- With the explicit export, soil Hg is near steady state after the spin-up jump (net storage +0.2 to +0.8), so the
  export is set by the inputs minus evasion; TVC export stays above Zolkos unless Kd >= 3.8.
- Kd_org has to be refit jointly in phase 5 after D3 (GPP uptake) raises Toolik uptake toward 4.1.

### Cost
Phase 2: 0.738 core-h (ledger 2.2372 -> 2.9755 of 7.537).

## 3. Phase 3: process soil re-emission, ground Hg(0) sink, stomatal foliar uptake (D3), melt routing, D8 (2026-11-24)

Commits (local feat/hg-t2, not pushed): e3f1c3c2 code, 4c0d6914 hg_gs_ingpp, 98aa09ee hg_gs_eff 0.75, then the docs commit.
Binary bins/binT2p3c (md5 362b542a..., GIT_SHA 98aa09ee). Runs: hg_t2/src/runP3.sh stages a-c (logs out/runP3_*.log).
All new switches default to 0 (hg_t2 behaviour of phase 2 unchanged); P3F = phase-2 defaults + hg_t2_evasion 1,
hg_t2_foliar 2, hg_t2_meltroute 1.

### What it does
- hg_t2_evasion 1: re-emitted fraction of decay-released Hg per layer = hg_fd (1 - f_snow) [T_layer > 0]
  exp(-z_mid / hg_ev_z) (replaces hg_f_evade; woody debris z = 0 with air T; geometry-change pot: decomposition-
  weighted). The hg0_soil column is replaced by a process ground Hg(0) sink v_ground x GEM x veg_scale x seconds, all
  months. New outputs HG_SOILEV (part of HG_EVASION), HG_SINK0 (part of HG_DEP).
- hg_t2_foliar 2 (D3, CLM5-Hg form): leaf Hg(0) uptake = GEM x integral of the Hg(0) stomatal conductance, with Medlyn
  g_sw = 1.6 m A / c_a from the PFT's monthly GPP (A), the month's CO2 and the molar air density (air T), Hg(0)/H2O
  diffusivity scaling sqrt(200.59/18.02) (Wesely 1989), x hg_gs_eff; plus the phase-1 LAI interception of dry Hg(II).
  hg_t2_foliar 1 (LAI scheme) stays available. hg_gs_ingpp 1 (default) uses the N-unlimited GPP (INGPP): dvmdostem's
  canopy conductance is not N-limited (N limitation scales GPP afterwards), and the N-limited GPP of the Toolik proxy
  run is 78 g C/m2/yr (INGPP 356) against several hundred observed for tussock tundra. Output HG_STOMATAL.
- hg_t2_meltroute 1 (needs hg_t2_diss): the month's snowmelt Hg (State::mpend, in the total and spin-up vector; 0 at
  month ends, folded into the top dissolved pool for the restart) is delivered day by day with the day's melt; the
  share qover / (qover + qinfl) leaves as overland flow (HG_MELTRUN, part of HG_RUNOFF), the rest enters diss[0].

### Parameters (model_settings; no parameters/ edits)
| name | value | source |
|---|---|---|
| hg_fd | 0.053 | ORCHIDEE-Hg fd_Hg (Smith-Downey et al. 2010 JGR 115 G03008) |
| hg_ev_z | 0.05 m | soil-air Hg(0) depleted below the top few cm, soils a net sink (Obrist et al. 2017 Nature 547:201) |
| hg_v_ground | 8.5e-5 m/s | Obrist 2017 winter net Hg(0) deposition 0.4 ng/m2/h (37 % of 6.5 ug/m2/yr) at GEM 1.3 |
| hg_gem | 1.3 ng/m3 | Arctic background GEM 1.2-1.5 (AMAP 2021), x the deposition series veg_scale |
| hg_gs_m | 4.0 | Medlyn 1 + g1/sqrt(D): g1 2.2-4.5 kPa^0.5 (CLM5 arctic grass / boreal shrub), D 0.5-1 kPa (Yuan et al. 2023 Environ Int 174:107904) |
| hg_gs_eff | 0.75 | fit here: Toolik stomatal uptake 5.32 at 1.0 -> 3.99 (target ~4.1, Obrist 2017) |
| hg_gs_ingpp | 1 | see above |
| D8 tissue ratios (sensitivity only) | stem 0.0044 / root 0.0067 ug/gC (shrubland medians 2 / 3 ng/g dw); 0.013 / 0.042 (means 6 / 19 ng/g) | Zhou & Obrist 2021 ES&T 55:14245, Table S2, / 0.45 g C per g dw |

### Sites (last 20 model years; ug/m2/yr; Hg:C ug/g). P3FP2 = phase-2 defaults (== P2F), P3F = phase 3 on
| site | run | R org / min | uptake total (stomatal) | export (overland / drain) | evasion (soil / snow) | net Hg(0) uptake | net storage | %MeHg JJAS | May-Jun export |
|---|---|---|---|---|---|---|---|---|---|
| TVC | P3FP2 | 0.43 / 1.28 | 5.30 (-) | 5.45 (1.53 / 3.92) | 4.51 | - | +0.65 | 0.98 | 25 % |
| TVC | P3F | 1.10 / 1.47 | 9.40 (4.83) | 13.66 (4.11 / 9.55) | 4.13 (1.98 / 2.15) | 4.2 | -1.14 | 0.97 | 26 % |
| HPC | P3FP2 | 0.33 / 1.53 | 2.79 (-) | 3.13 | 3.76 | - | +0.81 | 0.12 | 40 % |
| HPC | P3F | 0.51 / 1.80 | 5.13 (3.63) | 4.98 (1.10 / 3.88) | 4.78 (3.20 / 1.58) | 2.3 | +1.74 | 0.13 | 39 % |
| Toolik | P3FP2 | 0.35 / 1.29 | 0.53 (-) | 3.07 | 2.30 | - | +0.21 | 0.20 | 27 % |
| Toolik | P3F | 0.50 / 1.36 | 4.36 (3.99) | 4.79 (4.79 / 0) | 1.51 (0.09 / 1.42) | 6.0 (21 % Oct-mid May) | +4.20 | 0.20 | 33 % |

Net Hg(0) uptake = stomatal + ground sink - soil re-emission - snow re-emission (what a flux tower sees).
Targets: Toolik vegetation uptake ~4.1 (met, 3.99 stomatal / 4.36 total); Toolik net Hg(0) 6.5 +/- 0.7 (6.0, met);
winter share 37 % (21 %, low: snow re-emission Feb-May offsets the ground sink); organic Hg:C 0.29-0.41 (now high at all
sites: Toolik 0.50, HPC 0.51, TVC 1.10); export vs Zolkos 0.8-4.0 (HPC 4.98, Toolik 4.79, TVC 13.7: all above).

Single-process effects (gs_eff 1 runs, binT2p3/b): evasion alone (P3E) raises the input by the ground sink (+3.45) and
cuts the re-emission (Toolik 2.30 -> 1.49); stomatal alone with N-limited GPP (P3G): Toolik 1.17, TVC 6.2.
Melt routing (P3EG -> P3EGM): melt Hg 0.37 (Toolik) leaves almost entirely with the melt overland flow; the May-Jun
export share rises 27 -> 33 % (Toolik), TVC 25 -> 26 %, HPC unchanged; small, because 80-86 % of the snowpack Hg is
re-emitted before the melt (fig6).
D8 (P3D8a/b): Zhou & Obrist shrubland medians cut the stem/root uptake (T1 ratios 0.02 / 0.07) from 4.1 to 0.5 at TVC
(total uptake 9.40 -> 5.71, export 13.7 -> 11.7, R org 1.10 -> 0.94), 0.31 -> 0.03 at Toolik (4.36 -> 4.08), HPC 5.13 ->
4.26; the means give TVC 7.77. Stem/root uptake is 40 % of the TVC input with the T1 ratios.

### Checks
- Hg-off: sub_off PASS and full 8-case prod vs binM PASS (bit-identical), binT2p3c.
- hg_t2 = 0 (P3FT1) vs T1d VF4: all outputs bit-identical; only the restart hg_state length differs (as before).
- Phase-2 path (P3FP2, binT2p3c; and P3P2, binT2p3) vs P2F: every common output bit-identical (0 diffs).
- All phase-3 switches on (sub_t2e, 2 cases): C/N/water/energy identical to binM, ndiff 0; balance error max 1.3e-9 ug.
- Conservation: site balance errors <= 2e-9 ug (1e-14 relative).
- Unit tests: Hg suite 41/41 (2205 assertions, 6 new phase-3 cases); full unittest 523/523 cases.

### Findings and open items (for phase 5 / decisions)
- With the realistic Hg(0) inputs (stomatal + ground sink), every site now has more Hg input than export can carry under
  Zolkos-range yields at steady state: organic Hg:C and export both above target. Phase 5 needs a joint refit
  (Kd_org, hg_v_ground, hg_fd / hg_ev_z, D8 ratios); lowering Kd alone raises export further.
- TVC is the outlier (uptake 9.4, export 13.7): GPP 517, stem/root uptake 4.1, re-emission only 2.0. D8 medians remove
  most of the stem/root part.
- Winter net Hg(0) share 21 % vs 37 %: the snow re-emission (photoreduction in the photic layer) in Feb-May; the
  phase-1 hg_snow_kred is the candidate (10x lower changed little in phase 1 because the release form dominated; with
  D7 it is now the main term).
- Freshet: the melt routing works but moves little Hg; a larger spring pulse would need soil Hg mobilised in the
  first thaw days (top-layer desorption into the melt overland flow), not built.
- Default flips (hg_t2_evasion, hg_t2_foliar 2, hg_t2_meltroute under hg_t2) are left for phase 5 after the joint refit.

### Cost
Phase 3: 0.494 core-h (ledger 2.9755 -> 3.4692 of 7.537).

## 4. Phase 4: fire ash retention, RTS particulate export, canopy wash-off (2026-11-24)

Commit a94f93bb (code, local feat/hg-t2, not pushed) + docs commit. Binary bins/binT2p4 (md5 df11db94..., built at a94f93bb, clean).
Runs: hg_t2/src/runP4.sh (dist / dist2 / sites / ident / prod), case maker src/mk_p4.py, analysis src/an_p4.py, figure src/plots_p4.py.
All three items default off (hg_t2 phase-3 behaviour unchanged); they act with hg_t2 = 1 only.

### What it does
- hg_fire_ash (0-1, default 0): in the fire sub-step the emitted Hg is the pot's share of the fire C (T1). Of it, the part
  that came from burned soil (soil donors' share of the pot: removed layers + soil C losses) is kept at the fraction
  hg_fire_ash in the ash / char, added to the new top C layer (HG_ASH, internal); burned vegetation, standing dead and woody
  debris Hg is all emitted (Friedli et al. 2001 GRL 28:3223: essentially complete release from vegetation fuel). Linear in Hg.
- hg_t2_rts (0/1): RTS sub-step: the slumped material's Hg (pot share of RTS lateral C) is particulate export (HG_PARTIC, part of
  HG_EXPORT), and so is the Hg of the RTS-respired C share that is not re-emitted; T1 returned that Hg to the remaining column
  (receivers / donors). Mineral layers: dvmdostem's RTS removes C (organic and mineral) down to the removal depth; the
  tracer removes Hg with that C from the top of the old column, so the slumped mineral Hg leaves with the particulate flux
  at its layer Hg:C (no separate mineral term needed in this code; the design's "mineral Hg of removed layers" is covered).
- hg_wash (0-1, default 0), hg_wash_pref (mm/month, 10): the share hg_wash x P_rain / (P_rain + hg_wash_pref) of the month's
  canopy-intercepted dry Hg(II) (both canopy schemes) goes to the ground with the throughfall instead of the leaves
  (HG_WASH, part of HG_DEP); P_rain = the month's rainfall (m_a2l.rnfl). Design said "leaf Hg(II)": the tracer has no
  separate leaf Hg(II) pool, so the month's fresh surface deposit is washed (no new state, restart unchanged).

### Parameters (model_settings; no parameters/ edits)
| name | default | tested | source |
|---|---|---|---|
| hg_fire_ash | 0 | 0.1 (0.3 sensitivity) | Ku et al. 2018 ES&T 52:13366 (white ash >= 90 % of the Hg lost, black ash 34-83 % lost); Biswas et al. 2007 GBC 21 GB1005 (>95 % released above 300 C); Harden et al. 2004 GBC 18 GB3014 (FROSTFIRE forest floor ~70 % lost); range 0.03-0.3 |
| hg_t2_rts | 0 | 1 | St Pierre et al. 2018 ES&T 52:14099 (RTS Hg export is particulate) |
| hg_wash | 0 | 0.5 | Rea et al. 2000 ES&T 34:2418 (leaf-surface Hg(II) washes off); Wright et al. 2016 Environ Rev 24:205 (throughfall = wet + washed dry Hg(II)); assumption |
| hg_wash_pref | 10 mm/month | 10 | assumption (half-saturation rainfall) |

### Disturbance case (Toolik proxy CMT06, cell 0,0; sp 15 yr, tr 75 yr 1901-1975)
New scratch inputs hg_t2/inputs_p4: explicit fire 1926 (jday 196, severity 3) and RTS 1956 (jday 180; the repo file's RTS in
1902 moved), model_settings fire 1, thermokarst 1, tr dsb/fire true. dvmdostem: fire burns 3050 g C/m2 soil (0.113 m) and
20 g C/m2 veg; the RTS removes 0.129 m (TKLATC 6616 g C/m2 in 1956, then 1690, 1363, ... from the active slump; TKRESPC 625).

| run | settings | fire Hg emitted 1926 | ash | burned Hg:C | fire-year export | RTS 1956 export (particulate) | export 1956-75 | soil Hg 1975 (mg/m2) |
|---|---|---|---|---|---|---|---|---|
| P4DT1c | hg_t2 = 0 (T1) | 494.6 | - | 0.16 | 0.30 | 4208.0 | 18109 | 152.9 |
| P4DP3b | phase-3 settings (P3F) | 1876.3 | - | 0.61 | 24.0 | 3867.3 | 15946 | 132.1 |
| P4D1 | + ash 0.1, rts 1, wash 0.5 | 1690.0 | 186.6 | 0.61 | 82.7 | 4209.9 (4203.9) | 16278 (16215) | 131.9 |
| P4DA3 | + ash 0.3 | 1316.8 | 559.7 | 0.61 | 200.1 | 4361.6 (4355.3) | 16436 (16372) | 132.0 |

Units ug/m2 (per year where a year is given). Burned Hg:C = (emitted + ash) / burned C (ug/g).
- Ash: 99.4 % of the burned Hg is from soil, so ash = 0.0994 x the T1 fire Hg at 0.1. About 31 % of the ash Hg leaves in the
  fire year (export 24 -> 83 ug/m2): it is put in the thin new top layer, where the dissolved share is high. Ku et al. 2018
  find ash Hg largely recalcitrant, so this pulse is likely too large; a lower ash Hg mobility (own Kd) is a phase-5 option.
- RTS: the particulate flux is the T1 lateral export plus the respired-C Hg that T1 left in the column: +342 ug/m2 (+8.9 %)
  in 1956, +1.7 % over 1956-75. Particulate Hg:C 0.64 ug/g (soil + veg). The long tail is dvmdostem's ongoing lateral C of
  the active slump (also in T1).
- The T1 and T2 burned Hg:C differ (0.16 vs 0.61): T2 puts deposition and foliar litter Hg into the top layers.

### Wash-off at the sites (P4F vs P3F, last 20 years, ug/m2/yr)
| site | wash-off | foliar | dep | export | R org |
|---|---|---|---|---|---|
| TVC | 0.160 | 5.32 -> 5.16 | 7.24 -> 7.40 | 13.66 -> 13.66 | 1.100 -> 1.100 |
| HPC | 0.124 | 4.09 -> 3.96 | 6.37 -> 6.50 | 4.98 -> 4.98 | 0.514 -> 0.513 |
| Toolik | 0.017 | 4.05 -> 4.03 | 6.14 -> 6.16 | 4.79 -> 4.79 | 0.503 -> 0.503 |
Wash-off moves ~1/3 of the intercepted dry Hg(II) from the leaves to the ground; the budgets do not change.

### Checks
- Phase-4 switches off on binT2p4 vs binT2p3c, with fire + RTS (P4DP3b vs P4DP3a; P4DT1c vs P4DT1a for hg_t2 = 0): sp all 71
  variables, tr 327 of 330 bit-identical; the 3 others are the new HG_ASH / HG_PARTIC / HG_WASH files (the old binary has no
  values for them).
- Hg off: sub_off PASS (2 cases). The 8-case run (P4M1_prod) used by mistake the Hg-on case set hg_t1/cs/prod_hg (T1 path):
  RESULT FAIL only because the 10 Hg output files differ from binM (no tracer); all other outputs bit-identical in all 8 cases.
  The 8-case Hg-off run (regress/cs/prod) was not repeated (budget); it is left for the start of phase 5.
- All phase 3 + 4 switches on (sub_t2f, 2 cases): C/N/water/energy identical to binM, ndiff 0, balance error <= 1.3e-9 ug.
- Conservation: Toolik disturbance runs |HG_BALERR| <= 1.9e-9 ug/m2 (fire and RTS years included); sites <= 2e-9.
- Unit tests: Hg suite 46/46 (2303 assertions; 5 new: parse, ash retention + linearity, RTS particulate + linearity, wash-off in
  both canopy schemes, random steps with fire / RTS / daily transport and all phase 3 + 4 switches); full unittest 528/528 cases.

### Open items (phase 5)
- Ash Hg mobility (fire-year export pulse) and hg_fire_ash value; no observations at the three sites.
- RTS respiration uses hg_f_evade (0.02) for the re-emitted share (the fire / RTS sub-steps carry no layer temperatures).
- Defaults: hg_fire_ash, hg_t2_rts, hg_wash stay 0 until the phase-5 decision; they do not change the undisturbed site budgets
  except the small wash-off shift.

### Cost
Phase 4: 0.311 core-h (ledger 3.4692 -> 3.7805 of 7.537), of which 0.139 is the 8-case run on the wrong case set and 0.010 a
disturbance run cut short when I stopped the first chain (booked).

## 5. Phase 5a: TVC carbon spin-up (D4), joint refit, identifiability, k_meth (2026-11-24)

Phase 5 was split at a clean boundary. **5a** (this section) covers the D4 spin-up, the 64-member ensemble with the fit and
identifiability analysis, the verification runs, the k_meth scaling and all final checks. **5b** (optional, not run) would add a
refinement round with a wider hg_ev_z range and any structural option Elchin picks (see Open decisions). There are no code
changes, so the binary is still binT2p4 (md5 df11db94, a94f93bb). Scripts are in ~/tem-dev/hg_t2/src: d4_spin.py, mk_d4cs.py,
ens5.py, runP5.sh, col5.py, fit5.py, fit5b.py and plots_p5.py.

### D4: TVC carbon spin-up protocol (no carbon physics change)
- Method (d4_spin.py): chained 20-yr sc segments on the case's 1950-69 spin forcing, with Hg off (hg_tracer 0, hg_mehg 0).
- Each segment runs from a copy of the previous restart-sc.nc in a new directory, out/D4/TVC_Warming_nfeed_b/segNNN. Nothing is overwritten or reused.
- Stop rule: total soil C drift below 0.5 % per 100 yr.
- The first attempt, out/D4/TVC_Warming_nfeed (hg_mehg 1 without the tracer), failed. It is kept.
- Result: 34 segments (680 yr).

| | Start (restart0) | End (seg034) |
|---|---|---|
| Total soil C (kgC/m2) | ~67.8 | 59.2 |
| Organic-layer C (kgC/m2) | 16.8 | 10.4 |
| Mineral C (kgC/m2) | – | 48.8 |
| Organic thickness (m) | – | 0.28 |
| Drift (%/100 yr) | – | -0.46 |

- The vegc column in drift.json sums fill values and is not meaningful.
- cs/sites_d4 (mk_d4cs.py) copies hg_t1/cs/sites_hg with the paths rewritten. The TVC restart0.nc is a copy of the seg034 restart; its provenance is in case.json under restart0_src and d4. HPC is unchanged.
- D4 effect, P5B (phase-4 switches and defaults) against P4F:
  - TVC organic Hg:C fell from 1.10 to 0.81 and export from 13.7 to 10.4.
  - HPC and Toolik are bit-for-bit unchanged.
- Figure: figures_hg_t2/hg_t2_fig8_d4_spin.png.

### Joint refit and identifiability ensemble
- Switch set SW: mehg, t2, diss, evasion, foliar 2, meltroute, fire_ash 0.1, rts 1, wash 0.5.
- 64-member Latin hypercube (ens5.py, seed 11) over six parameters:

| Parameter | Range | Scale |
|---|---|---|
| hg_t2_log_kd_org | 3.0-4.6 | linear |
| hg_v_ground | 3e-5 to 3e-4 | log |
| hg_fd | 0.02-0.4 | log |
| hg_ev_z | 0.02-0.2 | log |
| D8 weight w | 0-1 | linear; r_stem = 0.0044 (0.02/0.0044)^w, r_root = 0.0067 (0.07/0.0067)^w |
| hg_snow_kred | 0.005-0.2 | log |

- Emulator (fit5.py): quadratic, one per observable, 18 observables (log for Hg:C, export, %MeHg and evasion).
  - Leave-one-out RMSE is at most 0.08 in log units against ensemble SDs of 0.03-1.0.
  - For Toolik net Hg(0) the LOO RMSE is 0.26 ug/m2/yr, against an SD of 3.2.
- Cost: the sum of squared normalised distances outside each target band (scale in brackets):
  - organic Hg:C 0.29-0.41 (0.10 in ln)
  - export 0.8-4.0 (Zolkos; 0.2 in ln)
  - Toolik uptake 4.1 (0.4)
  - Toolik net Hg(0) 6.5 (0.7)
  - Toolik winter share 0.37 (0.07)
  - mineral Hg:C 0.97-1.69 (0.15 in ln)
- Minimisation: random search over 2e5 points, then Nelder-Mead.
- Uncertainty (fit5b.py) comes from the near-optimal emulator set, cost <= best + 7.04 (the 68 % delta-chi2 for 6 parameters, n = 148 of 2e5). No point meets every band to within one scale, so there is no strict plausible set.

**Fitted set** (preset config/hg_t2_calibrated.json; the defaults are not changed):

| Parameter | Previous default | Fitted | 5-95 % (near-optimal) | Constrained by / source |
|---|---|---|---|---|
| hg_t2_log_kd_org | 3.5 | **3.46** | 3.36-3.55 | Organic Hg:C 0.29-0.41 (Obrist 2017 Toolik horizons) against Zolkos export |
| hg_v_ground (m/s) | 8.5e-5 | **1.1e-4** | 8.2e-5 to 1.25e-4 | Toolik net Hg(0) 6.5±0.7 (Obrist 2017) |
| hg_fd | 0.053 | **0.10** | 0.087-0.146 | Hg:C and evasion (trades off with ev_z) |
| hg_ev_z (m) | 0.05 | **0.20** | 0.10-0.20 (at the range bound) | Same |
| hg_r_stem / hg_r_root | 0.02 / 0.07 (T1) | **0.0044 / 0.0067** (D8a) | 0.0045-0.015 / 0.0069-0.045 | Toolik uptake 4.1 (weak); D8a = Zhou & Obrist 2021 wood/root |
| hg_snow_kred (1/h) | 0.090 | **0.090** (kept) | 0.007-0.12 | Not identifiable (cost profile flat, 7.3-10.8) |
| hg_k_meth_proxy | 87 (T1) | **71** | ±~5 % (linear scaling) | Thompson 2025: HPC 0.1 % (0.122 -> 0.100) and TVC 0.7 % (0.86 -> 0.70); both ask for ×0.81-0.82 |
| hg_k_meth_ch4 | 3250 (T1) | 3250 (kept) | – | Toolik CH4 fen case 2.06 % against the 2 % class |

- Verification: P5V (fitted set, k_meth_proxy 87, with the Toolik CH4 case) and P5V2 (k_meth_proxy 71).
- Emulator against model:

| Site | Quantity | Emulator | Model |
|---|---|---|---|
| TVC | Hg:C | 0.40 | 0.45 |
| TVC | Export | 5.9 | 6.5 |
| Toolik | Net Hg(0) | 6.44 | 6.70 |

- The emulator optimum used snow_kred 0.036 and the verification used 0.090.
- P5V2 differs from P5V only in %MeHg, which confirms that %MeHg is linear in k_meth and has no feedback on total Hg.

### Sites against targets (P5V2, last 20 model years; ug/m2/yr; Hg:C ug/g)

| | TVC (D4) | HPC | Toolik proxy | Toolik CH4 | Target | Met? |
|---|---|---|---|---|---|---|
| Organic Hg:C | 0.447 (P4F 1.10) | 0.276 (0.51) | 0.467 (0.50) | 0.467 | 0.29-0.41 | TVC/Toolik slightly high, HPC slightly low |
| Mineral Hg:C | 1.27 | 1.45 | 1.35 | 1.35 | 0.97-1.69 | yes |
| Export | 6.54 (13.7) | 2.98 (4.98) | 4.88 (4.79) | 4.88 | 0.8-4.0 (Zolkos) | HPC yes; TVC and Toolik high |
| Vegetation uptake | 5.42 | 4.14 | 4.06 | 4.08 | Toolik 4.1 | yes |
| Net Hg(0) | 2.17 | 0.69 | 6.70 | – | Toolik 6.5±0.7 | yes |
| Winter share of net Hg(0) | – | – | 0.28 | – | 0.37 | no (structural) |
| %MeHg (Jun-Sep active layer) | 0.70 | 0.100 | 0.165 | 2.06 | bog 0.7 / plateau 0.1 / fen 2 | yes |
| Evasion | 7.03 | 7.42 | 1.81 | 1.81 | – | – |
| Max Hg balance error (ug) | ≤ 3e-9 | | | | roundoff | yes |

TVC and HPC winter shares are not meaningful because net Hg(0) there is small. Figure: hg_t2_fig10_targets.png (ensemble spread,
target bands, and the P4F / P5B / P5V2 runs).

### Identifiability (hg_t2_fig9_identifiability.png: SRRC and near-optimal correlations)
1. **Kd_org** is well constrained (3.36-3.55). It raises Hg:C and lowers export at every site (SRRC +0.8-1.0 and -0.74 to -0.79), so Hg:C and export trade off along it.
   - The best compromise leaves Hg:C at TVC and Toolik about 10-15 % above the band, Hg:C at HPC 5 % below it, and export at TVC (6.5) and Toolik (4.9) above 4.0.
   - The site-to-site spread (HPC/Toolik Hg:C ratio about 0.6) is not controlled by any of the six parameters. That is a structural misfit, not a parameter one.
2. **v_ground** controls both Toolik net Hg(0) (SRRC 0.96) and the winter share (0.97). Reaching the 37 % winter share needs about 1.6e-4, which gives net Hg(0) of about 9 (above 7.2).
   - Along v_ground with the other parameters at their best values: 7.5e-5 gives 5.1 / 21 %, 1.2e-4 gives 7.1 / 31 %, 1.5e-4 gives 8.4 / 36 %.
   - snow_kred barely moves the winter share (SRRC -0.12), so it cannot carry it. The 37 % share needs a structural winter term, for example the sink under snow or D7.
3. **fd and ev_z** trade off (near-set correlation -0.48). Both set evasion and the Hg:C reduction at TVC and HPC. ev_z sits at the 0.2 m upper bound, so its upper uncertainty is open (5b).
4. **D8 stem/root ratios** act only on Toolik uptake (SRRC 1.0), and only weakly (4.06-4.34 over the whole range). The fit prefers the low D8a end, but anything up to w ≈ 0.8 stays near-optimal.
5. Other correlations: v_ground-fd +0.35, Kd-v_ground -0.26, ev_z-w +0.38 (from shared Hg:C and uptake). The Laplace (Hessian) correlations agree in sign. Its weakest direction is snow_kred (sd 0.38 in unit-cube units).

### Recommended switch defaults (proposals only; nothing changed in code)
- **Keep as is:** the hg_t2 master switch stays at 0 (D5, Elchin's decision). The T1 and Hg-off paths are untouched.
- **When hg_t2 = 1, proposed defaults** for the T2-only switches and parameters (none affect T1 or Hg off):

| Setting | Proposed default |
|---|---|
| hg_t2_diss | 1 |
| hg_t2_evasion | 1 |
| hg_t2_foliar | 2 |
| hg_t2_meltroute | 1 |
| hg_t2_rts | 1 |
| hg_fire_ash | 0.1 |
| hg_wash | 0.5 (hg_wash_pref 10) |
| hg_t2_log_kd_org | 3.46 |
| hg_v_ground | 1.1e-4 |
| hg_fd | 0.10 |
| hg_ev_z | 0.20 |
| hg_snow_kred | 0.090 (kept) |

- **Shared with T1, preset only:** r_stem 0.0044, r_root 0.0067 and k_meth_proxy 71. Either add T2-specific copies used only when hg_t2 = 1, or change the T1 values (a T1 recalibration). This is Elchin's decision.

### Checks (all on binT2p4, which is unchanged in phase 5)

| Check | Run / output | Result |
|---|---|---|
| Hg off, 8 cases (binT2p4 vs binM, regress/cs/prod) | P5I0 | PASS. Not repeated at the end of 5a because the binary and cases are unchanged (md5 re-checked). |
| Hg off, 2-case set | P5I1 sub_off | PASS |
| hg_t2 = 0 vs T1 (P5T1 vs VF4) | an2/cmp_P5T1_VF4.json | All common variables bit-identical. The only difference is the restart hg_state length (117 vs 94, the T2 slots). |
| hg_t2 = 0 vs P3FT1 | – | 0 differences |
| C/N/water/energy, all T2 on with the fitted set | P5I2 cs/sub_t2g + cmp_on | Identical (True), PASS |
| Conservation | 64 members + P5B/P5V/P5V2 | Max Hg balance error 2.5e-9 ug |
| Unit tests | dvmdostem_unittest | 528/528 cases (8 known may-fail assertions); Hg subset 46/46 (2303 assertions) |
| Core use | – | ≤4 processes, except for an overlap of about 7 s when the unit-test binary ran next to a 4-process identity run |

### Cost (core-h)

| Item | Core-h |
|---|---|
| P5I0 | 0.136 |
| D4 | 0.057 |
| P5B | 0.029 |
| E5 (64 members) | 1.878 |
| P5V (with CH4) | 0.045 |
| P5T1 | 0.029 |
| P5V2 | 0.029 |
| P5I1 + P5I2 | 0.068 |
| **Phase 5a total** | **2.270** |

The ledger stands at 6.048 of 7.537. That leaves 1.49, or 1.19 above the 0.3 reserve.

### Open decisions
1. Residual misfits:
   - TVC export 6.5 and Toolik export 4.9 (target ≤4.0).
   - Organic Hg:C at TVC 0.45 and Toolik 0.47 (≤0.41) and at HPC 0.28 (≥0.29).
   - Toolik winter share 28 % (37 %).

   Options are to accept them, or to make a structural change: a site- or moss-dependent Kd, or an export-path control (overland fraction / tau_sorb); and a winter ground Hg(0) sink under snow, or D7.
2. D5: whether hg_t2 becomes default-on with the tracer.
3. Shared T1 parameters (r_stem, r_root, k_meth_proxy): T2-specific copies, or a T1 recalibration.
4. 5b: a refinement round with hg_ev_z extended to 0.5 m, about 16 members ≈ 0.47 core-h, within the remaining budget.
5. D4: whether the 680-yr chained TVC restart (organic C 16.8 -> 10.4 kgC/m2) should become the TVC case restart for other work. So far it is used only in cs/sites_d4.

## 6. Phase 5b: structural options, preset for hg_t2 = 1, T2-only copies, D4 standard restart, refit (2026-11-24)

Elchin's decision (verbatim): "Go with your recommendations and do the structural fix within the remaining budget".
Commits on local feat/hg-t2 (not pushed): d787f0fd (code + tests, defaults off), 4339a0fc (fitted preset values +
config/hg_t2_calibrated.json), plus this docs commit. Binaries: binT2p5b (md5 34cb6d5a, d787f0fd), binT2p5c
(md5 b2e878b3, 4339a0fc, final). Scripts in ~/tem-dev/hg_t2/src: ens6.py, runP5b.sh, runP5b_chain1.sh, fit6.py, hg0_lin.py,
fit6c.py, runP5b_final.sh, plots_p6.py (runP5.sh gained the ENSGEN option; old copy src/runP5.sh.p5a).

### Structural changes (include/HgTracer.h; all inside hg_t2, each behind its own switch)
1. **Humic-layer Kd** (hg_t2_kdsplit 0/1, hg_t2_dlog_kd_hum in [-3, 3]): in diss_transport_, layers of type 2
   (deep/humic) use log10 Kd = hg_t2_log_kd_org + dlog; moss (0) and fibric (1) layers keep hg_t2_log_kd_org. This is the
   moss/fibric-vs-humic Kd tied to the organic-layer type (needs hg_t2_diss). dlog = 0 is bit-identical to off (unit test).
2. **Explicit export-path control** (hg_t2_xpath 0/1, hg_f_conn in [0, 1]): only the share f_conn of the organic-layer
   overland + lateral dissolved flux (including the melt overland share) leaves the column; the rest goes to the first
   mineral layer's dissolved pool (booked in the downward flux, so the column balance is unchanged). f_conn = 1 is
   identical to off (unit test).
3. **Ground Hg(0) sink under snow** (hg_t2_snowsink 0/1, hg_v_snow m/s): with snow fraction fs = min(1, SWE / swe_full),
   sink = ((1 - fs) v_ground + fs v_snow) x GEM x veg_scale x dt; off: v_ground in all months (phase 3; see the correction in section 8). Needs hg_t2_evasion.
4. **hg_t2 = 1 takes the calibrated preset** (hg_t2_preset, default 1; T2-only). t2_preset_() sets the fitted switches and
   values before the explicit keys are read, so any explicit hg_* key wins; hg_t2_preset 0 gives the old component defaults
   (used by the older unit tests). Dependent switches follow an explicit parent switch-off (e.g. hg_t2_diss 0 also turns
   meltroute/kdsplit/xpath off unless set). The master switch hg_t2 stays default 0 (D5 not flipped).
5. **T2-only copies** hg_t2_r_stem, hg_t2_r_root, hg_t2_k_meth_proxy: with hg_t2 = 1 they replace r_stem / r_root /
   k_meth_proxy; setting the T1 names together with hg_t2 = 1 is a parse error (no silent ambiguity). hg_t2 = 0 never reads
   them, so the T1 path is untouched (checked bit-identical below). k_meth_ch4 (3250) remains shared (fen class reproduces).
6. Unit tests (tests/unit/test_hg_tracer.cpp): 5 new cases (preset/parse, snow sink, kdsplit, xpath incl. melt share and
   retention, random all-on balance); standalone 51/51 cases (2450 assertions); full suite 533/533 cases (8 known may-fail).

### D4 restart as the standard TVC starting point
restarts/TVC_Warming_nfeed_D4_680yr_restart-sc.nc (copy of out/D4/TVC_Warming_nfeed_b/seg034/output/restart-sc.nc,
md5 450a15da, restarts/README.md) is the TVC restart0.nc of the standard case set cs/sites_d4 (cs/STANDARD.md; default
SRC_CS of runP5.sh). HPC is unchanged. Old restarts and T1 / identity case sets are untouched; every run copies into new dirs.

### Refit (E6 + exact Hg(0) solve)
- P5R (binT2p5b, new switches off, 5a values) reproduced P5V2 bit-identically.
- E6: 26-member LHS (seed 23, ~0.76 core-h, budget-limited from the planned 32) on binT2p5b with kdsplit, xpath and
  snowsink all on, 7 parameters: kd_org 3.0-4.2, dlog_kd_hum -0.5-1.5, f_conn 0.2-1, v_ground 3e-5-2e-4 (log),
  r_snow = v_snow / v_ground 0.5-8 (log), fd 0.03-0.3 (log), ev_z 0.05-0.5 m (log; the extended range). Max balance error 2.6e-9.
- fit6: quadratic emulator on E5 (64, structural terms 0) + E6 (26) = 90 points. Ro / export / Rm are well emulated
  (LOO RMSE in log 0.04-0.18 vs SD 0.5-1.0); Toolik net Hg(0) is not (LOO 3.5 vs SD 6).
- hg0_lin: the Toolik ground sink is exactly linear, sink0_m = a_m v_ground + b_m v_snow (residual 5e-16; snow covers
  Oct-May). fit6c therefore solves v_ground and v_snow exactly for net Hg(0) 6.5 and winter share 0.37 at every candidate and
  fits kd_org, dlog, f_conn, fd, ev_z on the emulator (cost bands as fit5).
- Best (fit6c): kd_org 3.453, dlog_kd_hum 0.242, f_conn 1.0, fd 0.083, ev_z 0.50 (upper bound), v_ground 4.88e-5,
  v_snow 1.283e-4 (r_snow 2.63). Emulator cost 0.32 (5a structure on the same emulator: 3.05).
- Near-optimal set (n = 15; 5-50-95 %): kd_org 3.04-3.26-3.46, dlog -0.21-0.17-0.84, f_conn 0.35-0.73-0.98 (corr with
  kd_org 0.94: export path and Kd substitute), fd 0.058-0.088-0.103, ev_z 0.31-0.44-0.49 (fd~ev_z -0.55). Fixed ev_z:
  0.1 -> cost 3.87, 0.2 -> 4.0, 0.3 -> 2.2; at shallow ev_z the fit uses xpath (f_conn 0.36) and a larger dlog.
- v_ground / v_snow are set by the two Toolik Obrist constraints: the solve maps net 6.5 with winter 0.30-0.44 to
  v_ground ~9e-5 to ~2e-5 and v_snow ~1.1e-4 to ~1.5e-4, so the winter-share uncertainty dominates v_ground.
- Verification P5G (explicit best set, --ch4) came out above the emulator for Ro / export (extrapolation at ev_z 0.5).
  P5H = P5G + xpath f_conn 0.6 lowered export (TVC 4.49, HPC 1.74, Toolik 3.99) but raised organic Hg:C (TVC 0.524, Toolik
  0.655): the retained Hg goes to the first mineral layer and is mixed back by cryo/bioturbation. Cost: P5H 28.4, P5G 6.86,
  5a P5V2 11.37 -> the preset is P5G (xpath available but off). Figures 11-13.

### Results (P5G = preset; last 20 yr; ug/m2/yr; Ro organic Hg:C ratio to T0 basis, Rm mineral)

| Quantity | TVC | HPC | Toolik proxy | Target | 5a (P5V2) |
|---|---|---|---|---|---|
| Ro | 0.409 | 0.268 | 0.448 | 0.29-0.41 | 0.447 / 0.276 / 0.467 |
| Rm | 1.254 | 1.430 | 1.351 | 0.97-1.69 | 1.268 / 1.445 / 1.354 |
| Export | 6.15 | 2.71 | 4.77 | 0.8-4.0 | 6.54 / 2.98 / 4.88 |
| Vegetation uptake | 5.42 | 4.14 | 4.06 | Toolik 4.1 | same |
| Net Hg(0) | 1.65 | 0.46 | 6.40 | Toolik 6.5 +- 0.7 | 2.17 / 0.69 / 6.70 |
| Winter share of net Hg(0) | - | - | 0.367 | 0.37 | 0.283 |
| Soil re-emission | 4.98 | 5.36 | 0.40 | - | |
| Ground sink | 4.05 | 3.77 | 4.23 | - | |
| %MeHg (JJAS active) | 0.709 | 0.100 | fen case 2.05 | 0.7 / 0.1 / 2 | 0.704 / 0.100 / 2.06 |

Fixed by 5b: Toolik winter share (28 -> 37 %) with net Hg(0) kept in band; TVC Ro now at the band edge (0.409); export down
4-9 %. Still outside: TVC export 6.15 and Toolik export 4.77 (> 4.0), Toolik Ro 0.448 (> 0.41), HPC Ro 0.268 (< 0.29).

### Checks (all on binT2p5c unless noted)

| Check | Result |
|---|---|
| hg_t2 = 1 preset only (P5F) vs explicit P5G | 0 differing variables (all 6 site/stage dirs) + CH4 fen case 0/379 |
| hg_t2 = 0 (P5T1b) vs T1 final (VF4) | only restart hg_state length (117 vs 94, the T2 restart extension) |
| P5T1b vs P5T1 (5a, binT2p4) | 0 differences |
| P5R (binT2p5b, new switches off) vs P5V2 | 0 differences |
| P5I4 Hg off 2-case (sub_off) | RESULT: PASS |
| P5I5 all on via preset (sub_t2h) | RESULT: PASS; C/N/water/energy identical with tracer on: True |
| P5I6 Hg off 8-case prod | RESULT: PASS |
| Balance error | <= 1e-9 ug/m2 in every run (E6 max 2.6e-9) |
| Unit tests | 533/533 (8 known may-fail); Hg standalone 51/51 (2450 assertions); run only inside the builds, no overlap |
| Cores | <= 4 processes; segments <= 75 yr; nohup; no parameters/ edits; no deletions |

### Cost (core-h)

| Item | Core-h |
|---|---|
| P5R | 0.029 |
| E6 (26) | 0.760 |
| P5G (with CH4) | 0.044 |
| P5H | 0.029 |
| P5F (with CH4) + P5T1b | 0.074 |
| P5I4 + P5I5 + P5I6 | 0.201 |
| **Phase 5b total** | **1.136** |

Ledger 7.184 of 7.537 (0.353 left, 0.053 above the 0.3 reserve). Phase 5 total 3.41; T2 total 5.125.

### Open items after 5b
1. Export at TVC (6.2) and Toolik (4.8) and the Ro edges remain; neither the humic Kd split nor flow-path retention fixes
   export without raising organic Hg:C, because retained Hg is mixed back. Next structural option: route retained Hg below
   the cryo/bioturbation zone or into a non-mixing (sorbed-mineral) pool, or the DOC-linked export (D6).
2. ev_z sits at the 0.5 m bound (deep re-emission is physically doubtful; Obrist soil-gas depletion is near-surface).
   ev_z 0.3 costs 2.2 (emulator) and is the conservative alternative.
3. dlog_kd_hum is weakly identified (near-optimal range spans 0); kdsplit could be dropped without much loss.
4. v_ground hinges on the winter-share target (0.30-0.44 -> v_ground 9e-5 to 2e-5).
5. D5: hg_t2 default stays 0. k_meth_ch4 is still shared with T1. Configs that set hg_r_stem / hg_r_root /
   hg_k_meth_proxy together with hg_t2 = 1 now fail at parse (use the hg_t2_* names).
6. The emulator uses 26 new members only (budget); a further ensemble would sharpen dlog / f_conn.

## 7. T2 summary (phases 1-5b, 2026-11-24)
T2 replaces T1's surface and transport closures with process pieces, all inside include/HgTracer.h and behind hg_t2
(default 0; needs hg_tracer = 1). With hg_t2 = 0 the T1 step is unchanged (bit-identical to T1 VF4 except the restart
vector length), and Hg off stays bit-identical to the base model (8-case prod). Nothing in T2 feeds back on C/N/water/energy
(checked identical with everything on). The tracer is affine and exactly conservative (balance at roundoff, spin-up jump).

| Phase | Content | Key commits |
|---|---|---|
| 1 | Foliar uptake (LAI), snowpack Hg with photoreduction / sublimation (D7), wet deposition timing, restart v3, outputs | (phase-1 commits, section 1) |
| 2 | Dissolved pool, daily sorption + overland / lateral / vertical export, cryo/bioturbation, moss sorbent; Kd_org refit | c32cd335, ecc13285, ae1a200b |
| 3 | Process soil re-emission (fd, temperature, snow, ev_z), ground Hg(0) sink, stomatal foliar uptake (D3), melt routing, D8 | e3f1c3c2, 4c0d6914, 98aa09ee |
| 4 | Fire ash retention, RTS particulate export, canopy wash-off | a94f93bb |
| 5a | D4 TVC carbon spin-up (680 yr), 64-member ensemble, joint refit, k_meth classes | cb018895, b7dbc42a |
| 5b | Humic Kd, export-path control, winter snow sink, hg_t2 = 1 preset, T2-only r_stem / r_root / k_meth_proxy, D4 standard | d787f0fd, 4339a0fc, docs |

Against T1d: TVC organic Hg:C 2.63 -> 0.41, HPC 1.06 -> 0.27; TVC export 4.26 -> 6.15 (now from explicit daily
pathways), HPC 5.80 -> 2.71, Toolik 0.24 -> 4.77; Toolik uptake 4.06 and net Hg(0) 6.40 with a 37 % winter share match
Obrist 2017; %MeHg classes match Thompson 2025. Main remaining gap: export at TVC / Toolik vs Zolkos-type yields and the
Ro edges (open item 1 above). Cost: T2 5.125 core-h (ledger 7.184 / 7.537). Use: hg_tracer 1, hg_mehg 1, hg_t2 1 (preset),
case set cs/sites_d4 (D4 TVC restart). Figures: docs_src/figures_hg_t2/hg_t2_fig1-13.

## 8. Phase 5c: flow-path retention routed below the mixing zone (option A), refit (2026-11-24)

Elchin's decision (verbatim): "Go with option A, below-mixing-zone routing, up to 1.2 core-h". Commits on local
feat/hg-t2 (not pushed): c19c2b9c (code + tests, default off), 2972fc49 (preset values + config/hg_t2_calibrated.json), plus
this docs commit. Binaries: binT2p5d (md5 7da12750, c19c2b9c), binT2p5e (md5 54a45cff, 2972fc49, final). Scripts in
~/tem-dev/hg_t2/src: runP5c.sh (phase-5c cap), ledger5c.py, runP5c_probe.sh, runP5c_fd.sh, gn7.py, runP5c_v1.sh,
runP5c_v2.sh, runP5c_final.sh, plots_p7.py. Backups of the 5b sources: p5c_src/*.p5b.
Correction to section 6: with hg_t2_snowsink = 0 the phase-3 ground sink uses v_ground in every month regardless of snow
(equivalent to v_snow = v_ground), not (1 - f_snow) v_ground.

### Why the 5b retention raised Hg:C, and the change
In 5b (P5H) the retained share 1 - f_conn of the organic-layer overland / lateral export went to the first mineral layer,
just below the organic horizon and inside the thawed zone. Cryo/bioturbation (diffusive in Hg per unit depth) then carried it
back up, so organic Hg:C rose (TVC Ro 0.41 -> 0.52, Toolik 0.45 -> 0.66) while export fell. A probe with the receiving depth
at 0.5 m (P5K, still in the active layer, ALD 0.6-1.0 m) behaved the same (Ro 0.52 / 0.29 / 0.65).

New switch (include/HgTracer.h, inside hg_t2, needs hg_t2_xpath = 1; component default 0):
- hg_t2_deeproute 0/1, hg_f_deep in [0, 1] (default 1), hg_z_deep (m, default 1.0): the share f_deep of the retained Hg goes
  to the dissolved pool of the layer containing depth z_deep (never above the first mineral layer; the deepest layer if the
  column is shallower); the rest stays with the first mineral layer as in 5b. Booked in HG_DOWN like the 5b retention.
- With z_deep below the active layer the Hg sorbs in frozen mineral soil, outside the monthly cryo/bioturbation of the thawed
  zone and below the daily drainage, so it does not return to the organic layers. It is a 1-D stand-in for retention at the
  base of the active layer / transient layer and its incorporation into permafrost.
- Geometry-only coefficients: the step stays affine and exactly conservative. Off, f_deep = 0 or z_deep above the first
  mineral layer are bit-identical to 5b (unit tests).
- Parse: needs hg_t2 = 1 and hg_t2_xpath = 1; under the preset an explicit hg_t2_xpath 0 also turns it off.
- Unit tests (tests/unit/test_hg_tracer.cpp, 4 new cases): parse / dependencies / ranges; routing to the layer at z_deep, the
  f_deep split and the identities above; below a frozen layer with strong mixing the organic layers end exactly as with no
  retention and all retained Hg is in the frozen layer; a random all-on balance test with geometry changes.
  Standalone Hg tests 55/55 (2579 assertions); full suite 537/537 (8 known may-fail assertions).

### Refit (local Gauss-Newton / Broyden on site runs)
The 5b emulator overpredicted near the optimum, so 5c fits directly on site runs:
- Probe P5J = 5b preset + xpath f_conn 0.6 + deep route z_deep 1.5 m: Ro unchanged at all sites (0.409 / 0.268 / 0.448, as P5G),
  export 3.69 / 1.62 / 2.86 (from 6.15 / 2.71 / 4.77), Rm 1.47 / 1.58 / 1.46 (from 1.25 / 1.43 / 1.35). The mechanism works as
  intended. z_deep only needs to be below the active layer, so it is fixed at 1.5 m (TVC / HPC / Toolik ALD <= 0.78 / 0.99 / 0.76 m).
- Jacobian: 7 one-at-a-time runs around P5J (kd_org -0.15, dlog_kd_hum +0.3, f_conn 0.45, fd x1.25, v_ground x1.3, v_snow x1.2,
  ev_z 0.3). In log units per step: kd_org moves Toolik Ro most (+0.30 per +0.15) and HPC Ro least (+0.07); f_conn moves only
  export (+0.29 per +0.15) and Rm; fd lowers Ro and Rm at TVC / HPC; v_ground / v_snow set Toolik net Hg(0) and winter share;
  dlog_kd_hum has no effect on any target (|effect| <= 0.05), so it is not fitted (kept at 0.242).
- Proposal: minimise the 12-band cost (fit5 bands) of the linear model inside a trust radius of 1.5 steps, bands shrunk by
  3-5 % of their width to aim inside; two verification rounds with Broyden updates (V7a/b, V8a/b), each with ev_z free (<= 0.5)
  and ev_z fixed at 0.3 (side comparison). Run costs (12 bands): P5G 6.86, P5J 1.46, V7a 0.013, V7b 0.039, V8a 0.008, V8b 0.014.
- Fitted set (V8a = preset): t2_log_kd_org 3.377, t2_dlog_kd_hum 0.242 (kdsplit on, not identified), hg_t2_xpath 1,
  hg_f_conn 0.53582, hg_t2_deeproute 1, hg_f_deep 1, hg_z_deep 1.5, hg_fd 0.068765, hg_ev_z 0.50, hg_v_ground 4.3731e-5,
  hg_v_snow 1.3038e-4; unchanged: snow_kred 0.090, t2_r_stem / r_root 0.0044 / 0.0067, t2_k_meth_proxy 71, k_meth_ch4 3250.
- ev_z: with the deep route ev_z is no longer constrained. V8b (ev_z 0.3, kd 3.355, f_conn 0.513, fd 0.0724, v_ground 4.19e-5,
  v_snow 1.303e-4) fits as well (cost 0.014; HPC Rm 1.71, 1 % above the band). It is the physically shallower option.
- The binding pair is HPC organic Ro (>= 0.29) vs HPC mineral Rm (<= 1.69): fd moves both the same way and the deep route adds
  to Rm. Both sit on their band edges in V8a.

### Results (last 20 yr; ug/m2/yr)

| Quantity | Target | TVC 5b -> 5c | HPC 5b -> 5c | Toolik 5b -> 5c |
|---|---|---|---|---|
| Ro (organic Hg:C) | 0.29-0.41 | 0.409 -> 0.407 | 0.268 -> 0.291 | 0.448 -> 0.391 |
| Rm (mineral Hg:C) | 0.97-1.69 | 1.254 -> 1.545 | 1.430 -> 1.689 | 1.351 -> 1.477 |
| Export | 0.8-4.0 | 6.15 -> 3.83 | 2.71 -> 1.85 | 4.77 -> 2.64 |
| Vegetation uptake | Toolik 4.1 | 5.42 -> 5.42 | 4.14 -> 4.14 | 4.06 -> 4.06 |
| Net Hg(0) | Toolik 6.5 +- 0.7 | 1.65 -> 2.52 | 0.46 -> 0.97 | 6.40 -> 6.50 |
| Winter share of net Hg(0) | Toolik 0.37 | - | - | 0.367 -> 0.370 |
| %MeHg (JJAS active layer) | 0.7 / 0.1 / fen 2 | 0.709 -> 0.709 | 0.100 -> 0.100 | fen 2.05 -> 2.03 |
| Evasion | - | 7.13 -> 6.24 | 6.94 -> 6.39 | 1.82 -> 1.71 |
| Column storage change | - | 0.13 -> 3.32 | 1.30 -> 2.66 | 4.36 -> 6.59 |

All 12 calibration bands are met (V8a cost 0.008, against 6.86 for the 5b preset). The rise in column storage is the routed
Hg accumulating in frozen mineral soil around 1.2-1.7 m (figure 16): a permafrost sink of about 2.3 (TVC), 0.9 (HPC) and
2.1 (Toolik) ug/m2/yr that the 5000-yr spin-up jump does not bring to equilibrium.

### Checks (binT2p5e unless noted)

| Check | Result |
|---|---|
| hg_t2 = 1 preset only (P5F2) vs explicit V8a (binT2p5d) | 0 differing variables (all 6 site/stage dirs) |
| P5R2 (binT2p5d, 5b preset) vs P5G | 0 differences (deep route off == 5b) |
| hg_t2 = 0 (P5T1c) vs T1 final (VF4) | only the restart hg_state length (117 vs 94) |
| P5T1c vs P5T1 (5a) | 0 differences |
| P5I7 Hg off 2-case (sub_off) | RESULT: PASS |
| P5I8 all on via the preset (new case set cs/sub_t2i) | RESULT: PASS; C/N/water/energy identical with the tracer on: True |
| P5I9 Hg off 8-case prod | RESULT: PASS |
| Balance error | <= 2e-9 ug/m2 in every run |
| Unit tests | 537/537 (8 known may-fail) in both builds; Hg standalone 55/55 (2579 assertions); run only inside the builds or with nothing else running |
| Rules | <= 4 processes; -j2 builds; segments <= 75 yr; nohup; no parameters/ edits; no deletions; restarts copied, never reused |

### Cost (core-h, ledger section "PHASE 5c" in LEDGER.txt)

| Item | Core-h |
|---|---|
| Probes P5R2, P5J, P5K | 0.088 |
| One-at-a-time runs F7k/h/c/f/g/s/e | 0.205 |
| Verification V7a/b, V8a/b | 0.117 |
| Final chain (P5F2 with CH4, P5T1c, P5I7-9) | 0.273 |
| **Phase 5c total** | **0.683 of 1.2** |

Ledger total 7.867 (5c start 7.184).

### Open items after 5c
1. The routed Hg is a permafrost sink that keeps growing: column storage rises by 2-3 ug/m2/yr at TVC / HPC. It needs a check
   against deep permafrost Hg profiles (Schuster 2018, Olson 2018) and a decision on whether the share of retention reaching
   the permafrost (f_deep < 1) should be lower. Under thaw (transient warming) this Hg becomes mobile, which is a real
   feedback, but it is untested here.
2. HPC organic Ro and mineral Rm are both on their band edges: the 1-D column cannot raise HPC organic Hg:C without raising
   mineral Hg:C.
3. ev_z is unconstrained (0.3 and 0.5 fit equally); the preset keeps 0.5 (best cost), V8b (0.3) is the physically shallower
   alternative.
4. t2_dlog_kd_hum has no effect: kdsplit could be dropped (structure simplification).
5. f_conn 0.54 means 46 % of the organic-layer export is retained, a single value for all sites. There is no site or hydrology
   dependence (riparian / hillslope connectivity), and Toolik export (2.6) now sits below the Kuparuk-type yields used earlier.
6. Unchanged: D5 (hg_t2 stays default 0), k_meth_ch4 shared with T1, D6 DOC export not implemented.

## 9. Phase 5d: cleanup (ev_z 0.3, humic offset 0, T2 copy of k_meth_ch4, refit) (2026-11-24)

Elchin agreed ("I agree") to the 5c recommendations. Commits on local feat/hg-t2 (not pushed), base 81a6fc5d:
91b456de (hg_t2_k_meth_ch4 copy + test), f60c7170 (preset + config/hg_t2_calibrated.json), plus this docs commit. Binary
binT2p5f (md5 5fb01b66, f60c7170). Scripts: src/gn8.py (gn7 + forced expansion point, band weights, f_deep column),
runP5d.sh (5d cap, Toolik-only vtlk), runT2tlk.sh, ledger5d.py, runP5d_final.sh, plots_p8.py. Backups: p5d_src/*.p5c.

### Changes
1. Preset ev_z 0.5 -> 0.3 (the physically shallower value; the targets do not constrain it).
2. Preset hg_t2_dlog_kd_hum 0.242 -> 0. The switch hg_t2_kdsplit and its code stay; with offset 0 it is bit-identical to off.
3. T2-only copy hg_t2_k_meth_ch4 (default 3250, the shared value, so all results are unchanged). With hg_t2 = 1 it replaces
   hg_k_meth_ch4; giving hg_k_meth_ch4 together with hg_t2 = 1 is a parse error, like the other copies. hg_t2 = 0 never reads it.
   config/hg_t2_calibrated.json now uses the T2 name. New unit test; standalone Hg tests 56/56 (2592 assertions).
4. Refit without new exploration runs (5d budget). gn8 reuses the 5c Jacobian (P5J + 7 one-at-a-time runs). It takes P5H
   (retention to the first mineral layer) as the f_deep = 0 column, Broyden-updates with V7a/b and V8a/b, and expands at V8b
   (ev_z 0.3). ev_z is fixed at 0.3, dlog at 0, z_deep at 1.5.
   - f_deep: lowering it sends retained Hg into the mixed zone. Per unit f_deep that raises Ro by 0.25 / 0.09 / 0.38
     (log, TVC / HPC / Toolik) but lowers Rm only by 0.12 / 0.07 / 0.04. Rm is not reduced selectively, so the fit keeps f_deep 1.
   - With ev_z 0.3 the linear model has no point that meets all 12 bands. The best cost is 0.0105, and one of TVC Ro, TVC
     export, HPC Ro or HPC Rm stays out by about 1 % (weighting tests). The preset keeps 11 in band with HPC Rm about 1.70.
   - Fitted set: t2_log_kd_org 3.356, dlog_kd_hum 0, f_conn 0.53419, f_deep 1.0, z_deep 1.5, fd 0.071857, ev_z 0.30,
     v_ground 4.1389e-5, v_snow 1.3046e-4; unchanged: snow_kred 0.090, t2_r_stem / r_root 0.0044 / 0.0067,
     t2_k_meth_proxy 71, t2_k_meth_ch4 3250.

### Results (P5F3 = hg_t2 = 1 preset only, binT2p5f; last 20 yr; ug/m2/yr)

| Quantity | Target | TVC 5c -> 5d | HPC 5c -> 5d | Toolik 5c -> 5d |
|---|---|---|---|---|
| Ro (organic Hg:C) | 0.29-0.41 | 0.407 -> 0.410 | 0.291 -> 0.288 (low) | 0.391 -> 0.377 |
| Rm (mineral Hg:C) | 0.97-1.69 | 1.545 -> 1.562 | 1.689 -> 1.706 (high) | 1.477 -> 1.477 |
| Export | 0.8-4.0 | 3.83 -> 3.995 | 1.85 -> 2.08 | 2.64 -> 2.65 |
| Vegetation uptake | Toolik 4.1 | 5.42 | 4.14 | 4.06 |
| Net Hg(0) | Toolik 6.5 +- 0.7 | 2.52 -> 2.83 | 0.97 -> 1.06 | 6.50 -> 6.50 |
| Winter share | Toolik 0.37 | - | - | 0.370 -> 0.370 |
| %MeHg (JJAS active) | 0.7 / 0.1 / fen 2 | 0.709 -> 0.708 | 0.100 -> 0.0997 | fen 2.03 -> 2.02 |
| Evasion | - | 6.24 -> 5.90 | 6.39 -> 6.25 | 1.71 -> 1.68 |
| Hg below 1 m: accumulation | - | 3.32 -> 3.48 | 1.60 -> 1.81 | - |
| Hg below 1 m: stock (mg/m2; 5b 37.4 / 14.8) | - | 51.5 -> 52.2 | 19.3 -> 19.7 | - |

10 of 12 bands are met. HPC Ro is 0.8 % below its band (0.2878) and HPC Rm 1.0 % above (1.706). The gn8 prediction was
HPC Ro 0.290 and HPC Rm 1.705 (error about 1 %). 12-band cost 0.018 (5c 0.008, 5b 6.86). With ev_z 0.3 the HPC
organic / mineral pair cannot both be met in this structure (section 8, open item 2). Keeping ev_z 0.5 (5c, V8a) is the only
tested way to meet all 12.

Deep storage: Hg below 1 m (all layers permanently frozen: max ALD 0.78 m at TVC, 0.99 m at HPC) grows by 3.5 ug/m2/yr at TVC
and 1.8 at HPC over the last 20 yr (5b: 0). The 5000-yr spin-up leaves +14.8 (TVC) and +4.9 (HPC) mg/m2 below 1 m compared
with 5b.

### Checks (binT2p5f)

| Check | Result |
|---|---|
| hg_t2 = 1 preset only (P5F3) vs explicit values with hg_t2_preset 0 (X5d, Toolik column) | 0 differing variables (Toolik sp + tr) |
| Routing off + 5b values (R5d, Toolik) vs P5G (5b, binT2p5c) | 0 differences |
| hg_t2 = 0 (T5d, Toolik) vs T1 final VF4 | only the restart hg_state length (117 vs 94) |
| T5d vs P5T1 (5a) | 0 differences |
| P5I10 Hg off 2-case (sub_off) | RESULT: PASS |
| P5I11 all on via the preset (new case set cs/sub_t2j) | RESULT: PASS; C/N/water/energy identical with the tracer on: True |
| P5I12 Hg off 8-case prod | RESULT: PASS |
| Balance error | <= 2e-9 ug/m2 |
| Unit tests | 538/538 (8 known may-fail); Hg standalone 56/56 (2592 assertions) |
| Rules | <= 4 processes; -j2 build; segments <= 75 yr; nohup; no parameters/ edits; no deletions; no GCS; restarts copied only |

The 0.3 core-h budget did not cover 3-site runs for the preset-vs-explicit, routing-off and hg_t2 = 0 checks. They ran on
the Toolik column only, which has every T2 switch active, and the unit tests cover the parse. The 3-site preset run, the
identity suite and all-on C/N/water/energy are full.

### Cost (core-h, LEDGER.txt section "PHASE 5d")

| Item | Core-h |
|---|---|
| P5F3 (3 sites + fen) | 0.044 |
| X5d, R5d, T5d (Toolik only) | 0.031 |
| P5I10 + P5I11 + P5I12 | 0.204 |
| **Phase 5d total** | **0.278 of 0.3** |

Ledger total 8.145.

### Open items after 5d
1. HPC Ro (0.288) and Rm (1.706) are just outside their bands with ev_z 0.3. Options: accept (about 1 %), return to
   ev_z 0.5 (5c preset met all 12), or make a structural change that separates HPC organic from mineral Hg:C.
2. The deep (permafrost) Hg store grows by 1.8-3.5 ug/m2/yr. It needs checking against deep profiles, and it would be
   released on thaw.
3. Unchanged: D5 (hg_t2 default 0), D6, a single f_conn for all sites.
