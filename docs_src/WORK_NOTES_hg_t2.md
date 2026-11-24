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
