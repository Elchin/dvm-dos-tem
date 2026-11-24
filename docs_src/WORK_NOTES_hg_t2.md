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
