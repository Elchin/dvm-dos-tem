# Work notes: `moss_dyn`, a dynamic live-moss layer (default off)

Branch `moss/dyn-moss` (local only, not pushed), branched from 7b6a629a (`wetland/rts-merge-rmf`).
Worktree ~/tem-dev/moss-wt on calib-elchin. Scratch test area ~/tem-dev/moss-work.
Requested by Elchin Jafarov on 2026-11-18: a switch that alternates between the static 1 cm moss and dynamic moss modules,
with the dynamic one coupled to carbon. Background review: ~/tem-dev/notes/moss_regrowth_review.md.

## 1. What the legacy code does (7b6a629a)
- The live moss layer is one layer, hard-coded at 0.01 m (`Ground.cpp` initSnowSoilLayers, "This is wrong - FIX",
  and `redivideMossLayers`). It holds no C. `cmt_dimground` maxdmossthick is read but never used.
- Fire and RTS remove it (`adjustSoilAfterburn`: zero-C layers are dropped), and the same call re-creates a 0.01 m layer.
- The moss PFTs (`cmt_dimvegetation` nonvasculartype > 0) already carry live C in `bd[ip].m_vegs.c[]`, grow with their own
  NPP, drop litter into the first non-moss (fibric) layer (`Soil_Bgc`), and are burned (WildFire) or removed (Thermokarst,
  `thermokarst_veg_remove_frac`) with the existing C bookkeeping. Nothing links that C to the moss layer thickness.

## 2. The switch
`model_settings` keys (all optional; absent = legacy). Parsed in `include/MossDyn.h` (`mossdyn::parse`), FATAL on bad values.

| key | default | meaning |
|---|---|---|
| `moss_dyn` | 0 | 0 legacy static 0.01 m (bit-identical); 1 prescribed regrowth; 2 moss-PFT-carbon coupled |
| `moss_dmax` | -1 | maximum thickness (m); < 0 = `cmt_dimground` maxdmossthick (CMT05 0.30, CMT70 0.035; see caveat 5.3) |
| `moss_min_thick` | 0.001 | floor (m), in (0, 0.05]. The layer is never removed: `fstmossl` is dereferenced without guards (Soil_Bgc.cpp:508, Soil_Env.cpp:966/1046) |
| `moss_t_half` | 5 | mode 1 half-recovery time (yr), Yi et al. 2009; 0 = always dmax (static prescribed thickness) |
| `moss_bulkden` | 25 | mode 2 live-moss dry bulk density (kg/m3) |
| `moss_cfrac` | 0.45 | mode 2 C fraction of moss dry mass (gC/g) |

- Mode 1: `d = clamp(dmax * t / (t + t_half), dmin, dmax)`, `t = cd.yrsdist` (years since fire/RTS; reset to 0 by both, already in
  the restart file). DOS-TEM form (Yi et al. 2009, doi:10.1029/2008JG000841; Jafarov et al. 2013 Eq. 1, doi:10.1088/1748-9326/8/3/035030).
  No carbon.
- Mode 2: `d = clamp(C_moss / (moss_bulkden * 1000 * moss_cfrac), dmin, dmax)`. `C_moss` = live C (all parts) of the PFTs with
  nonvasculartype > 0 and vegcov > 0, in gC/m2 ground. With the defaults 112.5 gC/m2 = 1 cm. Same idea as ORC-HL-VEG (Druel et al. 2017,
  doi:10.5194/gmd-10-4693-2017: h = biomass / density); moss bulk density range from O'Donnell et al. 2009 (doi:10.1097/SS.0b013e3181c4a7f8).
  A CMT with no nonvascular PFT (e.g. CMT70) keeps the legacy layer, with one stderr line per cohort.
- When: every January in `Cohort::updateMonthly_DIMgrd` (independent of the DSL switch), and right after a fire
  (`updateMonthly_Fir`, after `adjustSoilAfterburn`) or an RTS removal (`updateMonthly_Thermokarst`, before the RTS water
  accounting). Mode 1 uses t = 0 at a removal; mode 2 uses the post-removal moss PFT C.
- How: `Ground::setMossThickness(dz)` changes the dz of the single moss layer, keeps its volumetric liquid/ice and temperature,
  shifts the fronts, re-derives the layer properties and the z of the layers below. Skipped when |change| < 1e-12 m.
  Enthalpy changes are booked by the existing solver-1 windows (`enth_bk_layer` for January, `enth_bk_fire` for fire/RTS).
  The moss water change is summed in `Ground::moss_dz_water_mm` (diagnostic only; the legacy daily copy of fibric water into
  the moss layer, Soil_Env.cpp:1049-1056, overwrites moss liquid anyway).

### Why this coupling (mode 2)
I looked for the cleanest carbon link. The moss PFTs already are a closed C pool with NPP, litterfall to the fibric layer and
disturbance losses. Taking the thickness from them adds **no new C pool and no new C flux**, so the C budget cannot change by
construction. The alternatives were worse: (a) bring back the dead-moss C pool (`dmossc`, removed in 9c3ec516 because dead moss is
not separable from fibric) — a second home for the same C; (b) a separate live-moss C pool fed from moss NPP — double counts the
moss PFT C. The cost: mode 2 only works where a CMT has nonvascular PFTs (CMT05, 06, 07, 13, 31, ... yes; CMT70, 71, 65 no).

### Restart
No new restart fields are needed. The moss dz is saved in `DZsoil` (TYPEsoil 0) and rebuilt from it; mode 1 depends on `yrsdist`
and mode 2 on `vegc`, both already in the restart. Tested in 3.4.

## 3. Tests (all on calib-elchin, binary `dvmdostem_moss_v1`, md5 f9dffa63a0a6b040018f66083cbd38bc)
Test runs: one cell, two-stage RTS driver of the v4 bench (`run_rts_two_stage.py`), short spin-up p10 / e300 / s30, historical
(T 1901-2015, G 1950-2015), scenario 30 yr (2016-2045) with the bench RTS slump in Aug 2021. Settings = bench v4b (thermal_solver 1,
n_safetynet_conserve 2, thermokarst_regrow, rhmoistfrozen_override {"70": 0}). Sites: T = Toolik MAT (CMT05, moss PFTs 5-7),
G = Gydan c55 (CMT70). Parameters from a scratch copy of 7b6a629a parameters/. These are mechanics tests, not calibrated runs.

### 3.1 Unit tests (`tests/unit/test_moss_dyn.cpp`, 5 cases)
`make unittest`: 440/440 test cases pass (5 new, 1542 new assertions). The 8 failing assertions are the known pre-existing ones
(test_climate, test_richards_water_energy, test_seb, test_stefan_water; 2 each), unchanged from 7b6a629a.

### 3.2 Bit identity (moss_dyn absent / 0 vs the 7b6a629a binary, md5 1bb756b8…)
`tests/moss_dyn/identdir.py`: every variable of every NetCDF file (outputs and restarts of all stages; wall-clock `total_runtime` excluded).

| case | vs | result |
|---|---|---|
| T, keys absent | T, 7b6a629a binary | 486/486 variables, 62 files identical |
| T, `moss_dyn` 0 | same | 486/486 identical |
| T, `moss_dyn` 1, dmax 0.01, t_half 0 (equivalence) | same | 486/486 identical |
| G, keys absent | G, 7b6a629a binary | 486/486 identical |
| G, `moss_dyn` 2 (CMT70 fallback) | same | 486/486 identical; stderr "CMT70 has no nonvascular PFT ..." |

FATAL checks: `moss_dyn` 3, `moss_dyn` "1", `moss_min_thick` 0 all exit 1 with a FATAL message.

### 3.3 Carbon balance through the removal and regrowth
Monthly residual r = Δ(VEGC + SOC + DEADC + DWDC) + TKNETC over 1901/1950-2045 (TKNETC = NEE + RTS lateral export + RTS respiration
pulse − RTS regrow input), including the Aug 2021 slump (moss PFT C ~160 → 1.6 gC/m2 at T) and the regrowth.

| run | max abs r (gC/m2/month) | mean r |
|---|---|---|
| T mode 0 | 6.2e-4 | 1e-8 |
| T mode 1 | 7.8e-4 | 1e-6 |
| T mode 2 | 7.3e-4 | -1e-6 |
| G mode 0 | 7.8e-4 | -8e-7 |
| G mode 1 | 9.0e-4 | 1e-5 |

All at float round-off of the outputs. Mode 2 closes as well as legacy because it moves no carbon.

### 3.4 Restart round trip (`tests/moss_dyn/rcarry.py`), modes 1 and 2
For T1, T2 and G1, at the end of the historical run (restart-tr.nc, long undisturbed) and at the end of the scenario
(restart-sc.nc, 24 yr into the post-slump regrowth, yrsdist 25):
1. write: moss dz in memory at the restart time (MOSSDZ, Dec) == restart DZsoil[0];
2. read: a 2-yr run from that restart with moss_dyn 0 (nothing changes the thickness) has exactly that dz in every month;
3. state: a 2-yr run from that restart in the run's own mode gives, in its first January, exactly the mode formula evaluated on the
   restored yrsdist (mode 1) or moss PFT vegc (mode 2).

All 6 × 3 checks pass bit for bit (e.g. T2 restart-sc: 0.014260421413845486 m written and read; January target 0.014296891615125868 m
= formula). A split spin-up test (sp 60 vs sp 30 → restart → sp 30) is not usable as a bitwise check: legacy moss_dyn 0 already
differs in 29 of 63 restart variables (lai, vegc, soil T, toptA, prvltrfcnA, ...). That is pre-existing and not caused by this switch.

### 3.5 What the test runs show (Fig. 1)
![moss_dyn tests](figures/moss_dyn_tests.png)

- Mode 2 at T: 1.2-1.4 cm before the slump (moss PFT C 110-160 gC/m2), floor 0.1 cm for two years after it, back to ~1 cm by
  2027-2028 as the moss PFTs regrow. The thickness is set each January from the moss C of that moment, so it lags the moss C by up to a year. ALD after the slump is deeper than legacy in the first years (thinner moss).
- Mode 1 at T with the CMT05 parameter maxdmossthick (0.30 m) gives a 30 cm moss layer: it suppresses moss PFT C
  (16-28 vs ~140 gC/m2) and after the slump the regrowing moss (5 cm after 1 yr) holds ALD at 0.05-0.2 m. That value is the old
  dead-moss maximum, not a live-moss thickness (caveat 5.3).
- ALD in dvmdostem is measured from the top of the moss layer (Soil_Env.cpp:733-740), so a thicker moss adds its own thickness to ALD.

## 4. Gydan sensitivity: ALD vs moss thickness (8 cases + 1 control, `tests/moss_dyn/plot_sens.py`)
Gydan c55 (CMT70), mode 1, same short test settings as section 3 (not the calibrated A0/A2 set in ~/cal-gydan, so absolute
ALD is not comparable to the 0.78 m CALM target). Two thermal options: legacy ln-Kersten (L) and `org_tc_scheme` 1 with the
S21H values (O: kunf 0.27, kfrz 0.5134, dry 0.08571, thre_f 0.1795, thre_h 0.5638). Moss: static 0.1 / 3 / 10 cm (`moss_t_half` 0),
and regrowth to 5 cm with `moss_t_half` 5 (thinned to 0.1 cm by the Aug 2021 slump). Controls: L = G0 (legacy 1 cm), O = 1 cm.
ALD is from the moss top.

| thermal | moss | ALD 1986-2015 (m) | ALD below moss | ALD 2022-26 | ALD 2036-45 | max abs C resid |
|---|---|---|---|---|---|---|
| L | 1 cm (legacy) | 0.144 | 0.134 | 2.05 | 4.54 | 7.8e-4 |
| L | 0.1 cm | 0.127 | 0.126 | 2.16 | 4.71 | 6.8e-4 |
| L | 3 cm | 0.151 | 0.121 | 1.51 | 3.68 | 7.5e-4 |
| L | 10 cm | 0.104 | 0.004 | 0.32 | 0.34 | 6.6e-4 |
| L | regrow 5 cm | 0.180 | 0.131 | 1.61 | 3.42 | 7.2e-4 |
| O | 1 cm | 0.802 | 0.792 | 3.90 | 5.26 | 9.1e-4 |
| O | 0.1 cm | 0.763 | 0.762 | 3.54 | 5.25 | 6.5e-4 |
| O | 3 cm | 1.064 | 1.034 | 4.28 | 5.28 | 7.1e-4 |
| O | 10 cm | 3.748 | 3.648 | 5.35 | 5.36 | 6.5e-4 |
| O | regrow 5 cm | 1.306 | 1.257 | 4.76 | 5.30 | 7.1e-4 |

![Gydan moss sensitivity](figures/moss_dyn_gydan_sens.png)

What it shows (test settings, so read the signs, not the numbers):
- **The sign of the moss effect depends on the moss thermal scheme.** Legacy: more moss = shallower thaw below the moss
  (0.134 → 0.121 → 0.004 m for 1 / 3 / 10 cm), the insulation Elchin expected; with 10 cm the thaw stays at ~0.3 m even after the slump.
  S21H: more moss = much deeper thaw (0.79 → 1.03 → 3.65 m). In Feb 2000 the 10 cm S21H moss is −2.8 °C at its centre and −0.9 °C
  just below (1 cm case: −2.3 / −2.1 °C), and the 1995-2004 mean soil T is +1.4 °C just below the moss and still +0.4 °C at 0.5 m, versus +0.2 °C and −1.4 °C at the same depths with 1 cm.
  So the S21H moss insulates in winter more than in summer, which warms the soil. The likely reason is the moss hydrology
  shortcut: the moss copies only the fibric layer's *liquid* water each day, so it is wet when thawed (summer Sr ≈ 0.28, above
  thre_f, k well above dry) and has little ice when frozen (k near the dry value). Real moss is the opposite (dry in summer, icy
  in winter: the thermal offset). The legacy ln-Kersten puts summer moss at k = tcdry because Sr < 1/e, so it insulates in summer.
- Removal and regrowth (regrow 5 cm): with legacy thermal, ALD 2022-26 is 1.6 m vs 2.0 m with the static 1 cm moss (the two also
  differ before the slump: 5 vs 1 cm). The slump itself (organic and mineral removal, excess ice) dominates both.
- C balance stays at round-off in every case.

## 5. Caveats and open points
1. Mode 2 needs nonvascular PFTs. CMT70 (Gydan) has none: its "Omoss" PFT is flagged nonvasculartype 0. Using mode 2 at Gydan
   would need a parameter change (flag the PFT as moss), which is a physics/parameter decision for Elchin.
2. Removal is only by fire or RTS. There is no imposed moss-only removal (e.g. a list of years), because in mode 2 that would have
   to move moss PFT C somewhere (litter or emission) and in mode 1 needs a new restart field. Can be added if wanted.
3. `cmt_dimground` maxdmossthick values were set for the removed dead-moss layer (CMT05 0.30 m). Set `moss_dmax` explicitly
   (DOS-TEM used 0.04/0.06 m upland/lowland).
4. Moss thermal properties: with the legacy ln-Kersten at the summer saturations of these sites the moss conductivity is
   ~tcdry 0.025 W/m/K, so each extra cm of moss is a strong summer insulator. With `org_tc_scheme` 1 the sign reverses (section 4).
5. Moss hydrology is unchanged: the moss copies the fibric layer's liquid water daily (no own wetting/drying, little winter ice).
   Before using dynamic moss for ALD work, the moss moisture (or a moss-specific conductivity, e.g. Chadburn 2015 linear form) needs a decision.
6. The moss water change on a thickness change is not routed to runoff; it is summed in `moss_dz_water_mm` only.
