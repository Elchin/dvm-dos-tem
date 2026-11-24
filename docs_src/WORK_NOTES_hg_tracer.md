# Hg tracer (T1) + MeHg index (T3-lite): work notes

Branch `feat/hg-tracer` (local only), from `4f69f63b` (merge/upstream-develop-20261121). Work was done on calib-elchin
in worktree `~/tem-dev/hg_t1/wt`. Run outputs and scripts are in `~/tem-dev/hg_t1/{out,sites,an,src}`.
Date: 2026-11-23 (EST).

## 1. What it does

`model_settings.hg_tracer = 1` turns on a total-Hg (THg) tracer that follows carbon, in the style of Schaefer et al. 2020
(SiBCASA-Hg). `hg_mehg = 1` also turns on a diagnostic MeHg index.

- **Diagnostic only.** The tracer reads the C pools and fluxes and never writes them. C, N, water and energy
  are bit-identical whether it is on or off (verified, §3).
- **Default off.** Nothing runs and nothing is written unless it is switched on.

**Pools (µg Hg m⁻²), parallel to C:**
- live vegetation per PFT × {leaf, stem, root}
- standing dead per PFT
- woody debris
- one pool per soil layer

**One generic step (`hgtr::step`)** runs around every `Cohort::updateMonthly` sub-step that moves C:
env (daily CH4), DIMveg+DIMgrd, BGC, fire, RTS, end of month (after `soil_updatePools`), plus a month-start
gap step. It takes C snapshots before and after the sub-step and the sub-step's external C fluxes:

1. **Inputs.** NPP > 0 of each PFT part carries `r_leaf / r_stem / r_root` (Hg(0) uptake via plants). Optional
   `hg_dep_soil` (µg m⁻² yr⁻¹) goes to the top soil layer.
2. **Decomposition.** Releases C × (pool Hg:C). Of this, `f_evade` evades, `f_export` is exported and the rest
   re-binds in place. The decomposition C is the per-layer Σ RH (rawc, soma, sompr, somcr) plus CH4 production
   when the CH4 module is on, and woody-debris RH. `Soil_Bgc::deltac` confirms that the net layer loss is
   RH + CH4; the SOM-to-SOM products stay in the layer.
3. **Internal transfers** (litterfall, dead fall, burn residues, slump):
   - Donors put Hg into one pot at their Hg:C.
   - Aggregate sinks take their C share of the pot: fire emission (`fire_v2a + fire_soi2a`), RTS lateral
     (`soi2lat + vegc`) and RTS respiration.
   - Receivers get the rest in proportion to their C gain.
4. **Geometry change** (split, combine, burn, slump). The soil column is remapped bottom-up on cumulative C. A net
   loss is taken from the top; a net gain is added on top.

The step is exactly conservative: ΔHg = uptake + deposition − evasion − export − fire. Any C that cannot be matched
is reported as "unexplained C". In the 75-year TVC smoke test it totalled 0.06 g C m⁻², all in the BGC step.

**Initialization** (first month with no `hg_state` in the restart): Hg = C × `r_org` for moss, fibric and humic
layers and the vegetation/dead pools, and C × `r_min` for mineral layers.

**Restart.** An optional variable `hg_state(Y,X,hgstate)`, layout_version 1, is written to restart files only when the
tracer is on. It is read only when the tracer is on and present. Otherwise the tracer re-initializes from C. Under
MPI it is not written; a warning is printed.

### MeHg index (T3-lite)

Per soil layer, each month:

- Production: P = k_meth · Hg · A · q10_meth^((T−t_ref)/10)
- Demethylation: kd = k_demeth · q10_demeth^((T−t_ref)/10)
- Update: MeHg ← Meq + (MeHg − Meq)·e^(−kd), with Meq = P/kd, capped at Hg.

The driver A is the anaerobic C turnover fraction of the month:

- **CH4 driver** (used automatically where the CH4 module is on, `cd.ch4_on`): A = CH4 production / layer C.
  Methanogens are among the key putative methylators (Thompson 2025), so this is the mechanistically preferred driver.
- **Proxy** (CH4 module off, or `hg_mehg_driver = 1`): A = saturated fraction (from the water table) × RH / layer C.

Each driver has its own `k_meth`. MeHg does not feed back on THg or C.

## 2. Files touched (279 lines changed in 9 files, plus 3 new files)

| File | Change |
|---|---|
| `include/HgTracer.h` (new, 415 lines) | Params, parse, Snap/Flux/State, `initialize`, `step`, `mehg_month`, `pack`/`unpack` (namespace `hgtr`, header-only) |
| `include/ModelData.h`, `src/ModelData.cpp` | `hgtr::Params hg`; parse after the mossdyn block; FATAL on a bad or unknown `hg_*` key |
| `include/RestartData.h`, `src/RestartData.cpp` | optional `hg_state` read/write (serial NetCDF), reset in `reinitValue` |
| `include/Cohort.h`, `src/Cohort.cpp` | `hg_snap`, `hg_track`, `hg_month_begin`, `hg_mehg_month`; hooks in `updateMonthly`, all behind `if (hg_on)`; restart pack/unpack |
| `include/OutputHolder.h`, `src/Runner.cpp` | output maps; the block after DEEPC runs only when `md.hg.tracer == 1` |
| `config/output_spec.csv` | 10 new rows, all off: HGSOIL, MEHGSOIL (layered), HGVEG, HGDEAD, HG_UPTAKE, HG_EVASION, HG_EXPORT, HG_FIRE, MEHGPROD, HG_BALERR |
| `tests/unit/test_hg_tracer.cpp` (new) | 8 test cases (§3) |
| `docs_src/WORK_NOTES_hg_tracer.md`, `docs_src/figures_hg_t1/*.png` | this note |

`parameters/` is not touched. All new parameters are optional flat `model_settings` keys with defaults, so existing
config and parameter files load unchanged. A misspelled `hg_*` key is a FATAL error.

### Parameters (defaults)

| key | default | source / note |
|---|---|---|
| `hg_tracer` | 0 | switch |
| `hg_r_org` | 0.35 µg Hg/g C | T0: Toolik-region organic-horizon median (Obrist / Olson / Schuster; 0.29–0.41) |
| `hg_r_min` | 1.47 µg Hg/g C | T0: mineral median of the synthesis (0.97–1.69); Schuster 2018 overall median is 1.6 |
| `hg_r_leaf`, `hg_r_root`, `hg_r_stem` | 0.07, 0.07, 0.02 | tundra vegetation ~20–60 ng/g dw (Obrist 2017, Olson 2018) at 0.45 g C/g; wood ~5–10 ng/g |
| `hg_f_evade` | 0.02 (T1b; was 0.16) | see §8.1: Obrist 2017 net Hg(0) sink at Toolik; 0.16 = Schaefer 2020 (GTMM) |
| `hg_f_export` | 0.015 | Schaefer 2020 |
| `hg_dep_soil` | 0 µg m⁻² yr⁻¹ | direct Hg(II) deposition is off by default (Schaefer's steady state is ~7.9 µg m⁻² yr⁻¹ total) |
| `hg_mehg`, `hg_mehg_driver` | 0, 0 | switch; 0 = automatic (CH4 if on), 1 = always use the proxy |
| `hg_k_meth_ch4` | **3200** | calibrated: Toolik CMT06 with CH4 on → 2.0 % (fen), Thompson 2025 |
| `hg_k_meth_proxy` | **77** | calibrated: HPC S21H → 0.1 % (permafrost plateau), Thompson 2025 |
| `hg_k_demeth` | 3.0 month⁻¹ at t_ref | ~0.1 d⁻¹, from peat incubations; order-of-magnitude placeholder |
| `hg_q10_meth`, `hg_q10_demeth`, `hg_t_ref` | 1.0, 2.0, 10 °C | the C turnover already carries its own Q10 |

## 3. Tests, identity and conservation

- **Build:** `-j2`, candidate binary `bins/binH2` (md5 1693a49d…). The reference A is `upmerge/bins/binM`
  (md5 e9fe5eff…, the 4f69f63b build), used read-only.
- **Unit tests:** 490/490 test cases pass. 3,173,272 assertions, of which **8 fail. These are the same 8 known ones**
  (test_climate:47/52, test_richards_water_energy:198/215, test_seb:177/178, test_stefan_water:142/271).
  The 8 new Hg cases are:
  - parse: defaults, typo guard, ranges, mehg needs tracer
  - initialize
  - fixed-geometry BGC step (evasion/export fractions)
  - split/combine remap
  - fire
  - 3000 random steps with conservation to roundoff
  - restart pack/unpack round trip
  - MeHg steady state with Q10 and cap
- **Switch off, bit-identity vs 4f69f63b:**
  - J1, prod case set (8 cases, spin20 + n75): A==B True for 8/8; harness PASS.
  - J2, cs_full n75 (8 cases): see §3a.
  - The same checks on the first build (binH1, which differs only in the two k_meth defaults) passed as well:
    I1 1088 files / 4136 vars with 0 diffs; I2 1168 files / 4304 vars with 0 diffs.
- **Switch on** (J3, prod with `hg_tracer = 1, hg_mehg = 1`, A = binM vs B = binH2):
  - every non-Hg variable is identical in every case: 249 (spin) and 268 (main) vars per case, **0 diffs**;
  - the only extras are the 10 Hg output files and the restart variable `hg_state`;
  - conservation: |HG_BALERR| ≤ 4.8e-11 µg m⁻² (against stocks of 3.6–8.2e4 µg m⁻²), and the annual budget
    residual is ≤ 3.1e-11. That is roundoff.

### 3a. J2 (cs_full n75, off)
A==B is bit-identical (True) for all 8 cases (TVC/HPC × warming/no-warming × nfeed/nonfeed).

The harness result reads "FAIL" for one reason only: its closure gate reports "N closure outputs missing", because the
cs_full output spec has no N outputs. The upmerge check failed the same way on the same spec, so this is not an
identity failure.

## 4. Site results

Runs:
- TVC S22 and HPC S21H: spin20 + 1950–2024, each case's own restart.
- Toolik CMT06 (R0_X3_fl03 setup): tr 1901–1975 from a copy of its restart-sp. One run with the CH4 module off
  (proxy driver) and one with it on (wetland, CH4 driver).

Hg is initialized from C × R at the first month, so these are 75-year transient changes from that state, not an Hg
spin-up.

| | TVC | HPC | Toolik (proxy) | obs / T0 |
|---|---|---|---|---|
| total soil Hg first → last (mg m⁻²) | 82.0 → 82.2 | 35.8 → 35.7 | 161.9 → 161.5 | — |
| active layer Hg (mg m⁻²), ALD (m) | 16.4 → 21.8 (0.53 → 0.60) | 14.0 → 17.3 (0.64 → 0.84) | 28.8 (0.61) | — |
| Toolik 0–0.49 m | | | 22.2 | Olson 2018 pits 41.9 (16.6–136.5) |
| Toolik 0.32–1.47 m Hg (mg m⁻²) / SOC (kg m⁻²) | | | 54.3 / 36.5 | Schuster core 73.1 / 56.9 → model 0.74× Hg, 0.64× SOC |
| organic Hg:C first → last (µg/g C) | 0.388 → 0.463 | 0.362 → 0.367 | 0.350 → 0.355 | Toolik organic 0.29–0.41 |
| mineral Hg:C last | 1.49 | 1.47 | 1.49 | Schuster Toolik core ~1.29 (Hg/SOC) |
| uptake / evasion / export (µg m⁻² yr⁻¹, last 20 yr) | 21.9 / 22.9 / 2.1 | 12.8 / 12.5 / 1.2 | 3.5 / 7.6 / 0.7 | — |

- **Toolik.** The Hg stocks are within the observed range but low. The low bias follows the SOC bias: the Hg:C ratio
  of the Toolik deep core is reproduced to within 15%. As in T0, the gap is the model's SOC, not R. The model's
  organic/mineral transition (~0.3 m) is shallower than in the Toolik core (~0.5 m), Fig. 1.
- **Retention.** Re-binding concentrates Hg in thawing organic layers: TVC organic Hg:C rises 0.39 → 0.46 in 75
  years, and at TVC the Hg:C peaks at ~1.8 just above the ALD.
- **Fluxes.**
  - With f_evade = 0.16 and no deposition, the tundra is close to Hg-neutral at TVC and HPC.
  - At Toolik it is a net source: evasion 7.6 against uptake 3.5 µg m⁻² yr⁻¹. Toolik NPP in CMT06 is low, so
    uptake is small.
  - Obrist et al. 2017 measured Toolik tundra as a net Hg(0) sink, so evasion and uptake here are not yet
    constrained (open decision).
- **Fire.** No fire occurred in these runs (HG_FIRE = 0). The fire path is covered by the unit test only.

### MeHg (Fig. 2)

%MeHg of THg in the active layer, June–September mean, last 20 years. %MeHg is linear in k_meth (Meq ∝ k_meth), so
the calibration is one scaling step. It was confirmed by re-running with the new defaults:

| site / driver | k placeholder (30 / 0.6) | calibrated (3200 / 77), JJAS | annual mean | Thompson 2025 class |
|---|---|---|---|---|
| HPC, proxy | 0.0008 % | **0.101 %** (calibration point) | 0.31 % | plateau 0.1 % ✓ |
| TVC, proxy | 0.0064 % | 0.82 % | 1.07 % | upland tundra (no direct class; between bog and fen) |
| Toolik CMT06, proxy | 0.0009 % | 0.11 % | 0.22 % | wet sedge, closest to fen 2 % ✗ |
| Toolik CMT06, CH4 module | 0.0189 % | **2.01 %** (calibration point) | 1.07 % | fen 2 % ✓ |

- **The proxy does not discriminate well.** With the CH4 module off, the modelled Toolik water table leaves little
  saturated, actively decomposing C, so Toolik looks like a plateau. TVC looks wetter than the Toolik wet sedge.
- **CH4 vs proxy driver.** The CH4 driver gives the fen class at Toolik. We recommend it wherever wetland/CH4 is on,
  and treating the proxy as a placeholder (open decision).
- **Gross MeHg production is large relative to the pool** (MEHGPROD 0.9–12.6 mg m⁻² yr⁻¹ with kd = 3 month⁻¹).
  This is fast internal recycling within THg with no mass effect. Only %MeHg is the calibrated quantity. kd itself
  is not constrained.

## 5. Figures (`docs_src/figures_hg_t1/`)

1. `hg_t1_fig1_profiles.png`: Hg density and Hg:C profiles, first vs last year, at Toolik (with the Schuster 2018
   Toolik Lake core samples), TVC and HPC.
2. `hg_t1_fig2_mehg.png`: %MeHg by site and driver (placeholder vs calibrated k) against the Thompson 2025
   classes, and the JJAS %MeHg profiles.
3. `hg_t1_fig3_budget.png`: annual uptake, evasion, export and fire; |cumulative balance error| (all < 1e-9 µg m⁻²);
   total Hg.

## 6. Limitations and open decisions

1. **Default R values.** Horizon-specific 0.35 / 1.47 (T0), or Schuster's 1.6 everywhere? Either can be set in config.
2. **Vegetation Hg:C** (0.07 / 0.02 / 0.07) is a literature placeholder; HGVEG is 6–106 µg m⁻².
3. **No Hg spin-up.** Hg is initialized from C × R; Schaefer used 5000 years. Should a long Hg-only spin-up be added?
4. **Deposition** is 0 by default. A scalar is available; a site time series (`hg_dep.nc`) is not implemented.
5. **The evasion fraction vs Obrist 2017's net Hg(0) sink at Toolik.** The current setup makes Toolik a net source.
6. **MeHg.** Use the CH4 driver only, or keep the proxy? Is a single k per driver enough? Also kd/Q10 values, and
   whether to report JJAS or annual values.
7. **Single pot per sub-step.** Inter-layer SOM transfers (s2d/d2m) are mixed with litter in one pot, which slightly
   smooths Hg:C between layers. A two-pot routing would remove this.
8. `doclost` is always 0 in the model, so there is no DOC Hg path beyond f_export.
9. **The MPI build does not write `hg_state`** (it warns); serial is fine.
10. **Should a parameter file be added?** At present all values are config keys with defaults.

## 7. Core-hours

The ledger `~/tem-dev/hg_t1/LEDGER.txt` totals **0.974 core-h** (budget 1.5). All runs used ≤ 4 cores.

| runs | core-h |
|---|---|
| I1 + I2 + I3 (binH1: prod off, cs_full off, prod on) | 0.443 |
| J1 + J2 + J3 (the same set re-run with the final binH2) | 0.447 |
| TVC/HPC/Toolik site runs, twice (placeholder k and calibrated k), plus the smoke run | ~0.08 |

Builds were `-j2` and are not counted.

## 8. T1b (round 2): evasion revision, MeHg driver, spin-up, deposition series

Elchin's answers (2026-11-23): 1 keep separate organic and mineral R; 2 yes, revise evasion so Toolik is
not a net source; 3 keep both MeHg drivers and check the literature; 4 add spin-up and a "decomposition time
series". **Assumption:** item 4 was read as a *deposition* time series (Hg input history). Decomposition is
already driven month by month by the C model, so a separate series for it would not add anything. If Elchin
meant something else, it is a small follow-up.

Commit `a30a9264` (code, on `c832cc6b`) plus a docs commit. Binary `bins/binH4`, md5 725a95cd, built from
`a30a9264`. Every new feature is off by default. With them off the output is bit-identical to `binM`, and
with the tracer on the C output is identical.

### 8.1 New and changed parameters

| key | default | source / note |
|---|---|---|
| `hg_f_evade` | **0.02** (was 0.16) | Obrist 2017 (Nature 547, 201): Toolik is a net Hg(0) sink, 6.5±0.7 µg m⁻² yr⁻¹, with a soil Hg(0) sink in winter. Olson 2018 found almost no Hg re-emission at Toolik. Lim 2020 (BG 17, 3083) attributes net-sink tundra (Toolik) vs net-source peatland (Degerö) to climate. 0.16 (Schaefer 2020 / GTMM) makes Toolik a net Hg(0) source (−1.5). |
| `hg_f_export` | 0.015 (unchanged) | Schaefer 2020. Model export of 1.1–4.9 µg m⁻² yr⁻¹ is inside or near the pan-Arctic river THg yields of 0.8–4.0 g km⁻² yr⁻¹ (Zolkos et al. 2020, six great rivers; Sonke 2018 gives 44 Mg yr⁻¹). V5 (0.05) gives 3–10. |
| `hg_f_leach` | 0 | optional: fraction of decay-released Hg moved one layer down (exported from the bottom layer). Sensitivity only. |
| `hg_dep_file` | "" | optional CSV: `year[,month],{dep \| hg2_wet,hg2_dry,hg0_soil}[,veg_scale]`, µg m⁻² yr⁻¹. Annual rows are interpolated linearly and held constant beyond the ends. Monthly rows must match exactly. Without a file the constant `hg_dep_soil` is used. `veg_scale` multiplies the NPP-based vegetation Hg(0) uptake. |
| `hg_tr_year0`, `hg_sc_year0`, `hg_dep_spin_year` | 1901, 2016, 1850 | calendar year used to look up the file: tr = tr_year0 + yr, sc = sc_year0 + yr, pre-run/eq/sp = spin_year. |
| `hg_spinup`, `hg_spinup_years`, `hg_spinup_cycle` | 0, 5000, 1 | spin-up switch; years of equilibration (Schaefer 2020 used 5000); recorded cycle length in years. |

New outputs (off in `output_spec.csv`): `HG_DEP` (deposition to soil) and `HG_SPINADJ` (cumulative Hg added
by the spin-up jump). `HG_UPTAKE` is now vegetation uptake only.

**`config/hg_deposition_toolik_obrist2017.csv`.** The 2015 values are Obrist 2017: Hg(II) wet 0.2, Hg(II) dry
2.5 and soil Hg(0) 2.4. The soil Hg(0) is the 37 % snow-season part of the net Hg(0) deposition. The
growing-season 63 % is the model's vegetation uptake. A single factor F scales all columns for the history:

- 0.36 in 1850, rising to 1.15 in the 1970s–1990s;
- 1.00 in 2015 and 0.95 in 2025.

The rise is the ~3× increase since pre-industrial times in northern-Alaska lake cores (Fitzgerald 2005,
ES&T 39, 557). The post-1990 decline is ~−0.6 % yr⁻¹, a damped version of the −30 % 1990–2010 fall in
anthropogenic Hg(0) emissions (Zhang 2016, PNAS). The curve is stylised, not a model output. All sources are
open; no logins were used. **TVC and HPC use the same Toolik file, which is an assumption.**

### 8.2 Spin-up method

The tracer step is affine in the Hg state: s' = M s + b. M and b depend only on the recorded C snapshots
and fluxes. The spin-up works as follows:

1. The C × R initialization is done.
2. `hg_spinup_cycle` full years (January to January) of sub-step records are captured.
3. M and b are built by replaying the cycle on the zero vector and on each of the 63 unit vectors (41 pools
   plus 22 layers).
4. The map is raised to the power n = years / cycle by repeated squaring and applied once.

The jump is exact. It equals n explicit repetitions of the cycle (unit test: 37 years, to 1e-9 relative). The
difference is booked in `cum.spin`, so the balance still closes, and the state is marked `spun` in the restart.

- The jump runs once per cohort, in whichever stage first has the tracer on. It is not repeated in later
  stages.
- If the soil layer count changes across the recorded year, it retries up to 10 times and then gives up
  (`spun = 2`, logged).
- Restart layout v2 adds `cum.spin` and `spun`; v1 restarts are still read.
- Deposition and unassigned Hg go to the top soil layer that holds SOM C. The moss layer has no SOM C;
  depositing there made Hg pile up in it during the jump.
- With `HG_DIAG=1` in the environment, a replay check is printed (outputs ≈ inputs after the jump).

### 8.3 MeHg driver: literature check (both drivers kept)

- **Yang 2016** (Environ Pollut 214, 504; Barrow anoxic incubations): MeHg production rose >10× with warming
  from −2 to 8 °C. It correlated with **both** CH4 production and Fe(III) reduction. Organic soil made ~10×
  more than mineral soil.
- **Tarbier 2021** (ES&T 55, 6710; Stordalen thaw gradient): MeHg is 10× higher and %MeHg 13× higher in
  collapse fens than on plateaus (0.48 → 6.1 % in the top 50 cm). hgcA is mostly in **methanogens** and
  Firmicutes.
- **Fahnestock 2019** (GPL): potential methylators rise from 0.5 % (palsa) to 7.8 % (fen). They include
  Deltaproteobacteria and Clostridia (SRB/FeRB) and Methanomicrobia.
- **Thompson 2025:** methanogens are among the key putative methylators across permafrost thaw.

**Recommendation.** Keep `hg_mehg_driver = 0` as the default. It uses CH4 production where the CH4 module is
on, and the saturated-Rh proxy otherwise. The CH4 driver follows the strongest single control in the
literature, the wet and fen methanogenic niche, and it gives the observed fen/plateau contrast. The proxy
stands for anaerobic respiration in general (SRB/FeRB, which matter in drier plateau soils). It is the right
fallback but has a weak site signal because the model water table is deep. A blended driver (a·CH4 + b·anaerobic Rh) is
the natural next step once hgcA-resolved data are used to fit a and b.

### 8.4 Results (binH4; spin-up on; Toolik dep file; last 20 years; µg m⁻² yr⁻¹)

| site | veg uptake | dep→soil | evasion | export | **net Hg(0)** | net total | soil Hg mg m⁻² (start→end) | R org / R min | %MeHg JJAS active |
|---|---|---|---|---|---|---|---|---|---|
| Toolik CMT06 (proxy) | 3.50 | 5.11 | 1.44 | 1.08 | **+4.5 sink** (Obrist 6.5±0.7) | +6.1 | 167.0→167.3 | 0.93 / 1.48 | 0.17 |
| Toolik CMT06 (CH4) | 3.60 | 5.11 | 1.51 | 1.13 | **+4.5 sink** | +6.1 | 166.3→166.6 | 0.93 / 1.48 | 2.13 |
| TVC S22 | 21.9 | 5.10 | 5.39 | 4.04 | **+18.9 sink** | +17.5 | 72.2→73.5 | 0.65 / 1.27 | 0.56 |
| HPC S21H | 12.8 | 5.10 | 6.55 | 4.92 | **+8.6 sink** | +6.4 | 39.7→40.2 | 0.63 / 1.49 | 0.28 |

net Hg(0) = vegetation uptake + Hg(0) share of soil deposition (2.4/5.1) − evasion. Spin-up jump (5000 yr):
TVC 82.2 → 72.2 mg m⁻² (−10.0), HPC +4.0, Toolik +5.1 (proxy) / +4.4 (CH4) mg m⁻² (from `cum.spin` in the sp
restart). Toolik 0–0.49 m 28.4 mg m⁻² (Olson 41.9); 0.32–1.47 m 53.9 (Schuster 73.1).

Sensitivity (same setup; Toolik / TVC / HPC net Hg(0)): f_ev 0.16 → −1.5 (**source**) / +11.9 / −0.1;
f_ev 0.06 → +2.2 / +15.8 / +4.2; f_ev 0.02 + f_export 0.05 → +4.7 / +21.3 / +11.3 (R org 0.75 / 0.36 / 0.41,
but export 3.1 / 7.5 / 9.6, above river yields at TVC/HPC); f_leach 0.1 → R org does not fall (Toolik 0.98).
Toolik organic R: 0.46 (f_ev 0.16), 0.73 (0.06), 0.93 (0.02) vs obs 0.29–0.41, so evasion and organic R
trade off; see open decision 1. Figs 4–6 (`hg_t1b_fig4_profiles.png`, `_fig5_fluxes.png`, `_fig6_dep_spin.png`).
VF (binH4, final) reproduces V1b (binH3c) to all printed digits.

### 8.5 Checks

- Unit tests (binH4 build): 495/495 test cases; 8 assertion failures, the same known ones in
  test_climate/richards/seb/stefan (unrelated). New cases: dep-file parse/interp/errors, f_leach + veg_scale,
  spin-up jump = explicit 37-yr repetition + 200 000-yr convergence + geometry refusal, restart v2 round trip
  and v1 read, deposition into the top C layer.
- K1 prod, all Hg off, binM vs binH4: **PASS** (bit-identical, 8 cases).
- K2 cs_full n75, off: bit-identical. RESULT FAIL only on "N closure outputs missing" (same as J2 in T1;
  the case set does not write N closure outputs).
- K3 prod_hg (tracer + MeHg on): PASS; `cmp_on.py`: C/N/water/energy identical to binM: **True**.
- K4 prod_hg_full (tracer + MeHg + spin-up + dep file): PASS; `cmp_on.py`: **True**.
- Hg balance residual ≤ 2e-10 µg m⁻² in every site run, including the spin-up jump (Fig 6).

### 8.6 Core-hours

T1b runs (ledger `LEDGER2.txt`): **0.852 core-h** (sensitivity V1–V5 0.229, identity K1–K4 0.585, final VF 0.038). T1 total was 0.974.

### 8.7 Open decisions

1. Organic-horizon R overshoots with low evasion (0.6–0.9 vs obs 0.29–0.41). Options:
   - accept it;
   - raise `hg_f_export` (V5 lowers R but pushes export above the river yields);
   - add an explicit downward DOM-Hg transport.
2. `hg_k_meth_*` were calibrated in T1 without spin-up or deposition. With them, %MeHg shifts (see 8.4), so
   they may need recalibrating for the spin-up configuration.
3. Whether to turn spin-up on by default in production runs. It is cheap: one map build per cohort.
4. Gridded deposition (GEOS-Chem) and site-specific series for TVC and HPC.
5. Model vegetation uptake at Toolik (3.5) is low vs Obrist (~4.1 growing-season uptake; 29 µg m⁻² aboveground
   stock).
6. Leaching is applied only in same-geometry sub-steps.
