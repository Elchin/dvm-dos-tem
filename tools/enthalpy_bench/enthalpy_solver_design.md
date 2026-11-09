# dvmdostem: fully implicit, globally convergent enthalpy solver for soil/snow heat (Option C)

*Stage A design document, updated with Elchin's decisions (§0a) and Stage B results (§15a). TEM-dev, 2026-11-09 (EST). No model code was written. The code review was read-only against `wetland_deep_C` at `1dca209e` (box worktree `/workspace/tem-dev/wdc`). Line references are to that commit. The runtime profile in §15 comes from a separate `-pg` build of the same commit, made and run on the box (`/workspace/tem-dev/prof`, not the repo).*

---

## 0. Summary

| Item | Decision |
|---|---|
| Method | **Implicit (backward-Euler) cell-centred finite-volume enthalpy scheme solved with the nested Newton method of Casulli & Zanolli (NCZ)**, following Tubini et al. (2021). The h(T) interface is designed so that a **Katzenelson** solver (Hötten et al. 2026) can be added later as a cross-check for piecewise-affine curves. |
| State | Prognostic **layer enthalpy H_i [J m⁻²]**. It is computed from and written back to the existing `Layer::tem`, `liq` and `ice`, so the rest of the model keeps reading T, liq and ice. |
| Enthalpy reference | Liquid water at 0 °C has H = 0. This is the definition `Soil_Env::layerEnthalpy` already uses for `energy_diag`. Meltwater leaving at 0 °C therefore carries zero enthalpy, which makes the snow-melt hand-off exact. |
| Freezing curve | Per layer, selected by `thermal_sfcc`. **Production: 1 = Schaefer & Jafarov (2016) power law** (decision D1), taken from `Soil_Bgc::sjLayerPhi` (texture and organic-fraction dependent, b = −0.3/−0.5/−0.9/−1.0), with a linear ramp over [−ε, 0] (ε = 0.01 °C, legacy `zerodegc`) to make it continuous at 0 °C, and a floor at `minliq`. **0 = free water** with the same ramp, kept as an option and for verification (Katzenelson cross-check, legacy-like comparisons). Snow and rock always use free water (rock has no water). |
| Conductivity | Cell λ is a geometric blend of the existing frozen and unfrozen layer conductivities, weighted by ice fraction. Faces use the **harmonic (series) mean** (D2); Tubini's max rule is an optional, documented fallback, default off. λ is lagged at the start of the step, with optional Picard updates. |
| BCs | Top: the same Robin conductance the legacy solver and the SEB already use, Λ_top = 1/(1/100 + ½dz₁/λ₁) to T_s (n-factor or SEB skin temperature). Bottom: zero flux at the base of the rock column (≈ 50 m of rock below the soil), with a geothermal flux option: zero (default), one value, or a spatial map (D4). |
| Time step | One daily backward-Euler step by default (`thermal_substeps` = 1). There is no adaptive step cutting, because NCZ converges for any Δt. Tolerance is \|r_i\| < 1e-3 J m⁻² per cell, so column closure is ≤ 4e-7 W m⁻². |
| Fronts and ALD | **Diagnostic only.** They are reconstructed from the enthalpy/ice profile into the existing `frontsz`/`frontstype`, `frozen`, `frozenfrac` and `ald`/`alc` fields, so downstream code reads the same variables. |
| Remap | Layer split, merge, burn, carbon-driven thickness change and snow construction all go through one `EnthalpyRemap` helper. It **conserves H** for rearrangements and books ΔH explicitly for mass that enters or leaves. |
| Switch | `model_settings.thermal_solver`: 0 = legacy (default), 1 = enthalpy. With 0 the legacy code path is untouched, verified bit-identical by md5 in the harness. |
| Restart | No new restart variables. Both modes read and write the same T, LIQ, ICE, FROZEN, FROZENFRAC and front fields. |
| Cost | Measured (§15, §15a): the legacy thermal step is **79 % of CPU** (demo profile) and does about **200 tridiagonal solves per cell-day**. Stage B measurements on calib8 TVC/HPC columns are in §15a; the core-hour budget for the 4×4 km domain is deferred until Elchin asks (D6). |
| Stages | B: standalone library + tests + Katzenelson cross-check (D7), 9–11 d (done 2026-11-09, §15a). C: integration behind the switch, 12–16 d. D: TVC/HPC site evaluation, 6–8 d. E: calibration plan, 2 d, with the recalibration itself (15–25 d) shared with deep C. |

## 0a. Decisions 2026-11-09 (Elchin)

| # | Question | Decision | Sections changed |
|---|---|---|---|
| D1 | Production freezing curve | **`thermal_sfcc = 1` (Schaefer & Jafarov 2016, identical to `sjLayerPhi` below −ε, plus the 0 °C ramp)**. `thermal_sfcc = 0` (free water) stays as an option and as the verification case. | §0, §5, §12 |
| D2 | Interface conductivity | **Series (harmonic) mean.** Tubini's max rule only as an optional, documented fallback (`thermal_face_mode = 1`), default off. | §0, §3, §12 |
| D3 | Snow melt | **Yes:** with `thermal_solver = 1` only, radiation/SEB melt becomes a heat input to the snowpack (top snow cell source), and the "surface fixed at 0 °C while melting" rule is dropped. Legacy path unchanged. | §6 |
| D4 | Geothermal bottom flux | **Config option:** zero (default, legacy), a single value, or a spatial map input file. | §6, §12 |
| D5 | Full-day energy balance | **Yes:** full-day column energy balance ≤ 1e-6 W m⁻² becomes a pass/fail harness test once Stage C integrates. | §10, §14 |
| D6 | Core-hour budget | Deferred until Elchin asks. §15 and §15a give per-10,000-cell numbers only. | §15, §18 |
| D7 | Katzenelson cross-check | **Yes**, in Stage B, for the free-water case. | §2, §17 |

---

## 1. Why the legacy scheme cannot close, and why the front prototype failed

From the code, the Phase 2b/2c attribution (merge plan §8.1–8.2) and the calib8 runs with the corrected diagnostic:

1. **The Stefan front and the conduction solver are not coupled (M1).** `Stefan::updateFronts` (Stefan.cpp:38) moves fronts using surface degree-seconds and the resistance above the front. `TemperatureUpdator` (TemperatureUpdator.cpp:36–130) then solves the regions above, between and below the fronts separately, with each front layer held at ±`zerodegc` (Dirichlet). The latent heat the front uses and the heat the solver delivers to it are never reconciled. On front days the Stefan latent heat exceeds the conducted heat by **+2.0 to +2.9 W m⁻² (HPC) and +4.5 to +5.3 W m⁻² (TVC)**.
2. **Latent energy is added as capacity.** `hca[ind] = pce + hcap` adds a J m⁻² quantity to a J m⁻³ K⁻¹ one (TemperatureUpdator.cpp:163 and similar lines). This changes thermal inertia without being a source term.
3. **The discretisation does not match the enthalpy definition (about −0.7 W m⁻² on no-front days).** `CrankNicholson::geBackward` (CrankNicholson.cpp:14–50) uses `rc = (cap[il] + cap[iu])/dt`, so the node capacity is half of each of two neighbours. The face conductance is `cn[il] = λ_il/dz_il`, which uses one side only, not the series value. `energy_diag` meanwhile measures the full-layer enthalpy. No discrete energy is conserved exactly by this stencil, so the "conserved" quantity and the measured one differ by −0.65 to −0.86 W m⁻².
4. **Snow melt crosses the boundary of the thermal step.** `Stefan::meltingSnowLayer` (Stefan.cpp:165) converts snow ice to `liq` with a degree-second budget inside the thermal step. `Snow_Env::meltSnowLayersAfterT` (Snow_Env.cpp:185) removes that liquid later, and `rad_melt`/SEB melt remove more mass outside the thermal step. The energy taken from the column and the mass removed are booked in different places.
5. **Remaps do not conserve energy.** `Ground::splitOneSoilLayer` interpolates T from the gradient, and `combineTwoSoilLayersU2L/L2U` average T by thickness. Neither is capacity-weighted, and neither includes latent heat (Ground.cpp:1650–1850).

The **energy-limited front prototype** (`proto/front-energy`, adf4be54) addressed item 1 only. It closed most days but failed in multi-front years (annual up to 38 W m⁻²). It made ALD 0.07–0.46 m shallower, cooled soil T by 1–1.7 °C, failed water closure (2–5 % of P) and increased integrator safety-net events from 4 to 76. The lesson for this design is that **front bookkeeping is the problem**. The new solver removes fronts from the physics entirely, uses one stencil whose discrete energy is exactly the diagnostic enthalpy (items 1–3), and defines explicit, conservative hand-offs for melt and remap (items 4–5).

---

## 2. Choice of method: NCZ vs Katzenelson

Both papers were re-read: Tubini et al. (2021) in full, including §3.2 and Appendices B, C and E, and Hötten et al. (2026) in full from the ETNA PDF.

| | NCZ (Casulli & Zanolli 2010, 2012; Tubini et al. 2021) | Katzenelson (Hötten et al. 2026) |
|---|---|---|
| Admissible h(T) | Any h whose C_a = dh/dT is non-negative, of bounded variation, non-decreasing up to some T\* and non-increasing above it (conditions C1/C2). This covers free water with a smoothing ramp and power-law or van Genuchten SFCCs. | **Piecewise affine h(T) and piecewise constant λ only** (their §5–6). A power-law SFCC would have to be tabulated as piecewise linear. |
| Guarantee | Proven convergence for any Δt, given the prescribed initial guess. Exact conservation once converged. | Exact discrete solution in finitely many steps from almost every initial guess, for any Δt. Exponential worst case not ruled out (their §5). |
| Measured cost | Tubini Table G2 reports a few iterations per step at 1 h; Newton-Raphson and globally convergent Newton hit the iteration cap. | **1.48 (BE) / 1.93 (CN) linear solves per daily step** in LPJmL with n = 24, within a factor of 2 of implicit DECP. |
| Discretisation | Finite volume (fits dvmdostem's layer = control volume). | Weak FE form (P1 nodes), which would need re-deriving for layer-centred cells. |
| Conductivity | Can depend on T, either lagged or by Picard (Tubini App. C: conservation does not depend on λ convergence). | Piecewise constant in each phase. |
| Reference code | Pseudocode in Tubini App. B; FreeThaw1D (Java, open source). | Algorithm in the paper; LPJmL implementation, not checked. |

**Decision: NCZ is the primary solver.**

- It supports the texture- and organic-dependent unfrozen-water curve that dvmdostem already uses for frozen-soil respiration (`rh_frozen_mode`). Katzenelson would need a piecewise-linear approximation of every layer's curve, rebuilt whenever water content or organic fraction changes.
- It maps directly onto layer-centred finite volumes with a tridiagonal (Thomas) solve.
- It has a convergence proof and published pseudocode.

The `EnthalpyCurve` interface (§13) exposes piecewise-affine breakpoints when `thermal_sfcc = 0`. Per D7 the Katzenelson solver is part of Stage B (`ColumnSolver::stepKatz`), to cross-verify NCZ: both must produce the same discrete solution to round-off on free-water cases (verified, §15a). To make the free-water system exactly piecewise affine, SFCC 0 defines H as linear in T on the ramp [−ε, 0], and the cell's liquid water there is diagnosed from H, so (T, liq, ice, H) stay mutually consistent (§5).

---

## 3. Governing equations and discretisation

**Cells.** One control volume per `Layer` from `ground->toplayer` (snow, moss, fibric, humic, mineral) down to `botlayer` (rock). That is at most `MAX_GRN_LAY` = 6 + 22 + 5 = 33 cells. Cell centre z_i + ½dz_i, thickness dz_i.

**Discrete equation** (Tubini et al. 2021 Eq. 12, 1-D, per unit area):

    H_i(T_i^{n+1}) − H_i^n = Δt [ K_{i−½}(T_{i−1}^{n+1} − T_i^{n+1}) − K_{i+½}(T_i^{n+1} − T_{i+1}^{n+1}) ] + Δt S_i

with matrix form H(T) + A T = b (Tubini Eq. 16). A is symmetric, tridiagonal and an M-matrix. The nonlinearity is only on the diagonal.

- **Face conductance:** K_{i+½} = 1 / (½dz_i/λ_i + ½dz_{i+1}/λ_{i+1}), in W m⁻² K⁻¹. Tubini uses max(λ_i, λ_P) (their Eq. 13). We use the series value because dvmdostem stacks very different materials: moss and peat at 0.05–0.3 W m⁻¹ K⁻¹ over mineral soil at 1–2.5. A max would overstate flux through the insulating organic mat and bias ALD deeper. (Per D2 the max rule exists only as `thermal_face_mode = 1`, a documented fallback, default off; the Stage B unit test shows it misses the exact two-layer steady flux that the harmonic mean reproduces.) NCZ convergence only needs A to be a symmetric positive-definite M-matrix, which holds for any positive face conductances.
- **Top face:** K_{½} = 1/(1/cn_top + ½dz₁/λ₁) with cn_top = 100 W m⁻² K⁻¹. This is the legacy virtual-layer conductance (TemperatureUpdator.cpp:149) and the same Λ the SEB uses (`G = Lam*(Ts-T1)`, Soil_Env.cpp:395 comment), so SEB and soil stay consistent. Boundary value T_s (§6).
- **Bottom face:** flux = −q_geo (positive upward into the column), default 0.
- **Discrete energy = diagnostic energy.** H_i is defined by the same formula `energy_diag` measures, `Soil_Env::layerEnthalpy` (Soil_Env.cpp:361). The scheme is conservative in exactly the quantity we monitor, which removes the −0.7 W m⁻² no-front bias (item 3 in §1).

**Time integration.** Backward Euler (θ = 1). Crank–Nicolson would impose a monotonicity time-step limit proportional to Δz² (Tubini §5) and was the slower option for Katzenelson. Daily steps with backward Euler gave 1 h vs 1 d differences of ≤ 0.04 °C in mean T over 100 years (Tubini Fig. 12).

---

## 4. State variable and enthalpy closure

Per cell i, with reference state liquid water at 0 °C:

    H_i(T) = C_s,i T + c_w ρ_w θ_l,i(T) dz T + c_i ρ_i θ_ice,i(T) dz T − L ρ_i θ_ice,i(T) dz         [J m⁻²]

This is implemented through the masses (kg m⁻²) that `Layer` already stores. With total water W_i = liq + ice:

    liq_i(T) = W_i · φ_i(T)        ice_i(T) = W_i − liq_i(T)
    H_i(T)   = C_s,i T + liq_i(T) c_w T + ice_i(T) (c_i T − L)

- C_s,i = `vhcsolid·(1−poro)·dz` for soil (SoilLayer.cpp:33 without the water terms), `vhcsolid·dz` for rock and 0 for snow. Water constants come from physicalconst.h: SHCLIQ 4188, SHCICE 2117.27, LHFUS 3.337e5.
- **Prognostic H, diagnostic T.** Between thermal calls the model keeps `tem/liq/ice` on each Layer, so all existing readers still work. At entry, H_i is computed from (tem, liq, ice) with the general formula above, which is valid even if the triple is off the curve. The cell is then **equilibrated**: T is solved from H on the curve by scalar monotone Newton and bisection, and liq/ice are reset. Equilibration conserves H exactly and is the only way the thermal code changes a layer's phase partition.
- W_i, φ_i(·) and C_s,i are frozen for the duration of a step, as Tubini §6 requires for parameter changes. Water redistribution (Richards, infiltration, drainage) happens between thermal calls through the remap rules in §10.

---

## 5. Freezing characteristic (`thermal_sfcc`)

**Option 0: free water with ramp** (verification and legacy-comparison option; not the production default, D1). Implementation (Stage B): H is linear in T on [−ε, 0] between H(−ε) and H(0) = 0, so H is exactly piecewise affine (needed by Katzenelson), and liq on the ramp is solved from H = C_s·T + liq·c_w·T + (W − liq)(c_i·T − L), keeping liq and H consistent. φ(T) = 1 for T ≥ 0, linear from φ_min to 1 over [−ε, 0], and φ_min below −ε. Here φ_min = min(1, minliq/W) keeps the legacy residual liquid (`minliq = 0.05·poro·ρ_w·dz`, SoilLayer.cpp:169) and ε = 0.01 °C (legacy `zerodegc`). This is Tubini Eq. E1–E2 with a residual floor. A cell holding a front is at T ∈ [−ε, 0], with a frozen fraction set by its latent content, which is the fixed-grid equivalent of the legacy front layer.

**Option 1: Schaefer & Jafarov 2016 power law** (**production default, D1**), matching the curve dvmdostem already uses. `Soil_Bgc::sjLayerPhi` blends clay/silt/sand exponents (−0.3/−0.5/−0.9) with organic (−1.0) by the organic fraction from `sjOrgFraction`, with T_ref = 0.1 °C and T\* = 0.01 °C (Soil_Bgc.cpp:1590–1612). As coded, φ jumps at T = 0 from 1 to (0.1/0.01)^b, which is about 0.5 for clay and 0.1 for organic. The thermal curve needs h to be continuous (C1). It therefore uses:

    φ(T) = 1                                    T ≥ 0
    φ(T) = linear ramp from φ_SJ(−ε) to 1         −ε ≤ T < 0
    φ(T) = max(φ_SJ(T), φ_min)                  T < −ε

For T < −ε this is **identical to the BGC curve**, so thermal liquid water and frozen-soil respiration see the same unfrozen water. The ramp holds the part of the latent heat that the coded curve releases as a jump at 0 °C. **C2 check:** C_a is non-decreasing up to 0 °C, because the floor kink and the power-law derivative both increase toward 0 and the ramp slope (dφ/dT ≈ 51 °C⁻¹ for clay at ε = 0.01 °C) is larger than the power-law slope at −ε (≈ 1.3 °C⁻¹). Above 0 it drops to the sensible value. So T\* = 0⁻ for every layer, which is exactly the shape NCZ is built for (Tubini Fig. 1).

- Moss, fibric and humic layers: forg = 1, b = −1.0. Mineral layers: texture from the existing `rhf_fclay/fsilt/fsand` arrays plus forg from SOM. Rock: no water. Snow: always option 0 with φ_min = 0.
- Unit-test the curve against `sjLayerPhi` for T < −ε to 1e-12.
- **Caveat for Elchin:** option 1 changes physics as well as numerics. It produces a zero curtain and liquid water in frozen soil. Expect deeper heat penetration in winter and some ALD change, so it is evaluated separately from option 0 in Stage D.

---

## 6. Boundary conditions and snow

**Top.** T_s is the value `updateDailyGroundT` already computes: `tdrv·nfactor`, or the SEB skin temperature when `seb_scheme > 0` (Soil_Env.cpp:267–289). Two legacy special cases change:

- `meltsnow → T_top = 0` (TemperatureUpdator.cpp:141) is **dropped** under the enthalpy solver (D3; legacy unchanged). Melting is resolved by the snow cells' own enthalpy.
- The `nstep = 24` rule for an abrupt surface-state change (Soil_Env.cpp:293–300) is replaced by `thermal_substeps` (default 1). The legacy rule existed only to keep the front logic stable.

**Bottom.** Prescribed flux at the rock base: `ROCKTHICK` totals 50 m, so about 55 m in all. Per D4, `model_settings.geothermal_flux_mode`: 0 = zero flux (default, the legacy behaviour, `gflx = 0.0`, TemperatureUpdator.cpp:171); 1 = one value `geothermal_flux` [W m⁻²]; 2 = spatial map, read per cell from a NetCDF input `IO.geothermal_flux_file` (variable `geothermal_flux`, W m⁻², same grid as the run mask; FATAL if missing or not finite). The solver takes it as `BC::q_bot`.

**Snow cells.** Each snow `Layer` is a cell with H = ice(c_i T − L) + liq·c_w·T. The free-water curve gives T ≤ 0 automatically until all ice is melted. λ comes from `SnowLayer::getThermCond()` (density-based, unchanged).

**Snow-melt hand-off** (item 4 in §1):

1. During the step, melt appears in a snow cell as `liq > 0` at T ∈ [−ε, 0]. No degree-second budget is used, and `Stefan::meltingSnowLayer` is not called.
2. After the step, if a snow cell has H > 0 (fully melted and above 0 °C), its excess sensible heat H − 0 is moved to the next cell down. The cell is marked `frozen = −1`, which tells `constructSnowLayers` to remove it, as now. The move is conservative and booked as `edg.melt_excess`.
3. `Snow_Env::meltSnowLayersAfterT` then removes `liq` as melt, unchanged. Under the 0 °C liquid reference this removes c_w·T·liq ≈ 0 J (|T| ≤ ε). The small remainder is booked as `edg.melt_out`, so the snow module receives mass and the thermal module keeps all the energy.
4. **SEB and radiation melt.** When the SEB is on over snow, the surplus energy that now becomes `seb_melt`/`rad_melt` mass removal in `Snow_Env` (Snow_Env.cpp:57–66) is applied instead as a **flux source into the top snow cell** in the same implicit step (S_1 = Q_melt). This avoids melting the same snow twice. The mass-removal path is skipped when `thermal_solver = 1` (D3). The solver already accepts a per-cell source S_i (W m⁻²), tested in Stage B.
5. New snowfall (`constructSnowLayers(+dsmass, tsurface)`) enters at T = min(T_s, 0) with H = m(c_i T − L). This is booked as `edg.snowfall_in`. Sublimation removes ice at the layer's T and is booked as `edg.sublim_out`.

---

## 7. Time step and convergence

- Δt = 86400 s / `thermal_substeps` (1, 2, 4, 8 or 24; default 1). The nonlinear system is solved once per substep.
- **NCZ loop** (Tubini Eqs. 21–26, App. B), per substep:
  1. Assemble A and b from λ(T^n) and the BCs.
  2. Initial guess T^0 = min(T^n, T\*) per cell. Tubini §3.2 prescribes an initial guess on the convex side of T\*, which is 0⁻ here.
  3. Outer loop k: linearise h₂ at T^{k−1}.
  4. Inner loop l: linearise h₁ at T^{k,l−1} and solve the tridiagonal system (Thomas) for T^{k,l}. Stop when ‖r^{k,l}‖∞ < ε_H.
  5. Stop the outer loop when ‖r^k‖∞ < ε_H.
- **Tolerance:** ε_H = 1e-3 J m⁻² per cell. Cell enthalpies are 1e6–1e8 J m⁻², so this is 1e-9 to 1e-11 relative, well above the double-precision floor of about 2e-8 J m⁻². With 33 cells the column residual is ≤ 0.033 J m⁻² per day, about **4e-7 W m⁻²**, inside the 1e-6 W m⁻² target.
- **Caps and fallback:** 50 outer and 50 inner iterations. If a cap is ever reached (theory says it should not be), the substep is halved, a counter `thermal_ncz_fallbacks` is written to run.log in the same style as `soil_safetynet_events`, and the harness fails any case with a non-zero count.
- **Conductivity:** λ is lagged at T^n. `thermal_k_picard` = 0 by default; setting n > 0 re-evaluates λ at T^{n+1} and re-solves n times. Per Tubini App. C, conservation holds at every Picard pass, and only accuracy changes. Stage D measures whether one pass changes ALD by more than 1 cm.

---

## 8. Conductivity (frozen/unfrozen mixing)

- **Cell λ** = λ_u^{1−f} · λ_f^{f}, with ice fraction f = ice/(liq + ice) at T^n. λ_u and λ_f are the existing endpoint functions `getUnfThermCond()` and `getFrzThermCond()`, including `org_tc_scheme` and `kersten_form`, evaluated with current liq and ice. This reproduces legacy λ exactly for fully frozen or fully thawed layers. It replaces the legacy binary switch (`frozen == 1 ? frz : unf`, Layer.cpp:92), which puts partly frozen layers at the unfrozen value, with a continuous blend consistent with the geometric-mean saturated conductivity already in `MineralLayer` (`tcsatunf/tcsatfrz`, MineralLayer.cpp:62–63).
- **Snow and rock:** unchanged.
- **Faces:** harmonic, as in §3.

---

## 9. Fronts, ALD and downstream consumers

Fronts no longer drive physics. They are **reconstructed after each day** by `EnthalpyDiagnostics::rebuildFronts()`, which fills the same `Ground::frontsz/frontstype` deques and therefore `ed->d_sois.frontsz/frontstype` through `retrieveDailyFronts`.

- **Per-layer state.** `frozenfrac_i` = ice_i / (W_i − liq_min,i), clipped to [0, 1]. This is the fraction of freezable water that is frozen. For free water with a sharp front it equals the legacy geometric frozen thickness fraction. `frozen_i` = +1 if T_i < −ε (below the ramp), −1 if ice_i = 0, and 0 otherwise (cell in the ramp, i.e. it holds a front). `updateWholeFrozenStatus` and `ststate` keep working.
- **Front position in a ramp cell.** z = z_top + frozenfrac·dz if the frozen part is above, i.e. the colder neighbour is above (front type 1, freezing from the top). z = z_bot − frozenfrac·dz if the colder neighbour is below (type −1, thawing from the top).
- **Fronts between cells.** A 0 °C crossing between two non-ramp neighbours (possible with large Δt, or for SFCC option 1 where frozen cells are at T < −ε) is linearly interpolated between the cell centres, as the existing `ald` code does (Soil_Env.cpp:690–712). Fronts are sorted and merged within `MAX_NUM_FNT` with the existing `combineExtraFronts` rule.
- **ALD.** `ed->d_soid.ald` keeps its current definition (deepest 0 °C crossing, Soil_Env.cpp:681–712) for output continuity. `aldflag` (deepest thawing front) and `alc` are computed from the reconstructed fronts. Both should agree to within one ramp width, and a unit test checks this.

**Downstream consumers:**

| Consumer | Today | Under `thermal_solver = 1` |
|---|---|---|
| `Soil_Bgc::thawFraction`, `thawfrac_decomp`, `thawfrac_moist_mode` (Soil_Bgc.cpp:808–898, 1543) | `m_sois.frozenfrac`, `ts` | Same fields, filled from the enthalpy frozenfrac. Semantics are identical under SFCC 0. Under SFCC 1 a layer at −2 °C has frozenfrac < 1, i.e. some thawed share, which is physically right but changes decomposition. Flag it in calibration. |
| `rh_frozen_mode` (SJ2016 φ) | Own φ(T) from `ts` | Unchanged. With SFCC 1, thermal liquid water and BGC φ come from the same curve. An optional later step could pass thermal liq/W directly. |
| `stefan_front_resist_fix` (Stefan.cpp:4–17) | Front-layer resistance inside the Stefan solve | **Not applicable.** It is ignored with a startup warning, recorded in `config_record.js`. Setting both to 1 is allowed, so production configs keep working. |
| `humic_max_thick` (Organic.cpp:125) | Caps humic thickness and splits into more layers | Unchanged. It is layer geometry, and the remap rules (§10) apply. Thin layers are no problem for an implicit FV solve. |
| Richards / wetland (Richards.cpp:333 `frac_unfrozen = 1 − frozenfrac`; ice in hydraulic properties) | Front-based frozenfrac; ice from the Stefan partition | Same fields. Ice and liquid come from the curve. After Richards changes liq, the cell is re-equilibrated with advected enthalpy (§10). |
| `Cohort` non-vascular uptake (Cohort.cpp:707) | `minliq·(1−frozenfrac)` | Same field. |
| ALD/ALC outputs, `daily_frontsdepth/type`, restart `frontZ/frontFT` | Stefan fronts | Reconstructed fronts. |
| Frozen-top infiltration (`getInflFrozen`, `frozenInfiltration`) | `frozen` flags | Same flags (definition above). |

---

## 10. Layer restructuring and water redistribution: `EnthalpyRemap`

One helper owns every change to cell content outside the thermal solve. Each call (a) keeps H exact or (b) books ΔH explicitly, then equilibrates the affected cells.

| Event (code) | Rule | Booking |
|---|---|---|
| Split (`splitOneSoilLayer`) | H_child = H_parent·(dz_child/dz_parent) with matching liq/ice/solid fractions, so both children share the parent's T and phase. This replaces the gradient interpolation. | none (rearrangement) |
| Merge (`combineTwoSoilLayersU2L/L2U`) | H = H_u + H_l, then equilibrate with the combined properties. This replaces thickness-weighted T. | none |
| Snow combine/divide (`combineSnowLayers`, `divideSnowLayers`) | Same as merge/split | none |
| Organic thickness change from carbon (`updateOslThickness5Carbon`, deep C transfer, litter) | Solid mass enters or leaves at the cell's T: ΔH = c_s·Δm·T | `edg.solid_mass` |
| Fire (`adjustSoilAfterburn`) | H of the burned layers leaves the column. Any partly burned layer is split first. | `edg.fire_out` |
| Richards, infiltration, drainage, ET (`updateDailySM`) | Moved water carries its source-cell specific enthalpy (liquid at source T: c_w·T per kg). Precipitation and infiltration enter at max(T_air, 0) for rain or 0 °C for melt. After the moves, receiving cells are equilibrated, so water arriving in a frozen cell refreezes and releases latent heat conservatively. | `edg.water_adv` |
| `checkWaterValidity` clipping | Any liq/ice clipped by `maxliq/maxice` leaves at the cell's T | `edg.clip` (measured as 0 in legacy, maxice clipping 0) |

`energy_diag` is extended from a thermal-step window to a **full-day budget**: H_end − H_start = G_top + q_geo + Σ booked terms. The thermal-only residual must be ≤ 1e-6 W m⁻². Per D5, the full-day residual ≤ 1e-6 W m⁻² is a **pass/fail harness test** once Stage C has integrated the solver (unit tests enforce it from Stage C onward; until then it is reported).

---

## 11. Restart compatibility

- **No new restart variables.** Both solvers read and write the same fields: `TSsoil, LIQsoil, ICEsoil, FROZENsoil, FROZENFRACsoil, frontZ, frontFT` (RestartData.h:121–135) and the snow equivalents.
- **Legacy restart → enthalpy run:** the first `updateDailyGroundT` computes H from the stored (T, liq, ice), using the general formula valid off the curve, and equilibrates. Energy is conserved, and T may shift slightly on day 1. The maximum |ΔT| is logged once.
- **Enthalpy restart → legacy run:** the states are physically consistent. Legacy rebuilds its fronts from `frontZ/frontFT` or `rebuildFrontsFromLayerStates()` (Stefan.cpp:715). Supported, but not bit-reproducible across solvers.
- A global attribute `thermal_solver` is added to restart files for provenance. The restart reader only warns on a mismatch.
- Deep C already invalidates production restarts physically (merge plan §1.8), so the enthalpy solver adds no further restart burden for the 4×4 km runs.

---

## 12. Switch and bit-identity

```jsonc
"model_settings": {
  "thermal_solver": 0,       // 0 legacy Stefan + CN (default) | 1 enthalpy NCZ
  "thermal_sfcc": 1,         // 1 only: 1 SJ2016 power law (+ramp, production) | 0 free water+ramp
  "thermal_substeps": 1,     // 1 only: 1,2,4,8,24
  "thermal_k_picard": 0,     // 1 only: extra conductivity passes
  "thermal_nl_solver": 0,    // 1 only: 0 NCZ | 1 Katzenelson (sfcc 0 only; verification)
  "thermal_face_mode": 0,    // 1 only: 0 harmonic (default) | 1 Tubini max rule (fallback)
  "geothermal_flux_mode": 0, // 0 zero (legacy) | 1 single value | 2 spatial map (IO.geothermal_flux_file)
  "geothermal_flux": 0.0     // W m-2, used by mode 1
}
```

- Parsing follows the existing `ModelData.cpp` pattern: absent key means the default, out-of-range values are FATAL, and dependent keys under `thermal_solver = 0` are FATAL if set to non-defaults (as with `thawfrac_moist_mode`).
- `geothermal_flux_mode` is the one shared key. With 0 the legacy code takes literally the same path. Modes 1 and 2 under `thermal_solver = 0` set the legacy `gflx` (a physics change, so off by default).
- **Bit-identity:** all new code is behind `if (thermal_solver == 1)` at a single dispatch in `Soil_Env::updateDailyGroundT` plus the remap hooks. The legacy path is not refactored. The harness checks md5 identity on the calib8 and prod case sets in three configurations: key absent, key = 0, and the pre-change binary.

---

## 13. Data structures, classes and files

New files go under `include/` and `src/`, plain C++11 to match the Makefile's `-std=c++11` (Reburn moves to C++20 later, and the code stays compatible). No new dependencies.

| File | Content |
|---|---|
| `EnthalpyCurve.h/.cpp` | `struct CellProps {double Cs, W, liqmin, eps; int sfcc; double b_mix[4], forg;}`. Functions `phi(T)`, `H(T)`, `dHdT(T)` (C_a), `Tstar()`, `H1/H2` (Jordan split), `invert(H) → T` (safeguarded Newton + bisection), `breakpoints()` (SFCC 0, for Katzenelson). Pure functions with no Layer dependency. |
| `EnthalpySolver.h/.cpp` | `class EnthalpySolver { solve(const Column&, double dt, Result&) }`. Column-level arrays sized `MAX_GRN_LAY+2`, assembly of A and b, NCZ nested loop, Thomas algorithm, residuals and counters. Optional `KatzenelsonSolver` with the same interface. |
| `EnthalpyColumn.h/.cpp` | Gather and scatter between the `Ground` layer list and solver arrays: T, liq, ice, dz, props and λ, with `equilibrate(Layer*)`. |
| `EnthalpyRemap.h/.cpp` | Split, merge, mass-in/out and water-move helpers with booking (§10). Called from `Ground` and `Soil_Env` only when `thermal_solver == 1`. |
| `EnthalpyDiagnostics.h/.cpp` | Front reconstruction, frozen/frozenfrac, ALD/ALC, extended `EnergyDiag` terms, and `energy_daily.csv` columns added at the end of each line so existing parsers keep working. |
| Changes | `ModelData` (keys), `Cohort` (pass settings), `Soil_Env::updateDailyGroundT` (dispatch), `Ground` split/merge/burn (hook calls), `Snow_Env` (skip mass-removal melt, booking), `Soil_Bgc` (expose `sjLayerPhi` as a static free function so `EnthalpyCurve` uses the same code). |

`Stefan`, `TemperatureUpdator` and `CrankNicholson` are untouched.

---

## 14. Test plan

**Unit tests** (doctest, `tests/unit/test_enthalpy_*.cpp`, run with `make unittest`):

1. **Curve:** H is continuous and monotone; C1/C2 hold (C_a ≥ 0, non-decreasing to T\*, non-increasing after); `invert(H(T)) = T` to 1e-10 °C over [−40, 30] °C for every layer type and both SFCC options; SFCC 1 matches `sjLayerPhi` to 1e-12 for T < −ε.
2. **Neumann (Stefan) analytical solution:** water column, parameters of Tubini Table E1, dz = 5 mm and 1 mm, Δt = 1 h. Front error ≤ Tubini Table 2 values (to be read off at implementation). Repeated with dvmdostem-like 0.1 m cells and Δt = 1 d, with a documented, bounded error.
3. **Lunardini (InterFrost T1):** three T_m cases, Δt = 1 h, error in the range of Tubini Tables 3–4.
4. **Conservation:** random columns (snow + moss + organic + mineral + rock, 5–33 cells, random T, W and BCs) over 3650 daily steps. Per-step |ΔH − Δt·ΣG| ≤ 1e-3 J m⁻² per cell, i.e. **≤ 1e-6 W m⁻² column mean**. Remap and booking round-trip: split → merge returns H to 1e-9 relative.
5. **Large Δt:** the same problem with Δt = 1 h, 1 d and 10 d. Zero NCZ fallbacks. Front depth after one season differs from the 1 h result by less than one cell thickness at 1 d (Tubini Fig. 11 behaviour).
6. **Monotone front:** constant cold or warm surface. The reconstructed freezing or thaw front depth is non-decreasing in time, with no oscillation on cell crossings.
7. **No-freeze limit:** column at T > 1 °C everywhere, compared with the legacy `processColumnNofront` path. Agreement is set by CN vs BE accuracy (Tubini §5), so the check is daily-mean |ΔT| < 0.05 °C. The test also confirms that the legacy stencil does not conserve energy while the new one does, documenting the −0.7 W m⁻² term.
8. **Snow:** melt of a snow column under a warm boundary. Melt mass equals ΣQ/L. The excess heat of fully melted layers goes to soil; the `meltSnowLayersAfterT` removal books ≈ 0 J.
9. **Switch parsing:** absent, 0 and 1 keys; FATAL on bad combinations.
10. **Katzenelson cross-check** (if built): same discrete solution as NCZ to 1e-8 °C on SFCC 0 problems.

**Regression and `energy_diag` harness** (`tests/regression`, ≤ 75 yr and ≤ 75 cells, frozen baselines):

- `thermal_solver = 0`: bit-identical to the current baseline (md5) on calib8 and prod. This is a gate.
- `thermal_solver = 1` (SFCC 0, then SFCC 1), `energy_diag = 1`: thermal residual mean and max |res| ≤ **1e-6 W m⁻²** (new gate under `--etol-enth`, separate from the 3 W m⁻² legacy gate). Also gated (D5, from Stage C): **full-day column energy balance ≤ 1e-6 W m⁻²**, zero NCZ fallbacks, water closure PASS, safety-net events not above the legacy count + 10 %, and no run_status failures.
- Reported, not gated: ALD, ALC and soil T at 0.2/0.5/1/2 m against legacy. Reported per front-count class (0, 1, 2, ≥3 fronts), because multi-front years are where the prototype broke.

---

## 15. Performance estimate for 4×4 km HPC runs (Stage A estimate; superseded by the measurements in §15a. Core-hour budget for the full domain deferred, D6)

**Measured** (box, 2026-11-09 ~3:08 PM EST): `-pg` build of `1dca209e` (current default flags `-g`, no `-O`), demo Toolik 2-cell mask, `-p 10 -e 200`, 153,300 cell-days, 28.7 s wall. That is **187 µs per cell-day, or 68 ms per cell-year**, on a single Xeon core.

| Routine (inclusive) | % of sampled CPU |
|---|---|
| `Soil_Env::updateDailyGroundT` | **79.2** |
| `TemperatureUpdator::iterate` → `updateOneTimeStep` → `geBackward` + `cnForward` | 76.5 (geBackward alone 47.1 self) |
| `EnvData::grnd_endOfDay` | 6.0 |
| `Cohort::updateMonthly_Bgc` (all BGC) | 3.8 |
| `Soil_Env::updateDailySM` (+ Richards) | 1.9 |
| `Stefan::updateFronts` | 0.9 |

The legacy solver made **30.6 M tridiagonal solves for 153 k cell-days, about 200 per cell-day**. That is about 100 per sub-column solve, caused by the 0.01 °C Picard tolerance with step halving from 0.5 day down to 1e-5. There were also 441 M `warn_bad_tld` calls and 31 M Boost named-scope constructions.

**Cost model for NCZ** (N ≤ 33 cells). Per inner iteration: one Thomas solve (≈ 8N flops) and N evaluations of φ/H/C_a (one `pow` each for SFCC 1, or a table lookup). That is about 0.3 µs plus 33 × ~0.08 µs, roughly **3 µs at -O0**. Iterations per daily step:

- **Days with no phase-state change** (most of the year): 1 outer × 1–2 inner. Katzenelson in LPJmL needed 1.48 solves per day on average.
- **Freeze/thaw days:** 3–6 outer × 2–4 inner, at most 25 solves.

Estimated annual mean: 3–8 solves per day, or **10–25 µs per cell-day for the thermal step**, versus about 148 µs measured for legacy (0.79 × 187).

| Scenario | Thermal µs/cell-day | Whole model µs/cell-day | Factor vs legacy |
|---|---|---|---|
| Legacy (measured) | ≈148 | 187 | 1.0 |
| NCZ typical (5 solves/day, SFCC 0) | ≈15 | ≈55 | **≈0.3–0.6** |
| NCZ pessimistic (25 solves/day, SFCC 1, `pow`) | ≈75 | ≈115 | ≈0.6 |
| NCZ, `thermal_substeps = 4` and `k_picard = 1`, pessimistic | ≈300 | ≈340 | ≈1.8 |

Conservatively, quote **0.6×–1.3×** until Stage B measures it. These are estimates, not measurements. TR/SC runs with NetCDF output and wetland/CH4 enabled have a larger non-thermal share, which moves the factor toward 1.

**Core-hours.** Per century = cells × 100 yr × 68 ms / 3600 (legacy, -O0, this CPU):

- Per **10,000 cells: ≈ 19 core-hours per century** (legacy).
- Per **100,000 cells: ≈ 190 core-hours per century** (legacy), or about 60–250 with NCZ.

The 4×4 km cell count was not available to me (Q6). Scale linearly with it. An all-Alaska domain is about 1.5 M km², roughly 10⁵ cells at 16 km².

**Independent quick wins** (any solver, recommended for the HPC build): `-O2`, typically 2–4× for this kind of code (to be measured, not assumed); compiling out `warn_bad_tld` and the per-iteration `BOOST_LOG_NAMED_SCOPE` in the hot loop. Both are bit-identity-sensitive, so they go behind a build flag and get their own regression check.

---

## 15a. Compute cost (measured, Stage B)

Measured on calib-elchin, 2026-11-09 16:34–16:43 EST, ≤ 4 concurrent jobs (load average ≤ 4 on 16 cores). Branch `enthalpy/stageB`, tools in `tools/enthalpy_bench/` (README there). Everything in **A** is a measurement; everything in **B** is a projection derived from A.

**Method.**
- The legacy model is built from a throw-away `git archive` copy, instrumented by `bench_patch.py`: `std::chrono` around `Soil_Env::updateDailyGroundT` in Cohort.cpp, plus a counter of tridiagonal solves in `CrankNicholson::geBackward`. The branch sources stay unchanged.
- Cases: `calib8/TVC_NoWarming_nfeed` and `calib8/HPC_NoWarming_nfeed`, each run in full (spin000 7300 d plus main 27375 d).
- For the replay, each run dumps the column state before and after the legacy step (layer kind, dz, T, liq, ice, Cs, W, minliq, texture, organic fraction, λ_unfrozen, λ_frozen, T_s, legacy solves, legacy µs) for the first 3650 d of spin000 and of main. That gives 14,595 column-days, ten full years per run, covering every season including freeze-up and thaw days.
- `bench_enthalpy` replays every dumped day with the Stage B solver: NCZ with SFCC 0, NCZ with SFCC 1, and Katzenelson with SFCC 0. It uses the same Δt = 1 d, the same columns (snow + soil + rock, N ≈ 25–30 cells), Dirichlet T_s, and tolerance 1e-3 J m⁻². Each day is timed by `std::chrono` the same way as the legacy step, minimum of 3 repeats.
- "phase" days are days on which any cell crosses 0 °C or holds latent heat in the ramp. "nophase" days are all the others.
- Both builds use the model's own flags (`-g -ansi -std=c++11`), at the default optimisation (-O0) and with -O2.

### A. Measured

**A1. Legacy thermal step, in the model (timing-only runs, no dump), µs per column-day.**

| build | TVC spin000 | TVC main | HPC spin000 | HPC main | legacy solves/day (main) |
|---|---|---|---|---|---|
| -O0 | 440.5 | 442.7 | 491.6 | 506.0 | 739 (TVC), 754 (HPC) |
| -O2 | 206.5 | 207.9 | 234.0 | 240.5 | same (identical solve counts) |

The dump runs agree within 1–2 % (440.2–506.2 µs at -O0), so the dump overhead is outside the timed region.

**A2. Replay, all 14,595 column-days (TVC + HPC). µs per column-day and tridiagonal solves per day: mean / p95 / max.**

| solver | build | µs/day | solves/day | Newton outer iterations/day | converged | max daily closure |
|---|---|---|---|---|---|---|
| legacy (same days, from dump) | -O0 | 461 / 1428 / 4355 | 734 / 1887 / 5534 | n/a | n/a | not closed (§1) |
| legacy (same days, from dump) | -O2 | 217 / 682 / 2068 | 734 / 1887 / 5534 | n/a | n/a | n/a |
| NCZ, SFCC 0 | -O0 | 24.7 / 38.6 / 78.0 | 1.84 / 4 / 12 | 1.63 / 3 / 5 | 14595/14595 | 2.0e-11 W m⁻² |
| NCZ, SFCC 0 | -O2 | **4.6** / 7.0 / 14.0 | 1.84 / 4 / 12 | 1.63 / 3 / 5 | 100 % | 2.0e-11 |
| NCZ, SFCC 1 (production) | -O0 | 62.0 / 119.4 / 186.1 | 3.42 / 8 / 14 | 1.63 / 3 / 5 | 100 % | 3.4e-8 |
| NCZ, SFCC 1 (production) | -O2 | **27.9** / 56.5 / 91.2 | 3.42 / 8 / 14 | 1.63 / 3 / 5 | 100 % | 3.4e-8 |
| Katzenelson, SFCC 0 | -O0 | 17.1 / 25.4 / 60.5 | 1.51 / 4 / 13 | (pivots) 1.51 / 4 / 13 | 100 % | 2.5e-12 |
| Katzenelson, SFCC 0 | -O2 | **3.5** / 5.1 / 12.1 | 1.51 / 4 / 13 | 1.51 / 4 / 13 | 100 % | 2.5e-12 |

**A3. Split by day class at -O2, µs/day mean, with solves/day mean in brackets.**

| | legacy | NCZ sfcc0 | NCZ sfcc1 | Katz sfcc0 |
|---|---|---|---|---|
| phase days (6424) | 112 (633) | 5.3 (2.49) | 33.3 (4.59) | 3.9 (2.16) |
| no-phase days (8171) | 299 (813) | 4.1 (1.32) | 23.6 (2.50) | 3.2 (1.00) |

Per site at -O2, as mean µs/day for legacy / NCZ sfcc1: TVC 204 / 32.3 and HPC 230 / 23.5. The p95 of NCZ sfcc1 solves/day is 9 at TVC and 7 at HPC.

**A4. Thermal share of whole-run wall time on calib8.** This uses the main stage, the in-model thermal seconds over the run's integer wall seconds, and includes I/O and output. It is about **52 % at -O0** (TVC 12.1/24 s, HPC 13.9/26 s) and about **34 % at -O2** (5.7/18 s, 6.6/18 s). The 79 % figure (§15) is a gprof CPU share on the demo column at -O0 without I/O. The two are different measures, so both projections are given below.

**A5. Agreement with the legacy solution.** Per day, the maximum |T_new − T_legacy| over snow and soil cells has mean 1.6 °C, p95 4.8 °C, and max 16 °C, the max on a snow day at HPC. This is the physics difference between a closed-energy implicit step and the legacy split front scheme, starting from the same state with no feedback. It is not an accuracy measure: the legacy step does not conserve energy (§1). Stage D evaluates against observations.

### B. Projected (derived from A; not measured)

Here r = new µs/day ÷ legacy µs/day, from A2, same build, same days. The whole-model factor is (1 − s) + s·r. A factor below 1 means the model is faster.

| solver | build | r | thermal speed-up | factor, s = 0.79 (demo gprof share, as requested) | factor, s = calib8 wall share (A4) |
|---|---|---|---|---|---|
| NCZ sfcc0 | -O0 | 0.053 | 18.7× | 0.25 | 0.51 |
| NCZ sfcc0 | -O2 | 0.021 | 46.9× | 0.23 | 0.67 |
| **NCZ sfcc1** | -O0 | 0.134 | 7.4× | **0.32** | 0.55 |
| **NCZ sfcc1** | -O2 | 0.129 | 7.8× | **0.31** | 0.70 |
| Katz sfcc0 | -O0 | 0.037 | 27.0× | 0.24 | 0.50 |
| Katz sfcc0 | -O2 | 0.016 | 62.3× | 0.22 | 0.67 |

**Stage C overhead allowance.** The replay leaves out gather and scatter of layer state, the λ update, front and ALD reconstruction, `energy_diag`, and remap. Even if all of these triple the solver cost (3r), NCZ sfcc1 gives a factor of 0.53 (O0) or 0.52 (O2) at s = 0.79. Whatever the share, the model is not projected to get slower.

**Core-hours per century per 10,000 cells.** One century is 36,525 cell-days.

| | -O0 | -O2 |
|---|---|---|
| legacy thermal step alone (measured µs × days) | 46.8 | 22.0 |
| legacy whole model, via s = 0.79 | 59.2 | 27.9 |
| legacy whole model, calib8 wall time (incl. I/O) | 92.7 | 66.7 |
| enthalpy thermal step alone, NCZ sfcc1 (replay) | 6.3 | 2.8 |
| whole model with NCZ sfcc1, via s = 0.79 | **18.7** | **8.7** |
| whole model with NCZ sfcc1, calib8 wall share | 51 | 47 |
| whole model with NCZ sfcc0 / Katz, via s = 0.79 | 14.9 / 14.2 | 6.3 / 6.2 |

**Notes.**
- (i) SFCC 1 costs about 6× SFCC 0 at -O2 even though it needs only 1.9× the solves. The difference is the `pow()` calls inside the SJ2016 curve and its derivative, evaluated per cell per inner iteration. Caching the per-cell exponent or tabulating φ could remove most of it. That is a Stage C option, and it must keep the result bit-compatible with `sjLayerPhi` on its own grid.
- (ii) The -O2 model build only compiles with `-Wno-error=aggressive-loop-optimizations`. A pre-existing out-of-bounds read `cd->m_soil.z[il+1]` at `il = MAX_SOI_LAY−1` (src/EnvData.cpp:952 and :997) trips `-Werror` at -O2. This is not caused by Stage B. It is not fixed, because it is model-path code, but it should be fixed before -O2 HPC production builds.
- (iii) Solves per day are close to Hötten et al. (2026), 1.48 per backward-Euler daily step, for SFCC 0. SFCC 1 needs about twice as many because of its curvature below −ε.

### D. Stage B unit tests (`tests/unit/test_enthalpy_solver.cpp`; on calib-elchin `make unittest`: 290/290 test cases pass, the 281 existing tests plus 9 new; the 6 failed assertions are the pre-existing `[KNOWN BUG: xfail]` cases)

- SFCC 1 is bit-identical to `sjLayerPhi` over the texture and organic grid. The curve properties hold: H is monotone and continuous, H(0) ≈ 0, and the liq/ice bounds are respected.
- Neumann/Stefan analytic solution, 400 cells, dz 5 mm, Δt 1 h, 10 d: λ-ratio γ = 0.1704, front at 0.3245 m, maximum front error 7.9 mm (< 2 dz).
- Energy conservation over 365 days, with sources and q_bot = 0.06 W m⁻²: maximum daily closure 9.1e-12 W m⁻² for SFCC 0 and 2.3e-8 W m⁻² for SFCC 1 (gate 1e-6).
- Large Δt: converges at 1 h, 1 d and 10 d with no fallbacks. Thaw depth after 120 d is 0.26998 m (1 h), 0.26999 m (1 d) and 0.363 m (10 d).
- Monotone front propagation, with no spurious oscillation.
- Two-layer steady flux: the harmonic mean is exact for both Dirichlet and Robin. The max rule deviates, as documented.
- Katzenelson vs NCZ on 200 random free-water columns: agree to round-off. Katzenelson refuses SFCC 1.
- Microbenchmark (synthetic column, -O0): NCZ sfcc0 1.43 solves/day, 15 µs; NCZ sfcc1 2.82, 36 µs; Katzenelson 1.04, 10 µs.
- Model path: the library is compiled into dvmdostem but never called. calib8 regression, 8 cases, -j 4: A==B bit-identical in 8/8, determinism True, RESULT PASS.

### C. Stage B deviations from this design

1. SFCC 0 defines H as exactly piecewise affine on the ramp, with liq diagnosed from H (§5), instead of a linear φ ramp. This keeps Katzenelson exact. SFCC 1 keeps the linear-φ ramp.
2. The SJ2016 curve is a textual copy of `Soil_Bgc::sjLayerPhi` (bit-equality enforced by unit test) rather than a shared free function, because Stage B must not touch model sources. Stage C refactors it into one shared function.
3. One library (`EnthalpySolver.h/.cpp`: curve, NCZ, Thomas, Katzenelson) replaces the separate `EnthalpyCurve` and `EnthalpySolver` files of §13. It uses `std::vector`, so N is arbitrary. `EnthalpyRemap` moves to Stage C.
4. NCZ: because h₂ is linear in the Jordan split, the outer right-hand side stays b, and the outer loop only updates the linearisation point. Initial guess is T₀ = min(Tⁿ, 0). The convergence test is |residual| ≤ 1e-3 J m⁻² plus a round-off floor of 8 ε_mach·max(|b|, |AT| + |diag·T|). The floor is needed for Δt ≥ 1e6 s.
5. Replay states are taken at the legacy (T, liq, ice), which can be off the curve. The benchmark puts them on the curve via H(T) at the legacy T. Stage C must define the on-curve start state at switch-on and at restart (§11).

## 16. Risks and calibration impact

| Risk | Mitigation |
|---|---|
| **Calibration invalidated** (certain with `thermal_solver = 1`). ALD, soil T and thaw timing will change. The prototype changed ALD by 0.07–0.46 m, and the direction under the new solver is unknown. | Stage E: recalibrate with deep C (same campaign, merge plan §1.8). Keep T125/S23 with `thermal_solver = 0` for production until then. |
| SFCC 1 shifts decomposition (thawfrac, rh_frozen_mode) and Richards (liquid in frozen soil). | SFCC 1 is the production curve (D1). D still runs SFCC 0 first (numerics only), then SFCC 1 (physics), so the two effects can be attributed separately. |
| Snow melt rate changes: conducted melt replaces degree-second melt, SEB surplus moves to a flux source, and the meltsnow T_top = 0 rule is dropped. | Approved (D3). Stage D snow-off-date check vs TVC/HPC observations. |
| Daily backward Euler damps fast surface transients. Daily forcing already removes the diurnal signal. | `thermal_substeps` test 1/4/24. Tubini 1 h vs 1 d ≤ 0.04 °C mean. |
| Pre-existing -O2 `-Werror` failure (EnvData.cpp:952/997 out-of-bounds read, §15a note ii) blocks optimised HPC builds. | Fix in the model path as a separate, reviewed change (not part of Stage B). |
| Lagged λ across a freeze transition (first-order error). | `thermal_k_picard` test in D. |
| Remap hooks are missed in some code path (e.g. a new Reburn layer path), breaking closure silently. | Full-day budget diag plus a unit test per hook. Reburn `FireBurned` layer remap reviewed when Reburn is merged (Phase 3). |
| Integration with wetland/Richards ice could destabilise water as the prototype did. | Re-equilibration only between modules, water mass never changed by the thermal step; harness gates on water closure and safety-net events. |
| SFCC 1 `pow` cost: measured at 6× SFCC 0 at -O2 (§15a), still 7.8× faster than legacy. | Tabulated φ, or a per-cell cached exponent, in Stage C if needed. |
| C++11 vs Reburn C++20. | Code written to both. |

---

## 17. Staged plan

| Stage | Content | Exit criteria | Effort |
|---|---|---|---|
| **B: standalone library + tests** | `EnthalpySolver` (curve, NCZ, Thomas, Katzenelson), analytic tests, benchmark replay of calib8 TVC/HPC columns against the legacy thermal step (-O0 and -O2). No model wiring. `EnthalpyRemap` moved to C. | All unit tests pass. Conservation ≤ 1e-6 W m⁻². Zero fallbacks at Δt up to 10 d. Measured cost in §15a. | **9–11 days**: done 2026-11-09 (branch `enthalpy/stageB`, tag `stageB-solver`) |
| **C: integration behind switch** | ModelData keys (incl. `geothermal_flux_mode` zero/value/map, D4), dispatch in `updateDailyGroundT`, gather/scatter, snow hand-off, SEB melt source, remap hooks in Ground/Snow_Env/updateDailySM/fire, front/ALD reconstruction, restart attribute, extended energy_diag, harness `--etol-enth`. | `thermal_solver = 0` md5-identical. `= 1` passes all gates on calib8 and prod (SFCC 1 production, SFCC 0 verification), **including full-day balance ≤ 1e-6 W m⁻² as pass/fail (D5)**. Restart round-trip in both directions. | **12–16 days** |
| **D: site evaluation (TVC, HPC)** | 1950–2024 runs, legacy vs enthalpy (SFCC 0/1, substeps 1/4, k_picard 0/1). Compare with observed ALD, soil T profiles and snow-off dates. Multi-front-year analysis. Timing on a ≤ 75-cell 4×4 km tile with MPI on the HPC. | Closure ≤ 1e-6 W m⁻². ALD/soil T bias vs obs at least as good as legacy, or explained. Measured core-hours per century per 10,000 cells (full-domain budget only when Elchin asks, D6). | **6–8 days** (+ compute) |
| **E: calibration plan** | Which parameters move (n-factors, organic tc, `tcsolid`, snow density, possibly SJ exponents); targets; combine with the deep C recalibration; config for production 4×4 km. | Plan approved by Elchin. Recalibration itself shared with deep C (15–25 person-days in the merge plan). | **2 days** (plan) |

Total to a validated, switchable solver: **about 29–37 working days** (Katzenelson added to B), excluding recalibration; B is done.

---

## 18. Open questions for Elchin (answered 2026-11-09, see §0a; kept for the record)

1. **Default freezing curve for production:** SFCC 1 (SJ2016, physically consistent with rh_frozen_mode, changes physics) or SFCC 0 (free water, numerics-only change)? My recommendation is to evaluate both in D and decide then. → **D1: SFCC 1 for production; SFCC 0 kept as option/verification.**
2. **Face conductance:** harmonic (recommended for layered organic/mineral) or Tubini's max? Only harmonic would be implemented, unless you want both as an option for comparison. → **D2: harmonic; max rule as optional documented fallback, default off.**
3. **Snow melt:** OK to replace `rad_melt`/SEB mass removal with a top-flux source when `thermal_solver = 1`, and to drop the `meltsnow → T_top = 0` rule? → **D3: yes.**
4. **Geothermal flux:** keep 0 (legacy), or use a value or map for the 4×4 km domain? This matters for century-scale deep permafrost in 55 m columns. → **D4: config option: zero (default), single value, or spatial map.**
5. **Full-day energy budget:** should it become a gated harness criterion (≤ 1e-6 W m⁻²) once Stage C is done, or stay reported only? → **D5: yes, gated once Stage C integrates.**
6. **4×4 km domain size** (cell count) and target HPC node type, to turn §15 into a core-hour budget. → **D6: deferred until Elchin asks.**
7. **Katzenelson cross-check:** worth the extra ~2 days in Stage B? → **D7: yes (done in Stage B).**

---

## 19. References (checked against doi.org or the publisher page on 2026-11-09; Tubini 2021 and Hötten 2026 full texts re-read this session)

- Tubini, N., Gruber, S., and Rigon, R. (2021). A method for solving heat transfer with phase change in ice or soil that allows for large time steps while guaranteeing energy conservation. *The Cryosphere*, 15, 2541–2568. https://doi.org/10.5194/tc-15-2541-2021
- Hötten, D., Niebsch, J., Ramlau, R., and Zulehner, W. (2026). An exact piecewise affine solver for the discretized enthalpy formulation of the Stefan problem in permafrost modeling. *Electron. Trans. Numer. Anal.*, 65, 553–578. https://doi.org/10.1553/etna_vol65s553
- Casulli, V., and Zanolli, P. (2010). A nested Newton-type algorithm for finite volume methods solving Richards' equation in mixed form. *SIAM J. Sci. Comput.*, 32, 2255–2273. https://doi.org/10.1137/100786320
- Casulli, V., and Zanolli, P. (2012). Iterative solutions of mildly nonlinear systems. *J. Comput. Appl. Math.*, 236, 3937–3947. https://doi.org/10.1016/j.cam.2012.02.042
- Dall'Amico, M., Endrizzi, S., Gruber, S., and Rigon, R. (2011). A robust and energy-conserving model of freezing variably-saturated soil. *The Cryosphere*, 5, 469–484. https://doi.org/10.5194/tc-5-469-2011
- Nicolsky, D. J., Romanovsky, V. E., Alexeev, V. A., and Lawrence, D. M. (2007). Improved modeling of permafrost dynamics in a GCM land-surface scheme. *Geophys. Res. Lett.*, 34, L08501. https://doi.org/10.1029/2007GL029525
- Voller, V. R., Swaminathan, C. R., and Thomas, B. G. (1990). Fixed grid techniques for phase change problems: a review. *Int. J. Numer. Methods Eng.*, 30, 875–898. https://doi.org/10.1002/nme.1620300419
- Kurylyk, B. L., and Watanabe, K. (2013). The mathematical representation of freezing and thawing processes in variably-saturated, non-deformable soils. *Adv. Water Resour.*, 60, 160–177. https://doi.org/10.1016/j.advwatres.2013.07.016
- Goodrich, L. E. (1978). Efficient numerical technique for one-dimensional thermal problems with phase change. *Int. J. Heat Mass Transfer*, 21, 615–621. https://doi.org/10.1016/0017-9310(78)90058-3
- Schaefer, K., and Jafarov, E. (2016). A parameterization of respiration in frozen soils based on substrate availability. *Biogeosciences*, 13, 1991–2001. https://doi.org/10.5194/bg-13-1991-2016 (unfrozen-water curve used in `Soil_Bgc::sjPhiPure/sjLayerPhi`, Eqs. 4–6, Table 1)
- Johansen, O. (1975) and Lawrence & Slater (2008), conductivity: as cited in `SoilLayer.cpp`.
