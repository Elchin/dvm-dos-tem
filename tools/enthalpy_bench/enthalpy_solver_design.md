# dvmdostem: fully implicit, globally convergent enthalpy solver for soil/snow heat (Option C)

*Stage A design document, updated with Elchin's decisions (§0a), Stage B results (§15a) and Stage C integration results (§15b). TEM-dev, 2026-11-09 (EST). Stage A wrote no model code; Stages B and C are on local branches only. The Stage A code review was read-only against `wetland_deep_C` at `1dca209e` (box worktree `/workspace/tem-dev/wdc`). Line references are to that commit. The runtime profile in §15 comes from a separate `-pg` build of the same commit, made and run on the box (`/workspace/tem-dev/prof`, not the repo).*

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
| Time step | Four backward-Euler steps per day by default when `thermal_solver` = 1 (`thermal_substeps` = 4; Elchin, 2026-11-10). There is no adaptive step cutting, because NCZ converges for any Δt. Tolerance is a per-day budget (see §15d). |
| Fronts and ALD | **Diagnostic only.** They are reconstructed from the enthalpy/ice profile into the existing `frontsz`/`frontstype`, `frozen`, `frozenfrac` and `ald`/`alc` fields, so downstream code reads the same variables. |
| Remap | Layer split, merge, burn, carbon-driven thickness change and snow construction all go through one `EnthalpyRemap` helper. It **conserves H** for rearrangements and books ΔH explicitly for mass that enters or leaves. |
| Switch | `model_settings.thermal_solver`: 0 = legacy (default), 1 = enthalpy. With 0 the legacy code path is untouched, verified bit-identical by md5 in the harness. |
| Restart | No new restart variables. Both modes read and write the same T, LIQ, ICE, FROZEN, FROZENFRAC and front fields. |
| Cost | Measured (§15, §15a): the legacy thermal step is **79 % of CPU** (demo profile) and does about **200 tridiagonal solves per cell-day**. Stage B measurements on calib8 TVC/HPC columns are in §15a; the core-hour budget for the 4×4 km domain is deferred until Elchin asks (D6). |
| Stages | B: standalone library + tests + Katzenelson cross-check (D7), 9–11 d (done 2026-11-09, §15a). C: integration behind the switch, 12–16 d (done 2026-11-09, §15b). D: TVC/HPC site evaluation, 6–8 d. E: calibration plan, 2 d, with the recalibration itself (15–25 d) shared with deep C. |

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
- The `nstep = 24` rule for an abrupt surface-state change (Soil_Env.cpp:293–300) is replaced by `thermal_substeps` (default 4 when `thermal_solver` = 1; Elchin, 2026-11-10). The legacy rule existed only to keep the front logic stable.

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

- Δt = 86400 s / `thermal_substeps` (1, 2, 4, 8, 12, 24 or 48; default 4 when `thermal_solver` = 1; Elchin, 2026-11-10). The nonlinear system is solved once per substep.
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
  "thermal_substeps": 4,     // 1 only: 1,2,4,8,12,24,48 (default 4 when thermal_solver=1; Elchin, 2026-11-10)
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

## 15b. Stage C: integration behind `thermal_solver` (done 2026-11-09, branch `enthalpy/stageC`, tag `stageC-integration`; local only, nothing pushed)

Commits on calib-elchin `~/tem-dev/esc`: `accddbfc` (integration), `5d783561` merged via `9bfb8941`/`6826f0f3` (EnvData one-past-end fix, bit-identical, -O2 -Werror builds), `d5b68b3b` (FATAL layer dump on non-finite column H; closure gate fails on NaN), `e44c5433` (rock heat capacity from `ParentLayer`: `Layer::vhcsolid` is shadowed and stays MISSING for rock), plus the docs/tools checkpoint. Builds in `~/tem-dev/esc-work`: `ref` = `wetland_deep_C` 5d783561 (legacy reference), `O0`/`O2` = Stage C. All runs 2026-11-09 6:20–6:57 PM EST, ≤ 4 concurrent jobs.

### What is implemented
| Item | Where | Notes |
|---|---|---|
| Keys `thermal_solver` (0 default), `thermal_sfcc` (1 SJ2016 default, 0 free water), `thermal_face_mode` (0 harmonic, 1 max), `thermal_nl_solver` (0 NCZ, 1 Katzenelson, sfcc 0 only), `thermal_substeps`, `geothermal_flux_mode` (0 zero, 1 `geothermal_flux`, 2 map `IO.geothermal_flux_file`, NetCDF var `geothermal_flux`, W m⁻²) | ModelData.cpp | FATAL on out-of-range or on non-default dependent keys with `thermal_solver = 0` |
| One SJ2016 curve | `SJCurve.h/.cpp` | `Soil_Bgc::sjLayerPhi`/`sjOrgFraction` are wrappers; legacy bit-identical |
| Daily enthalpy step | `Soil_Env_enthalpy.cpp`, dispatched from `Soil_Env::updateDailyGroundT` | Robin top (cn 100), no fixed-0 °C-while-melting rule, SEB/radiation melt as a heat source in the top snow cell (Snow_Env mass removal skipped), geothermal bottom flux, dt halving fallback (counted), heat of fully melted snow cells passed down, pore-volume ice cap as legacy |
| Fronts, `frozen`, `frozenfrac`, ALD/ALC | reconstructed daily from the profile (`enthalpy::reconstructFronts`) | downstream readers unchanged |
| Switch-on / restart | fixed-T projection onto the curve on the first step after init or restart, booked as `init` | no new restart variables |
| Remap / non-thermal changes | `Ground::enthalpy_remap` (split/merge conserve H); Cohort booking windows for snow, soil water, layer/BGC thickness, fire (`bsnow`, `bwater`, `blayer`, `bfire`) | see deviation 2 |
| `stefan_front_resist_fix` | ignored with a WARNING when `thermal_solver = 1` | verified: outputs md5-identical with 0 or 1 |
| Harness | `closure.py --etol-enth` (default 1e-6 W m⁻²): max daily solve residual and full-day residual, zero non-converged steps, non-finite = FAIL | |
| Unit tests | `tests/unit/test_enthalpy_integration.cpp` | 300/300 test cases; 6 failed assertions are the known `[KNOWN BUG: xfail]` cases |

### Results
**C1. `thermal_solver = 0` bit-identity** (run_regression.sh, A = ref, B = Stage C -O0, 8 cases each, md5 of spin-up and main outputs): dc0_calib8, dc1_calib8, dc0_prod, dc1_prod all **A==B True, determinism True, PASS**.

**C2. `thermal_solver = 1` closure** (O0, 20-yr spin-up + 75-yr main, energy_diag 1). All PASS.

| set | cases | max daily solve residual | max full-day residual | non-converged | water % P | C max annual | N mean |
|---|---|---|---|---|---|---|---|
| calib8 dc0, SFCC 1 | 8/8 | 4.1e-8 W m⁻² | 4.1e-8 | 0 | −0.27…+0.20 | ≤ 0.0037 | ≤ 0.035 |
| calib8 dc1, SFCC 1 | 8/8 | 3.6e-8 | 3.6e-8 | 0 | −0.37…−0.10 | ≤ 0.0049 | ≤ 0.25 |
| prod dc1, SFCC 1 | 8/8 | 4.6e-8 | 4.6e-8 | 0 | −0.64…−0.06 | ≤ 0.0066 | ≤ 4e-5 |
| calib8 dc0, SFCC 0 | 8/8 | 6.4e-11 | 6.4e-11 | 0 | −0.21…+0.26 | ≤ 0.0037 | ≤ 0.023 |
| single-case options (HPC_Warming_nfeed): geothermal mode 1 (0.06), mode 2 (map 0.06), face max, substeps 4, Katzenelson | 5/5 | ≤ 4.0e-8 | ≤ 4.0e-8 | 0 | | | |
| restart legacy → enthalpy (TVC, HPC Warming_nfeed) | 2/2 | ≤ 4.2e-8 | ≤ 4.2e-8 | 0 | | | init booked 0.016–0.021 W m⁻² |
| restart enthalpy → legacy | 2/2 | legacy gate (3 W m⁻²) | | | | | |

Booked between-step terms, run means (calib8 SFCC 1): `bsnow` −0.89…−1.13 W m⁻² (snowfall enters as ice, −L per kg), `bwater` +0.03…+0.06, `blayer` ≤ 0.014, `clip` ≤ 0.004, `bfire` 0. Mode 2 with a constant map is md5-identical to mode 1. Mean NCZ solves per day 2.4–3.0 (SFCC 1).

**C3. First look, enthalpy (SFCC 1) − legacy, last 30 years of the 75-yr main stage.** ALD in m; soil T at 0.2/0.5/1/2 m in °C.

| set | ΔALD (m) | ΔT annual, range over depths and cases (°C) | ΔGPP | ΔNEE (gC m⁻² yr⁻¹) |
|---|---|---|---|---|
| calib8 dc0 TVC | −0.05…−0.14 | −1.0…−1.9 | −1.7…+1.2 % | −4…−26 (less source) |
| calib8 dc0 HPC | −0.26…−0.37 | 0…−1.9 | −2.5…+1.9 % | −3.3…+0.9 |
| calib8 dc1 TVC | −0.05…−0.14 | −1.1…−1.8 | −2.4…+1.1 % | −4…−27 |
| calib8 dc1 HPC | −0.38…−0.49 | −0.4…−1.8 | −1.6…−14.5 % (nfeed largest) | −4.3…+9.3 |
| prod dc1 TVC | −0.14…−0.20 | −2.6…−3.6 | −0.5…−2.8 % | −7.5…−34 |
| prod dc1 HPC | −0.16…−0.34 | −1.6…−2.3 | −0.8…−24 % (nfeed largest) | −2.1…+1.4 |

Attribution:
- **The cooling is not the freezing curve.** SFCC 0 (free water) gives nearly the same ALD and T as SFCC 1 (ΔALD within 0.04 m, ΔT within 0.25 °C).
- **It matches the legacy energy leak.** In the legacy runs the column loses −0.005…−1.49 W m⁻² through the top on average, balanced by a spurious +0.16…+1.65 W m⁻² residual. The enthalpy runs sit near zero mean top flux (−0.11…+0.05). The cooling is largest where the legacy residual is largest (TVC > HPC; prod, with about 2 W m⁻², > calib8). So the legacy parameterisation was tuned with an extra 0.2–2.5 W m⁻² of heat.
- **Option sensitivities** (HPC_Warming_nfeed, ALD vs the SFCC 1 base): geothermal 0.06 W m⁻² +0.05 m; face max rule +0.12 m; `thermal_substeps = 4` +0.04 m (soil T −0.1 °C), a daily-step time error larger than the Tubini estimate; Katzenelson vs NCZ (SFCC 0) +0.003 m, ≤ 0.14 °C after 75 yr of BGC feedback.

**C4. Whole-model runtime** (sequential, one process, load ≈ 1, calib8 dc0, 8 cases, 34,675 cell-days per case; wall time includes harness start-up, input and NetCDF output).

| build | legacy s/case | enthalpy SFCC 1 s/case | µs per cell-day legacy / enthalpy | factor | core-h per century per 10,000 cells legacy / enthalpy |
|---|---|---|---|---|---|
| -O0 | 31.6 | 18.7 | 912 / 541 | **0.59** | 92.5 / 54.9 |
| -O2 | 22.1 | 16.3 | 639 / 469 | **0.73** | 64.8 / 47.6 |

This matches the Stage B projection at the calib8 wall share (0.55 / 0.70). At -O2 the non-thermal part (BGC, I/O) now dominates.

### Stage C deviations from this design
1. `geothermal_flux_mode` 1/2 require `thermal_solver = 1`. The legacy solver keeps its zero-flux bottom and FATALs on a non-zero mode, so the legacy path stays untouched.
2. Remap: split/merge conserve H (`Ground::enthalpy_remap`). Carbon-driven thickness change, snow construction and the Richards/infiltration water moves are not rewritten per §10. Their enthalpy change is **measured and booked** per module window (`blayer`, `bsnow`, `bwater`, `bfire`) and the cells are re-equilibrated conservatively before the next solve. The full-day gate therefore proves that nothing changes H outside the solve or a booked window (`unattr` = 0); it does not independently check the physics of each booked term, e.g. the rain temperature. The booked magnitudes are physically plausible (`bsnow` ≈ −L × snowfall). **Resolved in C2 (§15c).**
3. No `thermal_solver` provenance attribute in restart files yet (§11). Restarts round-trip correctly in both directions without it. **Resolved in C2 (§15c).**
4. The SFCC 1 `pow` caching (§15a note i) was not needed: the measured whole-model factor is below 1.
5. The integrator safety net fires more in two cases, against the §14 "≤ legacy + 10 %" criterion: calib8 dc0 HPC_Warming_nfeed 0 → 4 (main), prod dc1 HPC_Warming_nfeed 113 → 127 (+12 %) and HPC_NoWarming_nfeed 105 → 114 (+9 %). These are N/C soil integrator events in colder soils. C/N closure still passes. Reported, not gated. Root cause in §15c D3.

## 15c. Stage C2 and Stage D (done 2026-11-09, local only, nothing pushed)

### C2: energy-conserving remap and restart label (tag `stageC2` = `15fcca85`, branch `enthalpy/stageC`)
Commits `f8f518b9` (remap of the snow and soil-water windows: every moved kg carries its enthalpy, mass entering or leaving carries a physically defined enthalpy — rain at max(T_air, 0 °C), melt and ponded water at 0 °C — and cell T is set from the resulting H at the module's liq/ice, exact because H = C·T − L·ice is linear in T at fixed phase; layer-window water check; restart attribute `thermal_solver`, written only when 1, absent = legacy, mismatch only warns) and `15fcca85` (C2b: without snow cells the direct-melt latent heat is shared by the top 0.10 m of soil by heat capacity instead of the first cell; a thin dry moss cell had reached −2850 °C in C2a). Carbon-driven thickness change (`blayer`, solid mass at cell T) and fire (`bfire`, burned H leaves) were already booked per the §10 rules. This closes Stage C deviations 2 and 3.

Verification (`tools/enthalpy_bench/stagec/c2_vm.sh 15fcca85`, log `esc-work/c2b_verify.log`, 2026-11-09 7:54–8:14 PM EST):

| check | result |
|---|---|
| build -O0 and -O2 (-Werror) | OK |
| unit tests | 304/304 test cases; the 6 failed assertions are the known `[KNOWN BUG: xfail]` cases |
| `thermal_solver = 0` md5 vs ref (`wetland_deep_C` 5d783561) | dc0_calib8, dc1_calib8, dc0_prod, dc1_prod: A==B bit-identical, deterministic, **PASS** |
| `thermal_solver = 1` full-day closure (calib8 dc0/dc1, prod dc1, calib8 dc0 SFCC 0; 8 cases each) | 32/32 PASS; max full-day residual 3.7e-8 W m⁻² (SFCC 0: 1.1e-10); unattributed 0; non-converged 0; remap fallbacks 0 |
| water / C / N closure (ts1) | water −0.42…+0.59 % of P; C max annual ≤ 0.007 g m⁻²; N mean ≤ 0.18 g m⁻² yr⁻¹ |
| remap terms, run means | `rsnow` −0.012…−0.045 W m⁻², `rwater` −0.04…+0.004 (calib8), −0.19…−0.37 (prod) |
| restart label | ts1 restart carries `:thermal_solver = 1`; ts0 restart has no attribute (unchanged format) |
| restart round-trip | legacy → enthalpy 2/2 PASS (≤ 3.5e-8 W m⁻²), enthalpy → legacy 2/2 OK; both warn once on the mismatch |

Integrator safety-net events (main stage, HPC nfeed cases; TVC and nonfeed cases 0 in all runs):

| set | legacy | Stage C | C2b |
|---|---|---|---|
| calib8 dc0 HPC_Warming_nfeed | 0 | 4 | 3 |
| calib8 dc1 HPC_NoWarming / Warming_nfeed | 6 / 4 | 103 / 92 | 64 / 60 |
| prod dc1 HPC_NoWarming / Warming_nfeed | 105 / 113 | 114 / 127 | 110 / 119 |

### D1. Model vs observations (TVC, HPC)
Tool `tools/enthalpy_bench/staged/stageD_eval.py` (obs in `staged/obs_stageD.json`, results `staged/eval_c2b.json`). Window 2010–2024 of the 1950–2024 main stage; forcing is bias-corrected CanRCM4, so the comparison is climatological. Values are means over the two NoWarming cases (nfeed, nonfeed). Obs: TVC soil T AWI PANGAEA.962726 (2–20 cm), ALD GSC thaw tubes 91TT11/93TT02 (TVC 0.66 m) and HPC thaw tube (0.75 m), HPC 3–5 cm surface T and the GSC 01TC02 envelope (3.5 km E of HPC), peak SWE TVC EOW survey 2005–24 (125 mm), HPC Pomeroy & Marsh 1997 via Krogh 2017 (152 mm), snow-off ECCC Trail Valley / Inuvik (first of 5 snow-free days after 1 Apr; TVC DOY 162, HPC 141). leg = legacy, ts1 = enthalpy SFCC 1, s0 = enthalpy SFCC 0.

TVC (bias = model − obs; T in °C, RMSE over the 12 monthly climatology values):

| run | ALD bias (m) | T2cm bias / RMSE | T20cm bias / RMSE | T10cm DJF / MAM / JJA / SON bias | peak SWE bias (mm) | snow-off bias (d) |
|---|---|---|---|---|---|---|
| calib8 dc0 leg | −0.23 | −2.9 / 4.2 | −3.0 / 4.2 | −6.1 / +0.7 / −1.4 / −4.8 | −36 | −38 |
| calib8 dc0 ts1 | −0.28 | −4.5 / 6.2 | −4.5 / 6.0 | −9.2 / −0.2 / −0.9 / −6.9 | −34 | −35 |
| calib8 dc0 s0 | −0.30 | −4.7 / 6.5 | −4.6 / 6.2 | −9.5 / +0.1 / −0.8 / −7.6 | −34 | −35 |
| calib8 dc1 leg | −0.23 | −2.9 / 4.2 | −3.0 / 4.1 | −6.0 / +0.7 / −1.5 / −4.8 | −36 | −38 |
| calib8 dc1 ts1 | −0.21 | −4.8 / 6.5 | −4.6 / 6.1 | −10.0 / −0.8 / −0.8 / −6.4 | −34 | −35 |
| calib8 dc1 s0 | −0.23 | −4.8 / 6.5 | −4.7 / 6.3 | −10.0 / −0.6 / −0.8 / −7.0 | −34 | −35 |
| prod dc1 leg | −0.14 | −1.5 / 2.9 | −1.4 / 2.4 | −2.5 / +1.5 / −0.8 / −3.6 | −15 | −33 |
| prod dc1 ts1 | −0.29 | −4.4 / 5.5 | −4.5 / 5.6 | −8.5 / −1.5 / −0.9 / −6.2 | −14 | −30 |

HPC:

| run | ALD bias (m) | T4cm bias / RMSE | T4cm DJF / MAM / JJA / SON bias | 1 m max / min bias | 2 m max / min bias | peak SWE bias (mm) | snow-off bias (d) |
|---|---|---|---|---|---|---|---|
| calib8 dc0 leg | −0.02 | −2.4 / 3.9 | −5.8 / +0.4 / −0.5 / −3.6 | −1.0 / −5.3 | −2.1 / −3.0 | −89 | −17 |
| calib8 dc0 ts1 | −0.31 | −3.4 / 4.8 | −7.5 / −1.2 / −0.6 / −4.2 | −2.2 / −6.8 | −3.8 / −4.4 | −87 | −5 |
| calib8 dc0 s0 | −0.31 | −3.4 / 4.9 | −7.5 / −1.2 / −0.6 / −4.2 | −2.2 / −6.8 | −3.8 / −4.4 | −87 | −5 |
| calib8 dc1 leg | +0.10 | −2.2 / 3.7 | −5.5 / +0.5 / −0.5 / −3.4 | −0.7 / −5.3 | −1.7 / −3.0 | −89 | −18 |
| calib8 dc1 ts1 | −0.16 | −3.4 / 4.8 | −7.7 / −1.7 / −0.6 / −3.4 | −1.8 / −7.9 | −3.7 / −5.4 | −87 | −5 |
| calib8 dc1 s0 | −0.21 | −3.2 / 4.6 | −7.4 / −1.5 / −0.6 / −3.4 | −1.7 / −8.0 | −3.8 / −5.1 | −88 | −5 |
| prod dc1 leg | −0.07 | −2.0 / 2.9 | −3.8 / −0.5 / −0.5 / −2.9 | −0.9 / −7.1 | −1.9 / −4.0 | −84 | −9 |
| prod dc1 ts1 | −0.31 | −3.9 / 5.3 | −8.2 / −3.2 / −0.5 / −3.9 | −3.3 / −7.5 | −4.9 / −4.9 | −83 | −5 |

Findings:
- Legacy is already too cold in winter at both sites (DJF −2.5…−6 °C), summer is close (JJA −0.5…−1.5 °C at 5–10 cm). The enthalpy solver removes the legacy spurious heat (§15b C3) and roughly doubles the cold-season bias (DJF −7.5…−10 °C, SON −3.4…−7.6 °C); JJA is unchanged or slightly better. Soil T RMSE rises from 2.4–4.3 to 4.6–6.5 °C. ALD bias worsens at HPC (−0.02…+0.10 → −0.16…−0.31 m) and is similar or slightly worse at TVC (−0.14…−0.23 → −0.21…−0.30 m).
- The cold-season error points at the snowpack, which the new solver no longer compensates for: modelled peak SWE is 14–36 mm low at TVC and 83–89 mm (55–59 %) low at HPC, and snow goes off 30–38 d early at TVC. The enthalpy solver's snowmelt-as-heat (D3) delays snow-off and improves it (HPC −17 → −5 d; TVC +3 d), the only metric that improves clearly.
- SFCC 0 and SFCC 1 differ by ≤ 0.15 °C in annual bias (≤ 0.7 °C in SON) and ≤ 0.05 m ALD: the bias change is the energy closure, not the freezing curve.
- Deep C (dc1) changes the comparison little with the new solver (TVC ALD +0.07 m, HPC +0.15 m vs dc0).
- Calibration implication (for Stage E): the likely candidates are the parameters controlling winter insulation and the snowpack (snow density `snwdennew`/`snwdenmax`, snowfall/precip partition or undercatch, moss/organic dry conductivity `org_tc_dry`, `org_tc_kfrz`), then ALD via `mineral_tcsolid`/organic thickness. Retuning only summer conductivity would not fix the winter bias.

### D2. Sub-step study (`thermal_substeps`, -O2 build, TVC and HPC Warming_nfeed, calib8 dc0, reference = 24 sub-steps)
Tool `staged/subtab.py`. Monthly soil T statistics over main years 1–40 (before the TVC late-run deepening); wall time per case (20-yr spin-up + 75-yr main, 1 s resolution, 2 concurrent jobs).

| run | wall s (TVC / HPC) | solves/day | ALD diff vs 24, yr 1–40 mean (TVC / HPC, m) | ALD rms yr 1–40 | mean T 1 m diff vs 24 (TVC / HPC, °C) | monthly T RMSE vs 24 at 20 cm (TVC / HPC) |
|---|---|---|---|---|---|---|
| legacy | 23 / 23 | – | +0.05 / +0.23 | 0.08 / 0.25 | +2.4 / +0.9 | 3.1 / 1.4 |
| 1 | 18 / 18 | 2.5 | −0.04 / −0.00 | 0.06 / 0.05 | +0.74 / −0.08 | 1.3 / 0.8 |
| 2 | 19 / 19 | 4.2–4.5 | −0.01 / +0.04 | 0.03 / 0.06 | −0.38 / −0.23 | 1.0 / 0.7 |
| 4 | 21 / 21 | 7.8–8.7 | −0.01 / +0.05 | 0.03 / 0.07 | +0.01 / −0.24 | 0.9 / 0.7 |
| 8 | 23 / 25 | 15–16 | +0.01 / +0.02 | 0.05 / 0.05 | +0.01 / −0.29 | 1.2 / 0.7 |
| 24 | 37 / 38 | 42–43 | 0 | 0 | 0 | 0 |

- The daily-step time error is small: ALD within 0.05 m of the 24-step reference for every sub-step count, comparable to the year-to-year scatter between runs (rms 0.03–0.07 m), and far smaller than the legacy − enthalpy difference. Monthly T RMSE does not converge below ≈ 0.7–1 °C with more sub-steps because 75 years of coupled snow/water/BGC feedback make the runs drift apart (chaotic, not discretisation, error).
- The one systematic signal is TVC with 1 step: 0.4–0.75 °C warmer at 0.2–1 m than the reference. 4 steps is within 0.07 °C of the reference at TVC (HPC: all step counts within 0.3 °C).
- Late in the TVC Warming run (years 50–75) ALD deepens to 0.9–1.3 m in the 2-, 8- and 24-step runs and less (0.5–0.9 m) in the 1- and 4-step runs. This is a threshold response of the warming case to small perturbations, not a convergence trend.
- Cost: +6 % (2), +17 % (4), +28–39 % (8), +106–111 % (24) wall time per case; with 4 sub-steps the -O2 model is still faster than legacy (21 vs 23 s).
- Recommendation: `thermal_substeps = 4` for calibration and production runs with `thermal_solver = 1`. Elchin approved this as the code default on 2026-11-10 (his words: “Yes, make 4 sub-steps the default”). An explicit config value still overrides it. See §15e.

### D3. Root cause of the integrator safety-net increase
Diagnostics with `soil_giveup_diag = 1` (v5n, default off) on HPC_Warming_nfeed (calib8 dc0, calib8 dc1, prod dc1), plus `soil_giveup_retry = 1` (v5n G2, default off). Tools `staged/gd.py`, `staged/gd2.py`.

| run | events (main) | pool / reason | layers (z top, median) | months | layer T median | deepest thawed layer (median index) | Σ root N uptake / Σ AVLN at month start (g N m⁻², median) | with `soil_giveup_retry = 1` |
|---|---|---|---|---|---|---|---|---|
| calib8 dc0 legacy | 0 | – | – | – | – | – | – | – |
| calib8 dc0 ts1 | 3 | AVLN negative | 3–4 (0.16 m) | Jun–Jul | 2.5 °C | 6 | 0.56 / 0.42 | 0 |
| calib8 dc1 legacy | 4 | AVLN negative | 5 (0.17 m) | Jul–Aug | 1.9 °C | 9 | 0.22 / 0.19 | – |
| calib8 dc1 ts1 | 60 | AVLN negative | 3–6 (0.17 m) | Jun–Jul | 2.4 °C | 7 | 0.30 / 0.26 | 0 |
| prod dc1 legacy | 113 | AVLN negative | 2–5 (0.17 m) | Jun–Aug | 1.9 °C | 8 | 0.37 / 0.31 | – |
| prod dc1 ts1 | 119 | AVLN negative | 3–5 (0.16 m) | Jun–Aug | 2.2 °C | 7 | 0.32 / 0.27 | 1 (unattributed) |

Root cause: every event, legacy and enthalpy, is the same N-cycle failure, not a thermal one. In early summer the root N uptake demand of a thawed 0.1–0.2 m organic layer exceeds that layer's available N at month start plus its net mineralisation, so AVLN is driven negative within the first 10–30 % of the month and the monthly soil integrator gives up and falls back to the safety net. Events sit at about 2 °C, not at the thaw front (0–3 % in the deepest thawed layer). The enthalpy solver makes this more frequent because its thaw is shallower and colder: the deepest thawed layer is 1–2 layers higher, so the same uptake is spread over fewer, colder layers with less mineralisation. That is why the increase is largest in calib8 dc1 HPC (deep C plus N feedback; legacy ALD 0.85 → 0.59 m). The G2 retry (`soil_giveup_retry = 1`) removes all events in calib8 and all but one in prod. The fix belongs in the N-uptake limiter (v5n G1 `n_extract_limiter`, or G2) and is a calibration-time decision, not a thermal-solver change. C and N closure pass in every run.

### Stage D scope not done (left for later or for Elchin)
- `thermal_k_picard` test, multi-front-year analysis and an HPC MPI tile timing (listed in §17 D) were not run; the approved Stage D scope was obs comparison, sub-steps and the safety-net root cause.
- HPC soil T obs below 5 cm are an envelope from a station 3.5 km away (01TC02); only monthly-mean extremes are compared.

---

## 15d. Sub-step anomaly (HPC_Warming_nfeed, 24 sub-steps): root cause and fix (2026-11-09, branch `enthalpy/substep-fix` off `stageD-eval`, local only)

**Finding (Stage 0, Calibration-expert).** With thermal_solver 1 on S23 C1 inputs, HPC_Warming_nfeed at 24 sub-steps/day jumped to ALD ~1.0 m from 1983 and up to 1.91 m (2010–24 mean 1.34 m), while 1 and 4 sub-steps stayed at 0.45–0.82 m. The worst full-day energy residual also grew with sub-steps (6e-8 → 1.7e-7 → 7.6e-7 W m⁻²).

**Verdict: a numerical/code artefact, not a physical tipping point.** There are two separate causes.

### A. ALD regime jump: unfrozen water of frozen cells deleted by the hydrology (bug, all sub-step counts)
- `Soil_Env::checkSoilLiquidWaterValidity` gives every `frozen == 1` layer effmaxliq = 0. With the new solver a fully frozen cell still holds SJ2016 unfrozen water (≈2 kg m⁻² per 0.1 m of mineral soil at −4 °C), and the routine dumps it to puddle/runoff.
- When no soil layer has `frozen == 0` (whole column frozen, `ststate = 1`, `drainl = NULL`), `last_active_layer` falls back to `lstsoill`, so the **whole column's** unfrozen water goes. That is about 118 kg m⁻² per day here (diag dump: liq → 0 in all 20 soil layers, P−G = 118 kg m⁻²). The routine is reached through the frozen-top branch's transpiration step, so this happens on whole-frozen days while vegetation transpires: spring cold snaps that fully refreeze a thin thawed top layer, and early autumn.
- The next thermal step puts each cell back on its curve at fixed H. It re-makes the liquid from ice, and the latent heat comes out of sensible heat. In 1981, DOY 150–153 and 262–265, the same thing happened on 8 days in a row. The mineral layers went from −5 to −24 °C **in June**, and the column lost its unfrozen water and part of its ice every day.
- Both budgets stay closed, so no gate caught it. Water leaves as runoff, and its enthalpy is booked in `bwater` (+3…+11 MJ m⁻² per day).
- Mineral VWC in layers 5–11 (March) fell from 0.24 (1975–80) to **0.05** (1985–2024). The dry, low-latent-heat column then thaws deep: ALD went 0.27 (1982) → 0.66 (1983) → 1.10 (1988).
- **Why 24 sub-steps:** whole-frozen transpiration days with a dump happen on 13 days at 24 sub-steps, against 1–3 at 1/4/8/12. Finer steps let a cold day refreeze the thin thawed top completely, which a daily backward-Euler step smears out. After the first dump the colder column makes the next day fully frozen too, a positive feedback. 48 sub-steps also diverge (8 dump days, ALD up to 1.02 m). So it is a threshold in a buggy path, not convergence.
- **The bug also affects 1 and 4 sub-steps:** s4 lost mineral water in 1975 (VWC 0.30 → 0.23). The same rule also dumps the liquid of `frozen == 1` layers above the first partly frozen layer (freeze-up from the top).
- **Fix (A):** with `thermal_solver == 1`, a `frozen == 1` layer is inert in that routine (effmin = effmax = its liq): no excess and no deficit. The legacy branch is untouched.

### B. Closure growing with sub-steps: the per-solve tolerance accumulates (numerics, not a bug)
- `opt.tol` = 1e-3 J m⁻² per cell is applied to every backward-Euler step, so the full-day residual bound grows about ×nsub.
- After fix A the worst day reached 1.6e-6 W m⁻² at 24 and 4.5e-6 at 48, **over the 1e-6 gate**.
- **Fix (B):** tol is a per-day budget, tol/nsub per step (restored after the loop; nsub = 1 is unchanged). Cost: +10 % solves at 24 sub-steps (52 vs 47 per day).

### Before/after: HPC_Warming_nfeed, S23 C1 inputs, calib-elchin -O0, 20-yr spin + 75-yr main (main 1950–2024)
| run | ALD 1980–89 | ALD 1990–2009 | ALD 2010–24 | max ALD | worst full-day residual (W m⁻²) | dump days | solves/day |
|---|---|---|---|---|---|---|---|
| stageD-eval s1 | 0.40 | 0.57 | 0.72 | 0.82 | 4.0e-08 | 1 | 3 |
| stageD-eval s4 | 0.38 | 0.55 | 0.70 | 0.82 | 1.1e-07 | 2 | 9 |
| stageD-eval s8 | 0.39 | 0.54 | 0.70 | 0.82 | 1.6e-07 | 3 | 18 |
| stageD-eval s12 | 0.40 | 0.56 | 0.70 | 0.82 | 4.2e-07 | 3 | 25 |
| **stageD-eval s24** | **0.65** | **1.16** | **1.34** | **1.91** | 5.7e-07 | **13** | 43 |
| stageD-eval s48 | 0.40 | 0.66 | 0.79 | 1.02 | **1.8e-06** | 8 | 79 |
| fix A only s24 | 0.38 | 0.51 | 0.62 | 0.80 | **1.6e-06** | 1 | 47 |
| fix A only s48 | 0.38 | 0.50 | 0.63 | 0.77 | **4.5e-06** | 1 | 87 |
| **fix A+B s1** | 0.44 | 0.52 | 0.66 | 0.82 | 4.7e-08 | 0 | 3 |
| **fix A+B s4** | 0.38 | 0.51 | 0.64 | 0.79 | 4.1e-08 | 1 | 10 |
| **fix A+B s8** | 0.39 | 0.52 | 0.63 | 0.80 | 3.1e-08 | 1 | 18 |
| **fix A+B s12** | 0.37 | 0.51 | 0.63 | 0.77 | 2.9e-08 | 1 | 27 |
| **fix A+B s24** | 0.38 | 0.51 | 0.61 | 0.76 | 2.0e-08 | 1 | 52 |
| **fix A+B s48** | 0.38 | 0.50 | 0.65 | 0.81 | 1.3e-08 | 1 | 103 |

How to read the table:
- "Dump days" are main-stage days with no front and `bwater` > 1e5 J m⁻². The one remaining fixed-run day (1978 DOY 145, 0.1 MJ m⁻²) is ordinary infiltration.
- The stageD-eval binary on the VM reproduces the box Stage 0 numbers (s4 / s24 means identical).
- The 12/48-allowing build (`base12`) is value-identical to stageD-eval at s24 (278 variables).
- Mineral VWC (layers 5–11, March 1985–2024): s4 0.233 and s24 0.047 before the fix, 0.274 after at both s4 and s24.
- 1981 TMINEB on DOY 165/270: s24 −21.2/−25.1 °C before, −3.9/−1.9 °C after.
- After the fix ALD converges across 1–48 sub-steps (2010–24: 0.61–0.66 m; s4 vs s24 within 0.03 m), and the residual falls with sub-steps instead of rising.

### Effect at 1 sub-step and on the legacy path
- `thermal_solver = 0`: run_regression.sh (ref = stageD-eval -O0, `--no-repeat`) gave A==B bit-identical for all 8 cases in each of dc0_calib8, dc1_prod and prod_t125_s23 (24/24).
  - dc0_calib8 and dc1_prod PASS.
  - prod_t125_s23 reports FAIL only for "N closure outputs missing" (the production output spec, the same in A and B, as in Stage 0).
- `thermal_solver = 1`, 1 sub-step: **not unchanged.** Fix B is a no-op at nsub = 1 (fixA s1 ≡ fixAB s1, 278 variables), but fix A changes every case from the spin-up on, because frozen-layer water is no longer lost.
  - ts1_dc0_calib8, last-30-yr ALD, REF → NEW: HPC NoW nfeed 0.43 → 0.31, NoW nonfeed 0.44 → 0.31, W nfeed 0.50 → 0.45, W nonfeed 0.51 → 0.45; TVC −0.01…−0.07 m.
  - S23 HPC_W_nfeed s1 2010–24: 0.72 → 0.66.
  - All 8 calib8 cases close: full-day energy ≤ 3.8e-8 W m⁻², water −0.26…+0.04 % P, C max annual ≤ 0.0034 g m⁻², 0 non-converged steps.
  - The wetter, more latent column gives shallower ALD, so the Stage 0 / Stage D thermal_solver 1 numbers (all sub-step counts) should be re-run on this branch before calibration.
- Unit tests: 307/307 test cases. The 6 failed assertions are the known xfails. There are 3 new cases in `tests/unit/test_enthalpy_substep_fix.cpp`:
  - (A) whole-column freeze keeps the curve water, and the next step needs no re-equilibration;
  - (A) the legacy rule is unchanged with thermal_solver 0;
  - (B) the tolerance is restored after the step, and the residual is ≤ cells × tol at 1/4/12/24/48 sub-steps.
- `thermal_substeps` now also accepts 12 and 48 (validation only).

## 15e. Default `thermal_substeps` is 4 (Elchin, 2026-11-10)

Elchin approved making 4 sub-steps the default (his words, 2026-11-10: “Yes, make 4 sub-steps the default”), following the Stage D / post-fix recommendation in §15c D2 and the rerun on this branch.

- When `thermal_solver = 1` and `thermal_substeps` is absent, the code now uses 4 (was 1). An explicit value in the config still overrides it (1, 2, 4, 8, 12, 24 or 48).
- When `thermal_solver = 0`, the key must stay absent or 1 (same FATAL rule as before). The constructor / unconfigured `ModelData` default remains 1, because that object is the legacy solver.
- The shipped `config/config.js` template documents the default in a comment. The active template leaves `thermal_solver` unset (legacy 0), so an active `"thermal_substeps": 4` would be FATAL.
- Checks (this change, local only, nothing pushed): `thermal_solver = 0` output on `dc0_calib8` is md5 / bit-identical to `7c5e4759`; `thermal_solver = 1` with the key absent is bit-identical to an explicit `thermal_substeps = 4` run on that old binary; an explicit 1 still matches the old 1-step result. Unit tests include a default-is-4 case; the known xfail assertions are unchanged.

## 15f. Combine default-4 sub-steps with `richards_dt_fix` (Elchin, 2026-11-10)

Elchin approved turning the Richards step-halving fix on for all new-solver work and Calibration-expert's Stage E baseline, and keeping production legacy off until recalibration (his words, 2026-11-10: “Yes: richards_dt_fix on for new-solver work and Stage E baseline, production stays off until recalibration”).

- Local branch `enthalpy/ssf-rdt` (worktree `~/tem-dev/srd` on calib-elchin) combines `enthalpy/substep-fix` @ `f8617974` (default 4 sub-steps, §15e) with `fix/richards-dt` @ `c12218cd` (switch, fix, counters, tests). The fix commits were cherry-picked onto `f8617974`: `bd652441` → `582cfc3b` (tests) and `c12218cd` → `8ffc55ed` (switch/fix). `include/ModelData.h` and `src/ModelData.cpp` auto-merged with no conflicts (the default-4 and `richards_dt_fix` edits are in different regions). Follow-up `05ef374c` writes `richards_dt_fix=1` in new-solver case-set tooling (`mk_variant.py`, `extra.sh`, `stageD_vm.sh`, commented example in `config/config.js`). Nothing pushed. The **code default of `richards_dt_fix` stays 0**.
- **Stage E baseline:** branch `enthalpy/ssf-rdt`, `thermal_solver = 1`, `thermal_substeps` default 4 (key absent), `thermal_sfcc = 1`, `richards_dt_fix = 1` set explicitly in the case-set config (not the code default).
- **Production legacy** (`thermal_solver = 0`) keeps `richards_dt_fix = 0` until the merged branch is recalibrated. The shipped `config/config.js` template leaves the key unset and only documents a commented new-solver example.
- Unit tests on `8ffc55ed` (−g, 4 cores): **322/322 test cases passed**. 8 expected failed assertions, marked xfail / `may_fail`: the same 6 pre-existing (climate 2, SEB albedo 2, Stefan `combineExtraFronts`, Richards early-return) plus the 2 documented cases in `test_richards_water_energy.cpp` (bug R1 at `dt_fix = 0`, and the frozen-lens case). `test_richards_dtfix` passed. The R1 `may_fail` case still documents the default-off behaviour.
- Bit-identity on `prod_t125_s23` (all 8 cases, 20-yr spin-up + 75-yr main, n = 75). Data-level comparison of every spin-up and main output (`.nc` variables and csv); `Git_SHA` and `run_status.nc:total_runtime` ignored. 32 runs, all OK.
  - **(a)** `thermal_solver = 0`, `richards_dt_fix` absent (= 0) vs the production legacy baseline (`richards_dtfix/runs/leg_f0`, already identical to `7c5e4759`): **8/8 identical**.
  - **(b)** `thermal_solver = 1`, `thermal_substeps` unset (default 4), `richards_dt_fix` absent (= 0) vs `f8617974` default-4 on the same configs: **8/8 identical**. The same outputs are also identical to explicit `thermal_substeps = 4`, `richards_dt_fix = 0` (`richards_dtfix/runs/ts1s4_f0`).
  - **(c)** `thermal_solver = 1`, `thermal_substeps = 4`, `richards_dt_fix = 1` vs the `c12218cd` fix-on results (`richards_dtfix/runs/ts1s4_f1`): **8/8 identical**.
- Closure on the (c) runs (hence also on `c12218cd` s4 fix-on): full-day energy residual ≤ 4.6×10⁻⁸ W m⁻² (well under the 1×10⁻⁶ gate). Water residual as a fraction of precipitation is −0.001…+0.059 % P (largest +0.059 % at HPC_Warming_nfeed); several cases are ~0.001 %. That leftover is the known non-Richards term, not bug R1. With `dt_fix = 0` the new-solver water residual is −0.10…−0.53 % P. The harness `closure_pass` is still false on this case set only because “N closure outputs missing” (the production output spec; same in legacy).

## 15g. Melt-out closure fix: heat of a melted snow cell lost in a zero-mass snow cell (Elchin, 2026-11-10, branch `fix/meltout-closure`)

Symptom: Calibration-expert's r1 case (TVC_Warming_nfeed, CMT51 `snwdennew` 75, `nfactor(w)` 0.3, `snow_scheme` 0, Stage E baseline
`thermal_solver` 1 / sfcc 1 / 4 sub-steps / `richards_dt_fix` 1) missed the daily gate on 1 Jun 2024: res = resday = -57.36 J m-2
(6.6e-4 W m-2), 0 fallbacks, no warning. Entering that day the pack was a remnant of a 27 May snowfall: top cell 0.16 mm,
0.009 kg m-2 at -0.01 C, and a second snow cell with **no water at all** (W = 0, dz = 0, kept by the snow module).

Root cause (instrumented, not guessed): the solve itself closed (sum of the four sub-step closures 1e-7 J m-2) and the scatter was exact.
The whole -57.36 J m-2 was in step 5, "heat of fully melted snow cells goes to the cell below": the top cell melted out and ended at
+1.52 C (H = 57.36 J m-2); that heat was given to the next cell, the empty snow cell. Its heat capacity is zero (Cs = 0, W = 0), so
`equilibrate()` -> `invertH()` returns T = 0 and the heat vanished.

Fix (Soil_Env_enthalpy.cpp, step 5; thermal_solver 1 only, the legacy path never runs this code): the receiver is the next cell below that
can hold heat (Cs > 0 or liq + ice > 0); zero-capacity cells are skipped. When the next cell has heat capacity (every other case) the
arithmetic is the same as before, so results only change on days where heat was being lost.

Tests: `tests/unit/test_enthalpy_meltout.cpp` (r1 1 Jun 2024 state on a small column, a cascade over two empty cells, an empty cell
between two snow cells, and the unchanged path); three of the four fail without the fix (res -57.9, -419, -27 J m-2).
Regression: `tests/regression/check_meltout_r1.sh` re-runs r1 and applies the 1e-6 W m-2 gate (closure.py).
Open (not changed here): the snow module keeps zero-mass snow layers; removing them belongs in Snow_Env and would also change legacy (done behind a switch in §15h).


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
| **C: integration behind switch** | ModelData keys (incl. `geothermal_flux_mode` zero/value/map, D4), dispatch in `updateDailyGroundT`, gather/scatter, snow hand-off, SEB melt source, remap hooks in Ground/Snow_Env/updateDailySM/fire, front/ALD reconstruction, restart attribute, extended energy_diag, harness `--etol-enth`. | `thermal_solver = 0` md5-identical. `= 1` passes all gates on calib8 and prod (SFCC 1 production, SFCC 0 verification), **including full-day balance ≤ 1e-6 W m⁻² as pass/fail (D5)**. Restart round-trip in both directions. | **12–16 days**: integration done 2026-11-09 (branch `enthalpy/stageC`, tag `stageC-integration`, §15b); C2 remap + restart attribute done (tag `stageC2`, §15c) |
| **D: site evaluation (TVC, HPC)** | 1950–2024 runs, legacy vs enthalpy (SFCC 0/1, substeps 1/4, k_picard 0/1). Compare with observed ALD, soil T profiles and snow-off dates. Multi-front-year analysis. Timing on a ≤ 75-cell 4×4 km tile with MPI on the HPC. | Closure ≤ 1e-6 W m⁻². ALD/soil T bias vs obs at least as good as legacy, or explained. Measured core-hours per century per 10,000 cells (full-domain budget only when Elchin asks, D6). | **6–8 days** (+ compute): obs comparison, sub-steps and safety-net root cause done 2026-11-09 (tag `stageD-eval`, §15c) |
| **E: calibration plan** | Which parameters move (n-factors, organic tc, `tcsolid`, snow density, possibly SJ exponents); targets; combine with the deep C recalibration; config for production 4×4 km. | Plan approved by Elchin. Recalibration itself shared with deep C (15–25 person-days in the merge plan). | **2 days** (plan): assigned to Calibration-expert (Elchin 2026-11-09); handoff `/workspace/tem-dev/stageE_handoff.md` |

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
