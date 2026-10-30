# daily-climate-mode: references for every switch and fix (v5-daily … v5q)

Standing rule (Elchin, 2026-10-22): every code change carries its scientific references in three places: code comments, commit message and README. This file is the full bibliography. Code comments cite entries by `[Key]`.

Every DOI was checked against Crossref or DataCite, and every URL was opened (HTTP 200), on 2026-10-22 (v5m, v5n and v5o entries: Crossref, 2026-10-27; v5p entries: Crossref, 2026-10-28). Author, year, title and journal were compared with the registry record.

Where a fix is a code or numerics bug with no literature basis, it says **no published source**, followed by the code location, the original commit and the principle or standard it restores.

Switch defaults: every switch below is 0/off by default, which is legacy behaviour and bit-identical.

## 1. Switch → references

### Forcing and climate

| switch (model_settings) | commit(s) | what it does | references |
|---|---|---|---|
| `climate_timestep` = "daily" (daily forcing mode) | 76bc0858, 3a7cc57f, bb93da53, 0e22205a, 17d53b56 | SC stage reads a daily (365_day) projected-climate file directly instead of interpolating monthly values with the stochastic precipitation split. | **No published source for the mode itself.** Rationale: drive the model with the observed daily meteorology rather than synthetic disaggregation of monthly means. The calendar follows the CF `365_day`/`noleap` calendar [CF18]. Daily extraterrestrial radiation is FAO-56 eq. 21 [Allen1998]. Model context: [Yi2009], [Genet2013]. |
| `vp_hpa_fix` | 6c78bd4c | The vapor_press input is in hPa, but svp is in Pa. 1 = convert to Pa before VPD and RH. | Input units: CRU TS `vap` is in hPa [Harris2020]. VPD = e_s − e_a in consistent units: FAO-56 eqs. 11–19 [Allen1998]. The VPD enters Penman–Monteith [Monteith1965]. Principle restored: dimensional consistency. |
| `svp_ice_fix` (B1) | 07f6c8d0 | t ≤ 0 °C: legacy 6.112·exp(17.27 t/(272.62+t)) mixes the Tetens water coefficient with the ice denominator. 1 = **WMO (2008) ice formula e_i = 6.112·exp(22.46 t/(272.62+t)) hPa** (valid −65…0 °C), also used in the SEB esat. | [WMO2008] (CIMO Guide, Annex 4.B; formula reproduced in [EOL-VP]). Close to [Buck1981] (22.452, 272.55). Reference values for the tests: [MurphyKoop2005]. The legacy water branch uses 17.63 where WMO-2008 has 17.62. That predates this work, is not changed, and amounts to about 0.06 % at 20 °C. |

### Soil thermal and n-factor

| switch | commit(s) | what it does | references |
|---|---|---|---|
| `nfactor_sameday` | 9b3a4840 (+ seam d718934d) | Picks the summer or winter n-factor from today's Tair. Legacy used the previous day's value, which was uninitialised on day 1 of each month. | n-factor definition and seasonal (thaw/freeze) application: [Klene2001]. The uninitialised selector is a code bug (Cohort::updateMonthly_Env, upstream 419a4eb3). Principle restored: no uninitialised state; forcing taken from the same day. |
| `kersten_form` | 4b1041fd | Unfrozen Kersten number: 1 = Ke = log10(Sr)+1 instead of ln(Sr)+1. | [Johansen1975] (fine-grained soils). Same form in CLM ([CLM5-SoilT] §2.6.3, unfrozen K_e) and [LawrenceSlater2008]. |
| `org_tc_scheme` (+ `org_tc_*`) | 79668759 | Moisture-dependent conductivity of moss and organic layers: Ke = κSr/(1+(κ−1)Sr), with Sr relative to the peat retention capacity. | [CoteKonrad2005]; peat retention capacity θ_re: [Letts2000]. Linear dry–saturated moss form as in JULES: [Chadburn2015]. |
| `mineral_tcsolid`, `mineral_vhcsolid`, `mineral_tcdry_form` | 3e225cec | Settable solid conductivity and heat capacity. Dry conductivity (0.135ρ_b+64.7)/(2700−0.947ρ_b). | [Johansen1975]; [Farouki1981]; CLM [CLM5-SoilT] eqs. 2.6.80 and 2.6.82 (and the solid heat capacity in §2.6.3); [LawrenceSlater2008]. |
| `seb_snow_tc` | d41c7f79 | Snow conductivity: 0 = Goodrich (1982) k = 2.9e-6 ρ²; 1 = Sturm et al. (1997) k(ρ) for ρ ≤ 600 kg m⁻³. | [Sturm1997]; legacy [Goodrich1982]. |

### Surface energy balance (`seb_scheme`, d41c7f79)

| switch / formulation | what it does | references |
|---|---|---|
| `seb_scheme` 0/1/2 | Daily skin-temperature energy balance under snow (2) or always (1, experimental). It solves SW_abs + LW_in − εσT_s⁴ − (1−ε)LW_in − H − LE − G = 0 for T_s; the surplus at T_s = 0 °C becomes melt. | General structure of a land-surface skin energy balance: [CLM5-Fluxes]; [Brutsaert1982]. |
| bulk transfer | C_HN = k²/[ln(z_u/z_0)·ln(z_T/z_0h)], with z_0h = z_0/10. H = ρc_p C_H U (T_s − T_a). Snow LE = L_s ρ C_H U (0.622/p)(e_sat(T_s) − e_a). | [Brutsaert1982] (bulk transfer formulation; z_0/z_0h = 10, i.e. kB⁻¹ = ln 10 ≈ 2.3, is a common simplifying choice, not a site-fitted value). e_sat from [WMO2008]. |
| stability | Stable: C_H·1/(1+b·Ri_b) with b = `seb_stab_b` (10). Unstable: ·min(3, (1−16·Ri_b)^½). | **Simplified forms; no paper uses exactly these expressions.** Stable: Louis-type bulk-Richardson-number reduction [Louis1979]. Unstable: Businger–Dyer φ_h = (1−16ζ)^−½ [Dyer1974], applied with Ri_b in place of ζ. |
| `seb_emis_snow` 0.99 | Snow emissivity. | [Warren1982] (snow is nearly black in the thermal infrared). |
| `seb_lw_source` | 0 = forcing `lwdown`; 1 = parameterised. Clear sky: ε_clr = 1 − (1+ξ)exp(−(1.2+3ξ)^½), ξ = 46.5 e/T [Prata1996]. Cloud: ε = c + (1−c)ε_clr, with c = 1 − SW/SW_clear [CrawfordDuchon1999]. | [Prata1996]; [CrawfordDuchon1999]. |
| skin-T solver | Bracketed Newton iteration with a bisection safeguard. | [Press2007] (safeguarded Newton, §9.4 "rtsafe"). |
| ground heat flux | G = λ(T_s − T_1), with λ = 1/(1/100 + 0.5 dz_1/k_1), the conductance of the top half-layer. | Standard finite-volume conduction as in [CLM5-SoilT]. The 100 W m⁻² K⁻¹ virtual-layer conductance has **no published source** (numerical coupling choice). |

### Soil hydrology

| switch | commit(s) | what it does | references |
|---|---|---|---|
| `frozen_infil_eice` (draft) | 939790f5 | Infiltration into frozen soil with ice impedance 10^(−E·F_ice)·k_sat (E = 6 in CLM5). | [Swenson2012]; [CLM5-Hydro] eq. 2.7.48. |
| `soil_evap_beta` (draft) | 8053cabc | β = 0.25[1 − cos(π θ/θ_fc)]² on soil evaporation; θ_fc at ψ = −3.4 m. | [LeePielke1992]; θ(ψ) from [ClappHornberger1978]. |
| `org_drain_rate` (+ `org_drain_*`) | 014dda0d | Lateral "fill and spill" drainage of thawed near-surface organic layers toward θ_re. | [QuintonMarsh1999]; θ_re from [Letts2000]; compare the CLM5 perched-water-table drainage [CLM5-Hydro] §2.7.4. |
| `richards_infil_fix` (B2) | 07f6c8d0 | When Richards::update returns early, the infiltration went to a private EnvData copy and was lost. 1 = puddle/runoff. | **No published source.** Code: Richards::update early return / Soil_Env::rebookUnappliedRichards; the bug is in upstream 419a4eb3. Principle restored: water mass conservation (continuity in Richards' equation [Richards1931]). |
| `nonvasc_trans_withdraw` | 9dae595b | Non-vascular transpiration is withdrawn from the first shallow organic layer; uptake that is not applied is removed from reported ET. | **No published source** for the scheme. Rationale: water mass conservation (reported ET must leave storage). Context: moss/understory water-vapour exchange is a measured ground-layer flux [Heijmans2004]. |

### Snow and soil fronts

| switch | commit(s) | what it does | references |
|---|---|---|---|
| `stefan_fronts_fix` (B3) | 07f6c8d0 | combineExtraFronts set `mindz` only once, so it hung with more than one surplus pair. 1 = reset on every pass. | **No published source** for the fix. Code: Stefan::combineExtraFronts (upstream 419a4eb3). Principle restored: the loop terminates (an invariant re-established each iteration). Algorithm context: the two-directional Stefan freeze/thaw algorithm [Woo2004], [Yi2009]. |
| `snow_albedo_fix` (B4) | 07f6c8d0 | Below −10 °C getAlbedoVis returned a hard-coded 0.9 (above snwalbmax, with a jump at −10 °C). 1 = snwalbmax. | **No published source.** Code: Snow_Env::getAlbedoVis (upstream 419a4eb3). Principle restored: consistency with the parameter-file maximum and continuity at −10 °C. Note: fresh-snow visible albedo can be ~0.9–0.98 [Warren1982], so this is a consistency fix, not a physics claim. |

### Soil biogeochemistry

| switch | commit(s) | what it does | references |
|---|---|---|---|
| `soil_nmin_tol_fix` (C1) | d427d649 | The error test of the signed NMIN accumulator used ptol·\|NMIN\|. 1 = ptol·max(\|NMIN\|, \|NIMMOB\|, \|AVLN\|). | **No published source for the fix.** Code: Integrator::boundcon (upstream 419a4eb3). Standard restored: adaptive step control should scale the error by a magnitude that cannot pass through zero [Press2007] §17.2; [HairerNorsettWanner1993] §II.4. |
| `n_rootless_fix` (N1; default 2 since v5o, 5812ff19) | 29d72af8 | N mineralised in rootless layers left ORGN without reaching AVLN. | **No published source.** Code: Soil_Bgc::deltastate (upstream). Principle restored: N mass conservation. |
| `n_som_move_fix` (N2; default 1 since v5o, 5812ff19) | 20b2723e | SOM N moved down the profile was overwritten ('=' instead of '+='), shared wrongly among receivers, and taken from cumulative donor C. | **No published source.** Code: Soil_Bgc::deltastate (upstream). Principle restored: N moves with C (fixed C:N of moved SOM), N mass conservation. |
| `n_nonvasc_extract` (N3; default 1 since v5o, 5812ff19) | 20b2723e, d3d3fbba | Non-vascular N uptake had no soil extraction. 1 = extract from moss/fibric layers, capped at 95 % of available N. | Mass conservation (**no published source** for the scheme). Ecological basis that bryophytes take N from the moss/peat layer and deposition: [Turetsky2003]. |
| `n_safetynet_conserve` (L1, v5m) | 23e8077b (merge bf05797e) | When adapt() gives up, the safety net (be8a62fa) delivered the remaining inputs but not the remaining (1 − t0) root N extraction, which the vegetation had already taken up: N was created at every event. 1 = debit it from soil ORGN (rooted layers by extraction share, then any layer), then AVLN, carry any rest to the next month; deposit the remaining N input where deltastate does. C unchanged (no C leak). | **No published source for the fix.** Code: Integrator::applySoilSafetyNet (be8a62fa). Principle restored: closed N mass balance of the TEM N cycle, with plant N uptake taken from the soil N pools [McGuire1992]; [Raich1991]; multi-PFT N competition in DVM-TEM [Euskirchen2009]; DOS-TEM soil N [Yuan2012]. Organic N as first donor: arctic/boreal plants take up organic N [Nasholm1998]; [SchimelBennett2004]. Numerics: resetting/clipping an ODE state is not mass-conserving unless the removed mass is booked [Shampine2005]. Why closure matters for spin-up: [ThorntonRosenbloom2005]. |
| `n_extract_limiter` (G1, v5n) + `n_extract_tau` | 904a4a7c, 5c597aec (merge 4bb7cd99) | Root N extraction is a zero-order sink on AVLN during the soil month; AVLN in rooted layers reaches 0 mid-month and adapt() gives up. 1 = AVLN part r·min(1, AVLN/(τ·r)), rest from ORGN of the same layer; total extraction unchanged. τ default 0.2 (month fraction), read only when G1 (or a G2 retry month) is on. | Saturating plant N uptake in TEM [Raich1991]; [McGuire1992]. Uptake down-regulated when mineral N cannot supply demand [Thornton2007]. Organic N uptake [Nasholm1998]; [SchimelBennett2004]. Numerics: a zero-order sink drives a non-negative state below zero in finite time; positivity requires a state-dependent sink [Shampine2005]; [Burchard2003]. |
| `soil_giveup_retry` (G2, v5n) | 904a4a7c (merge 4bb7cd99) | Only when adapt() gives up: re-integrate the month from its start with the G1 limiter. | as G1. |
| `n_netnmin_rescale_guard` (G3, v5n) | 5c597aec (merge 4bb7cd99) | deltan rescales net N mineralisation by nscale = nminadj/Σnetnmin when avlnlost < 0; a negative nscale (net immobilisation, positive adjustment) flips every layer's sign. 1 = skip a negative nscale. | **No published source for the fix.** Code: Soil_Bgc::deltan (upstream). Principle: a scale applied to a signed flux must be ≥ 0. Net mineralisation as gross mineralisation minus immobilisation [SchimelBennett2004]; immobilisation limited by mineral N [Thornton2007]. |
| `soil_giveup_diag` (v5n) | 904a4a7c, 5c597aec (merge 4bb7cd99) | Diagnostics only (give-up cause, pool, layer, budget terms; csv + stderr). | none needed. |

**v5o (5812ff19): N1-N3 on by default.** The fixes are corrections of code bugs; there is **no published source** for them. The principle they restore is a closed TEM N cycle in which plant N uptake is debited from soil inorganic N and no N is created or destroyed outside the booked fluxes: McGuire et al. 1992 [doi:10.1029/92GB00219](https://doi.org/10.1029/92GB00219); Raich et al. 1991 [doi:10.2307/1941899](https://doi.org/10.2307/1941899). N3's basis (bryophytes take up N from the moss/peat layer): Turetsky 2003 [doi:10.1639/05](https://doi.org/10.1639/05). Legacy values (0) stay selectable and reproduce v5n bit-for-bit. The N closure check in tests/regression/closure.py is enforced only with all three on.

### ET partition (v5l, 798388f5)

| switch | what it does | references |
|---|---|---|
| `et_fpc_once_fix` (S1) | Apply the foliar projected cover once in transpiration: radiation and PAR per covered area, then Cohort's single tran·FPC. | Original double scaling: commit 28167031 ("Water fix 4/4", Cohort.cpp `tran *= fpc` on top of the FPC-scaled radiation in Vegetation_Env::updateRadiation). Per-area fluxes scaled by cover **once**: LPJ ET = min(S, D)·f_V [Gerten2004]; LPJmL4 supply per FPC [Schaphoff2018]; CLM5 PFT-area weighting of transpiration [CLM5-Hydro] eq. 2.7.68. FPC = 1 − exp(−k·LAI) (Beer's law) [MonsiSaeki1953]. |
| `nv_snow_bury_fix` (S2), `nv_snow_bury_depth` | No snow interception by moss/lichen under the snowpack (default depth 0 m: the moss layer is the top of the soil column). | Snow interception is a canopy process above the snowpack [HedstromPomeroy1998]. Low vegetation buried by snow no longer interacts with snowfall [Sturm2001]. The moss layer as part of the soil column in DOS-TEM: [Yi2009]. |
| `rain_int_lai_fix` (S3) | Capacity = ircoef·LAI (mm per ground area), storage carried between days and evaporated first; legacy 0.21·P + 0.5 mm. | p_liq = 0.1 mm per unit (L+S) (CLM5 eq. 2.7.12, after [Dickinson1993]) [CLM5-Hydro]. Forest/crop S_max = 0.935 + 0.498·LAI − 0.00575·LAI² [vonHoyningenHuene1981] as implemented in LISFLOOD [vanderKnijff2010]. Observed TVC interception: 16–29 % of summer rain under dwarf birch [Zwieback2019]; 16.6 ± 3.2 % watershed-wide [Chang2022]. |
| `nv_rain_surface` (S4) | Moss/lichen rain "interception" is surface water: 1 = booked as EVAP_MOSS; 2 = reaches the ground. | Moss treated as a ground/soil layer: [Chadburn2015] (JULES), [Yi2009] (DOS-TEM). Moss/understory evaporation as a ground-layer flux: [Heijmans2004]. |

### Output-only diagnostics (v5p, branch v5p-daily-split)

| switch | what it does | references |
|---|---|---|
| `daily_c_split` | Output only, no feedback: each monthly C flux (per PFT, per soil layer) is split over the days of its calendar month with the model's own environmental scalar of that day as weight, X_d = M·w_d/Σw (GPP: the day-varying part of getGPP; RM: Raq10(Ta); RG: max(GPP_d−RM_d, 0); RH: rhmoist·rhq10 of the layer; NEE_d = RECO_d − GPP_d), giving `*_DS` daily outputs whose calendar-month sums equal the monthly values (no smoothing across months). All-zero weights → uniform, flagged in `DS_FALLBACK`. A statistical split, not daily physics. | **No published source for this exact method.** Spreading coarse model fluxes over finer time steps with driver-based weights: [OlsenRanderson2004] (3-hourly fluxes from monthly CASA fluxes), [Fisher2016] (3-hourly downscaling of monthly MsTMIP NEE). Exact preservation of the coarse (monthly) total, temporal disaggregation — cited for the principle only, their regression/quadratic-minimisation methods are not used: [Denton1971], [ChowLin1971]. |

### Biological N fixation (v5q, branch v5q-bnf)

| switch | what it does | references |
|---|---|---|
| `bnf_scheme` (0 off / 1 PRO / 2 AET / 3 CONST) + `bnf_pro_c`, `bnf_pro_d`, `bnf_aet_a`, `bnf_aet_b`, `bnf_const` | Annual BNF from the previous model year into AVLN of the top soil layer (`NINPUT`), monthly by the previous-year NPP shape; PRO B = 1.8·(1 − exp(−0.003·NPP)), AET B = max(0, 0.00234·ET − 0.0172) with ET = transpiration + canopy evaporation; needs `n_rootless_fix` = 2; output `BNF`. | PRO Eq. A2, AET Eq. A1, parameters Table A1 of [Meyerholt2016]; regressions from [Cleveland1999]; BNF into soil mineral N as in CLM-CN [Thornton2007]; caveat: no significant BNF–ET/NPP relation in field data, free-living ≥ 1/3 of the total [DaviesBarnard2020]. |
| `bnf_asym` | Constant free-living term (no climate modulation). | [Meyerholt2016] (asymbiotic BNF to soil); moss-associated cyanobacterial fixation: [DeLuca2002], [Stewart2011], [Rousk2013]. |
| `bnf_avln_cap` | Linear taper max(0, 1 − AVLN/cap). | **No published source for the linear form**; FOR scheme converges to 0 above 2 gN/m2 inorganic N [Meyerholt2016]. |
| `bnf_to_plant` (Phase 2) | 1 = scheme (symbiotic) part to plant labile N of active PFTs, split by previous-year PFT NPP; `bnf_asym` stays in the soil; output `BNF_VEG`, extra input in `closure.py`. | Symbiotic BNF to the plant, asymbiotic to the soil [Meyerholt2016]; **no published source for the NPP-proportional PFT split**; rationale: delivers fixed N where the previous year's growth was, without a per-PFT fixer parameter. |

### N leaching bookkeeping (v5q-ndrain, branch v5q-ndrain)

| switch | what it does | references |
|---|---|---|
| `n_leach_fix` (0 legacy / 1 fix) | 1 = the soil AVLN pools lose exactly the booked inorganic N leaching `avlnlost` (`NLOST`): overland part from layer 0, drainage part from the drainage-zone layers by their current AVLN. Legacy `deltastate` removed avlnlost·qdrain (qdrain in mm/month), not the overland part, and nothing when totdzliq ≤ 0.01 (same code in upstream dvm-dos-tem c0fa75cb). | **No published source** (code bug; N mass conservation). Observed arctic DIN export for plausibility of the resulting leaching: [McNamara2008], [Khosh2017]. |

### No negative soil N pools (v5q-negpool, branch v5q-negpool)

| switch | what it does | references |
|---|---|---|
| `soil_negpool_fix` (0 legacy / 1 fix) | 1 = the "shuhua 2007" clamp in `Integrator::y2cstate_soi` can no longer write a negative ORGN (empty layer: AVLN = 0, ORGN = max(0, ORGN + AVLN), shortfall to `n_debt`); after `c2ystate_soi` any negative layer AVLN/ORGN or WDEBRISN is healed to 0 with the deficit booked to `n_debt` and repaid by `payNDebt`. Removes the permanent give-up/safety-net lock (RH = 0) triggered by an empty moss layer (CAP125 TVC_Warming_nfeed 2080). | **No published source for the fix** (code bug; N mass conservation). Principle: clipping/resetting a non-negative ODE state is mass-conserving only if the removed mass is booked [Shampine2005]. |

### Deep->mineral SOM transfer gated on thaw (v5q-d2m-gate, branch v5q-d2m-gate)

| switch | what it does | references |
|---|---|---|
| `d2m_transfer_gate` (0 legacy / 1 thaw-gated / 2 off) | 1 = the deep(humic)->mineral transfer of `Soil_Bgc::deltastate` step (II)3)-4) goes only into the part of the top 0.20 m of mineral soil above max(this year's running max ALD, previous year's ALD), scaled by that thawed fraction; C and N together at the donor C:N, receiver N `+=`, column C and N conserved; no transfer when ALD = 0. 2 = no transfer (physics-only spin-up phase). Previous-year ALD in the restart (`d2m_ald_prevyr`, optional; gate 0 files unchanged). | Permafrost SOM stays where frozen, transport only in the active layer [Koven2009]; [Koven2013]; CLM5/CTSM gating on max(altmax, altmax_lastyear) [Lawrence2019]. Spin-up to a stable ALD before carbon [Melton2019] (context). **No published source** for the f-scaling of the dvmdostem transfer. |

### Layer N follows C in the dynamic soil-layer remap (v5q-layern-fix, branch v5q-layern-fix)

| switch | what it does | references |
|---|---|---|
| `layer_n_remap_fix` (0 legacy / 1 fix) | 1 = (a) `Cohort::updateMonthly_DIMgrd` and the post-burn path copy layer `orgn`/`avln` from/to `bdall` around `Ground::redivideSoilLayers` together with the C pools (`Soil_Bgc::assignNitrogen{Bd2Layer,Layer2Bd}Monthly`; legacy copies only C, so after a layer split/combine the bd N stays at the old layer index while C moves); (b) organic `splitOneSoilLayer`: N and AVLN go to the new lower layer in the same fraction as its C (legacy: all N to the lower part, the upper part keeps its C with N = 0); a new deep layer formed from the last fibric layer takes N in the fraction of the C moved (legacy: all of the fibric N, fibric N = 0); (c) `combineTwoSoilLayersL2U` sums AVLN (legacy typo `usl->avln =+lsl->avln` assigns). Column C and N conserved; per-layer C:N follows C. No calibrated parameter changed. | Layer structure and the organic C–thickness relation of the dynamic organic soil module: [Yi2009], [Yi2010]. **No published source** for the remap rule: conservation / stoichiometry bookkeeping bug fix (N follows C). |

### Dynamic soil-layer C/N bookkeeping fixes (v5q-layer-bookkeep, branch v5q-layer-bookkeep)

| switch | what it does | references |
|---|---|---|
| `layer_bookkeep_fix` (0 legacy / 1 fix; 1 requires `layer_n_remap_fix` 1) | 1 = (1) `Ground::redivideShlwLayers`: a new fibric layer (none left) takes its C, orgn and avln out of the layer below at that layer's C:N (legacy: rawc = the month's aboveground litterfall C, created, N = 0); (2) `Ground::redivideDeepLayers`, `deepnum == 0`: the C/N of the removed humic layers go to the next soil layer (legacy: lost); (3) `BgcData::prvltrfcnque` follows its layer through the January re-division (legacy: stays at the old layer index, so kd of a layer is set from another layer's litterfall C:N history); (4) organic `splitOneSoilLayer`: the lower part gets at most the whole layer's C (legacy: upper C can go negative). Column C and N conserved; per-layer C:N kept. No calibrated parameter changed. | Layer structure and the organic C–thickness relation of the dynamic organic soil module: [Yi2009], [Yi2010]. **No published source** for the fixes: these are mass-conservation bookkeeping bug fixes. |

### Decomposition weighted by the layer thawed fraction (v5q-thawfrac-decomp, branch v5q-thawfrac-decomp)

| switch | what it does | references |
|---|---|---|
| `thawfrac_decomp` (0 legacy / 1 blend) | 1 = in `Soil_Bgc::deltac` the decomposition modifier of each soil layer is E = f·rhmoist(sws_th)·Q10^(max(ts,0)/10) + (1−f)·rhmoistfrozen·Q10^(min(ts,0)/10) with f = 1 − frozenfrac (monthly mean of the Stefan-front layer frozen fraction; missing → legacy ts indicator) and sws_th = min(1, sws/f) (liquid in the unfrozen part); applied to all C pools and the wood debris; net N mineralization from the thawed-part RH only (frozen part: legacy closed gate). Legacy uses the indicator ts < 0 → frozen rate, so a partly thawed zero-curtain layer (ts ≈ −0.005 °C) decomposes at rhmoistfrozen (0 at TVC). No calibrated parameter changed. | Layer-resolved decomposition with frozen-soil limitation [Koven2013]; CLM5 [Lawrence2019]; layer frozen fraction from the Stefan freeze/thaw fronts of the DOS-TEM soil module [Yi2009]. **No published source** for the linear thawed-fraction blend and the sws/f choice. |

### N drainage leaching from the current drainage-zone AVLN (v5q-leach-curn, branch v5q-leach-curn)

| switch | what it does | references |
|---|---|---|
| `n_leach_curn_fix` (0 legacy / 1 fix; 1 requires `n_leach_fix` 1) | 1 = in `Soil_Bgc::deltan` the drainage part of the inorganic N leaching `avlnlost` uses the current (integrator-state) AVLN of the drainage zone (fdz as `prepareIntegration`) instead of the month-start AVLN, so it goes to 0 as the zone empties. Legacy: a constant-rate sink during the month that drove an emptied frozen moss layer 0 negative (adapt() give-up, safety net; TVC P2, V5Q_STEP8 §13). Cap, split and per-layer removal (`n_leach_fix`) unchanged. | **No published source for the fix** (code defect). Principle: a zero-order sink drives a non-negative state below zero in finite time; positivity requires a state-dependent sink [Shampine2005]; [Burchard2003]. Leaching context as `n_leach_fix`: [McNamara2008], [Khosh2017]. |

### Parameter test (not a switch)

| change | status | references |
|---|---|---|
| TVC (CMT51) Sedges klai 0.15 → 0.5–1.0 | Tested only, not adopted (V5L_REPORT §5) | Vascular/vegetation cover at TVC: [WallaceBaltzer2019] (paper and data [WallaceBaltzer2019data]); [GrunbergBoike2019]; [Yazbeck2025]. Shrub LAI and interception: [Chang2022], [Chang2021data], [Zwieback2019]. Beer-law cover: [MonsiSaeki1953]. |

## 2. Non-switch fixes (always on in daily-climate-mode)

Format: commit — code location — principle — references.

- **Snow:**
  - 48053527 — Ground::updateSnowHorizon keeps Snow::extramass across reset() — snow mass conservation — **no published source**.
  - a671ad77 — snow bookkeeping (extramass melt; melt not removed twice) — mass conservation — **no published source**.
  - e3fcc91c — remove emptied snow layers (Ground::constructSnowLayers) — layer-index invariant, prevents a hang — **no published source**.
  - 228b0f08 (draft) — radiation-driven melt at Tair ≥ −1 °C routed to melt water instead of sublimation — energy/mass partition (melt ≠ sublimation) — **no published source**.
  - 08a525ba (draft) — restore the snow-surface n-factor (removes a TESTING line) — [Klene2001].
- **Water balance:**
  - 333d7197 — out-of-bounds write of root_water_up (Soil_Env::updateDailySM) — memory safety — **no published source**.
  - 337020ac — supply-limited ET (floor top-up and unapplied evaporation removed from reported ET) — water mass conservation — **no published source**.
  - 28167031 — reported EET and soil uptake both use Σ tran·FPC — consistent water accounting — area weighting as in [CLM5-Hydro] eq. 2.7.68. Its double FPC scaling of transpiration is undone by S1 above.
  - d9c30fc1 — freeze-expansion clip and organic-shrink water booked as drainage; puddle kept across years — mass conservation — **no published source**.
  - 40529d67 — soil evaporation days-since-rain reset only on days with ≥ 2 mm rain — **no published source**; rationale: drizzle should not reset the dry-surface clock.
- **Richards and Stefan:**
  - 032e8f13 — Richards work arrays sized MAX_SNW_LAY+MAX_SOI_LAY+2; clearRichardsArrays off-by-one — memory safety — **no published source**.
  - 7a7fe61a — Richards fail-fast guards — robustness — **no published source**.
  - 6044df4c, b7e43750 — Stefan front rebuild when the deque is empty or all fronts are swept — keeps the front state consistent with the layer states — algorithm context [Woo2004].
- **Thermal diagnostics:** d2c7d35e — ALD from the 0 °C isotherm (old value kept as ALD_FLAG). Rationale: active-layer thickness is the maximum annual thaw depth, i.e. the depth of the 0 °C isotherm. Standard definition; no specific citation verified here.
- **Carbon/nitrogen:**
  - 503bb21b, de49c6d5, ed28e049, e86633fe — root N extraction made N-conserving (Soil_Bgc) — N mass conservation — **no published source**.
  - be8a62fa — soil integrator safety net delivers the remaining monthly inputs when adapt() gives up — C mass conservation — **no published source**. It does not apply the remaining root N extraction (N created per event); the switchable fix is `n_safetynet_conserve` (v5m, above).
  - fe93d93c — NLOST reports only the AVLN loss — diagnostic definition — **no published source**.
- **Forcing and CO2:** 0e22205a — scenario CO2 indexed by calendar year, fail if not covered — correct time alignment — **no published source**.
- **Output and runtime:**
  - d6eee2c6, 689d8b97, 3a7cc57f — daily/monthly output time indexing — CF time coordinates [CF18].
  - 17d53b56 — startup undefined behaviour (array sizes, NUL-safe read, dangling pointer) — C++ memory safety — **no published source**.
  - 0f76919c — exit code 2 when a cell fails — **no published source**.
  - v5p-transp-units — `TRANSPIRATION` units label `mm/day` → `mm/time` in config/output_spec.csv (values are step sums: EnvData::veg_endOfDay/veg_endOfMonth); metadata only — **no published source**; units attribute must describe the stored value, [CF18] §3.1.
  - v5p-daily-water-fix — daily `TRANSPIRATION`, `EET`, `PET` NetCDF values taken from the cohort daily terms that EnvData sums into the monthly outputs (EET/PET no longer summed over PFTs, which counted ground evaporation once per PFT; TRANSPIRATION one value per day instead of per year); output only — **no published source**; rationale: daily values of an accumulated flux must sum to the monthly value.

## 3. Bibliography

- [Allen1998] Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). Crop evapotranspiration — Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56. FAO, Rome. https://www.fao.org/4/x0490e/x0490e00.htm
- [Brutsaert1982] Brutsaert, W. (1982). Evaporation into the Atmosphere: Theory, History and Applications. D. Reidel/Springer. https://doi.org/10.1007/978-94-017-1497-6
- [Buck1981] Buck, A.L. (1981). New equations for computing vapor pressure and enhancement factor. J. Appl. Meteorol. 20, 1527–1532. https://doi.org/10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2
- [Burchard2003] Burchard, H., Deleersnijder, E., Meister, A. (2003). A high-order conservative Patankar-type discretisation for stiff systems of production–destruction equations. Appl. Numer. Math. 47, 1–30. https://doi.org/10.1016/S0168-9274(03)00101-6
- [CF18] CF Conventions Committee. NetCDF Climate and Forecast (CF) Metadata Conventions, v1.8, section "Calendar". https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html
- [Chadburn2015] Chadburn, S., Burke, E., Essery, R., Boike, J., Langer, M., Heikenfeld, M., Cox, P., Friedlingstein, P. (2015). An improved representation of physical permafrost dynamics in the JULES land-surface model. Geosci. Model Dev. 8, 1493–1508. https://doi.org/10.5194/gmd-8-1493-2015
- [Chang2022] Chang, Q., Zwieback, S., DeVries, B., Berg, A. (2022). Application of L-band SAR for mapping tundra shrub biomass, leaf area index, and rainfall interception. Remote Sens. Environ. 268, 112747. https://doi.org/10.1016/j.rse.2021.112747
- [Chang2021data] Chang, Q., Zwieback, S., Berg, A., DeVries, B. (2021). Shrub aboveground biomass and Leaf Area Index maps, Trail Valley Creek Watershed, Northwest Territories. Borealis (dataset). https://doi.org/10.5683/SP2/O7MCCR
- [ChowLin1971] Chow, G.C., Lin, A. (1971). Best linear unbiased interpolation, distribution, and extrapolation of time series by related series. Rev. Econ. Stat. 53, 372–375. https://doi.org/10.2307/1928739
- [ClappHornberger1978] Clapp, R.B., Hornberger, G.M. (1978). Empirical equations for some soil hydraulic properties. Water Resour. Res. 14, 601–604. https://doi.org/10.1029/WR014i004p00601
- [Cleveland1999] Cleveland, C.C., Townsend, A.R., Schimel, D.S., Fisher, H., Howarth, R.W., Hedin, L.O., Perakis, S.S., Latty, E.F., Von Fischer, J.C., Elseroad, A., Wasson, M.F. (1999). Global patterns of terrestrial biological nitrogen (N2) fixation in natural ecosystems. Global Biogeochem. Cycles 13, 623–645. https://doi.org/10.1029/1999GB900014
- [CLM5] Lawrence, D.M. et al. (2019). The Community Land Model Version 5: Description of new features, benchmarking, and impact of forcing uncertainty. J. Adv. Model. Earth Syst. 11, 4245–4287. https://doi.org/10.1029/2018MS001583
- [CLM5-Hydro] CLM5.0 Technical Note, ch. 2.7 Hydrology. https://escomp.github.io/CTSM/release-clm5.0/tech_note/Hydrology/CLM50_Tech_Note_Hydrology.html
- [CLM5-Fluxes] CLM5.0 Technical Note, ch. 2.5 Momentum, Sensible Heat, and Latent Heat Fluxes. https://escomp.github.io/CTSM/release-clm5.0/tech_note/Fluxes/CLM50_Tech_Note_Fluxes.html
- [CLM5-SoilT] CLM5.0 Technical Note, ch. 2.6 Soil and Snow Temperatures. https://escomp.github.io/CTSM/release-clm5.0/tech_note/Soil_Snow_Temperatures/CLM50_Tech_Note_Soil_Snow_Temperatures.html
- [CoteKonrad2005] Côté, J., Konrad, J.-M. (2005). A generalized thermal conductivity model for soils and construction materials. Can. Geotech. J. 42, 443–458. https://doi.org/10.1139/t04-106
- [CrawfordDuchon1999] Crawford, T.M., Duchon, C.E. (1999). An improved parameterization for estimating effective atmospheric emissivity for use in calculating daytime downwelling longwave radiation. J. Appl. Meteorol. 38, 474–480. https://doi.org/10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2
- [DaviesBarnard2020] Davies-Barnard, T., Friedlingstein, P. (2020). The global distribution of biological nitrogen fixation in terrestrial natural ecosystems. Global Biogeochem. Cycles 34, e2019GB006387. https://doi.org/10.1029/2019GB006387
- [DeLuca2002] DeLuca, T.H., Zackrisson, O., Nilsson, M.-C., Sellstedt, A. (2002). Quantifying nitrogen-fixation in feather moss carpets of boreal forests. Nature 419, 917–920. https://doi.org/10.1038/nature01051
- [Denton1971] Denton, F.T. (1971). Adjustment of monthly or quarterly series to annual totals: an approach based on quadratic minimization. J. Am. Stat. Assoc. 66, 99–102. https://doi.org/10.1080/01621459.1971.10482227
- [Dickinson1993] Dickinson, R.E., Henderson-Sellers, A., Kennedy, P.J. (1993). Biosphere-Atmosphere Transfer Scheme (BATS) Version 1e as coupled to the NCAR Community Climate Model. NCAR Tech. Note NCAR/TN-387+STR. https://doi.org/10.5065/D67W6959
- [Dyer1974] Dyer, A.J. (1974). A review of flux-profile relationships. Boundary-Layer Meteorol. 7, 363–372. https://doi.org/10.1007/BF00240838
- [EOL-VP] NCAR/UCAR Earth Observing Laboratory. Water vapor pressure formulations (reproduces the WMO 2008 CIMO ice formula). https://www.eol.ucar.edu/data-software/conventions-and-standards/water-vapor-pressure-formulations
- [Euskirchen2009] Euskirchen, E.S., McGuire, A.D., Chapin, F.S. III, Yi, S., Thompson, C.C. (2009). Changes in vegetation in northern Alaska under scenarios of climate change, 2003–2100: implications for climate feedbacks. Ecol. Appl. 19, 1022–1043. https://doi.org/10.1890/08-0806.1
- [Farouki1981] Farouki, O.T. (1981). Thermal Properties of Soils. CRREL Monograph 81-1, U.S. Army CRREL, Hanover NH. https://doi.org/10.21236/ADA111734
- [Fisher2016] Fisher, J.B., Sikka, M., Huntzinger, D.N., Schwalm, C., Liu, J. (2016). Technical note: 3-hourly temporal downscaling of monthly global terrestrial biosphere model net ecosystem exchange. Biogeosciences 13, 4271–4277. https://doi.org/10.5194/bg-13-4271-2016
- [Genet2013] Genet, H. et al. (2013). Modeling the effects of fire severity and climate warming on active layer thickness and soil carbon storage of black spruce forests across the landscape in interior Alaska. Environ. Res. Lett. 8, 045016. https://doi.org/10.1088/1748-9326/8/4/045016
- [Gerten2004] Gerten, D., Schaphoff, S., Haberlandt, U., Lucht, W., Sitch, S. (2004). Terrestrial vegetation and water balance — hydrological evaluation of a dynamic global vegetation model. J. Hydrol. 286, 249–270. https://doi.org/10.1016/j.jhydrol.2003.09.029
- [Goodrich1982] Goodrich, L.E. (1982). The influence of snow cover on the ground thermal regime. Can. Geotech. J. 19, 421–432. https://doi.org/10.1139/t82-047
- [GrunbergBoike2019] Grünberg, I., Boike, J. (2019). Vegetation map of Trail Valley Creek, Northwest Territories, Canada. PANGAEA (dataset). https://doi.org/10.1594/PANGAEA.904270 (metadata: https://epic.awi.de/id/eprint/50259/1/metadata_vegetation_map.pdf)
- [HairerNorsettWanner1993] Hairer, E., Nørsett, S.P., Wanner, G. (1993). Solving Ordinary Differential Equations I: Nonstiff Problems, 2nd ed. Springer. https://doi.org/10.1007/978-3-540-78862-1
- [Harris2020] Harris, I., Osborn, T.J., Jones, P., Lister, D. (2020). Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. Sci. Data 7, 109. https://doi.org/10.1038/s41597-020-0453-3
- [HedstromPomeroy1998] Hedstrom, N.R., Pomeroy, J.W. (1998). Measurements and modelling of snow interception in the boreal forest. Hydrol. Process. 12, 1611–1625. https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4
- [Heijmans2004] Heijmans, M.M.P.D., Arp, W.J., Chapin, F.S. III (2004). Carbon dioxide and water vapour exchange from understory species in boreal forest. Agric. For. Meteorol. 123, 135–147. https://doi.org/10.1016/j.agrformet.2003.12.006
- [vonHoyningenHuene1981] von Hoyningen-Huene, J. (1981). Die Interzeption des Niederschlags in landwirtschaftlichen Pflanzenbeständen. Arbeitsbericht, Deutscher Verband für Wasserwirtschaft und Kulturbau (DVWK), Braunschweig. No DOI. Cited as such in the LISFLOOD references (https://ec-jrc.github.io/lisflood-model/5_1_references/). Published as DVWK Schriften 57 (1983), Parey, Hamburg/Berlin, ISBN 3-490-05797-X; full text https://wiki.bluemodel.org/images/9/9e/DVWK_57_I.pdf
- [Johansen1975] Johansen, O. (1975). Thermal conductivity of soils. PhD thesis, Univ. of Trondheim. CRREL Draft English Translation 637 (1977). https://doi.org/10.21236/ADA044002
- [Khosh2017] Khosh, M.S., McClelland, J.W., Jacobson, A.D., Douglas, T.A., Barker, A.J., Lehn, G.O. (2017). Seasonality of dissolved nitrogen from spring melt to fall freezeup in Alaskan Arctic tundra and mountain streams. J. Geophys. Res. Biogeosci. 122, 1718–1737. https://doi.org/10.1002/2016jg003377
- [Klene2001] Klene, A.E., Nelson, F.E., Shiklomanov, N.I., Hinkel, K.M. (2001). The N-factor in natural landscapes: variability of air and soil-surface temperatures, Kuparuk River Basin, Alaska, U.S.A. Arct. Antarct. Alp. Res. 33, 140–148. https://doi.org/10.1080/15230430.2001.12003416
- [LawrenceSlater2008] Lawrence, D.M., Slater, A.G. (2008). Incorporating organic soil into a global climate model. Clim. Dyn. 30, 145–160. https://doi.org/10.1007/s00382-007-0278-1
- [LeePielke1992] Lee, T.J., Pielke, R.A. (1992). Estimating the soil surface specific humidity. J. Appl. Meteorol. 31, 480–484. https://doi.org/10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2
- [Letts2000] Letts, M.G., Roulet, N.T., Comer, N.T., Skarupa, M.R., Verseghy, D.L. (2000). Parametrization of peatland hydraulic properties for the Canadian Land Surface Scheme. Atmos.-Ocean 38, 141–160. https://doi.org/10.1080/07055900.2000.9649643
- [Louis1979] Louis, J.-F. (1979). A parametric model of vertical eddy fluxes in the atmosphere. Boundary-Layer Meteorol. 17, 187–202. https://doi.org/10.1007/BF00117978
- [McGuire1992] McGuire, A.D., Melillo, J.M., Joyce, L.A., Kicklighter, D.W., Grace, A.L., Moore, B. III, Vorosmarty, C.J. (1992). Interactions between carbon and nitrogen dynamics in estimating net primary productivity for potential vegetation in North America. Global Biogeochem. Cycles 6, 101–124. https://doi.org/10.1029/92GB00219
- [McNamara2008] McNamara, J.P., Kane, D.L., Hobbie, J.E., Kling, G.W. (2008). Hydrologic and biogeochemical controls on the spatial and temporal patterns of nitrogen and phosphorus in the Kuparuk River, arctic Alaska. Hydrol. Process. 22, 3294–3309. https://doi.org/10.1002/hyp.6920
- [Meyerholt2016] Meyerholt, J., Zaehle, S., Smith, M.J. (2016). Variability of projected terrestrial biosphere responses to elevated levels of atmospheric CO2 due to uncertainty in biological nitrogen fixation. Biogeosciences 13, 1491–1518. https://doi.org/10.5194/bg-13-1491-2016
- [MonsiSaeki1953] Monsi, M., Saeki, T. (1953). Über den Lichtfaktor in den Pflanzengesellschaften und seine Bedeutung für die Stoffproduktion. Jpn. J. Bot. 14, 22–52. English translation: On the factor light in plant communities and its importance for matter production. Ann. Bot. 95, 549–567 (2005). https://doi.org/10.1093/aob/mci052
- [Monteith1965] Monteith, J.L. (1965). Evaporation and environment. Symp. Soc. Exp. Biol. 19, 205–234. https://pubmed.ncbi.nlm.nih.gov/5321565/
- [MurphyKoop2005] Murphy, D.M., Koop, T. (2005). Review of the vapour pressures of ice and supercooled water for atmospheric applications. Q. J. R. Meteorol. Soc. 131, 1539–1565. https://doi.org/10.1256/qj.04.94
- [Nasholm1998] Näsholm, T., Ekblad, A., Nordin, A., Giesler, R., Högberg, M., Högberg, P. (1998). Boreal forest plants take up organic nitrogen. Nature 392, 914–916. https://doi.org/10.1038/31921
- [OlsenRanderson2004] Olsen, S.C., Randerson, J.T. (2004). Differences between surface and column atmospheric CO2 and implications for carbon cycle research. J. Geophys. Res. Atmos. 109(D2), 2003JD003968 (article number as registered). https://doi.org/10.1029/2003JD003968
- [Prata1996] Prata, A.J. (1996). A new long-wave formula for estimating downward clear-sky radiation at the surface. Q. J. R. Meteorol. Soc. 122, 1127–1151. https://doi.org/10.1002/qj.49712253306
- [Press2007] Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P. (2007). Numerical Recipes: The Art of Scientific Computing, 3rd ed. Cambridge Univ. Press. https://numerical.recipes/book.html
- [QuintonMarsh1999] Quinton, W.L., Marsh, P. (1999). A conceptual framework for runoff generation in a permafrost environment. Hydrol. Process. 13, 2563–2581. https://doi.org/10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D
- [Raich1991] Raich, J.W., Rastetter, E.B., Melillo, J.M., Kicklighter, D.W., Steudler, P.A., Peterson, B.J., Grace, A.L., Moore, B. III, Vorosmarty, C.J. (1991). Potential net primary productivity in South America: application of a global model. Ecol. Appl. 1, 399–429. https://doi.org/10.2307/1941899
- [Richards1931] Richards, L.A. (1931). Capillary conduction of liquids through porous mediums. Physics 1, 318–333. https://doi.org/10.1063/1.1745010
- [Rousk2013] Rousk, K., Jones, D.L., DeLuca, T.H. (2013). Moss-cyanobacteria associations as biogenic sources of nitrogen in boreal forest ecosystems. Front. Microbiol. 4, 150. https://doi.org/10.3389/fmicb.2013.00150
- [Schaphoff2018] Schaphoff, S., von Bloh, W., Rammig, A., et al. (2018). LPJmL4 — a dynamic global vegetation model with managed land — Part 1: Model description. Geosci. Model Dev. 11, 1343–1375. https://doi.org/10.5194/gmd-11-1343-2018
- [SchimelBennett2004] Schimel, J.P., Bennett, J. (2004). Nitrogen mineralization: challenges of a changing paradigm. Ecology 85, 591–602. https://doi.org/10.1890/03-8002
- [Koven2009] Koven, C., Friedlingstein, P., Ciais, P., Khvorostyanov, D., Krinner, G., Tarnocai, C. (2009). On the formation of high-latitude soil carbon stocks: Effects of cryoturbation and insulation by organic matter in a land surface model. Geophys. Res. Lett. 36, L21501. https://doi.org/10.1029/2009GL040150
- [Koven2013] Koven, C.D., Riley, W.J., Subin, Z.M., Tang, J.Y., Torn, M.S., Collins, W.D., Bonan, G.B., Lawrence, D.M., Swenson, S.C. (2013). The effect of vertically resolved soil biogeochemistry and alternate soil C and N models on C dynamics of CLM4. Biogeosciences 10, 7109–7131. https://doi.org/10.5194/bg-10-7109-2013
- [Lawrence2019] Lawrence, D.M., et al. (2019). The Community Land Model Version 5: Description of New Features, Benchmarking, and Impact of Forcing Uncertainty. J. Adv. Model. Earth Syst. 11, 4245–4287. https://doi.org/10.1029/2018MS001583
- [Melton2019] Melton, J.R., Verseghy, D.L., Sospedra-Alfonso, R., Gruber, S. (2019). Improving permafrost physics in the coupled Canadian Land Surface Scheme (v.3.6.2) and Canadian Terrestrial Ecosystem Model (v.2.1) (CLASS-CTEM). Geosci. Model Dev. 12, 4443–4467. https://doi.org/10.5194/gmd-12-4443-2019
- [Shampine2005] Shampine, L.F., Thompson, S., Kierzenka, J.A., Byrne, G.D. (2005). Non-negative solutions of ODEs. Appl. Math. Comput. 170, 556–569. https://doi.org/10.1016/j.amc.2004.12.011
- [Stewart2011] Stewart, K.J., Coxson, D., Grogan, P. (2011). Nitrogen inputs by associative cyanobacteria across a low Arctic tundra landscape. Arct. Antarct. Alp. Res. 43, 267–278. https://doi.org/10.1657/1938-4246-43.2.267
- [Sturm1997] Sturm, M., Holmgren, J., König, M., Morris, K. (1997). The thermal conductivity of seasonal snow. J. Glaciol. 43, 26–41. https://doi.org/10.3189/S0022143000002781
- [Sturm2001] Sturm, M., Holmgren, J., McFadden, J.P., Liston, G.E., Chapin, F.S. III, Racine, C.H. (2001). Snow–shrub interactions in Arctic tundra: a hypothesis with climatic implications. J. Climate 14, 336–344. https://doi.org/10.1175/1520-0442(2001)014<0336:SSIIAT>2.0.CO;2
- [Swenson2012] Swenson, S.C., Lawrence, D.M., Lee, H. (2012). Improved simulation of the terrestrial hydrological cycle in permafrost regions by the Community Land Model. J. Adv. Model. Earth Syst. 4. https://doi.org/10.1029/2012MS000165
- [Thornton2007] Thornton, P.E., Lamarque, J.-F., Rosenbloom, N.A., Mahowald, N.M. (2007). Influence of carbon-nitrogen cycle coupling on land model response to CO2 fertilization and climate variability. Global Biogeochem. Cycles 21, GB4018. https://doi.org/10.1029/2006GB002868
- [ThorntonRosenbloom2005] Thornton, P.E., Rosenbloom, N.A. (2005). Ecosystem model spin-up: estimating steady state conditions in a coupled terrestrial carbon and nitrogen cycle model. Ecol. Model. 189, 25–48. https://doi.org/10.1016/j.ecolmodel.2005.04.008
- [Turetsky2003] Turetsky, M.R. (2003). The role of bryophytes in carbon and nitrogen cycling. The Bryologist 106, 395–409. https://doi.org/10.1639/05
- [vanderKnijff2010] van der Knijff, J.M., Younis, J., de Roo, A.P.J. (2010). LISFLOOD: a GIS-based distributed model for river basin scale water balance and flood simulation. Int. J. Geogr. Inf. Sci. 24, 189–212. https://doi.org/10.1080/13658810802549154
- [WallaceBaltzer2019] Wallace, C.A., Baltzer, J.L. (2020; online 2019). Tall shrubs mediate abiotic conditions and plant communities at the taiga–tundra ecotone. Ecosystems 23, 828–841. https://doi.org/10.1007/s10021-019-00435-0
- [WallaceBaltzer2019data] Wallace, C.A., Baltzer, J.L. (2019). Abiotic and vegetation community data from Alnus alnobetula shrub and open tundra sites at Trail Valley Creek. Borealis (dataset). https://doi.org/10.5683/SP2/F4PQRY
- [Warren1982] Warren, S.G. (1982). Optical properties of snow. Rev. Geophys. 20, 67–89. https://doi.org/10.1029/RG020i001p00067
- [WMO2008] World Meteorological Organization (2008). Guide to Meteorological Instruments and Methods of Observation, WMO-No. 8, 7th ed., Part I ch. 4, Annex 4.B "Formulae for the computation of measures of humidity". WMO e-Library record (current edition): https://library.wmo.int/records/item/41650-guide-to-meteorological-instruments-and-methods-of-observation ; the 2008 formula is reproduced in [EOL-VP].
- [Woo2004] Woo, M.-K., Arain, M.A., Mollinga, M., Yi, S. (2004). A two-directional freeze and thaw algorithm for hydrologic and land surface modelling. Geophys. Res. Lett. 31, 2004GL019475 (article number as registered). https://doi.org/10.1029/2004GL019475
- [Yazbeck2025] Yazbeck, T., Bohrer, G., Sonnentag, O., et al. (2025). Shrub expansion can counteract carbon losses from warming tundra. J. Geophys. Res. Biogeosci. 130, e2024JG008721. https://doi.org/10.1029/2024JG008721
- [Yi2009] Yi, S., McGuire, A.D., Harden, J., Kasischke, E., et al. (2009). Interactions between soil thermal and hydrological dynamics in the response of Alaska ecosystems to fire disturbance. J. Geophys. Res. Biogeosci. 114. https://doi.org/10.1029/2008JG000841
- [Yi2010] Yi, S., McGuire, A.D., Kasischke, E., Harden, J., Manies, K., Mack, M., Turetsky, M. (2010). A dynamic organic soil biogeochemical model for simulating the effects of wildfire on soil environmental conditions and carbon dynamics of black spruce forests. J. Geophys. Res. Biogeosci. 115. https://doi.org/10.1029/2010JG001302
- [Yuan2012] Yuan, F.-M., Yi, S.-H., McGuire, A.D., Johnson, K.D., Liang, J., Harden, J.W., Kasischke, E.S., Kurz, W.A. (2012). Assessment of boreal forest historical C dynamics in the Yukon River Basin: relative roles of warming and fire regime change. Ecol. Appl. 22, 2091–2109. https://doi.org/10.1890/11-1957.1
- [Zwieback2019] Zwieback, S., Chang, Q., Marsh, P., Berg, A. (2019). Shrub tundra ecohydrology: rainfall interception is a major component of the water balance. Environ. Res. Lett. 14, 055005. https://doi.org/10.1088/1748-9326/ab1049
