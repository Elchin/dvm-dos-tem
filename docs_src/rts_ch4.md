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
  `pond_sat0 = 0.5` means that right after the slump half of the unsaturated part of the thawed column
  (scar floor/ponded mud) is anoxic. `pond_tau = 10 yr`: the slump floor drains, stabilises and revegetates on
  decadal time scales (Kokelj et al. 2013, J. Geophys. Res. Earth Surf.; Lantz & Kokelj 2008). These values are
  PLACEHOLDERS to be constrained by slump-floor chamber data.
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
| pond_sat0, pond_tau, lat_ch4_frac | 0.5, 10 yr, 0.02 | same | see above (PLACEHOLDERS) |

**Calibration targets:**
- FLUXNET-CH4 (Knox et al. 2019, BAMS; Delwiche et al. 2021, ESSD): most freshwater wetlands emit
  0-20 gC/m2/yr, and wet tundra is among the lowest-emitting wetland classes.
- Arctic wet tundra sites: typically ~1-10 gC/m2/yr (e.g. Samoylov, Wille et al. 2008; Utqiagvik/Atqasuk,
  Zona et al. 2016, where the cold season contributes >= 50%).
- BAWLD-CH4 (Kuhn et al. 2021, ESSD): tundra wetlands have the highest medians, while dry tundra
  (WT > 5 cm below the surface) is low and includes net uptake.

We use 5 gC/m2/yr for fully saturated wet sedge tundra. The dry case (CMT07 heath, model WT) then gives
~0.1-0.2 gC/m2/yr. Neither scheme has CH4 uptake by upland soils (no atmospheric CH4 oxidation sink).
