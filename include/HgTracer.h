#ifndef HGTRACER_H_
#define HGTRACER_H_
// Mercury (Hg) tracer that follows carbon (model_settings.hg_tracer, default 0 = off, nothing runs).
// Diagnostic only: it reads the C pools and C fluxes and never writes them, so C/N/water/energy are
// bit-identical with the switch on or off. See docs_src/WORK_NOTES_hg_tracer.md.
//
// Refs:
//   Schaefer et al. 2020, Nat. Commun. 11, 4650, doi:10.1038/s41467-020-18398-5 (SiBCASA-Hg: Hg moves
//     with every C transfer; of the Hg released by decay 16% evades as Hg(0), 1.5% is exported, the rest
//     re-binds to the soil organic matter of the same layer).
//   Schuster et al. 2018, GRL 45, 1463, doi:10.1002/2017GL075571 (RHgC median 1.6 ug Hg / g C).
//   T0 inventory (docs_src/HG_T0_inventory.md, Azaroff et al. 2026 synthesis, doi:10.17043/azaroff-2026-
//     permafrost-mercury-1): organic-horizon RHgC 0.32-0.35, mineral 0.97-1.69 (synthesis median 1.47).
//   Obrist et al. 2017, Nature 547, 201, doi:10.1038/nature22997; Olson et al. 2018, GBC 32, 1058 (tundra veg.
//     THg ~20-60 ng/g dry mass -> ~0.05-0.1 ug Hg / g C at 0.45 g C / g; doi:10.1029/2017GB005840).
//   Thompson et al. 2025, ES&T 59, 19457, doi:10.1021/acs.est.5c04510 (soil %MeHg 0.1 plateau, 0.7 bog,
//     2.0 fen) for the MeHg index (hg_mehg).
//
// T1b additions (2026-11-23): Hg spin-up (hg_spinup), deposition time series (hg_dep_file), optional
// downward leaching of decay-released Hg (hg_f_leach), and a lower default evasion fraction:
//   hg_f_evade 0.16 (Schaefer 2020, from GTMM / Smith-Downey et al. 2010, temperate soils) -> 0.02. At
//   Toolik, Obrist et al. 2017 (doi:10.1038/nature22997) measured a year-round net Hg(0) sink (6.5 +- 0.7
//   ug/m2/yr) with soil pore-air Hg(0) below the atmosphere; Olson et al. 2018 note "almost no Hg
//   re-emission" there; Lim et al. 2020 (BG 17, 3083) attribute the Toolik sink vs. a net Hg(0) source at a
//   warmer boreal mire (Degero, Osterwalder et al. 2018) to a climate control on re-emission.
//   Deposition (hg_dep_file, CSV): Hg(II) wet + Hg(II) dry + non-vegetation (soil/snow) Hg(0) uptake to the
//   top soil layer, plus veg_scale, a multiplier of the vegetation uptake (atmospheric Hg(0) relative to the
//   present day). See config/hg_deposition_toolik_obrist2017.csv.
//   Spin-up (hg_spinup = 1): after hg_spinup_cycle full years from the C x RHgC initialization, the affine
//   map of those years' tracer steps (Hg_end = M Hg_start + b; exact, the step is linear in Hg) is raised
//   to the power hg_spinup_years / hg_spinup_cycle and applied: the state after hg_spinup_years years of
//   repeating that C cycle and deposition (Schaefer 2020 used a 5000-yr spin-up). Booked as cum.spin.
//   Deposition and unassigned Hg go to the top soil layer that holds SOM C (the moss layer usually has
//   none), so Hg cannot pile up in a layer without C turnover. HG_DIAG (env) prints a replay check.
//
// T1c additions (2026-11-23): downward transport of dissolved Hg (hg_transport, default 1 with the tracer),
// spin-up on by default with the tracer (hg_spinup 1), GEOS-Chem site deposition files.
//   Once a month (a C-neutral sub-step, recorded for the spin-up like every other step; any geometry) the
//   Hg of each unfrozen layer partitions between the solid phase M (kg/m2: organic matter = C / 0.5;
//   mineral (1 - porosity) dz 2650 kg/m3) and the layer water liq (L/m2) with Kd (L/kg); the dissolved
//   (DOM-bound) part leaves with the model's lateral layer drainage (export, HG_DRAIN) and moves down with the
//   percolating water (see transport_; frozen layers stop the water). Hg reaching the bottom is not lost: the
//   vertical flux is zero at the frost table / bottom of the column (it can only leave by lateral drainage,
//   HG_DRAIN). Kd: mineral 10^3.8 L/kg (Allison & Allison 2005, EPA/600/R-05/074, soil Hg(II) median,
//   range 10^2.2-10^5.8); organic 10^3.3 calibrated to the Toolik organic-horizon Hg:C (0.29-0.41).
//   hg_spinup_cycle default 10 years (was 1): a single recorded year can have almost no transport (snowmelt
//   on frozen ground, dry summer), and the 5000-yr jump would then carry that year's state.
//
// T1d (2026-11-24): spin-up closure (hg_spinup_closure, default 0). The recorded C cycle is not periodic (layer
//   C drifts, organic layers thin or thicken, layers are re-divided, the layer count may differ at the two
//   ends), so the jump closes the record with a map from its end geometry back to its start geometry:
//   Hg (and MeHg) mass moved within each horizon, organic layers by relative depth in the organic horizon and
//   mineral layers by depth below the organic-mineral interface (depth_map_). The repeated map
//   start -> cycle -> end -> closure -> start is then periodic, the jump ends with one plain cycle in the
//   current geometry, and the first hg_spinup_cycle years are always usable (no layer-count condition).
//   Validated against explicit 300- and 1000-yr spin stages (WORK_NOTES 10). Defaults recalibrated:
//   hg_f_export 0.005 (was 0.015), hg_k_meth_ch4 3250, hg_k_meth_proxy 87.
//
// Units: C in g C/m2, Hg and MeHg in ug/m2, ratios in ug Hg / g C, fluxes per model month.
//
// Pools (parallel to C): live vegetation per PFT x {leaf, stem, root}; standing dead per PFT; woody
// debris; one Hg pool per soil layer (pool-weighted: the 4 SOM pools of a layer share its Hg:C ratio).
//
// One generic step (hgtr::step) is applied around every model sub-step that moves C (environment with
// daily CH4, ground re-division, BGC, fire, RTS) from C snapshots before (a) and after (b) plus the
// sub-step's external C fluxes (Flux):
//   1. inputs: NPP > 0 of each PFT part carries r_leaf/r_stem/r_root (atmospheric Hg(0) uptake);
//      optional direct soil deposition hg_dep_soil (ug/m2/yr) to the top soil layer in the BGC step.
//   2. decomposition (RH, CH4 production, woody debris RH) of a pool releases C x (Hg:C of the pool);
//      f_evade of it evades, f_export is exported, the rest stays in the pool (re-binding);
//      DOC loss exports its Hg.
//   3. internal transfers: dC_int = C(b) - C(a) - inputs + outputs for each pool. Pools that lost C
//      (donors: litterfall, standing dead fall, burned / slumped material, ...) put Hg at their Hg:C ratio
//      into one pot; aggregate sinks (fire emission, RTS lateral export, RTS respiration) take their
//      C share of the pot; the rest goes to the pools that gained C (receivers) in proportion to dC_int.
//   4. soil geometry change (layer split/combine, burn, slump): the soil column is re-mapped bottom-up
//      on cumulative C (uniform Hg:C within an old layer); a net soil C loss is taken from the top, a net
//      gain is added on top as a receiver.
// The step is exactly conservative: total Hg change = uptake + deposition - evasion - export - fire.
// C that leaves a donor with no receiver or sink leaves its Hg behind (concentration rises), C that
// appears in a receiver with no donor brings no Hg (dilution); both are reported as unexplained C.
//
// MeHg index (hg_mehg = 1, needs hg_tracer = 1): per soil layer, no feedback on Hg or C:
//   P = k_meth x Hg x A x q10_meth^((T - t_ref)/10),   kd = k_demeth x q10_demeth^((T - t_ref)/10)
//   A = anaerobic C turnover fraction of the month: CH4 production / layer C where the CH4 module is on
//       (cd.ch4_on), else saturated fraction (water table) x RH / layer C (proxy); separate k_meth for
//       the two drivers (k_meth_ch4, k_meth_proxy).
//   MeHg(t+1) = Meq + (MeHg(t) - Meq) exp(-kd), Meq = P / kd, capped at the layer Hg.
// T2 phase 1 (hg_t2 = 1, needs hg_tracer = 1; default 0 = T1 unchanged; docs_src/HG_T2_design.md):
// atmosphere-surface interface. Diagnostic like T1 (reads LAI, snow, water and radiation; writes nothing).
//   Foliar uptake (hg_t2_foliar, ORCHIDEE-Hg v1.0 stomate_Hg.f90 l.1046-1056): leaf Hg input of PFT p =
//     hg_fol_scale x f_fix(LAI_p) x vegcov_p x (1 - f_snow) x (dry Hg(0) + dry Hg(II)), f_fix = hg_fol_fmax x
//     min(1, LAI / hg_fol_lai) (0.25, 1.25), canopy total capped at 1; replaces the leaf NPP x hg_r_leaf uptake
//     (stem / root keep NPP x ratio). The intercepted dry Hg(II) leaves the ground deposition; the non-fixed dry
//     Hg(0) is re-emitted (GEOS-Chem dry Hg(0) is gross; not booked).
//   Wet deposition by precipitation (hg_t2_wetprec): month = annual x P_month / P_year.
//   Snowpack Hg (hg_t2_snow): Hg(II) on the snow-covered fraction f_snow = min(1, SWE / hg_swe_full) goes to the
//     snow pool; photoreduction of its photic part min(1, hg_snow_photic / z_snow) at hg_snow_kred x 24 <cos Z>
//     (ORCHIDEE l.845-849, 0.090 /h x cos zenith) -> evasion; release with melt (to the top C soil layer) and
//     sublimation (to the air) in proportion to melt, sublimation and SWE (ORCHIDEE l.1018-1020).
//   Photoreduction of ground Hg(II) (hg_t2_photored, default 0): hg_pr_a exp(-hg_pr_b SW) (ORCHIDEE l.852).
//   Booking: foliar -> uptake; Hg(II) to ground + snow and the soil Hg(0) sink -> dep; photoreduction and
//   snow sublimation -> evasion. The snow pool is part of the total, the restart (layout v3) and the spin-up
//   vector, so the step stays affine in Hg and exactly conservative.
// T2 phase 2 (hg_t2_diss = 1, needs hg_t2 = 1; default 0; docs_src/HG_T2_design.md section 3): explicit dissolved
// pool diss[l] (ug/m2, not part of soil[l], the sorbed Hg). Decay-released Hg that does not evade goes to diss of its
// layer (hg_f_export is not used: the export is explicit). Once a month (transport sub-step, C-neutral) the month's
// days are stepped with the model's daily water (Soil_Env qinfl, qover, Richards layer_drain) and daily layer liq / ice:
//   sorption  diss += (phi (soil + diss) - diss) (1 - exp(-1 day / hg_tau_sorb)), phi = liq / (liq + Kd M)
//             (Kd: hg_t2_log_kd_org / hg_log_kd_min, M as T1c; moss layers, which carry no soil C: M = dz x
//             hg_t2_rho_moss, the moss dry bulk density, phase-2 refit WORK_NOTES_hg_t2 section 2);
//   overland  diss of the layers starting above hg_z_runoff (ORCHIDEE 4.5 cm) x qover / (qover + their liq) -> export;
//   lateral   diss_l x layer_drain_l / liq_l -> export (drainage);
//   vertical  diss_l x hg_v_adv x q_l / liq_l to the layer below, q_l as T1c (infiltration not drained laterally,
//             removed linearly with depth; zero at the frost table); overland, lateral and vertical fractions act
//             together on the day's diss (no ordering bias), their sum capped at 0.9;
// then cryo/bioturbation of the sorbed Hg in the thawed zone (ORCHIDEE cryoturbate_Hg): D = hg_d_cryo (1 - z / hg_z_cryo)
//   + hg_d_bio (m2/yr), explicit, flux = D dt (c_l - c_l+1) / dz_mid with c = Hg / dz, coefficient capped at 0.25 dz.
// A geometry change (re-division, burn, slump) re-sorbs diss into soil first (diss re-equilibrates within days).
// All coefficients come from water, ice, C and geometry: the step stays affine in Hg (diss is in the spin-up vector).
// T2 phase 3 (docs_src/HG_T2_design.md section 3, WORK_NOTES_hg_t2 section 3), all with hg_t2 = 1, default off:
//   hg_t2_evasion = 1: the re-emitted fraction of decay-released Hg becomes per layer
//     hg_fd x (1 - f_snow) x [T_layer > 0] x exp(-z_mid / hg_ev_z) (woody debris: z = 0, air T; with a geometry change the
//     decomposition-weighted mean), replacing hg_f_evade; the soil Hg(0) sink (hg0_soil column) becomes a process
//     ground sink hg_v_ground x GEM x veg_scale x seconds, all months (Obrist 2017: winter net Hg(0) uptake under snow).
//   hg_t2_foliar = 2 (D3): leaf uptake = stomatal Hg(0) from the PFT's GPP (Medlyn g_sw = 1.6 m A / c_a, Hg(0)
//     diffusivity scaling, molar air density from air T, CO2 of the month; CLM5-Hg, Yuan et al. 2023) x GEM x veg_scale
//     x hg_gs_eff, plus the LAI interception of dry Hg(II) (phase 1); hg_t2_foliar = 1 keeps the LAI scheme.
//   hg_t2_meltroute = 1 (with hg_t2_diss): the month's snowmelt Hg (State::mpend, in the spin-up vector) is delivered by
//     the transport sub-step day by day with the day's melt; the share qover / (qover + qinfl) of the day leaves as
//     overland flow (export), the rest enters the top layer's dissolved pool.
//   Coefficients come from C, GPP, water, snow and temperature only: the step stays affine and exactly conservative.
// T2 phase 4 (docs_src/HG_T2_design.md section 3, WORK_NOTES_hg_t2 section 4), all with hg_t2 = 1, default off:
//   hg_fire_ash > 0: of the fire-emitted Hg that came from burned soil (the soil donors' share of the step's pot), the
//     fraction hg_fire_ash stays in the ash / char on the new top C layer (HG_ASH, internal); burned vegetation, standing
//     dead and woody debris Hg is all emitted (Friedli et al. 2001).
//   hg_t2_rts = 1: RTS particulate export: the slumped material's Hg leaves as particulate export (HG_PARTIC, part of
//     HG_EXPORT), including the Hg of the RTS-respired C share that is not re-emitted (T1 returned it to the column).
//   hg_wash > 0: canopy wash-off: the share hg_wash x P_rain / (P_rain + hg_wash_pref) of the month's canopy-intercepted dry
//     Hg(II) goes to the ground with the throughfall instead of the leaves (HG_WASH, part of HG_DEP).
// D7 (T2): hg_snow_sublim_loss = 0 (default): Hg of sublimated snow stays in the pack (concentrates) and leaves with
// the melt (release fraction melt / (melt + SWE)); 1 = ORCHIDEE (sublimation share lost to the air).
// T2 phase 5b (docs_src/HG_T2_design.md section 18, WORK_NOTES_hg_t2 section 6), all need hg_t2 = 1:
//   hg_t2_kdsplit = 1 (with hg_t2_diss): humic (deep organic, type 2) layers sorb with log10 Kd = hg_t2_log_kd_org +
//     hg_t2_dlog_kd_hum; moss and fibric (shallow) layers keep hg_t2_log_kd_org (decomposed humic peat binds Hg(II) to more
//     reduced-S sites per unit C than fibric moss litter: Skyllberg et al. 2003 Environ Sci Technol 37:2048, Grigal 2003).
//   hg_t2_xpath = 1 (with hg_t2_diss): explicit export path. Of the dissolved Hg leaving the organic layers above the first
//     mineral layer by overland flow (incl. the melt overland share) or lateral drainage, the fraction hg_f_conn reaches the
//     stream (export); the rest is re-sorbed along the hillslope / riparian flow path, booked into the dissolved pool of the
//     first mineral layer (HG_DOWN; a 1-D stand-in for flow-path retention, e.g. INCA-Hg / Futter et al. 2012).
//   hg_t2_snowsink = 1 (with hg_t2_evasion): the ground Hg(0) sink of the snow-covered fraction f_snow uses its own
//     deposition velocity hg_v_snow (net Hg(0) uptake through the snowpack into the soil, Obrist et al. 2017: 37 % of
//     the annual net uptake in Oct-May): sink = GEM x veg_scale x s x ((1 - f_snow) hg_v_ground + f_snow hg_v_snow).
//   hg_t2_preset (default 1, acts only with hg_t2 = 1): the absent hg_t2 settings take the phase-5 calibrated values
//     (Params::t2_preset_); hg_t2_preset = 0 keeps the component defaults (all T2 options off). The T1 / Hg-off paths
//     (hg_t2 = 0) never read these. T2-only copies hg_t2_r_stem, hg_t2_r_root, hg_t2_k_meth_proxy replace hg_r_stem,
//     hg_r_root, hg_k_meth_proxy when hg_t2 = 1 (giving the T1 names with hg_t2 = 1 is an error).
//   All new coefficients come from water, snow, geometry and type only: the step stays affine and exactly conservative.
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <json/json.h>
#include "cohortconst.h"
#include "layerconst.h"

namespace hgtr {

struct Params {
  int    tracer     = 0;      // hg_tracer: 0 off (default), 1 on
  double r_org      = 0.35;   // hg_r_org: initial RHgC of moss/fibric/humic layers (T0 Toolik-region organic median)
  double r_min      = 1.47;   // hg_r_min: initial RHgC of mineral layers (T0 synthesis mineral median)
  double r_leaf     = 0.07;   // hg_r_leaf: RHgC of new leaf tissue (NPP), ~30 ng/g dry mass
  double r_stem     = 0.02;   // hg_r_stem: RHgC of new stem tissue (wood ~5-10 ng/g)
  double r_root     = 0.07;   // hg_r_root: RHgC of new root tissue
  double f_evade    = 0.02;   // hg_f_evade: fraction of decay-released Hg that evades (0.16 Schaefer 2020; 0.02: Obrist 2017 Toolik net Hg(0) sink, see top)
  double f_export   = 0.005;  // hg_f_export: fraction of decay-released Hg exported (Schaefer 2020: 0.015). T1d: 0.005, since
                              // dissolved export now leaves explicitly with the drainage (hg_transport); keeps TVC export
                              // within the pan-Arctic river yields 0.8-4.0 ug m-2 yr-1 (Zolkos 2020), see WORK_NOTES 10
  double f_leach    = 0.0;    // hg_f_leach: fraction of decay-released Hg moved to the next soil layer down (DOC-bound), bottom -> export
  // T1c downward transport of dissolved (DOM-bound) Hg with the percolating soil water, once a month, any
  // geometry: Hg leaving layer l = q_l * Hg_l / (liq_l + Kd * M_l) (Kd partitioning between the solid phase
  // M_l (kg/m2) and the layer water liq_l (L/m2)); q = infiltration at the top -> drainage at the bottom of the
  // thawed zone (linear in depth); Hg leaving the deepest thawed layer leaves with the drainage (export).
  int    transport  = 1;      // hg_transport: 1 on (with hg_tracer), 0 off
  double log_kd_org = 3.3;    // hg_log_kd_org (log10 L/kg): organic layers; calibrated (T1c) to the Toolik organic-horizon Hg:C 0.29-0.41
  double log_kd_min = 3.8;    // hg_log_kd_min (log10 L/kg): mineral layers; Allison & Allison 2005 soil Hg(II) median 3.8 (2.2-5.8)
  double dep_soil   = 0.0;    // hg_dep_soil: constant deposition to the top soil layer (ug/m2/yr); fallback without hg_dep_file
  std::string dep_file;       // hg_dep_file: optional CSV deposition series (see dep_at); "" = constant hg_dep_soil
  int    tr_year0   = 1901;   // hg_tr_year0: calendar year of tr-run year 0 (deposition lookup)
  int    sc_year0   = 2016;   // hg_sc_year0: calendar year of sc-run year 0
  int    spin_year  = 1850;   // hg_dep_spin_year: deposition year used in pre-run / eq-run / sp-run (pre-industrial)
  int    spinup     = 1;      // hg_spinup: 1 (default; acts only with hg_tracer = 1), 0 off; 1 = analytical spin-up after the first hg_spinup_cycle years
  double spinup_years = 5000; // hg_spinup_years: equivalent spin-up length (yr; Schaefer 2020: 5000)
  int    spinup_closure = 0;  // hg_spinup_closure: T1d, how the recorded cycle is closed (end -> start geometry): 0 Hg mass
                              // by depth overlap (default), 1 C-neutral model step (bottom-up C re-map), 2 none (raw repetition)
  int    spinup_cycle = 10;   // hg_spinup_cycle: model years whose C + water cycle is recorded and repeated (1-50); T1c: 10 (was 1)
                              // so that the interannual variability of infiltration / thaw (transport) is averaged;
                              // T1d: the first full years after the initialization, closed by a C-neutral closure step
  // deposition table (from hg_dep_file): rows with month 0 = annual (interpolated in year), 1-12 = monthly
  std::vector<int> dep_y, dep_m; std::vector<double> dep_v, dep_vs;
  int    mehg       = 0;      // hg_mehg: 0 off (default), 1 MeHg index
  int    mehg_driver = 0;     // hg_mehg_driver: 0 auto (CH4 production if the CH4 module is on), 1 proxy always
  // k_meth: calibrated (T1, 2026-11-23) to Thompson et al. 2025 growing-season (Jun-Sep) active-layer %MeHg in
  // steady state: CH4 driver -> 2.0 % (fen) at Toolik CMT06 wet sedge with the CH4 module on (k 30 gave 0.0189 %);
  // proxy driver -> 0.1 % (permafrost plateau) at HPC S21H (k 0.6 gave 0.00078 %). %MeHg is linear in k_meth.
  // T1d recalibration with spin-up (10-yr cycle, closure 0), transport and deposition series (GEOS-Chem TVC/HPC, Obrist
  // Toolik): Toolik CH4 driver 1.97 % at 3200 -> 3250 (fen 2 %); HPC proxy 0.089 % at 77 -> 87 (plateau 0.1 %).
  double k_meth_ch4   = 3250.0; // hg_k_meth_ch4 (-): MeHg production per unit Hg x CH4-C turnover fraction
  double k_meth_proxy = 87.0;   // hg_k_meth_proxy (-): same, per saturated-fraction RH C turnover fraction
  double k_demeth   = 3.0;    // hg_k_demeth (1/month at t_ref): ~0.1 /day, peat incubations
  double q10_meth   = 1.0;    // hg_q10_meth: extra T response of methylation (the C turnover already has one)
  double q10_demeth = 2.0;    // hg_q10_demeth
  double t_ref      = 10.0;   // hg_t_ref (C)
  // T2 phase 1 (see top). Deposition components from hg_dep_file (ug/m2/yr): wet, dry Hg(II), gross dry Hg(0)
  // (column hg0_dry; if absent hg0_soil / 0.37, the factor of the T1 GEOS-Chem files), soil Hg(0) sink.
  int    t2          = 0;     // hg_t2: 0 off (default), 1 T2 phase-1 surface interface (needs hg_tracer = 1, hg_dep_file components)
  int    t2_foliar   = 1;     // hg_t2_foliar: LAI-based foliar uptake (with hg_t2)
  int    t2_snow     = 1;     // hg_t2_snow: snowpack Hg pool (with hg_t2)
  int    t2_wetprec  = 1;     // hg_t2_wetprec: wet deposition in proportion to the month's precipitation (with hg_t2)
  int    t2_photored = 0;     // hg_t2_photored: ORCHIDEE photoreduction of ground-deposited Hg(II) (decision D2; default off)
  double fol_fmax    = 0.25;  // hg_fol_fmax: maximum fraction of dry deposition fixed by leaves (ORCHIDEE l.1051)
  double fol_lai     = 1.25;  // hg_fol_lai: LAI at which the fixed fraction saturates (ORCHIDEE l.1049)
  double fol_scale   = 1.0;   // hg_fol_scale: calibration factor of the foliar uptake (decision D3)
  double snow_kred   = 0.090; // hg_snow_kred (1/h at cos Z = 1): snow Hg(II) photoreduction (ORCHIDEE l.845)
  double snow_photic = 0.10;  // hg_snow_photic (m): photic depth of snow (King & Simpson 2001, UV e-folding 5-20 cm)
  double snow_rho    = 250.0; // hg_snow_rho (kg/m3): snow density for the photic fraction
  double swe_full    = 10.0;  // hg_swe_full (mm): SWE of full snow cover
  double pr_a        = 0.67;  // hg_pr_a: photoreduced fraction of deposited Hg(II) at SW = 0 (ORCHIDEE l.852)
  double pr_b        = 0.016; // hg_pr_b (m2/W): its decrease with shortwave (ORCHIDEE l.852)
  int    snow_sublim_loss = 0; // hg_snow_sublim_loss: 0 Hg stays in the pack when snow sublimates (D7, default), 1 ORCHIDEE (lost)
  // T2 phase 2 (see top)
  int    t2_diss     = 0;     // hg_t2_diss: explicit dissolved pool + daily export + cryo/bioturbation (needs hg_t2 = 1)
  double t2_log_kd_org = 3.5; // hg_t2_log_kd_org (log10 L/kg): organic Kd with hg_t2_diss; phase-2 refit (WORK_NOTES_hg_t2
                              // section 2) to Toolik organic Hg:C 0.29-0.41 (Obrist 2017) and HPC 0.32-0.35 (T0): 3.5
                              // (Toolik 0.35, HPC 0.33, TVC 0.43); provisional until the phase-5 joint refit
  double tau_sorb    = 0.25;  // hg_tau_sorb (days): sorption relaxation time (Hg(II)-soil sorption near equilibrium within
                              // hours, Yin et al. 1996 J Environ Qual 25:837; phase-2 scan: tau <= 0.5 d all alike)
  double rho_moss    = 30.0;  // hg_t2_rho_moss (kg/m3): moss-layer dry bulk density as its sorbent mass (live + dead
                              // feather moss / Sphagnum 10-50 kg/m3, O'Donnell et al. 2009 Soil Sci 174:646)
  double v_adv       = 0.2;   // hg_v_adv: dissolved Hg advection as a fraction of the water flux (ORCHIDEE)
  double z_runoff    = 0.045; // hg_z_runoff (m): depth of the overland-flow mixing zone (ORCHIDEE runoff, top 4.5 cm)
  double d_cryo      = 1.0e-3; // hg_d_cryo (m2/yr): cryoturbation diffusivity at the surface (ORCHIDEE cryoturbate_Hg)
  double z_cryo      = 3.0;   // hg_z_cryo (m): cryoturbation decreases linearly to 0 at this depth (ORCHIDEE)
  double d_bio       = 1.0e-4; // hg_d_bio (m2/yr): bioturbation diffusivity (ORCHIDEE)
  // T2 phase 3 (see top)
  int    t2_evasion  = 0;     // hg_t2_evasion: process soil Hg(0) re-emission + ground Hg(0) sink (needs hg_t2 = 1)
  int    t2_meltroute = 0;    // hg_t2_meltroute: snowmelt Hg delivered daily with the meltwater (needs hg_t2_diss = 1)
  double fd          = 0.053; // hg_fd: Hg(0) re-emitted fraction of decay-released Hg at the surface, thawed, snow-free
                              // (ORCHIDEE-Hg fd_Hg, after Smith-Downey et al. 2010 JGR 115 G03008)
  double ev_z        = 0.05;  // hg_ev_z (m): e-folding depth of the re-emission (soil-air Hg(0) is depleted below the
                              // top few cm at Toolik: Obrist et al. 2017 Nature 547:201, soil-gas profiles)
  double v_ground    = 8.5e-5; // hg_v_ground (m/s): Hg(0) deposition velocity to the ground (litter / moss / soil), all
                              // year; 0.4 ng/m2/h under snow at GEM 1.3 (Obrist 2017: 37 % of 6.5 ug/m2/yr, Oct-May)
  double gem         = 1.3;   // hg_gem (ng/m3): atmospheric Hg(0) (Arctic background 1.2-1.5, AMAP 2021), x veg_scale
  double gs_m        = 4.0;   // hg_gs_m: Medlyn 1 + g1 / sqrt(D) (g1 2.2-4.5 kPa^0.5 CLM5 arctic grass / boreal shrub,
                              // D 0.5-1 kPa): stomatal conductance g_sw = 1.6 m A / c_a (Yuan et al. 2023 Environ Int 174:107904)
  double gs_eff      = 0.75;  // hg_gs_eff: calibration factor of the stomatal Hg(0) uptake (hg_t2_foliar = 2; decision D3):
                              // phase-3 fit to the Toolik vegetation Hg(0) uptake ~4.1 ug/m2/yr (Obrist 2017; 1.0 gave 5.3)
  int    gs_ingpp    = 1;     // hg_gs_ingpp: 1 the stomatal uptake follows the PFT's N-unlimited GPP (INGPP; dvmdostem's canopy
                              // conductance is not N-limited, its N limitation scales GPP afterwards), 0 the N-limited GPP
  // T2 phase 4 (see top)
  double fire_ash    = 0.0;   // hg_fire_ash: fraction of the burned soil Hg retained in ash / char (0 = T1: all emitted).
                              // Ku et al. 2018 ES&T 52:13366 (white ash: >= 90 % lost; black ash 34-83 % lost); Biswas et al. 2007
                              // GBC 21 GB1005 (>95 % released above 300 C); Harden et al. 2004 GBC 18 GB3014 (FROSTFIRE, ~70 % of
                              // the burned forest-floor Hg lost); phase-4 test value 0.1, range 0.03-0.3
  int    t2_rts      = 0;     // hg_t2_rts: RTS particulate export (St Pierre et al. 2018 ES&T 52:14099) (needs hg_t2 = 1)
  double wash        = 0.0;   // hg_wash: maximum wash-off fraction of the month's canopy-intercepted dry Hg(II) (0 = off;
                              // Rea et al. 2000 ES&T 34:2418, Wright et al. 2016 Environ Rev 24:205: throughfall carries the
                              // washed dry Hg(II)); phase-4 test value 0.5
  double wash_pref   = 10.0;  // hg_wash_pref (mm/month): half-saturation rainfall of the wash-off
  // T2 phase 5b (see top)
  int    t2_preset   = 1;     // hg_t2_preset: 1 absent hg_t2 settings take the calibrated values (t2_preset_), 0 component defaults
  int    t2_kdsplit  = 0;     // hg_t2_kdsplit: humic-layer Kd offset (needs hg_t2_diss = 1)
  double t2_dlog_kd_hum = 0.0; // hg_t2_dlog_kd_hum (log10 L/kg): humic minus moss/fibric organic log10 Kd
  int    t2_xpath    = 0;     // hg_t2_xpath: explicit export path / flow-path retention (needs hg_t2_diss = 1)
  double f_conn      = 1.0;   // hg_f_conn: share of the organic-layer dissolved export that reaches the stream
  int    t2_snowsink = 0;     // hg_t2_snowsink: separate ground Hg(0) sink velocity under snow (needs hg_t2_evasion = 1)
  double v_snow      = 8.5e-5; // hg_v_snow (m/s): Hg(0) deposition velocity to the ground through the snowpack
  double t2_r_stem   = 0.02;  // hg_t2_r_stem: T2 copy of hg_r_stem (used when hg_t2 = 1)
  double t2_r_root   = 0.07;  // hg_t2_r_root: T2 copy of hg_r_root (used when hg_t2 = 1)
  double t2_k_meth_proxy = 87.0; // hg_t2_k_meth_proxy: T2 copy of hg_k_meth_proxy (used when hg_t2 = 1)
  bool   dep_comp    = false; // the deposition file has components (hg2_wet / hg2_dry / hg0_soil columns)
  std::vector<double> dep_w, dep_d2, dep_d0, dep_s0;   // components per row (ug/m2/yr)
};

// T2 phase 5: calibrated values of the hg_t2 settings (applied with hg_t2 = 1 and hg_t2_preset = 1 before the explicit
// model_settings; WORK_NOTES_hg_t2 sections 5-6, config/hg_t2_calibrated.json). Never used with hg_t2 = 0.
inline void t2_preset_(Params& p) {
  // phase 5b fit (P5G; WORK_NOTES_hg_t2 section 6): Toolik net Hg(0) 6.4 / winter 37 % from v_ground + v_snow (exact linear
  // solve), Hg:C / export from Kd, fd, ev_z (emulator on E5 + E6, 90 runs); hg_t2_xpath left off (f_conn fit 1.0)
  p.t2_diss = 1; p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1; p.fire_ash = 0.1; p.t2_rts = 1; p.wash = 0.5;
  p.t2_log_kd_org = 3.453; p.t2_kdsplit = 1; p.t2_dlog_kd_hum = 0.242; p.t2_xpath = 0; p.f_conn = 1.0;
  p.t2_snowsink = 1; p.v_ground = 4.88e-5; p.v_snow = 1.283e-4; p.fd = 0.083; p.ev_z = 0.50; p.snow_kred = 0.090;
  p.t2_r_stem = 0.0044; p.t2_r_root = 0.0067; p.t2_k_meth_proxy = 71.0;
}

// Deposition CSV: '#' comments; header with columns year, [month], and either dep or any of hg2_wet,
// hg2_dry, hg0_soil (summed), optional veg_scale (default 1). Rates in ug/m2/yr (monthly rows too).
inline std::string load_dep(const std::string& fn, Params& p) {
  std::ifstream f(fn.c_str());
  if (!f) { return "hg_dep_file: cannot open " + fn; }
  std::string line; std::vector<std::string> hdr; int ln = 0;
  p.dep_y.clear(); p.dep_m.clear(); p.dep_v.clear(); p.dep_vs.clear();
  p.dep_w.clear(); p.dep_d2.clear(); p.dep_d0.clear(); p.dep_s0.clear(); p.dep_comp = false;
  int iw = -1, id2 = -1, is0 = -1, id0 = -1;   // T2 components (hg0_dry: gross dry Hg(0), not added to the T1 rate)
  auto split = [](const std::string& s) { std::vector<std::string> o; std::stringstream ss(s); std::string x;
    while (std::getline(ss, x, ',')) { size_t a = x.find_first_not_of(" \t\r"), b = x.find_last_not_of(" \t\r");
      o.push_back(a == std::string::npos ? "" : x.substr(a, b - a + 1)); } return o; };
  int iy = -1, im = -1, ivs = -1; std::vector<int> iv;
  while (std::getline(f, line)) {
    ln++;
    if (line.find_first_not_of(" \t\r") == std::string::npos || line[line.find_first_not_of(" \t\r")] == '#') { continue; }
    std::vector<std::string> c = split(line);
    if (hdr.empty()) {
      hdr = c;
      for (size_t i = 0; i < c.size(); i++) {
        if (c[i] == "year") iy = (int)i; else if (c[i] == "month") im = (int)i; else if (c[i] == "veg_scale") ivs = (int)i;
        else if (c[i] == "dep" || c[i] == "hg2_wet" || c[i] == "hg2_dry" || c[i] == "hg0_soil") iv.push_back((int)i);
        if (c[i] == "hg2_wet") iw = (int)i; else if (c[i] == "hg2_dry") id2 = (int)i; else if (c[i] == "hg0_soil") is0 = (int)i;
        else if (c[i] == "hg0_dry") id0 = (int)i;
      }
      p.dep_comp = (iw >= 0 || id2 >= 0 || is0 >= 0);
      if (iy < 0 || iv.empty()) { return "hg_dep_file " + fn + ": header needs 'year' and 'dep' or hg2_wet/hg2_dry/hg0_soil"; }
      continue;
    }
    if (c.size() < hdr.size()) { return "hg_dep_file " + fn + ": short row at line " + std::to_string(ln); }
    try {
      const int y = std::stoi(c[iy]); const int m = (im >= 0) ? std::stoi(c[im]) : 0;
      double v = 0.0; for (int k : iv) { v += std::stod(c[k]); }
      const double vs = (ivs >= 0) ? std::stod(c[ivs]) : 1.0;
      if (m < 0 || m > 12 || !(v >= 0.0) || !(vs >= 0.0) || !std::isfinite(v) || !std::isfinite(vs)) { throw std::invalid_argument("range"); }
      if (!p.dep_y.empty() && m == 0 && p.dep_m.back() == 0 && y <= p.dep_y.back()) { throw std::invalid_argument("annual years must increase"); }
      p.dep_y.push_back(y); p.dep_m.push_back(m); p.dep_v.push_back(v); p.dep_vs.push_back(vs);
      auto col = [&](int k) { const double x = (k >= 0) ? std::stod(c[k]) : 0.0; if (!(x >= 0.0) || !std::isfinite(x)) { throw std::invalid_argument("range"); } return x; };
      const double w = col(iw), d2 = col(id2), s0 = col(is0);
      p.dep_w.push_back(w); p.dep_d2.push_back(d2); p.dep_s0.push_back(s0); p.dep_d0.push_back((id0 >= 0) ? col(id0) : s0 / 0.37);
    } catch (const std::exception& e) {
      return "hg_dep_file " + fn + ": bad value at line " + std::to_string(ln) + " (" + e.what() + ")";
    }
  }
  if (p.dep_y.empty()) { return "hg_dep_file " + fn + ": no data rows"; }
  return "";
}

struct DepVal { double rate = 0.0; double veg_scale = 1.0; double wet = 0.0, dry2 = 0.0, dry0 = 0.0, soil0 = 0.0; };
// Deposition rate (ug/m2/yr) and vegetation-uptake scale for a calendar year and month (1-12): an exact
// monthly row if present, else linear interpolation of the annual rows (constant beyond the ends), else the
// constant hg_dep_soil with scale 1.
inline DepVal dep_at(const Params& p, int year, int month) {
  DepVal d; d.rate = p.dep_soil;
  if (p.dep_y.empty()) { return d; }
  for (size_t i = 0; i < p.dep_y.size(); i++) {
    if (p.dep_m[i] == month && p.dep_y[i] == year) {
      d.rate = p.dep_v[i]; d.veg_scale = p.dep_vs[i];
      if (i < p.dep_w.size()) { d.wet = p.dep_w[i]; d.dry2 = p.dep_d2[i]; d.dry0 = p.dep_d0[i]; d.soil0 = p.dep_s0[i]; }
      return d;
    }
  }
  int lo = -1, hi = -1;
  for (size_t i = 0; i < p.dep_y.size(); i++) {
    if (p.dep_m[i] != 0) { continue; }
    if (p.dep_y[i] <= year) { lo = (int)i; }
    if (p.dep_y[i] >= year && hi < 0) { hi = (int)i; }
  }
  if (lo < 0 && hi < 0) { return d; }
  if (lo < 0) { lo = hi; } if (hi < 0) { hi = lo; }
  const double w = (p.dep_y[hi] > p.dep_y[lo]) ? double(year - p.dep_y[lo]) / double(p.dep_y[hi] - p.dep_y[lo]) : 0.0;
  d.rate = p.dep_v[lo] + w * (p.dep_v[hi] - p.dep_v[lo]);
  d.veg_scale = p.dep_vs[lo] + w * (p.dep_vs[hi] - p.dep_vs[lo]);
  if ((size_t)hi < p.dep_w.size()) {
    d.wet = p.dep_w[lo] + w * (p.dep_w[hi] - p.dep_w[lo]); d.dry2 = p.dep_d2[lo] + w * (p.dep_d2[hi] - p.dep_d2[lo]);
    d.dry0 = p.dep_d0[lo] + w * (p.dep_d0[hi] - p.dep_d0[lo]); d.soil0 = p.dep_s0[lo] + w * (p.dep_s0[hi] - p.dep_s0[lo]);
  }
  return d;
}

inline bool num_(const Json::Value& v) { return !v.isBool() && v.isNumeric() && std::isfinite(v.asDouble()); }

// Parse model_settings; returns "" or the FATAL message. Absent keys = defaults; unknown hg_* keys are errors.
inline std::string parse(const Json::Value& ms, Params& p) {
  p = Params();
  if (!ms.isObject()) { return ""; }
  // T2 phase 5: with hg_t2 = 1 (and hg_t2_preset absent or 1) the calibrated values replace the component defaults
  // (validated below with the other keys); explicit settings still win
  const auto num1 = [&](const char* k) { return ms.isMember(k) && num_(ms[k]) && ms[k].asDouble() == 1.0; };
  const bool pre = num1("hg_t2") && (!ms.isMember("hg_t2_preset") || num1("hg_t2_preset"));
  if (pre) { t2_preset_(p); }
  struct D { const char* k; double* v; double lo, hi; };
  D ds[] = {{"hg_r_org", &p.r_org, 0, 1e3}, {"hg_r_min", &p.r_min, 0, 1e3}, {"hg_r_leaf", &p.r_leaf, 0, 1e3},
            {"hg_r_stem", &p.r_stem, 0, 1e3}, {"hg_r_root", &p.r_root, 0, 1e3}, {"hg_f_evade", &p.f_evade, 0, 1},
            {"hg_f_export", &p.f_export, 0, 1}, {"hg_f_leach", &p.f_leach, 0, 1}, {"hg_dep_soil", &p.dep_soil, 0, 1e4},
            {"hg_spinup_years", &p.spinup_years, 1, 1e6}, {"hg_log_kd_org", &p.log_kd_org, 0, 9}, {"hg_log_kd_min", &p.log_kd_min, 0, 9},
            {"hg_k_meth_ch4", &p.k_meth_ch4, 0, 1e9}, {"hg_k_meth_proxy", &p.k_meth_proxy, 0, 1e9},
            {"hg_k_demeth", &p.k_demeth, 1e-9, 1e6}, {"hg_q10_meth", &p.q10_meth, 1e-3, 100},
            {"hg_q10_demeth", &p.q10_demeth, 1e-3, 100}, {"hg_t_ref", &p.t_ref, -50, 50},
            {"hg_fol_fmax", &p.fol_fmax, 0, 1}, {"hg_fol_lai", &p.fol_lai, 1e-3, 100}, {"hg_fol_scale", &p.fol_scale, 0, 100},
            {"hg_snow_kred", &p.snow_kred, 0, 100}, {"hg_snow_photic", &p.snow_photic, 0, 100}, {"hg_snow_rho", &p.snow_rho, 10, 1000},
            {"hg_swe_full", &p.swe_full, 1e-3, 1e4}, {"hg_pr_a", &p.pr_a, 0, 1}, {"hg_pr_b", &p.pr_b, 0, 10},
            {"hg_t2_log_kd_org", &p.t2_log_kd_org, 0, 9}, {"hg_tau_sorb", &p.tau_sorb, 1e-3, 1e5}, {"hg_t2_rho_moss", &p.rho_moss, 0, 1000}, {"hg_v_adv", &p.v_adv, 0, 1},
            {"hg_z_runoff", &p.z_runoff, 0, 1}, {"hg_d_cryo", &p.d_cryo, 0, 1}, {"hg_z_cryo", &p.z_cryo, 1e-3, 100}, {"hg_d_bio", &p.d_bio, 0, 1},
            {"hg_fd", &p.fd, 0, 1}, {"hg_ev_z", &p.ev_z, 1e-4, 100}, {"hg_v_ground", &p.v_ground, 0, 1}, {"hg_gem", &p.gem, 0, 100},
            {"hg_gs_m", &p.gs_m, 1, 100}, {"hg_gs_eff", &p.gs_eff, 0, 100},
            {"hg_fire_ash", &p.fire_ash, 0, 1}, {"hg_wash", &p.wash, 0, 1}, {"hg_wash_pref", &p.wash_pref, 1e-3, 1e4},
            {"hg_t2_dlog_kd_hum", &p.t2_dlog_kd_hum, -3, 3}, {"hg_f_conn", &p.f_conn, 0, 1}, {"hg_v_snow", &p.v_snow, 0, 1},
            {"hg_t2_r_stem", &p.t2_r_stem, 0, 1e3}, {"hg_t2_r_root", &p.t2_r_root, 0, 1e3}, {"hg_t2_k_meth_proxy", &p.t2_k_meth_proxy, 0, 1e9}};
  struct I { const char* k; int* v; int lo, hi; };
  I is[] = {{"hg_tracer", &p.tracer, 0, 1}, {"hg_mehg", &p.mehg, 0, 1}, {"hg_mehg_driver", &p.mehg_driver, 0, 1},
            {"hg_spinup", &p.spinup, 0, 1}, {"hg_transport", &p.transport, 0, 1}, {"hg_spinup_cycle", &p.spinup_cycle, 1, 50}, {"hg_spinup_closure", &p.spinup_closure, 0, 2}, {"hg_tr_year0", &p.tr_year0, -100000, 100000},
            {"hg_sc_year0", &p.sc_year0, -100000, 100000}, {"hg_dep_spin_year", &p.spin_year, -100000, 100000},
            {"hg_t2", &p.t2, 0, 1}, {"hg_t2_foliar", &p.t2_foliar, 0, 2}, {"hg_t2_snow", &p.t2_snow, 0, 1},
            {"hg_t2_wetprec", &p.t2_wetprec, 0, 1}, {"hg_t2_photored", &p.t2_photored, 0, 1},
            {"hg_t2_diss", &p.t2_diss, 0, 1}, {"hg_snow_sublim_loss", &p.snow_sublim_loss, 0, 1},
            {"hg_t2_evasion", &p.t2_evasion, 0, 1}, {"hg_t2_meltroute", &p.t2_meltroute, 0, 1}, {"hg_gs_ingpp", &p.gs_ingpp, 0, 1},
            {"hg_t2_rts", &p.t2_rts, 0, 1}, {"hg_t2_preset", &p.t2_preset, 0, 1}, {"hg_t2_kdsplit", &p.t2_kdsplit, 0, 1},
            {"hg_t2_xpath", &p.t2_xpath, 0, 1}, {"hg_t2_snowsink", &p.t2_snowsink, 0, 1}};
  for (const std::string& key : ms.getMemberNames()) {
    if (key.rfind("hg_", 0) != 0) { continue; }
    bool known = false;
    for (D& d : ds) {
      if (key == d.k) {
        known = true; const Json::Value& v = ms[key];
        if (!num_(v) || v.asDouble() < d.lo || v.asDouble() > d.hi) {
          return "model_settings." + key + " must be a number in [" + std::to_string(d.lo) + ", " + std::to_string(d.hi) + "]";
        }
        *d.v = v.asDouble();
      }
    }
    for (I& d : is) {
      if (key == d.k) {
        known = true; const Json::Value& v = ms[key];
        if (!num_(v) || v.asDouble() != std::floor(v.asDouble()) || v.asInt() < d.lo || v.asInt() > d.hi) {
          return "model_settings." + key + " must be an integer in [" + std::to_string(d.lo) + ", " + std::to_string(d.hi) + "]";
        }
        *d.v = v.asInt();
      }
    }
    if (key == "hg_dep_file") {
      known = true;
      if (!ms[key].isString()) { return "model_settings.hg_dep_file must be a string (path of a CSV file)"; }
      p.dep_file = ms[key].asString();
    }
    if (!known) { return "model_settings." + key + ": unknown hg_* setting"; }
  }
  if (pre) {   // preset options whose prerequisite was switched off explicitly follow it (unless given explicitly too)
    const auto set = [&](const char* k) { return ms.isMember(k); };
    if (p.t2_diss != 1 && !set("hg_t2_meltroute")) { p.t2_meltroute = 0; }
    if (p.t2_diss != 1 && !set("hg_t2_kdsplit")) { p.t2_kdsplit = 0; }
    if (p.t2_diss != 1 && !set("hg_t2_xpath")) { p.t2_xpath = 0; }
    if (p.t2_evasion != 1 && !set("hg_t2_snowsink")) { p.t2_snowsink = 0; }
    if (p.t2_foliar < 1 && !set("hg_wash")) { p.wash = 0.0; }
  }
  if (p.t2 == 1) {   // T2-only copies of the parameters shared with T1 (the T1 path, hg_t2 = 0, keeps the T1 values)
    for (const char* k : {"hg_r_stem", "hg_r_root", "hg_k_meth_proxy"}) {
      if (ms.isMember(k)) { return std::string("model_settings.") + k + " is a T1 setting; with hg_t2 = 1 use hg_t2_" + (k + 3); }
    }
    p.r_stem = p.t2_r_stem; p.r_root = p.t2_r_root; p.k_meth_proxy = p.t2_k_meth_proxy;
  }
  if (p.f_evade + p.f_export + p.f_leach > 1.0) { return "model_settings.hg_f_evade + hg_f_export + hg_f_leach must be <= 1"; }
  // hg_spinup (default 1) and hg_transport (default 1) act only with hg_tracer = 1 (ignored otherwise)
  if (!p.dep_file.empty()) { std::string e = load_dep(p.dep_file, p); if (!e.empty()) { return e; } }
  if (p.mehg == 1 && p.tracer != 1) { return "model_settings.hg_mehg = 1 needs hg_tracer = 1"; }
  if (p.t2 == 1 && p.tracer != 1) { return "model_settings.hg_t2 = 1 needs hg_tracer = 1"; }
  if (p.t2_diss == 1 && p.t2 != 1) { return "model_settings.hg_t2_diss = 1 needs hg_t2 = 1"; }
  if (p.t2_evasion == 1 && p.t2 != 1) { return "model_settings.hg_t2_evasion = 1 needs hg_t2 = 1"; }
  if (p.t2_foliar == 2 && p.t2 != 1) { return "model_settings.hg_t2_foliar = 2 needs hg_t2 = 1"; }
  if (p.t2_meltroute == 1 && !(p.t2_diss == 1 && p.transport == 1)) { return "model_settings.hg_t2_meltroute = 1 needs hg_t2_diss = 1 and hg_transport = 1"; }
  if (p.fire_ash > 0.0 && p.t2 != 1) { return "model_settings.hg_fire_ash > 0 needs hg_t2 = 1"; }
  if (p.t2_rts == 1 && p.t2 != 1) { return "model_settings.hg_t2_rts = 1 needs hg_t2 = 1"; }
  if (p.wash > 0.0 && !(p.t2 == 1 && p.t2_foliar >= 1)) { return "model_settings.hg_wash > 0 needs hg_t2 = 1 and hg_t2_foliar >= 1"; }
  if (p.t2_kdsplit == 1 && !(p.t2 == 1 && p.t2_diss == 1)) { return "model_settings.hg_t2_kdsplit = 1 needs hg_t2 = 1 and hg_t2_diss = 1"; }
  if (p.t2_xpath == 1 && !(p.t2 == 1 && p.t2_diss == 1)) { return "model_settings.hg_t2_xpath = 1 needs hg_t2 = 1 and hg_t2_diss = 1"; }
  if (p.t2_snowsink == 1 && !(p.t2 == 1 && p.t2_evasion == 1)) { return "model_settings.hg_t2_snowsink = 1 needs hg_t2 = 1 and hg_t2_evasion = 1"; }
  if (p.t2 == 1 && !p.dep_comp) { return "model_settings.hg_t2 = 1 needs an hg_dep_file with hg2_wet / hg2_dry / hg0_soil columns"; }
  return "";
}

const int NPOOL = NUM_PFT * NUM_PFT_PART + NUM_PFT + 1;   // veg parts, standing dead, woody debris

// C snapshot of the column (g C/m2).
struct Snap {
  double veg[NUM_PFT][NUM_PFT_PART] = {};
  double dead[NUM_PFT] = {};
  double wd = 0.0;
  int nsl = 0;
  double soil[MAX_SOI_LAY] = {};
  double dz[MAX_SOI_LAY] = {};
  int type[MAX_SOI_LAY] = {};
};

// External C fluxes of one sub-step (g C/m2 over the sub-step).
struct Flux {
  double npp[NUM_PFT][NUM_PFT_PART] = {};   // net uptake (may be < 0: C loss without Hg)
  double decomp[MAX_SOI_LAY] = {};          // soil decomposition C (RH + CH4 production), by layer
  double doc[MAX_SOI_LAY] = {};             // DOC C exported, by layer
  double decomp_wd = 0.0;                   // woody debris RH
  double dep_soil = 0.0;                    // direct Hg deposition to the top soil layer (ug/m2), not C
  double veg_scale = 1.0;                   // multiplier of the vegetation uptake ratios (deposition series)
  double burn = 0.0;                        // fire C emission (aggregate, from the donor pot)
  double lateral = 0.0;                     // RTS lateral C export (aggregate, from the donor pot)
  double pot_decomp = 0.0;                  // RTS respiration (aggregate decomposition from the pot)
  // T1c transport sub-step (tr = 1; C unchanged): month's water (mm = L/m2) and layer state
  int    tr = 0;
  double q_in = 0.0, q_out = 0.0;           // infiltration, total drainage (mm/month; q_out: record only)
  double drain_l[MAX_SOI_LAY] = {};         // lateral drainage by soil layer (mm/month; Richards layer_drain)
  double liq[MAX_SOI_LAY] = {}, ice[MAX_SOI_LAY] = {}, por[MAX_SOI_LAY] = {};   // monthly mean liq / ice (mm), porosity (-)
  // T2 phase 1 surface sub-step (t2 = 1; BGC step only): month's deposition components (ug/m2), canopy, snow, radiation
  int    t2 = 0;
  double dep_wet = 0.0, dep_dry2 = 0.0, dep_dry0 = 0.0, dep_soil0 = 0.0;   // wet Hg(II), dry Hg(II), gross dry Hg(0), soil Hg(0) sink
  double lai[NUM_PFT] = {}, vegcov[NUM_PFT] = {};
  double swe = 0.0, melt = 0.0, sublim = 0.0;   // monthly-mean SWE (mm), melt and sublimation (mm/month)
  double sw = 0.0, cosz = 0.0, ndays = 30.0;    // mean shortwave (W/m2), 24-h mean cos zenith, days
  // T2 phase 2 (hg_t2_diss, transport sub-step): the month's days (mm/day; liq / ice mm at the end of the day)
  int    nd = 0;
  double dq_in[31] = {}, dq_over[31] = {};
  double d_drain[31][MAX_SOI_LAY] = {}, d_liq[31][MAX_SOI_LAY] = {}, d_ice[31][MAX_SOI_LAY] = {};
  // T2 phase 3: BGC step: GPP by PFT (g C/m2/month; INGPP with hg_gs_ingpp = 1), air T (C), CO2 (ppm), layer T (C); transport step: daily melt (mm/day)
  double gpp[NUM_PFT] = {}; double ta = 0.0, co2 = 0.0; double ts[MAX_SOI_LAY] = {};
  double d_melt[31] = {};
  double rain = 0.0;                        // T2 phase 4 (hg_wash, BGC step): the month's rainfall (mm)
};

struct Acc {
  double uptake = 0, dep = 0, evasion = 0, export_ = 0, fire = 0, unexpl_closs = 0, unexpl_cgain = 0, spin = 0;
  double drain = 0, down = 0;   // T1c diagnostics: drainage export (part of export_), organic -> mineral transfer
  double foliar = 0, snowdep = 0, snowred = 0, snowmelt = 0, photored = 0;   // T2 diagnostics (parts of uptake / dep / evasion; melt internal)
  double runoff = 0, mixdown = 0;   // T2 phase 2: overland-flow export (part of export_), net organic -> mineral mixing
  double ash = 0, partic = 0, wash = 0;   // T2 phase 4: fire Hg kept in ash (internal), RTS particulate export (part of export_),
                                         // canopy wash-off (part of dep)
  double soilev = 0, sink0 = 0, stomatal = 0, meltrun = 0;   // T2 phase 3: soil re-emission (part of evasion), ground Hg(0) sink (part of dep),
                                                             // stomatal uptake (part of foliar), melt Hg exported with the melt overland flow (part of runoff)
  void clear() { *this = Acc(); }
  void add(const Acc& o) {
    uptake += o.uptake; dep += o.dep; evasion += o.evasion; export_ += o.export_; fire += o.fire;
    unexpl_closs += o.unexpl_closs; unexpl_cgain += o.unexpl_cgain; spin += o.spin; drain += o.drain; down += o.down;
    foliar += o.foliar; snowdep += o.snowdep; snowred += o.snowred; snowmelt += o.snowmelt; photored += o.photored;
    runoff += o.runoff; mixdown += o.mixdown;
    soilev += o.soilev; sink0 += o.sink0; stomatal += o.stomatal; meltrun += o.meltrun;
    ash += o.ash; partic += o.partic; wash += o.wash;
  }
  double net() const { return uptake + dep - evasion - export_ - fire + spin; }   // spin: spin-up jump
};

struct State {
  bool init = false;
  double veg[NUM_PFT][NUM_PFT_PART] = {};
  double dead[NUM_PFT] = {};
  double wd = 0.0;
  int nsl = 0;
  double soil[MAX_SOI_LAY] = {};
  double mehg[MAX_SOI_LAY] = {};
  double snow = 0.0;      // T2: snowpack Hg (ug/m2)
  double diss[MAX_SOI_LAY] = {};   // T2 phase 2: dissolved Hg by layer (ug/m2; not part of soil)
  double mpend = 0.0;     // T2 phase 3 (hg_t2_meltroute): snowmelt Hg released this month, delivered by the transport sub-step
  Acc m, y;               // this month / this year (outputs)
  Acc cum;                // since initialization (conservation diagnostic)
  double h_init = 0.0;    // total Hg at initialization
  double mehg_prod_m = 0.0, mehg_prod_y = 0.0;   // MeHg gross production (ug/m2) month / year
  int spun = 0;           // spin-up: 0 not done, 1 done, 2 gave up (geometry never repeated)
  int spin_tries = 0;
  double vegtot() const { double s = 0; for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) s += veg[i][j]; return s; }
  double deadtot() const { double s = wd; for (int i = 0; i < NUM_PFT; i++) s += dead[i]; return s; }
  double soiltot() const { double s = 0; for (int i = 0; i < MAX_SOI_LAY; i++) s += soil[i]; return s; }
  double mehgtot() const { double s = 0; for (int i = 0; i < MAX_SOI_LAY; i++) s += mehg[i]; return s; }
  double disstot() const { double s = 0; for (int i = 0; i < MAX_SOI_LAY; i++) s += diss[i]; return s; }
  double total() const { return vegtot() + deadtot() + soiltot() + snow + disstot() + mpend; }
  // conservation residual since initialization (ug/m2); roundoff only
  double balance_error() const { return total() - h_init - cum.net(); }
};

inline bool is_mineral(int type) { return type >= 3; }   // CohortData m_soil.type: 0 moss, 1 shallow, 2 deep, 3 mineral
inline bool is_moss(int type) { return type == 0; }

// Initialize Hg from C x RHgC (soil by layer type; vegetation, standing dead and woody debris at the
// tissue ratios: leaf ratio for standing dead and woody debris).
inline void initialize(const Params& p, State& s, const Snap& c) {
  const double rv[NUM_PFT_PART] = {p.r_leaf, p.r_stem, p.r_root};
  for (int i = 0; i < NUM_PFT; i++) {
    for (int j = 0; j < NUM_PFT_PART; j++) { s.veg[i][j] = std::fmax(0.0, c.veg[i][j]) * rv[j]; }
    s.dead[i] = std::fmax(0.0, c.dead[i]) * p.r_leaf;
  }
  s.wd = std::fmax(0.0, c.wd) * p.r_leaf;
  s.nsl = c.nsl;
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    s.soil[l] = (l < c.nsl) ? std::fmax(0.0, c.soil[l]) * (is_mineral(c.type[l]) ? p.r_min : p.r_org) : 0.0;
    s.mehg[l] = 0.0;
  }
  s.snow = 0.0;
  for (int l = 0; l < MAX_SOI_LAY; l++) { s.diss[l] = 0.0; }
  s.init = true;
  s.cum.clear();
  s.h_init = s.total();
}

inline double ratio_(double h, double c) { return (c > 1.0e-12 && h > 0.0) ? h / c : 0.0; }

inline bool same_geometry(const Snap& a, const Snap& b) {
  if (a.nsl != b.nsl) { return false; }
  for (int l = 0; l < a.nsl; l++) { if (a.dz[l] != b.dz[l]) { return false; } }
  return true;
}

// T1c: transport of dissolved Hg (see Params). Contiguous unfrozen layers from the top (not ice-dominated:
// ice <= 0.5 mm or ice <= liq; a dry moss layer passes the water). Water of the month (mm = L/m2):
//   lateral drainage d_l = the model's own layer drainage (Richards layer_drain, saturated layers above the
//   drainage depth) -> its dissolved Hg is exported (HG_DRAIN);
//   vertical flux at the bottom of layer l: (q_in - D) (1 - z_l / Z) + (D - sum_{k<=l} d_k) >= 0, i.e. the
//   infiltration not drained laterally (ET, storage) is removed linearly with depth, the drained water where
//   the model drains it; zero at the bottom of the unfrozen zone (no flow into frozen ground).
// Dissolved concentration c_l = Hg_l / (liq_l + Kd M_l); M: organic layers C / 0.5 (organic matter),
// mineral (1 - por) dz 2650 kg/m3. MeHg moves with its layer's Hg (same fraction). Explicit upwind, outflow
// fraction capped at 0.9. Linear in Hg.
inline void transport_(const Params& p, State& s, const Snap& b, const Flux& f, Acc& acc) {
  const int nb = std::min(b.nsl, MAX_SOI_LAY);
  const double qi = std::fmax(0.0, f.q_in);
  int tb = -1; double Z = 0.0, D = 0.0;
  for (int l = 0; l < nb; l++) {
    if (f.ice[l] > 0.5 && f.ice[l] > f.liq[l]) { break; }   // frozen (ice-dominated) layer: the water stops
    tb = l; Z += std::fmax(0.0, b.dz[l]); D += std::fmax(0.0, f.drain_l[l]);
  }
  if (tb < 0 || !(Z > 0.0) || !(qi > 0.0 || D > 0.0)) { return; }
  const double kdo = std::pow(10.0, p.log_kd_org), kdm = std::pow(10.0, p.log_kd_min);
  double out[MAX_SOI_LAY] = {}, lat[MAX_SOI_LAY] = {}, mout[MAX_SOI_LAY] = {}, z = 0.0, Dc = 0.0;
  for (int l = 0; l <= tb; l++) {
    z += std::fmax(0.0, b.dz[l]); const double dl = std::fmax(0.0, f.drain_l[l]); Dc += dl;
    const double qv = (l < tb) ? std::fmax(0.0, (qi - D) * (1.0 - z / Z) + (D - Dc)) : 0.0;
    const bool mn = is_mineral(b.type[l]);
    const double M = mn ? (1.0 - std::fmin(1.0, std::fmax(0.0, f.por[l]))) * std::fmax(0.0, b.dz[l]) * 2650.0
                        : std::fmax(0.0, b.soil[l]) * 2.0e-3;
    const double den = f.liq[l] + (mn ? kdm : kdo) * M;
    if (!(den > 0.0)) { continue; }
    double fv = qv / den, fl = dl / den; const double ft = fv + fl;
    if (ft > 0.9) { fv *= 0.9 / ft; fl *= 0.9 / ft; }
    const double h = std::fmax(0.0, s.soil[l]), m = std::fmax(0.0, s.mehg[l]);
    out[l] = fv * h; lat[l] = fl * h; mout[l] = (fv + fl) * m;
  }
  for (int l = 0; l <= tb; l++) {
    s.soil[l] -= out[l] + lat[l];
    const double mfrac = (out[l] + lat[l] > 0.0) ? out[l] / (out[l] + lat[l]) : 0.0;
    s.mehg[l] -= mout[l];
    acc.export_ += lat[l]; acc.drain += lat[l];
    if (l < tb) {
      s.soil[l + 1] += out[l]; s.mehg[l + 1] += mout[l] * mfrac;
      if (!is_mineral(b.type[l]) && is_mineral(b.type[l + 1])) { acc.down += out[l]; }
    }
    if (s.mehg[l] > s.soil[l]) { s.mehg[l] = std::fmax(0.0, s.soil[l]); }
    if (s.mehg[l] < 0.0) { s.mehg[l] = 0.0; }
  }
}

// T2 phase 2: daily sorption, overland / lateral / vertical transport of the dissolved Hg, then monthly mixing (see top).
inline void diss_transport_(const Params& p, State& s, const Snap& b, const Flux& f, Acc& acc) {
  const int nb = std::min(b.nsl, MAX_SOI_LAY);
  if (nb <= 0) { return; }
  const double kdo = std::pow(10.0, p.t2_log_kd_org), kdm = std::pow(10.0, p.log_kd_min);
  const double kdh = (p.t2_kdsplit == 1) ? std::pow(10.0, p.t2_log_kd_org + p.t2_dlog_kd_hum) : kdo;   // phase 5b: humic Kd
  const double ks = 1.0 - std::exp(-1.0 / p.tau_sorb);
  // phase 5b (hg_t2_xpath): first mineral layer (receives the flow-path retention), share reaching the stream
  int lm = -1; for (int l = 0; l < nb; l++) { if (is_mineral(b.type[l])) { lm = l; break; } }
  const double fcx = (p.t2_xpath == 1 && lm > 0) ? p.f_conn : 1.0;
  double M[MAX_SOI_LAY] = {}, ztop[MAX_SOI_LAY] = {};
  { double z = 0.0;
    for (int l = 0; l < nb; l++) {
      const bool mn = is_mineral(b.type[l]);
      M[l] = mn ? (1.0 - std::fmin(1.0, std::fmax(0.0, f.por[l]))) * std::fmax(0.0, b.dz[l]) * 2650.0
                : (is_moss(b.type[l]) ? std::fmax(0.0, b.dz[l]) * p.rho_moss : std::fmax(0.0, b.soil[l]) * 2.0e-3);
      ztop[l] = z; z += std::fmax(0.0, b.dz[l]);
    } }
  const int nd = std::min(31, std::max(0, f.nd));
  // T2 phase 3 (hg_t2_meltroute): the month's snowmelt Hg arrives with the day's melt; the overland share of the
  // day's surface water (qover / (qover + qinfl)) is exported, the rest enters the top layer's dissolved pool
  const double mp0 = s.mpend; double msum = 0.0;
  for (int d = 0; d < nd; d++) { msum += std::fmax(0.0, f.d_melt[d]); }
  for (int d = 0; d < nd; d++) {
    if (mp0 > 0.0) {
      double x = (msum > 0.0) ? mp0 * std::fmax(0.0, f.d_melt[d]) / msum : ((d == 0) ? mp0 : 0.0);
      if (d == nd - 1 || x > s.mpend) { x = s.mpend; }
      s.mpend -= x;
      const double qo = std::fmax(0.0, f.dq_over[d]), qi = std::fmax(0.0, f.dq_in[d]);
      const double xo = (qo + qi > 0.0) ? x * qo / (qo + qi) : 0.0;
      const double xr = xo * (1.0 - fcx), xe = xo - xr;   // phase 5b: flow-path retention of the melt overland share
      s.diss[0] += x - xo; acc.export_ += xe; acc.runoff += xe; acc.meltrun += xe;
      if (xr > 0.0) { s.diss[lm] += xr; acc.down += xr; }
    }
    // sorption in every layer toward the Kd equilibrium of the day's liquid water
    for (int l = 0; l < nb; l++) {
      const double liq = std::fmax(0.0, f.d_liq[d][l]);
      const double den = liq + (is_mineral(b.type[l]) ? kdm : ((b.type[l] == 2) ? kdh : kdo)) * M[l];
      const double phi = (den > 0.0) ? liq / den : 0.0;
      const double dx = (phi * (s.soil[l] + s.diss[l]) - s.diss[l]) * ks;
      s.diss[l] += dx; s.soil[l] -= dx;
      if (s.soil[l] < 0.0) { s.diss[l] += s.soil[l]; s.soil[l] = 0.0; }
      if (s.diss[l] < 0.0) { s.soil[l] += s.diss[l]; s.diss[l] = 0.0; }
    }
    // thawed zone from the top
    int tb = -1; double Z = 0.0, D = 0.0;
    for (int l = 0; l < nb; l++) {
      if (f.d_ice[d][l] > 0.5 && f.d_ice[d][l] > f.d_liq[d][l]) { break; }
      tb = l; Z += std::fmax(0.0, b.dz[l]); D += std::fmax(0.0, f.d_drain[d][l]);
    }
    if (tb < 0 || !(Z > 0.0)) { continue; }
    // overland flow: mixing zone = thawed layers starting above z_runoff (at least the top layer)
    const double qo = std::fmax(0.0, f.dq_over[d]);
    double fo = 0.0; int nz = 0;
    if (qo > 0.0) {
      double lz = 0.0;
      for (int l = 0; l <= tb; l++) { if (l == 0 || ztop[l] < p.z_runoff) { lz += std::fmax(0.0, f.d_liq[d][l]); nz = l + 1; } }
      fo = qo / (qo + lz);
    }
    // lateral drainage and vertical advection (T1c water profile); all fractions act on the same diss
    const double qi = std::fmax(0.0, f.dq_in[d]);
    double ov[MAX_SOI_LAY] = {}, out[MAX_SOI_LAY] = {}, lat[MAX_SOI_LAY] = {}, z = 0.0, Dc = 0.0;
    for (int l = 0; l <= tb; l++) {
      z += std::fmax(0.0, b.dz[l]); const double dl = std::fmax(0.0, f.d_drain[d][l]); Dc += dl;
      const double qv = (l < tb) ? std::fmax(0.0, (qi - D) * (1.0 - z / Z) + (D - Dc)) : 0.0;
      const double liq = std::fmax(0.0, f.d_liq[d][l]);
      double fv = 0.0, fl = 0.0, fr = (l < nz) ? fo : 0.0;
      if (liq > 1.0e-9) { fv = p.v_adv * qv / liq; fl = dl / liq; }
      const double ft = fv + fl + fr;
      if (ft > 0.9) { fv *= 0.9 / ft; fl *= 0.9 / ft; fr *= 0.9 / ft; }
      ov[l] = fr * s.diss[l]; out[l] = fv * s.diss[l]; lat[l] = fl * s.diss[l];
    }
    for (int l = 0; l <= tb; l++) {
      s.diss[l] -= ov[l] + out[l] + lat[l];
      if (fcx < 1.0 && l < lm) {   // phase 5b: organic-layer export, the share 1 - f_conn retained on the flow path
        const double ro = ov[l] * (1.0 - fcx), rl = lat[l] * (1.0 - fcx);
        ov[l] -= ro; lat[l] -= rl; s.diss[lm] += ro + rl; acc.down += ro + rl;
      }
      acc.export_ += ov[l] + lat[l]; acc.runoff += ov[l]; acc.drain += lat[l];
      if (l < tb) { s.diss[l + 1] += out[l]; if (!is_mineral(b.type[l]) && is_mineral(b.type[l + 1])) { acc.down += out[l]; } }
    }
  }
  // cryo / bioturbation of the sorbed Hg in the thawed zone of the month (monthly mean liq / ice, as T1c)
  if (p.d_cryo > 0.0 || p.d_bio > 0.0) {
    int tb = -1;
    for (int l = 0; l < nb; l++) { if (f.ice[l] > 0.5 && f.ice[l] > f.liq[l]) { break; } tb = l; }
    double fl[MAX_SOI_LAY] = {};
    for (int l = 0; l < tb; l++) {
      const double d0 = std::fmax(0.0, b.dz[l]), d1 = std::fmax(0.0, b.dz[l + 1]);
      if (!(d0 > 0.0 && d1 > 0.0)) { continue; }
      const double zi = ztop[l] + d0;
      const double Dz = p.d_cryo * std::fmax(0.0, 1.0 - zi / p.z_cryo) + p.d_bio;
      double k = Dz / 12.0 / (0.5 * (d0 + d1));
      k = std::fmin(k, 0.25 * std::fmin(d0, d1));
      fl[l] = k * (s.soil[l] / d0 - s.soil[l + 1] / d1);   // > 0: downward
    }
    for (int l = 0; l < tb; l++) {
      s.soil[l] -= fl[l]; s.soil[l + 1] += fl[l];
      if (!is_mineral(b.type[l]) && is_mineral(b.type[l + 1])) { acc.mixdown += fl[l]; }
    }
  }
  if (s.mpend != 0.0) { s.diss[0] += s.mpend; s.mpend = 0.0; }   // no days (nd = 0): to the top layer
  for (int l = 0; l < nb; l++) { if (s.mehg[l] > s.soil[l] + s.diss[l]) { s.mehg[l] = std::fmax(0.0, s.soil[l] + s.diss[l]); } }
}

// T2: 24-h mean cosine of the solar zenith angle (night = 0) at latitude lat (deg) and solar declination of
// day-of-year doy (Cooper 1969): (h0 sin(phi) sin(delta) + cos(phi) cos(delta) sin(h0)) / pi, h0 = sunset hour angle.
inline double cosz24(double lat_deg, double doy) {
  const double pi = 3.14159265358979323846;
  const double phi = lat_deg * pi / 180.0, dl = 23.45 * pi / 180.0 * std::sin(2.0 * pi * (284.0 + doy) / 365.0);
  const double x = -std::tan(phi) * std::tan(dl);
  const double h0 = (x >= 1.0) ? 0.0 : (x <= -1.0) ? pi : std::acos(x);
  return std::fmax(0.0, (h0 * std::sin(phi) * std::sin(dl) + std::cos(phi) * std::cos(dl) * std::sin(h0)) / pi);
}
// mid-month day of year (month 0-11) and days in the month (no leap years, as the model)
inline double mid_doy(int m0) { static const int c[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334}; static const int d[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; const int m = std::min(11, std::max(0, m0)); return c[m] + 0.5 * d[m] + 0.5; }
inline double days_in(int m0) { static const int d[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; return d[std::min(11, std::max(0, m0))]; }
// ORCHIDEE-Hg fixed fraction of dry deposition taken up by leaves (stomate_Hg.f90 l.1046-1052)
inline double f_fix(const Params& p, double lai) { return (lai > 1.0e-12) ? p.fol_fmax * std::fmin(1.0, lai / p.fol_lai) : 0.0; }

// T2 phase 3 (D3, hg_t2_foliar = 2): stomatal Hg(0) uptake per g C of GPP (ug/gC). Medlyn: g_sw = 1.6 m A / c_a
// (mol m-2 s-1), integrated over the month with A = GPP / 12; Hg(0) conductance g_sw / sqrt(M_Hg / M_H2O)
// (Wesely 1989 diffusivity scaling, 3.34); m3 via the molar air density P / (R T); x GEM (ng/m3) x veg_scale x hg_gs_eff.
// CLM5-Hg (Yuan et al. 2023) uses the same stomatal resistance in its Wesely dry deposition; Ra, Rb and the
// cuticular path are lumped into hg_gs_eff.
inline double stomatal_k(const Params& p, const Flux& f) {
  const double ca = ((f.co2 > 1.0) ? f.co2 : 400.0) * 1.0e-6;
  const double cair = 101325.0 / (8.314 * (std::fmax(-40.0, std::fmin(45.0, f.ta)) + 273.15));
  return p.gs_eff * p.gem * f.veg_scale * 1.6 * p.gs_m / (12.0 * ca * 3.3365 * cair) * 1.0e-3;
}
// T2 phase 3: ground Hg(0) sink of the month (ug/m2) = v_ground x GEM x veg_scale x seconds
inline double ground_sink(const Params& p, const Flux& f) { return p.v_ground * p.gem * f.veg_scale * std::fmax(0.0, f.ndays) * 86.4; }
// T2 phase 5b (hg_t2_snowsink): the snow-covered fraction fs takes up Hg(0) at hg_v_snow
inline double ground_sink_snow(const Params& p, const Flux& f, double fs) {
  return ((1.0 - fs) * p.v_ground + fs * p.v_snow) * p.gem * f.veg_scale * std::fmax(0.0, f.ndays) * 86.4;
}
// T2 phase 3: re-emitted fraction of the Hg released by decay in a layer centred at depth zm (m) of temperature ts (C)
inline double fev_layer(const Params& p, double fsnow, double zm, double ts) {
  return (ts > 0.0) ? p.fd * (1.0 - fsnow) * std::exp(-std::fmax(0.0, zm) / p.ev_z) : 0.0;
}

// T2 phase-1 surface part of a BGC step (f.t2 = 1): foliar uptake per PFT (returned in fol, ug/m2), snowpack, ground
// deposition; returns the Hg to add to the top C soil layer (deposition + melt release). Linear in Hg.
inline double surface_t2_(const Params& p, State& s, const Flux& f, Acc& acc, double fol[NUM_PFT], bool have_soil) {
  for (int i = 0; i < NUM_PFT; i++) { fol[i] = 0.0; }
  const double fs = (p.t2_snow == 1) ? std::fmin(1.0, std::fmax(0.0, f.swe) / p.swe_full) : 0.0;
  const double wet = std::fmax(0.0, f.dep_wet), d2 = std::fmax(0.0, f.dep_dry2), d0 = std::fmax(0.0, f.dep_dry0);
  // foliar uptake
  double fc = 0.0;   // canopy fraction that fixes dry deposition
  // T2 phase 4 (hg_wash): wash-off share of the intercepted dry Hg(II) (to the ground with the throughfall)
  const double fw = (p.wash > 0.0 && f.rain > 0.0) ? p.wash * f.rain / (f.rain + p.wash_pref) : 0.0;
  double washed = 0.0;
  if (p.t2_foliar == 1) {
    double w[NUM_PFT] = {}; double tw = 0.0;
    for (int i = 0; i < NUM_PFT; i++) {
      if (f.vegcov[i] > 0.0) { w[i] = p.fol_scale * f_fix(p, f.lai[i]) * f.vegcov[i] * (1.0 - fs); tw += w[i]; }
    }
    const double sc = (tw > 1.0) ? 1.0 / tw : 1.0;
    for (int i = 0; i < NUM_PFT; i++) {
      const double wo = w[i] * sc * d2 * fw; washed += wo;
      fol[i] = w[i] * sc * (d0 + d2) - wo; fc += w[i] * sc; acc.foliar += fol[i];   // booked as uptake in the leaf loop
    }
  } else if (p.t2_foliar == 2) {   // T2 phase 3 (D3): stomatal Hg(0) uptake from GPP + LAI interception of dry Hg(II)
    double w[NUM_PFT] = {}; double tw = 0.0;
    for (int i = 0; i < NUM_PFT; i++) {
      if (f.vegcov[i] > 0.0) { w[i] = p.fol_scale * f_fix(p, f.lai[i]) * f.vegcov[i] * (1.0 - fs); tw += w[i]; }
    }
    const double sc = (tw > 1.0) ? 1.0 / tw : 1.0;
    const double k = stomatal_k(p, f);   // ug Hg per g C of GPP
    for (int i = 0; i < NUM_PFT; i++) {
      const double st = k * std::fmax(0.0, f.gpp[i]);
      const double wo = w[i] * sc * d2 * fw; washed += wo;
      fol[i] = w[i] * sc * d2 - wo + st; fc += w[i] * sc; acc.foliar += fol[i]; acc.stomatal += st;
    }
  }
  // Hg(II) to the ground / snow (dry Hg(II) intercepted by the canopy removed); soil Hg(0) sink to the soil
  const double g = wet + d2 * (1.0 - fc) + washed;
  acc.wash += have_soil ? washed : washed * fs;   // booked in dep with the rest of g
  const double gsnow = g * fs, gground = g - gsnow;
  const double fpr = (p.t2_photored == 1) ? std::fmin(1.0, std::fmax(0.0, p.pr_a * std::exp(-p.pr_b * std::fmax(0.0, f.sw)))) : 0.0;
  const double pr = gground * fpr;
  double to_soil = 0.0;
  const double s0 = (p.t2_evasion == 1) ? ((p.t2_snowsink == 1) ? ground_sink_snow(p, f, fs) : ground_sink(p, f))
                                        : std::fmax(0.0, f.dep_soil0);   // phase 3: process ground Hg(0) sink (5b: under snow)
  if (have_soil) {
    acc.dep += g + s0; if (p.t2_evasion == 1) { acc.sink0 += s0; }
    acc.photored += pr; acc.evasion += pr;
    to_soil = gground - pr + s0;
  } else {   // no soil column: only the snow part is kept (as T1, deposition needs a soil layer)
    acc.dep += gsnow;
  }
  // snowpack
  s.snow += gsnow; acc.snowdep += gsnow;
  if (s.snow > 0.0) {
    const double swe = std::fmax(0.0, f.swe);
    const double zs = swe / p.snow_rho;   // m (swe mm = kg/m2)
    const double fph = (zs > 0.0) ? std::fmin(1.0, p.snow_photic / zs) : 1.0;
    const double red = s.snow * fph * (1.0 - std::exp(-p.snow_kred * 24.0 * std::fmax(0.0, f.cosz) * std::fmax(0.0, f.ndays)));
    s.snow -= red; acc.snowred += red; acc.evasion += red;
    const double mlt = std::fmax(0.0, f.melt), sub = std::fmax(0.0, f.sublim);
    const double sub_out = (p.snow_sublim_loss == 1) ? sub : 0.0;   // D7: 0 = sublimation concentrates Hg in the pack
    double rel = (swe < 0.5) ? 1.0 : (mlt + sub_out) / (mlt + sub_out + swe);
    if (!have_soil) { rel = 0.0; }
    if (rel > 0.0) {
      const double out = s.snow * rel;
      const double fsub = (mlt + sub_out > 0.0) ? sub_out / (mlt + sub_out) : 0.0;
      const double xs = out * fsub, xm = out - xs;
      s.snow -= out; acc.evasion += xs; acc.snowred += xs; acc.snowmelt += xm;
      if (p.t2_meltroute == 1 && p.t2_diss == 1) { s.mpend += xm; } else { to_soil += xm; }   // phase 3: delivered daily with the melt
    }
    if (s.snow < 0.0) { s.snow = 0.0; }
  }
  return to_soil;
}

// One sub-step. a, b: C before / after; f: the sub-step's external C fluxes. Updates s (and s.m, s.cum).
inline void step(const Params& p, State& s, const Snap& a, const Snap& b, const Flux& f) {
  if (!s.init) { return; }
  Acc acc;
  const bool td = (p.t2 == 1 && p.t2_diss == 1);   // T2 phase 2: explicit dissolved pool (export explicit, no f_export)
  const double f_exp = td ? 0.0 : p.f_export;
  const double rv[NUM_PFT_PART] = {p.r_leaf, p.r_stem, p.r_root};
  const int na = std::min(a.nsl, MAX_SOI_LAY), nb = std::min(b.nsl, MAX_SOI_LAY);
  // T2 phase 3 (hg_t2_evasion, BGC step): re-emitted fraction of decay-released Hg by layer, woody debris, the pot
  const bool ev3 = (f.t2 == 1 && p.t2 == 1 && p.t2_evasion == 1);
  double fev[MAX_SOI_LAY]; double fev_wd = p.f_evade;
  for (int l = 0; l < MAX_SOI_LAY; l++) { fev[l] = p.f_evade; }
  if (ev3) {
    const double fsn = (p.t2_snow == 1) ? std::fmin(1.0, std::fmax(0.0, f.swe) / p.swe_full) : 0.0;
    double z = 0.0;
    for (int l = 0; l < na; l++) { const double d = std::fmax(0.0, a.dz[l]); fev[l] = fev_layer(p, fsn, z + 0.5 * d, f.ts[l]); z += d; }
    fev_wd = fev_layer(p, fsn, 0.0, f.ta);
  }
  const bool same = same_geometry(a, b);
  if (td && !same) {   // geometry change: the dissolved Hg re-sorbs (follows the C re-map)
    for (int l = 0; l < MAX_SOI_LAY; l++) { if (s.diss[l] != 0.0) { s.soil[std::min(l, std::max(0, na - 1))] += s.diss[l]; s.diss[l] = 0.0; } }
  }
  // top soil layer that holds C (the moss layer usually has no SOM C): target of deposition and of
  // unassigned Hg, so that Hg never piles up in a layer without C turnover
  int topc = 0;
  for (int l = 0; l < nb; l++) { if (b.soil[l] > 1.0e-6) { topc = l; break; } }
  // T2 phase 1 (hg_t2): surface interface (foliar uptake, snowpack, deposition components)
  const bool t2 = (f.t2 == 1 && p.t2 == 1);
  double fol[NUM_PFT] = {}; double t2_soil = 0.0;
  if (t2) { t2_soil = surface_t2_(p, s, f, acc, fol, nb > 0); }

  // ---- non-soil pools (flattened): veg parts, standing dead, woody debris
  double* hp[NPOOL]; double dint[NPOOL]; double r[NPOOL];
  int k = 0;
  for (int i = 0; i < NUM_PFT; i++) {
    for (int j = 0; j < NUM_PFT_PART; j++, k++) {
      hp[k] = &s.veg[i][j];
      const double in = f.npp[i][j];
      const double hin = (t2 && p.t2_foliar >= 1 && j == 0) ? fol[i] : std::fmax(0.0, in) * rv[j] * f.veg_scale;
      *hp[k] += hin; acc.uptake += hin;
      r[k] = ratio_(*hp[k], a.veg[i][j] + std::fmax(0.0, in));
      dint[k] = b.veg[i][j] - a.veg[i][j] - in;
    }
  }
  for (int i = 0; i < NUM_PFT; i++, k++) { hp[k] = &s.dead[i]; r[k] = ratio_(s.dead[i], a.dead[i]); dint[k] = b.dead[i] - a.dead[i]; }
  {
    hp[k] = &s.wd; r[k] = ratio_(s.wd, a.wd);
    const double rel = std::fmin(s.wd, std::fmax(0.0, f.decomp_wd) * r[k]);
    acc.evasion += rel * fev_wd; acc.export_ += rel * f_exp; s.wd -= rel * (fev_wd + f_exp);
    if (ev3) { acc.soilev += rel * fev_wd; }
    dint[k] = b.wd - a.wd + f.decomp_wd;
    k++;
  }

  // ---- soil
  double sd[MAX_SOI_LAY] = {};      // same geometry: dC_int by layer
  double rs[MAX_SOI_LAY] = {};      // Hg:C at the start of the step
  double pot_decomp = f.pot_decomp, lateral = f.lateral;
  double fev_pot = p.f_evade;       // T2 phase 3: re-emitted fraction of the pot's decomposition
  double soilnet = 0.0;             // different geometry: net soil C change (after the external fluxes)
  double h0soil[MAX_SOI_LAY];
  for (int l = 0; l < MAX_SOI_LAY; l++) { h0soil[l] = s.soil[l]; rs[l] = (l < na) ? ratio_(s.soil[l], a.soil[l]) : 0.0; }
  double leach_in[MAX_SOI_LAY] = {};   // hg_f_leach: decay-released Hg moved one layer down (added after the MeHg map)
  if (same) {
    for (int l = 0; l < na; l++) {
      double rel = std::fmin(s.soil[l], std::fmax(0.0, f.decomp[l]) * rs[l]);
      const double frl = fev[l] + f_exp;
      acc.evasion += rel * fev[l]; acc.export_ += rel * f_exp; s.soil[l] -= rel * frl;
      if (ev3) { acc.soilev += rel * fev[l]; }
      if (td) { const double dv = std::fmin(s.soil[l], rel * (1.0 - frl)); s.soil[l] -= dv; s.diss[l] += dv; }   // T2: released Hg -> dissolved
      if (p.f_leach > 0.0) {
        const double lx = std::fmin(s.soil[l], rel * p.f_leach); s.soil[l] -= lx;
        if (l + 1 < na) { leach_in[l + 1] += lx; } else { acc.export_ += lx; }
      }
      const double ex = std::fmin(s.soil[l], std::fmax(0.0, f.doc[l]) * rs[l]);
      acc.export_ += ex; s.soil[l] -= ex;
      sd[l] = b.soil[l] - a.soil[l] + f.decomp[l] + f.doc[l];
    }
  } else {
    double c0 = 0.0, c1 = 0.0, dec = 0.0, doc = 0.0;
    for (int l = 0; l < na; l++) { c0 += a.soil[l]; dec += f.decomp[l]; doc += f.doc[l]; }
    for (int l = 0; l < nb; l++) { c1 += b.soil[l]; }
    if (ev3 && dec + std::fmax(0.0, f.pot_decomp) > 0.0) {   // phase 3: decomposition-weighted re-emission of the pot
      double w = std::fmax(0.0, f.pot_decomp) * p.f_evade;
      for (int l = 0; l < na; l++) { w += std::fmax(0.0, f.decomp[l]) * fev[l]; }
      fev_pot = w / (dec + std::fmax(0.0, f.pot_decomp));
    }
    pot_decomp += dec; lateral += doc;          // per-layer fluxes with a geometry change: from the pot
    soilnet = c1 - c0 + dec + doc;
  }
  if (f.dep_soil > 0.0 && nb > 0) { acc.dep += f.dep_soil; }   // added to the top C layer after the remap
  const double dep_add = t2 ? t2_soil : acc.dep;                 // T2: deposition to the ground + snow melt release

  // ---- donors / receivers
  double D = 0.0, R = 0.0;
  for (int i = 0; i < NPOOL; i++) { if (dint[i] < 0) D -= dint[i]; else R += dint[i]; }
  if (same) { for (int l = 0; l < na; l++) { if (sd[l] < 0) D -= sd[l]; else R += sd[l]; } }
  else { if (soilnet < 0) D -= soilnet; else R += soilnet; }
  const double S = std::fmax(0.0, f.burn) + std::fmax(0.0, lateral) + std::fmax(0.0, pot_decomp);
  const double g = (D > 0.0) ? std::fmin(1.0, (R + S) / D) : 0.0;
  const double cflow = D * g;
  acc.unexpl_closs = D - cflow;

  double P = 0.0; double give[NPOOL] = {}; double give_soil[MAX_SOI_LAY] = {}; double give_col = 0.0;
  // T2 phase 4 switches of this step (fire step: f.burn > 0; RTS step: f.lateral / f.pot_decomp > 0)
  const bool ash4 = (p.t2 == 1 && p.fire_ash > 0.0 && f.burn > 0.0 && nb > 0);
  const bool rts4 = (p.t2 == 1 && p.t2_rts == 1 && (f.lateral > 0.0 || f.pot_decomp > 0.0));
  for (int i = 0; i < NPOOL; i++) {
    if (dint[i] < 0) { give[i] = std::fmin(*hp[i], -dint[i] * r[i] * g); *hp[i] -= give[i]; P += give[i]; }
  }
  // MeHg follows Hg losses of a layer proportionally (same geometry); with a geometry change it is mapped.
  double mehg_new[MAX_SOI_LAY] = {};
  if (same) {
    for (int l = 0; l < na; l++) {
      if (sd[l] < 0) { give_soil[l] = std::fmin(s.soil[l], -sd[l] * rs[l] * g); s.soil[l] -= give_soil[l]; P += give_soil[l]; }
      mehg_new[l] = (h0soil[l] > 0.0) ? s.mehg[l] * (s.soil[l] / h0soil[l]) : 0.0;
    }
  } else {
    // old column (top -> bottom), take the net loss from the top
    double cc[MAX_SOI_LAY], hh[MAX_SOI_LAY], mm[MAX_SOI_LAY];
    for (int l = 0; l < na; l++) { cc[l] = std::fmax(0.0, a.soil[l]); hh[l] = s.soil[l]; mm[l] = s.mehg[l]; }
    if (soilnet < 0) {
      double t = -soilnet * g;
      for (int l = 0; l < na && t > 0.0; l++) {
        const double x = std::fmin(cc[l], t);
        const double fx = (cc[l] > 0.0) ? x / cc[l] : 0.0;
        const double h = hh[l] * fx, m = mm[l] * fx;
        cc[l] -= x; hh[l] -= h; mm[l] -= m; t -= x; give_col += h; P += h;
      }
      if (t > 0.0) { acc.unexpl_closs += t; }
    }
    // bottom-up map onto the new layers
    double newh[MAX_SOI_LAY] = {}, newm[MAX_SOI_LAY] = {}, excess[MAX_SOI_LAY] = {};
    int i = na - 1;
    for (int j = nb - 1; j >= 0; j--) {
      double need = std::fmax(0.0, b.soil[j]);
      while (need > 0.0 && i >= 0) {
        if (cc[i] <= 0.0) { newh[j] += hh[i]; newm[j] += mm[i]; hh[i] = 0; mm[i] = 0; cc[i] = 0; i--; continue; }
        const double x = std::fmin(cc[i], need);
        const double fx = x / cc[i];
        const double h = (fx >= 1.0) ? hh[i] : hh[i] * fx, m = (fx >= 1.0) ? mm[i] : mm[i] * fx;
        newh[j] += h; newm[j] += m; hh[i] -= h; mm[i] -= m; cc[i] -= x; need -= x;
        if (fx >= 1.0) { cc[i] = 0; i--; }
      }
      excess[j] = need;
    }
    // leftover old material (C loss with no sink): Hg stays, in the top new layer
    double lh = 0.0, lm = 0.0;
    for (int l = 0; l <= i && l < na; l++) { lh += hh[l]; lm += mm[l]; }
    if (nb > 0) { newh[topc] += lh; newm[topc] += lm; } else { P += lh; give_col += lh; }
    for (int l = 0; l < MAX_SOI_LAY; l++) { s.soil[l] = (l < nb) ? newh[l] : 0.0; mehg_new[l] = (l < nb) ? newm[l] : 0.0; }
    // receivers: the top excess
    R = 0.0;
    for (int ii = 0; ii < NPOOL; ii++) { if (dint[ii] > 0) R += dint[ii]; }
    for (int j = 0; j < nb; j++) { sd[j] = excess[j]; R += excess[j]; }
  }
  for (int l = 0; l < MAX_SOI_LAY; l++) { s.mehg[l] = mehg_new[l]; }
  if (same) { for (int l = 0; l < na; l++) { s.soil[l] += leach_in[l]; } }

  // ---- sinks from the pot
  if (P > 0.0 && cflow > 0.0) {
    const double sc = (S > cflow) ? cflow / S : 1.0;
    const double hb = P * std::fmax(0.0, f.burn) * sc / cflow;
    const double hl = P * std::fmax(0.0, lateral) * sc / cflow;
    const double hd = P * std::fmax(0.0, pot_decomp) * sc / cflow;
    double ash = 0.0;   // phase 4: burned soil Hg kept in the ash (soil donors' share of the pot)
    if (ash4) {
      double gs = give_col, ga = give_col;
      for (int l = 0; l < MAX_SOI_LAY; l++) { gs += give_soil[l]; ga += give_soil[l]; }
      for (int i = 0; i < NPOOL; i++) { ga += give[i]; }
      if (ga > 0.0) { ash = hb * p.fire_ash * gs / ga; }
    }
    acc.fire += hb - ash; acc.ash += ash; acc.evasion += hd * fev_pot;
    if (ev3) { acc.soilev += hd * fev_pot; }
    if (rts4) {   // phase 4: the slumped material leaves as particulate export, with the Hg of its respired C not re-emitted
      const double xp = hl + hd * (1.0 - fev_pot);
      acc.export_ += xp; acc.partic += xp;
      P -= hb + hl + hd;
    } else {
      acc.export_ += hl + hd * f_exp;
      P -= hb + hl + hd * (fev_pot + f_exp);
    }
    if (P < 0.0) { P = 0.0; }
    if (ash > 0.0) { s.soil[topc] += ash; }
  }
  // ---- receivers
  if (P > 0.0) {
    if (R > 0.0) {
      double left = P;
      int last = -1;
      for (int i = 0; i < NPOOL; i++) { if (dint[i] > 0) last = i; }
      int lastl = -1;
      const int ns = same ? na : nb;
      for (int l = 0; l < ns; l++) { if (sd[l] > 0) lastl = l; }
      for (int i = 0; i < NPOOL; i++) {
        if (dint[i] > 0) { const double h = P * dint[i] / R; *hp[i] += h; left -= h; }
      }
      for (int l = 0; l < ns; l++) {
        if (sd[l] > 0) { const double h = P * sd[l] / R; s.soil[l] += h; left -= h; }
      }
      // roundoff remainder to the last receiver (exact conservation)
      if (lastl >= 0) { s.soil[lastl] += left; } else if (last >= 0) { *hp[last] += left; }
    } else {
      // nothing received C: the remaining pot goes back to the donors
      double tg = 0.0; for (int i = 0; i < NPOOL; i++) tg += give[i];
      for (int l = 0; l < MAX_SOI_LAY; l++) tg += give_soil[l];
      tg += give_col;
      double left = P;
      if (tg > 0.0) {
        for (int i = 0; i < NPOOL; i++) { if (give[i] > 0) { const double h = P * give[i] / tg; *hp[i] += h; left -= h; } }
        for (int l = 0; l < MAX_SOI_LAY; l++) { if (give_soil[l] > 0) { const double h = P * give_soil[l] / tg; s.soil[l] += h; left -= h; } }
      }
      if (nb > 0) { s.soil[topc] += left; } else { s.wd += left; }
    }
  }
  if (dep_add > 0.0) { s.soil[topc] += dep_add; }
  acc.unexpl_cgain = std::fmax(0.0, R - (cflow - std::fmin(S, cflow)));
  if (f.tr == 1 && p.transport == 1) { if (td) { diss_transport_(p, s, b, f, acc); } else { transport_(p, s, b, f, acc); } }
  s.nsl = nb;
  s.m.add(acc); s.y.add(acc); s.cum.add(acc);
}

// ---- Hg spin-up (hg_spinup = 1)
struct Rec { Snap a, b; Flux f; };   // one recorded tracer step
const int NVEC = NPOOL + 2 * MAX_SOI_LAY + 2;   // + T2 snowpack Hg, + T2 dissolved Hg by layer, + T2 pending melt Hg (0 when off)
inline void to_vec_(const State& s, double* x) {
  int k = 0;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) x[k++] = s.veg[i][j];
  for (int i = 0; i < NUM_PFT; i++) x[k++] = s.dead[i];
  x[k++] = s.wd;
  for (int l = 0; l < MAX_SOI_LAY; l++) x[k++] = s.soil[l];
  x[k++] = s.snow;
  for (int l = 0; l < MAX_SOI_LAY; l++) x[k++] = s.diss[l];
  x[k++] = s.mpend;
}
inline void from_vec_(State& s, const double* x) {
  int k = 0;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) s.veg[i][j] = x[k++];
  for (int i = 0; i < NUM_PFT; i++) s.dead[i] = x[k++];
  s.wd = x[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.soil[l] = x[k++];
  s.snow = x[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.diss[l] = x[k++];
  s.mpend = x[k++];
}
// Affine map x -> M x + b (row-major M, NVEC x NVEC).
struct Aff { std::vector<double> M, b; };
inline Aff compose_(const Aff& B, const Aff& A) {   // B o A
  Aff C; C.M.assign(NVEC * NVEC, 0.0); C.b.assign(NVEC, 0.0);
  for (int i = 0; i < NVEC; i++) {
    for (int k = 0; k < NVEC; k++) {
      const double bik = B.M[i * NVEC + k];
      if (bik == 0.0) continue;
      for (int j = 0; j < NVEC; j++) C.M[i * NVEC + j] += bik * A.M[k * NVEC + j];
      C.b[i] += bik * A.b[k];
    }
    C.b[i] += B.b[i];
  }
  return C;
}
// Apply hg_spinup_years of the recorded cycle to s (s = the state at the end of the cycle, geometry
// rec.back().b). T1d: geometry closure. The C model is not exactly periodic over the recorded years (layer
// C drifts, layers are re-divided, the layer count may differ at the two ends). Repeating the raw record
// would book the C drift of one cycle again in every repetition (layers that gained C over the cycle keep
// receiving Hg, the donors keep losing it), which biased the equilibrium (T1c: TVC +13 vs -12.5 mg m-2 for
// 10- vs 1-yr cycles). The closure is one C-neutral step from the end snapshot back to the start snapshot
// (the same step function: internal C transfers at the layers' Hg:C, bottom-up re-map if the geometry
// differs), so that the repeated map  start -> cycle -> end -> closure -> start  is periodic in C and
// geometry and its fixed point is the Hg equilibrium for the recorded C state. The jump is
//   s -> closure -> (cycle + closure)^(n-1) -> cycle,   n = round(spinup_years / cycle),
// i.e. it ends in the current geometry. If the record is periodic (end snapshot == start snapshot) the
// closure step changes nothing and the jump is n plain repetitions. Exactly conservative except for the
// closure's dilution / concentration of unexplained C (no Hg created or lost: Hg is only moved). Returns false
// (s unchanged) if s is not at the end geometry of the record or the result is not finite.
inline bool snap_equal_(const Snap& a, const Snap& b) {
  if (a.nsl != b.nsl) { return false; }
  for (int l = 0; l < MAX_SOI_LAY; l++) { if (a.soil[l] != b.soil[l] || a.dz[l] != b.dz[l] || a.type[l] != b.type[l]) { return false; } }
  for (int i = 0; i < NUM_PFT; i++) { for (int j = 0; j < NUM_PFT_PART; j++) { if (a.veg[i][j] != b.veg[i][j]) { return false; } } if (a.dead[i] != b.dead[i]) { return false; } }
  return a.wd == b.wd;
}
// Closure mode 0: soil Hg / MeHg of the layers of geometry `from` are moved to the layers of `to` within each soil
// horizon: organic layers (moss, fibric, humic) by overlap in relative depth z / Z_org of the organic horizon,
// mineral layers by overlap in depth below the organic-mineral interface (below the new column -> its bottom
// layer). A horizon missing in `to` sends its Hg to the nearest layer of the other horizon; zero-thickness layers
// keep their index. Linear, exactly conservative, the identity for the same geometry. T1d: plain depth overlap
// moved mineral Hg up into a thinning organic horizon in every repetition (HPC lost its mineral Hg).
inline void depth_map_(State& t, const Snap& from, const Snap& to) {
  const int nf = std::min(from.nsl, MAX_SOI_LAY), nt = std::min(to.nsl, MAX_SOI_LAY);
  if (nt <= 0) { return; }
  if (same_geometry(from, to)) {   // exact identity (same layer count and thicknesses)
    bool st = true; for (int l = 0; l < nt; l++) { if (from.type[l] != to.type[l]) { st = false; } }
    if (st) { t.nsl = to.nsl; return; }
  }
  for (int l = 0; l < MAX_SOI_LAY; l++) { if (t.diss[l] != 0.0) { t.soil[l] += t.diss[l]; t.diss[l] = 0.0; } }   // T2: re-sorb before the map
  auto coords = [](const Snap& g, int n, double* c0, double* c1, bool* org, int& no, int& nm) {
    double zo = 0.0; no = 0; nm = 0;
    for (int l = 0; l < n; l++) { org[l] = !is_mineral(g.type[l]); if (org[l]) { zo += std::fmax(0.0, g.dz[l]); no++; } else { nm++; } }
    double z = 0.0, zm = 0.0;
    for (int l = 0; l < n; l++) {
      const double d = std::fmax(0.0, g.dz[l]);
      if (org[l]) { c0[l] = (zo > 0.0) ? z / zo : 0.0; z += d; c1[l] = (zo > 0.0) ? z / zo : 0.0; }
      else { c0[l] = zm; zm += d; c1[l] = zm; }
    }
  };
  double f0[MAX_SOI_LAY], f1[MAX_SOI_LAY], t0[MAX_SOI_LAY], t1[MAX_SOI_LAY]; bool fo[MAX_SOI_LAY], to_[MAX_SOI_LAY];
  int nfo, nfm, nto, ntm;
  coords(from, nf, f0, f1, fo, nfo, nfm); coords(to, nt, t0, t1, to_, nto, ntm);
  int topo = -1, boto = -1, topm = -1, botm = -1;
  for (int l = 0; l < nt; l++) { if (to_[l]) { if (topo < 0) topo = l; boto = l; } else { if (topm < 0) topm = l; botm = l; } }
  double h[MAX_SOI_LAY] = {}, m[MAX_SOI_LAY] = {};
  for (int k = 0; k < MAX_SOI_LAY; k++) {
    const double hk = t.soil[k], mk = t.mehg[k];
    if (hk == 0.0 && mk == 0.0) { continue; }
    if (k >= nf) { h[nt - 1] += hk; m[nt - 1] += mk; continue; }
    const bool og = fo[k];
    const int bot = og ? (boto >= 0 ? boto : topm) : (botm >= 0 ? botm : boto);
    if ((og && topo < 0) || (!og && topm < 0)) { h[bot] += hk; m[bot] += mk; continue; }   // horizon missing in `to`
    const double dk = f1[k] - f0[k];
    if (!(dk > 0.0)) { const int l = std::min(k, nt - 1); h[l] += hk; m[l] += mk; continue; }
    double left = 1.0;
    for (int l = 0; l < nt && left > 0.0; l++) {
      if (to_[l] != og) { continue; }
      const double ov = std::fmin(f1[k], t1[l]) - std::fmax(f0[k], t0[l]);
      if (ov <= 0.0) { continue; }
      const double fr = std::fmin(left, ov / dk); h[l] += hk * fr; m[l] += mk * fr; left -= fr;
    }
    if (left > 0.0) { h[bot] += hk * left; m[bot] += mk * left; }
  }
  for (int l = 0; l < MAX_SOI_LAY; l++) { t.soil[l] = h[l]; t.mehg[l] = m[l]; }
  t.nsl = to.nsl;
}
inline bool spinup_jump(const Params& p, State& s, const std::vector<Rec>& rec, int /*nsl_start*/) {
  if (rec.empty() || s.nsl != rec.back().b.nsl) { return false; }
  const Snap& a0 = rec.front().a; const Snap& b1 = rec.back().b;
  if (p.spinup_closure == 2 && a0.nsl != b1.nsl) { return false; }
  const bool periodic = snap_equal_(a0, b1);
  const Flux none;
  auto closure = [&](State& t) {
    if (periodic || p.spinup_closure == 2) { t.nsl = a0.nsl; return; }
    if (p.spinup_closure == 1) { step(p, t, b1, a0, none); } else { depth_map_(t, b1, a0); }
  };
  auto pools = [&](const State& src, State& dst) { std::vector<double> x(NVEC); to_vec_(src, &x[0]); from_vec_(dst, &x[0]); dst.nsl = src.nsl; };
  auto run = [&](const double* x0, double* x1) {   // start geometry -> cycle -> closure -> start geometry
    State t; t.init = true; t.nsl = a0.nsl; from_vec_(t, x0);
    for (const Rec& r : rec) { step(p, t, r.a, r.b, r.f); }
    closure(t);
    to_vec_(t, x1);
  };
  long n = std::lround(p.spinup_years / std::max(1, p.spinup_cycle));
  if (n < 1) n = 1;
  // s (end geometry) -> start geometry
  State t0s; t0s.init = true; pools(s, t0s); closure(t0s); t0s.nsl = a0.nsl;
  std::vector<double> x(NVEC), xn(NVEC, 0.0);
  to_vec_(t0s, &x[0]);
  if (n > 1) {
    Aff A; A.M.assign(NVEC * NVEC, 0.0); A.b.assign(NVEC, 0.0);
    std::vector<double> z(NVEC, 0.0), y(NVEC);
    run(&z[0], &A.b[0]);
    for (int j = 0; j < NVEC; j++) {
      std::vector<double> e(NVEC, 0.0); e[j] = 1.0;
      run(&e[0], &y[0]);
      for (int i = 0; i < NVEC; i++) A.M[i * NVEC + j] = y[i] - A.b[i];
    }
    long m = n - 1;
    Aff R; R.M.assign(NVEC * NVEC, 0.0); R.b.assign(NVEC, 0.0);
    for (int i = 0; i < NVEC; i++) R.M[i * NVEC + i] = 1.0;
    Aff B = A;
    while (m > 0) { if (m & 1) R = compose_(B, R); m >>= 1; if (m > 0) B = compose_(B, B); }
    for (int i = 0; i < NVEC; i++) {
      double v = R.b[i];
      for (int j = 0; j < NVEC; j++) v += R.M[i * NVEC + j] * x[j];
      if (!std::isfinite(v)) { return false; }
      xn[i] = std::fmax(0.0, v);
    }
  } else { xn = x; }
  // last cycle (no closure): back to the end geometry
  State t; t.init = true; t.nsl = a0.nsl; from_vec_(t, &xn[0]);
  for (const Rec& r : rec) { step(p, t, r.a, r.b, r.f); }
  to_vec_(t, &xn[0]);
  for (int i = 0; i < NVEC; i++) { if (!std::isfinite(xn[i])) { return false; } xn[i] = std::fmax(0.0, xn[i]); }
  double h0[MAX_SOI_LAY]; for (int l = 0; l < MAX_SOI_LAY; l++) h0[l] = s.soil[l];
  const double tt0 = s.total();
  from_vec_(s, &xn[0]);
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    s.mehg[l] = (h0[l] > 0.0) ? std::fmin(s.soil[l], s.mehg[l] * s.soil[l] / h0[l]) : 0.0;
  }
  const double d = s.total() - tt0;
  s.m.spin += d; s.y.spin += d; s.cum.spin += d;
  return true;
}
// Relative soil C drift of a record (|C_end - C_start| / C_start; diagnostic of how far the C cycle is from periodic).
inline double rec_c_drift(const std::vector<Rec>& rec) {
  if (rec.empty()) { return 0.0; }
  double c0 = 0.0, c1 = 0.0;
  for (int l = 0; l < MAX_SOI_LAY; l++) { c0 += std::fmax(0.0, rec.front().a.soil[l]); c1 += std::fmax(0.0, rec.back().b.soil[l]); }
  return (c0 > 0.0) ? std::fabs(c1 - c0) / c0 : 0.0;
}

// MeHg index for one month. a_turn[l]: anaerobic C turnover fraction of the month (-); ts[l]: layer T (C).
inline void mehg_month(const Params& p, State& s, const double* a_turn, const double* ts, int nsl, bool ch4_driver) {
  const double km = ch4_driver ? p.k_meth_ch4 : p.k_meth_proxy;
  double prod = 0.0;
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    const double hl = s.soil[l] + s.diss[l];   // T2: layer Hg = sorbed + dissolved (diss = 0 without hg_t2_diss)
    if (l >= nsl || !(hl > 0.0)) { s.mehg[l] = 0.0; continue; }
    const double t = std::isfinite(ts[l]) ? ts[l] : p.t_ref;
    const double pr = km * hl * std::fmax(0.0, a_turn[l]) * std::pow(p.q10_meth, (t - p.t_ref) / 10.0);
    const double kd = p.k_demeth * std::pow(p.q10_demeth, (t - p.t_ref) / 10.0);
    const double meq = pr / kd;
    s.mehg[l] = meq + (s.mehg[l] - meq) * std::exp(-kd);
    s.mehg[l] = std::fmin(std::fmax(s.mehg[l], 0.0), hl);
    prod += pr;
  }
  s.mehg_prod_m += prod; s.mehg_prod_y += prod;
}

// Restart layout (version 4): [0] version, veg (NUM_PFT x 3), dead (NUM_PFT), wd, soil (MAX_SOI_LAY),
// mehg (MAX_SOI_LAY), h_init, cum uptake, dep, evasion, export, fire, spin, spun (diagnostic continuity), snow (T2),
// diss (MAX_SOI_LAY, T2 phase 2). Still read: v3 (T2 phase 1, NSTATE - MAX_SOI_LAY values), v2 (T1b-T1d, NSTATE -
// MAX_SOI_LAY - 1) and v1 (T1, NSTATE - MAX_SOI_LAY - 3).
const int NSTATE = 1 + NUM_PFT * NUM_PFT_PART + NUM_PFT + 1 + 3 * MAX_SOI_LAY + 9;
inline void pack(const State& s, double* v) {
  int k = 0; v[k++] = 4.0;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) v[k++] = s.veg[i][j];
  for (int i = 0; i < NUM_PFT; i++) v[k++] = s.dead[i];
  v[k++] = s.wd;
  for (int l = 0; l < MAX_SOI_LAY; l++) v[k++] = s.soil[l];
  for (int l = 0; l < MAX_SOI_LAY; l++) v[k++] = s.mehg[l];
  v[k++] = s.h_init; v[k++] = s.cum.uptake; v[k++] = s.cum.dep; v[k++] = s.cum.evasion; v[k++] = s.cum.export_; v[k++] = s.cum.fire;
  v[k++] = s.cum.spin; v[k++] = (double)s.spun;
  v[k++] = s.snow;
  for (int l = 0; l < MAX_SOI_LAY; l++) v[k++] = s.diss[l] + ((l == 0) ? s.mpend : 0.0);   // mpend is 0 at month ends
}
inline bool unpack(State& s, const double* v, int nsl) {
  if (!(v[0] == 1.0 || v[0] == 2.0 || v[0] == 3.0 || v[0] == 4.0)) { return false; }
  int k = 1;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) s.veg[i][j] = v[k++];
  for (int i = 0; i < NUM_PFT; i++) s.dead[i] = v[k++];
  s.wd = v[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.soil[l] = v[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.mehg[l] = v[k++];
  s.h_init = v[k++]; s.cum.clear(); s.cum.uptake = v[k++]; s.cum.dep = v[k++]; s.cum.evasion = v[k++]; s.cum.export_ = v[k++]; s.cum.fire = v[k++];
  s.spun = 0; s.spin_tries = 0;
  if (v[0] >= 2.0) { s.cum.spin = v[k++]; s.spun = (int)std::lround(v[k++]); }
  s.snow = (v[0] >= 3.0) ? v[k++] : 0.0;
  for (int l = 0; l < MAX_SOI_LAY; l++) { s.diss[l] = (v[0] >= 4.0) ? v[k++] : 0.0; }
  s.nsl = nsl; s.init = true; s.m.clear(); s.y.clear();
  return true;
}

} // namespace hgtr
#endif /*HGTRACER_H_*/
