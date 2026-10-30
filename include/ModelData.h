#ifndef MODELDATA_H_
#define MODELDATA_H_
#include "Bnf.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <json/value.h>

#include "ArgHandler.h"
#include "util_structs.h"
#include "layerconst.h"
#include "errorcode.h"

#include <netcdf.h>
#ifdef WITHMPI
#include <mpi.h>
#include <netcdf_par.h>
#endif

using namespace std;

class ModelData {
public:

  ModelData(Json::Value controldata);

  ModelData();
  ~ModelData();

  void update(ArgHandler const * arghandler);
  std::string describe_module_settings();

  void create_netCDF_output_files(int ysize, int xsize, const std::string & stage, int stage_year_count, bool copy_grid_mapping);

  string loop_order; // time-major or space-major

  int force_cmt; // used to override the veg map (calibration mode only)

  int eq_yrs;
  int pr_yrs;
  int sp_yrs;
  int tr_yrs;
  int sc_yrs;

  //General config settings
  std::string run_name;
  std::string run_description;

  //Config Stage Settings
  bool inter_stage_pause; // Controls pauses between EQ, SP, TR, SC
  //PR modules
  bool pr_env, pr_bgc, pr_nfeed, pr_avln;
  bool pr_baseline, pr_dsb, pr_dsl, pr_dyn_lai;
  //EQ modules
  bool eq_env, eq_bgc, eq_nfeed, eq_avln;
  bool eq_baseline, eq_dsb, eq_dsl, eq_dyn_lai;
  //SP modules
  bool sp_env, sp_bgc, sp_nfeed, sp_avln;
  bool sp_baseline, sp_dsb, sp_dsl, sp_dyn_lai;
  //TR modules
  bool tr_env, tr_bgc, tr_nfeed, tr_avln;
  bool tr_baseline, tr_dsb, tr_dsl, tr_dyn_lai;
  //SC modules
  bool sc_env, sc_bgc, sc_nfeed, sc_avln;
  bool sc_baseline, sc_dsb, sc_dsl, sc_dyn_lai;

  //Config IO settings
  string parameter_dir;
  string hist_climate_file;
  string proj_climate_file;
  string veg_class_file;
  string topo_file;
  string fri_fire_file;
  string hist_exp_fire_file;
  string proj_exp_fire_file;
  string soil_texture_file;
  string drainage_file;
  string co2_file;
  string proj_co2_file;
  string runmask_file;
  string output_dir;
  string restart_from;      // Restart from a previous run
  string output_spec_file;
  bool output_monthly;
  bool nc_eq; // NetCDF output flags for each stage
  bool nc_sp;
  bool nc_tr;
  bool nc_sc;
  int output_interval; //How many years to store for output
  //The following two config values are temporarily stored in
  // ModelData, to be transferred to Climate.
  int baseline_start;//Start year for baseline EQ climate
  int baseline_end;//End year for baseline EQ climate

  // Opt-in daily climate mode (model_settings.climate_timestep).
  // "monthly" (default, key absent) -> original behavior.
  // "daily" -> the projected climate file (used by the SC stage) is a
  //            daily series (365_day calendar, 365 steps per model year)
  //            that is consumed directly, without monthly->daily
  //            interpolation or the stochastic precipitation split.
  std::string climate_timestep;
  bool daily_proj_climate;
  // Refs (daily mode): no published source for the mode itself; calendar = CF 365_day/noleap (CF Conventions 1.8); daily extraterrestrial radiation FAO-56 eq. 21 (Allen et al. 1998); see docs/daily-mode-REFERENCES.md.

  // DRAFT (daily-structural): infiltration into a frozen top soil layer.
  // model_settings.frozen_infil_eice: < 0 (default, key absent) = off (legacy:
  // no infiltration while the top shallow/moss layer is frozen); >= 0 = ice
  // impedance exponent E in 10^(-E*F_ice) (CLM5/Swenson et al. 2012 use 6).
  double frozen_infil_eice;
  // Refs: Swenson et al. 2012, doi:10.1029/2012MS000165; CLM5 Tech Note eq. 2.7.48 (Omega = 6); see docs/daily-mode-REFERENCES.md.

  // DRAFT (tvc-summer): model_settings.kersten_form, 0 = legacy ln(Sr)+1
  // (default), 1 = log10(Sr)+1 for unfrozen soil thermal conductivity.
  int kersten_form;
  // Refs: Johansen 1975 (CRREL transl. 637), doi:10.21236/ADA044002; Lawrence & Slater 2008, doi:10.1007/s00382-007-0278-1; see docs/daily-mode-REFERENCES.md.
  // DRAFT (daily-structural 4/4): model_settings.soil_evap_beta: 0 (default) =
  // off (legacy, no soil-moisture limit on soil evaporation); 1 = scale soil
  // evaporation by beta = 0.25*(1-cos(pi*min(1,theta_liq/theta_fc)))^2 of the
  // first shallow organic layer (Lee & Pielke 1992; CLM4.0), theta_fc at -3.4 m.
  int soil_evap_beta;
  // Refs: Lee & Pielke 1992, doi:10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2; theta_fc from Clapp & Hornberger 1978, doi:10.1029/WR014i004p00601.

  // Option 1 (v5h): moss/organic thermal conductivity
  // (Cote & Konrad 2005 / CLASS); see SoilLayer.h. org_tc_scheme 0 = off.
  int org_tc_scheme;
  double org_tc_kunf, org_tc_kfrz, org_tc_dry, org_tc_thre_f, org_tc_thre_h,
         org_tc_depth;
         // Refs: Cote & Konrad 2005, doi:10.1139/t04-106; Letts et al. 2000, doi:10.1080/07055900.2000.9649643; Chadburn et al. 2015, doi:10.5194/gmd-8-1493-2015.

  // Option 2 (v5h): lateral ("fill and spill") drainage of
  // thawed near-surface organic layers (Quinton & Marsh 1999; Quinton et al.
  // 2000) towards the peat retention capacity theta_re (Letts et al. 2000, as in
  // CLASS/CLASSIC organic soils), cf. the CLM5 perched-water-table lateral
  // drainage (Lawrence et al. 2018/2019). Per day, for every fully thawed
  // fibric/humic layer with top < org_drain_depth:
  //   drain = org_drain_rate * max(0, liq - theta_re*dz*1000)   [mm]
  // booked as subsurface drainage (qdrain). org_drain_rate <= 0 (default) = off.
  // Keys (defaults): org_drain_rate 0 [1/day], org_drain_thre_f 0.275,
  // org_drain_thre_h 0.62, org_drain_depth 0.3 [m].
  double org_drain_rate, org_drain_thre_f, org_drain_thre_h, org_drain_depth;
  // Refs: Quinton & Marsh 1999, doi:10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D; Letts et al. 2000, doi:10.1080/07055900.2000.9649643; cf. CLM5 Tech Note sec. 2.7.4.

  // Option 3 (v5h): settable mineral thermal properties; see
  // SoilLayer.h. Keys: mineral_tcsolid (-1 = texture), mineral_vhcsolid
  // (-1 = texture), mineral_tcdry_form (0 = legacy geometric mean).
  double mineral_tcsolid, mineral_vhcsolid;
  int mineral_tcdry_form;
  // Refs: Johansen 1975, doi:10.21236/ADA044002; Farouki 1981 (CRREL Monograph 81-1), doi:10.21236/ADA111734; CLM5 Tech Note eqs. 2.6.80, 2.6.82.

  // Fix (a) (v5i): vapor pressure units, see Climate.h. vp_hpa_fix 0 (default)
  // = legacy (hPa vapor_press used as if Pa); 1 = convert to Pa.
  int vp_hpa_fix;
  // Refs: input vapour pressure in hPa (CRU TS vap: Harris et al. 2020, doi:10.1038/s41597-020-0453-3); VPD = es - ea in consistent units (FAO-56, Allen et al. 1998); Penman-Monteith (Monteith 1965, PMID 5321565).
  // Fix (b) (v5i): n-factor selection day, see Cohort::updateMonthly_Env.
  // nfactor_sameday 0 (default) = legacy previous-day Tair; 1 = same-day Tair.
  int nfactor_sameday;
  // Refs: n-factor (Klene et al. 2001, doi:10.1080/15230430.2001.12003416); the fix removes an uninitialised selector (code bug, upstream 419a4eb3, commit 9b3a4840).
  // Fix (c) (v5j), OFF by default: model_settings.nonvasc_trans_withdraw.
  // 0 (default) = legacy: non-vascular (moss/lichen) transpiration is reported
  //   in EET but never removed from soil water (no r_e_ij for those PFTs), and
  //   root uptake that the soil-water step cannot apply (frozen top branch,
  //   layers outside the Richards domain) is silently dropped.
  // 1 = non-vascular transpiration is withdrawn from the first shallow organic
  //   layer (the layer the hydrologically excluded moss layer is slaved to),
  //   capped by that layer's liquid above its minimum, and zero when that layer
  //   is below 0.01 C (same criterion as vascular roots); uptake that could not
  //   be applied is removed from the reported transpiration (as for evaporation).
  int nonvasc_trans_withdraw;
  // Refs: no published source for the scheme; rationale: water mass conservation (reported ET must leave storage); context: Heijmans et al. 2004, doi:10.1016/j.agrformet.2003.12.006.
  // Fixes B1-B4 (v5k), all OFF by default (0 = legacy), see Climate.h
  // (svp_ice_fix), Soil_Env.h (richards_infil_fix), Stefan.h (stefan_fronts_fix)
  // and Snow_Env.h (snow_albedo_fix).
  int svp_ice_fix, richards_infil_fix, stefan_fronts_fix, snow_albedo_fix;
  // Refs: at each declaration named above (Climate.h, Soil_Env.h, Stefan.h, Snow_Env.h) and in docs/daily-mode-REFERENCES.md.
  // Fix C1 (v5k), OFF by default (0 = legacy): soil integrator NMIN error
  // tolerance, see Integrator.h (soil_nmin_tol_fix).
  int soil_nmin_tol_fix;
  // Fix L1 (v5m), OFF by default (0 = legacy): the soil integrator safety net
  // debits the root N extraction it skips, see Integrator.h (n_safetynet_conserve).
  int n_safetynet_conserve;
  // v5n, all OFF by default (0 = legacy): G1 n_extract_limiter (+ n_extract_tau),
  // see Soil_Bgc.h; G2 soil_giveup_retry and soil_giveup_diag, see Integrator.h.
  int n_extract_limiter, soil_giveup_retry, soil_giveup_diag;
  int n_netnmin_rescale_guard;  // G3, see Soil_Bgc.h
  double n_extract_tau;
  // Refs (L1): no published source for the fix; N mass conservation of the TEM N cycle (McGuire et al. 1992, doi:10.1029/92GB00219; Yuan et al. 2012, doi:10.1890/11-1957.1); organic N as donor (Nasholm et al. 1998, doi:10.1038/31921); see docs/daily-mode-REFERENCES.md.
  // Refs (C1): no published source for the fix; standard restored: step-size error scaled by a magnitude that cannot pass through zero (Press et al. 2007, Numerical Recipes 3rd ed., sec. 17.2; Hairer, Norsett & Wanner 1993, doi:10.1007/978-3-540-78862-1).
  // Fix N1 (v5k), default 2 since v5o (0 = legacy, selectable): see Soil_Bgc.h (n_rootless_fix 0/1/2).
  int n_rootless_fix;
  // Refs (N1): no published source; code bug in Soil_Bgc::deltastate (upstream 419a4eb3, fix commit 29d72af8); restores N mass conservation.
  // Fixes N2/N3 (v5k), default 1 since v5o (0 = legacy, selectable): see Soil_Bgc.h (n_som_move_fix) and Vegetation_Bgc.h (n_nonvasc_extract).
  int n_som_move_fix, n_nonvasc_extract;
  // v5p: statistical daily split of the monthly C fluxes, output only
  // (0 = off, default; 1 = on; anything else FATAL). With 1 the *_DS daily
  // outputs (and DS_FALLBACK) may be requested in output_spec.csv; requesting
  // them with 0 is FATAL. No effect on the model state. See DailyCSplit.h.
  // Refs: Olsen & Randerson 2004, doi:10.1029/2003JD003968; Fisher et al. 2016,
  //   doi:10.5194/bg-13-4271-2016; principle of exact temporal disaggregation:
  //   Denton 1971, doi:10.1080/01621459.1971.10482227; Chow & Lin 1971, doi:10.2307/1928739.
  //   No published source for the exact method.
  int daily_c_split;
  // v5q: biological N fixation, model_settings.bnf_scheme (0 = off, default;
  // 1 PRO, 2 AET, 3 CONST) and the optional bnf_* parameters; see Bnf.h.
  // Refs: Meyerholt et al. 2016, doi:10.5194/bg-13-1491-2016; Cleveland et al. 1999,
  //   doi:10.1029/1999GB900014; Thornton et al. 2007, doi:10.1029/2006GB002868.
  bnf::Params bnf;
  // v5q-ndrain: N leaching bookkeeping fix, model_settings.n_leach_fix (0 = legacy, default;
  // 1 = the soil AVLN pools lose exactly the booked NLOST); see Soil_Bgc.h. Integers 0/1 only.
  // Refs: no published source (code bug, N mass conservation); see Soil_Bgc.h.
  int n_leach_fix;
  // v5q-negpool: model_settings.soil_negpool_fix (0 = legacy, default; 1 = no negative soil N
  // pools from the legacy AVLN clamp, negatives healed at month start, removed N booked as debt);
  // see Integrator.h. Integers 0/1 only.
  int soil_negpool_fix;
  // v5q-d2m-gate: model_settings.d2m_transfer_gate (0 = legacy, default; 1 = deep->mineral SOM
  // transfer only into thawed receiver thickness, max(running max ALD, previous-year ALD);
  // 2 = transfer off). Integers 0/1/2 only; see Soil_Bgc.h for the scheme and references.
  int d2m_transfer_gate;
  // v5q-layern-fix: model_settings.layer_n_remap_fix (0 = legacy, default, bit-identical;
  // 1 = layer organic N (orgn) and inorganic N (avln) follow the C pools through the
  // dynamic-soil-layer split/combine/remap, so per-layer C:N stays aligned and column N is
  // conserved). Integers 0/1 only; see Ground.h / Soil_Bgc.h.
  // Refs: layer structure and the organic-layer C<->thickness relation of the dynamic organic
  //   soil (DOS) module: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010,
  //   doi:10.1029/2010JG001302. The fix itself is a mass/stoichiometry bookkeeping
  //   (conservation) bug fix: no published source for the remap rule (N follows C).
  int layer_n_remap_fix;
  // v5q-layer-bookkeep: model_settings.layer_bookkeep_fix (0 = legacy, default, bit-identical;
  // 1 = four dynamic-soil-layer C/N bookkeeping (conservation) bug fixes, see Ground.h;
  // requires layer_n_remap_fix = 1). Integers 0/1 only, anything else FATAL.
  // Refs: layer structure / organic C<->thickness relation: Yi et al. 2009,
  //   doi:10.1029/2008JG000841; Yi et al. 2010, doi:10.1029/2010JG001302. The fixes are
  //   mass-conservation bookkeeping bug fixes: no published source.
  int layer_bookkeep_fix;
  // Refs (N2): no published source; code bug in Soil_Bgc::deltastate (upstream 419a4eb3, fix commit 20b2723e); restores N mass conservation (N moves with C).
  // Refs (N3): no published source for the scheme; restores N mass conservation; bryophyte N uptake from the moss/peat layer: Turetsky 2003, doi:10.1639/05.
  // ET-split fixes (v5l), all OFF by default (0 = legacy, bit-identical); see
  // Vegetation_Env.h for the physics of each switch:
  //   et_fpc_once_fix    0/1  S1: apply the foliar projected cover (FPC) once in transpiration
  //   nv_snow_bury_fix   0/1  S2: no snow interception by moss/lichen buried by snow
  //     nv_snow_bury_depth (m, default 0.0): snowpack depth above which they count as buried
  //   rain_int_lai_fix   0/1  S3: rain interception capacity = ircoef*LAI with canopy storage
  //   nv_rain_surface    0/1/2 S4: moss/lichen rain interception is moss-surface water
  //     (1: evaporated as EVAP_MOSS, not canopy evaporation; 2: goes to the ground surface)
  int et_fpc_once_fix, nv_snow_bury_fix, rain_int_lai_fix, nv_rain_surface;
  double nv_snow_bury_depth;
  // Refs (S1-S4): see Vegetation_Env.h and docs/daily-mode-REFERENCES.md.
  static bool v5l_switches_valid(int fpc1, int bury, double bury_depth, int rint, int nvsurf) {
    return (fpc1 == 0 || fpc1 == 1) && (bury == 0 || bury == 1) && (bury_depth >= 0.0)
           && (rint == 0 || rint == 1) && (nvsurf >= 0 && nvsurf <= 2);
  }
  // Option 4 (v5j), OFF by default: daily surface energy balance (SEB) top
  // boundary (ported from the seb-explore exploration branch).
  // model_settings.seb_scheme: 0 = off (n-factor top boundary, unchanged);
  //   2 = SEB over snow only (calibration option): on days with snow on the
  //       ground the snow-surface temperature, melt and sublimation come from a
  //       daily skin energy balance (SW, forcing LW, H, sublimation, G); on
  //       snow-free days the n-factor boundary (nfactor(s)/nfactor(w)) is used
  //       as before (nfactor(w) then only acts on snow-free days with
  //       Tair <= 0, e.g. shoulder seasons; tested with the site values).
  //   1 = EXPERIMENTAL (exploration only, not recommended): SEB on all days,
  //       n-factors unused; snow-free LE = the model's own ET.
  //   Daily climate mode only. LW down: forcing variable 'lwdown' (W m-2) in the
  //   daily projected-climate file, wind: 'wind' (m s-1); see seb_lw_source.
  // Other keys (defaults): seb_z0_snow 0.001 m, seb_z0_soil 0.01 m (momentum
  //   roughness; z0h = z0/10), seb_zwind 10 m, seb_ztemp 2 m (forcing heights),
  //   seb_emis_snow 0.99, seb_emis_soil 0.97, seb_wind -1 (<0 = use forcing
  //   'wind', else constant m/s; 3.0 if no forcing), seb_minwind 0.5 m/s,
  //   seb_stab_b 10 (stable: CH *= 1/(1+b*Ri); 0 = neutral), seb_tran_frac 1
  //   (fraction of transpiration latent heat taken from the ground skin),
  //   seb_pres 100000 Pa, seb_lw_source 0 (0 = forcing 'lwdown' if present,
  //   1 = always parameterized: Prata 1996 clear sky + Crawford & Duchon 1999
  //   cloud from the model's SW cloudiness), seb_diag 0 (1 = write
  //   output_dir/seb_daily.csv).
  int seb_scheme;
  double seb_z0_snow, seb_z0_soil, seb_zwind, seb_ztemp, seb_emis_snow,
         seb_emis_soil, seb_wind, seb_minwind, seb_stab_b, seb_tran_frac,
         seb_pres;
  int seb_lw_source, seb_diag;
  // Refs (SEB): bulk transfer, Brutsaert 1982, doi:10.1007/978-94-017-1497-6; simplified stability functions after Louis 1979, doi:10.1007/BF00117978
  //   (stable) and Dyer 1974, doi:10.1007/BF00240838 (unstable); LW: Prata 1996, doi:10.1002/qj.49712253306; Crawford & Duchon 1999,
  //   doi:10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2; snow emissivity Warren 1982, doi:10.1029/RG020i001p00067; esat WMO-No. 8 (2008);
  //   safeguarded Newton, Press et al. 2007 (Numerical Recipes 3rd ed., sec. 9.4); see docs/daily-mode-REFERENCES.md.
  // seb_snow_tc: 0 = legacy Goodrich 1982 snow conductivity (default); 1 = Sturm et al. 1997
  //   k(rho) for rho <= 600 kg m-3 (floor 0.023). Independent of seb_scheme.
  int seb_snow_tc;
  // Refs: Sturm et al. 1997, doi:10.3189/S0022143000002781; legacy Goodrich 1982, doi:10.1139/t82-047.

  // Maps holding data about variables to be output at specific timesteps
  // C++11 would allow the use of unordered_maps, which have a faster
  // by-key access time.
  std::map<std::string, OutputSpec> daily_netcdf_outputs;
  std::map<std::string, OutputSpec> monthly_netcdf_outputs;
  std::map<std::string, OutputSpec> yearly_netcdf_outputs;

  std::string pid_tag;
  std::string caldata_tree_loc;
  int last_n_json_files;
  bool archive_all_json;
  bool tar_caljson;

  int changeclimate; // 0: default (up to run stage); 1: dynamical; -1: static
  int changeco2; // 0: default (up to run stage); 1: dynamical; -1: static
  bool dynamic_LAI; // True: calculate LAI as a function of vegc, False: use static_lai from CohortLookup 
  bool useseverity; // using fire severity inputs

  int cell_timelimit; //Time limit in seconds for cell computation time
  time_t cell_stime; //Start time per cell. Move to Runner?

  bool outSiteDay;

  bool get_envmodule();
  void set_envmodule(const std::string &s);
  void set_envmodule(const bool v);

  bool get_bgcmodule();
  void set_bgcmodule(const std::string &s);
  void set_bgcmodule(const bool v);

  bool get_dynamic_lai_module();
  void set_dynamic_lai_module(const std::string &s);
  void set_dynamic_lai_module(const bool v);

  bool get_dslmodule();
  void set_dslmodule(const std::string &s);
  void set_dslmodule(const bool v);

  bool get_dsbmodule();
  void set_dsbmodule(const std::string &s);
  void set_dsbmodule(const bool v);

  bool get_nfeed();
  void set_nfeed(const std::string &s);
  void set_nfeed(const bool v);

  bool get_avlnflg();
  void set_avlnflg(const std::string &s);
  void set_avlnflg(const bool v);

  bool get_baseline();
  void set_baseline(const std::string &s);
  void set_baseline(const bool v);


private:

  // the following 3 switches will control N modules in BGC
  bool avlnflg; // inorganic N in/out module - partial open N cycle

  bool nfeed; // when true, allows N uptake to be limited by soil conditions
  // which then controls plant growth. Basically its a switch
  // for soil-plant N process module.

  bool baseline;   // When true allowing ninput and nlost to be used for
  // adjusting c/n of soil - partial open N cycle


  bool envmodule;             // (Bio?)physical module on/off
  bool bgcmodule;             // BGC module on/off
  bool dynamic_lai_module;    // dynamic lai module on/off

  bool dslmodule;  // dynamic soil layer module on/off
  bool dsbmodule;  // disturbance module on/off

  // Note: the dsl module ON lets the thickness, number and type of layers
  //       change. With dsl module OFF, the C and N content will change, but
  //       the thickness and number of layers should not change.

};

#endif /*MODELDATA_H_*/
