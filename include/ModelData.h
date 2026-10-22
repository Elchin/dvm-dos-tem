#ifndef MODELDATA_H_
#define MODELDATA_H_
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

  // DRAFT (daily-structural): infiltration into a frozen top soil layer.
  // model_settings.frozen_infil_eice: < 0 (default, key absent) = off (legacy:
  // no infiltration while the top shallow/moss layer is frozen); >= 0 = ice
  // impedance exponent E in 10^(-E*F_ice) (CLM5/Swenson et al. 2012 use 6).
  double frozen_infil_eice;

  // DRAFT (tvc-summer): model_settings.kersten_form, 0 = legacy ln(Sr)+1
  // (default), 1 = log10(Sr)+1 for unfrozen soil thermal conductivity.
  int kersten_form;
  // DRAFT (daily-structural 4/4): model_settings.soil_evap_beta: 0 (default) =
  // off (legacy, no soil-moisture limit on soil evaporation); 1 = scale soil
  // evaporation by beta = 0.25*(1-cos(pi*min(1,theta_liq/theta_fc)))^2 of the
  // first shallow organic layer (Lee & Pielke 1992; CLM4.0), theta_fc at -3.4 m.
  int soil_evap_beta;

  // Option 1 (v5h): moss/organic thermal conductivity
  // (Cote & Konrad 2005 / CLASS); see SoilLayer.h. org_tc_scheme 0 = off.
  int org_tc_scheme;
  double org_tc_kunf, org_tc_kfrz, org_tc_dry, org_tc_thre_f, org_tc_thre_h,
         org_tc_depth;

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

  // Option 3 (v5h): settable mineral thermal properties; see
  // SoilLayer.h. Keys: mineral_tcsolid (-1 = texture), mineral_vhcsolid
  // (-1 = texture), mineral_tcdry_form (0 = legacy geometric mean).
  double mineral_tcsolid, mineral_vhcsolid;
  int mineral_tcdry_form;

  // Fix (a) (v5i): vapor pressure units, see Climate.h. vp_hpa_fix 0 (default)
  // = legacy (hPa vapor_press used as if Pa); 1 = convert to Pa.
  int vp_hpa_fix;
  // Fix (b) (v5i): n-factor selection day, see Cohort::updateMonthly_Env.
  // nfactor_sameday 0 (default) = legacy previous-day Tair; 1 = same-day Tair.
  int nfactor_sameday;
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
  // Fixes B1-B4 (v5k), all OFF by default (0 = legacy), see Climate.h
  // (svp_ice_fix), Soil_Env.h (richards_infil_fix), Stefan.h (stefan_fronts_fix)
  // and Snow_Env.h (snow_albedo_fix).
  int svp_ice_fix, richards_infil_fix, stefan_fronts_fix, snow_albedo_fix;
  // Fix N1 (v5k), default 0: see Soil_Bgc.h (n_rootless_fix 0/1/2).
  int n_rootless_fix;
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
  // seb_snow_tc: 0 = legacy Goodrich 1982 snow conductivity (default); 1 = Sturm et al. 1997
  //   k(rho) for rho <= 600 kg m-3 (floor 0.023). Independent of seb_scheme.
  int seb_snow_tc;

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
