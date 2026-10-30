#include <exception>
#include <iostream>
#include <cstdlib>
#include <json/value.h>
#include <boost/filesystem.hpp>
#include "../include/ModelData.h"
#include "../include/DailyCSplit.h"
#include "../include/timeconst.h"

#include "../include/TEMLogger.h"
#include "../include/TEMUtilityFunctions.h"

extern src::severity_logger< severity_level > glg;


/** Returns a string with first colum r justified and
 * of with 'w'. Can be used to build tables likle this:
 *
 *       somestr: 0
 *       somestr: 0
 *       somestr: 0
 *
 *  with a newline, for use in a table.
 */
std::string table_row(int w, std::string d, bool v) {
  std::stringstream s;
  s << std::setw(w) << std::setfill(' ') << d << ": " << v << "\n";
  return s.str();
}

ModelData::~ModelData() {}

ModelData::ModelData(Json::Value controldata):force_cmt(-1) {

  BOOST_LOG_SEV(glg, debug) << "Creating a ModelData. New style constructor with injected controldata...";

  //General config settings
  run_name = controldata["general"]["output_global_attributes"]["run_name"].asString();
  run_description = controldata["general"]["output_global_attributes"]["description"].asString();

  //Config Stage Settings  
//  std::string stgstr(controldata["stage_settings"]["run_stage"].asString());

  inter_stage_pause = controldata["stage_settings"]["inter_stage_pause"].asBool();
  tr_yrs = controldata["stage_settings"]["tr_yrs"].asInt();
  sc_yrs = controldata["stage_settings"]["sc_yrs"].asInt();

  //PR stage module settings
  pr_env = controldata["stage_settings"]["pr"]["env"].asBool();
  pr_bgc = controldata["stage_settings"]["pr"]["bgc"].asBool();
  pr_nfeed = controldata["stage_settings"]["pr"]["nfeed"].asBool();
  pr_avln = controldata["stage_settings"]["pr"]["avlnflg"].asBool();
  pr_baseline = controldata["stage_settings"]["pr"]["baseline"].asBool();
  pr_dsb = controldata["stage_settings"]["pr"]["dsb"].asBool();
  pr_dsl = controldata["stage_settings"]["pr"]["dsl"].asBool();
  pr_dyn_lai = controldata["stage_settings"]["pr"]["dyn_lai"].asBool();
  //EQ stage module settings
  eq_env = controldata["stage_settings"]["eq"]["env"].asBool();
  eq_bgc = controldata["stage_settings"]["eq"]["bgc"].asBool();
  eq_nfeed = controldata["stage_settings"]["eq"]["nfeed"].asBool();
  eq_avln = controldata["stage_settings"]["eq"]["avlnflg"].asBool();
  eq_baseline = controldata["stage_settings"]["eq"]["baseline"].asBool();
  eq_dsb = controldata["stage_settings"]["eq"]["dsb"].asBool();
  eq_dsl = controldata["stage_settings"]["eq"]["dsl"].asBool();
  eq_dyn_lai = controldata["stage_settings"]["eq"]["dyn_lai"].asBool();
  //SP stage module settings
  sp_env = controldata["stage_settings"]["sp"]["env"].asBool();
  sp_bgc = controldata["stage_settings"]["sp"]["bgc"].asBool();
  sp_nfeed = controldata["stage_settings"]["sp"]["nfeed"].asBool();
  sp_avln = controldata["stage_settings"]["sp"]["avlnflg"].asBool();
  sp_baseline = controldata["stage_settings"]["sp"]["baseline"].asBool();
  sp_dsb = controldata["stage_settings"]["sp"]["dsb"].asBool();
  sp_dsl = controldata["stage_settings"]["sp"]["dsl"].asBool();
  sp_dyn_lai = controldata["stage_settings"]["sp"]["dyn_lai"].asBool();
  //TR stage module settings
  tr_env = controldata["stage_settings"]["tr"]["env"].asBool();
  tr_bgc = controldata["stage_settings"]["tr"]["bgc"].asBool();
  tr_nfeed = controldata["stage_settings"]["tr"]["nfeed"].asBool();
  tr_avln = controldata["stage_settings"]["tr"]["avlnflg"].asBool();
  tr_baseline = controldata["stage_settings"]["tr"]["baseline"].asBool();
  tr_dsb = controldata["stage_settings"]["tr"]["dsb"].asBool();
  tr_dsl = controldata["stage_settings"]["tr"]["dsl"].asBool();
  tr_dyn_lai = controldata["stage_settings"]["tr"]["dyn_lai"].asBool();
  //SC stage module settings
  sc_env = controldata["stage_settings"]["sc"]["env"].asBool();
  sc_bgc = controldata["stage_settings"]["sc"]["bgc"].asBool();
  sc_nfeed = controldata["stage_settings"]["sc"]["nfeed"].asBool();
  sc_avln = controldata["stage_settings"]["sc"]["avlnflg"].asBool();
  sc_baseline = controldata["stage_settings"]["sc"]["baseline"].asBool();
  sc_dsb = controldata["stage_settings"]["sc"]["dsb"].asBool();
  sc_dsl = controldata["stage_settings"]["sc"]["dsl"].asBool();
  sc_dyn_lai = controldata["stage_settings"]["sc"]["dyn_lai"].asBool();

  //Config IO Settings
  parameter_dir     = controldata["IO"]["parameter_dir"].asString();
  hist_climate_file = controldata["IO"]["hist_climate_file"].asString();
  proj_climate_file = controldata["IO"]["proj_climate_file"].asString();
  veg_class_file    = controldata["IO"]["veg_class_file"].asString();
  fri_fire_file     = controldata["IO"]["fri_fire_file"].asString();
  hist_exp_fire_file= controldata["IO"]["hist_exp_fire_file"].asString();
  proj_exp_fire_file= controldata["IO"]["proj_exp_fire_file"].asString();
  topo_file         = controldata["IO"]["topo_file"].asString();
  drainage_file     = controldata["IO"]["drainage_file"].asString();
  soil_texture_file = controldata["IO"]["soil_texture_file"].asString();
  co2_file          = controldata["IO"]["co2_file"].asString();
  proj_co2_file     = controldata["IO"]["proj_co2_file"].asString();
  runmask_file      = controldata["IO"]["runmask_file"].asString();
  output_dir        = controldata["IO"]["output_dir"].asString();
  restart_from      = controldata["IO"]["restart_from"].asString();
  output_spec_file  = controldata["IO"]["output_spec_file"].asString();
  output_monthly    = controldata["IO"]["output_monthly"].asInt();
  nc_eq             = controldata["IO"]["output_nc_eq"].asBool();
  nc_sp             = controldata["IO"]["output_nc_sp"].asBool();
  nc_tr             = controldata["IO"]["output_nc_tr"].asBool();
  nc_sc             = controldata["IO"]["output_nc_sc"].asBool();
  output_interval   = controldata["IO"]["output_interval"].asInt();

  //Config Calibration IO Settings
  pid_tag           = controldata["calibration-IO"]["pid_tag"].asString();
  caldata_tree_loc  = controldata["calibration-IO"]["caldata_tree_loc"].asString();

  cell_timelimit = controldata["model_settings"]["cell_timelimit"].asInt();
  dynamic_LAI       = controldata["model_settings"]["dynamic_lai"].asInt(); // checked in Cohort::updateMonthly_DIMVeg
  baseline_start = controldata["model_settings"]["baseline_start"].asInt();
  baseline_end   = controldata["model_settings"]["baseline_end"].asInt();

  // Opt-in daily climate mode. Absent/empty key keeps original monthly
  // behavior. Any value other than "monthly" or "daily" is a hard error.
  climate_timestep = controldata["model_settings"]["climate_timestep"].asString();
  frozen_infil_eice = controldata["model_settings"].get("frozen_infil_eice", -1.0).asDouble();
  kersten_form = controldata["model_settings"].get("kersten_form", 0).asInt();
  soil_evap_beta = controldata["model_settings"].get("soil_evap_beta", 0).asInt();
  // Option 1 (v5h), defaults = off / CLASS organic values
  org_tc_scheme = controldata["model_settings"].get("org_tc_scheme", 0).asInt();
  org_tc_kunf   = controldata["model_settings"].get("org_tc_kunf", 0.6).asDouble();
  org_tc_kfrz   = controldata["model_settings"].get("org_tc_kfrz", 0.25).asDouble();
  org_tc_dry    = controldata["model_settings"].get("org_tc_dry", -1.0).asDouble();
  org_tc_thre_f = controldata["model_settings"].get("org_tc_thre_f", 0.275).asDouble();
  org_tc_thre_h = controldata["model_settings"].get("org_tc_thre_h", 0.62).asDouble();
  org_tc_depth  = controldata["model_settings"].get("org_tc_depth", 10.0).asDouble();
  // Option 2 (v5h), default rate 0 = off
  org_drain_rate   = controldata["model_settings"].get("org_drain_rate", 0.0).asDouble();
  org_drain_thre_f = controldata["model_settings"].get("org_drain_thre_f", 0.275).asDouble();
  org_drain_thre_h = controldata["model_settings"].get("org_drain_thre_h", 0.62).asDouble();
  org_drain_depth  = controldata["model_settings"].get("org_drain_depth", 0.3).asDouble();
  // Option 3 (v5h), defaults = legacy
  mineral_tcsolid    = controldata["model_settings"].get("mineral_tcsolid", -1.0).asDouble();
  mineral_vhcsolid   = controldata["model_settings"].get("mineral_vhcsolid", -1.0).asDouble();
  mineral_tcdry_form = controldata["model_settings"].get("mineral_tcdry_form", 0).asInt();
  // Fix (a) (v5i), default 0 = off (legacy)
  vp_hpa_fix = controldata["model_settings"].get("vp_hpa_fix", 0).asInt();
  // Fix (b) (v5i), default 0 = off (legacy)
  nfactor_sameday = controldata["model_settings"].get("nfactor_sameday", 0).asInt();
  // Fix (c) (v5j), default 0 = off (legacy)
  nonvasc_trans_withdraw = controldata["model_settings"].get("nonvasc_trans_withdraw", 0).asInt();
  // Fixes B1-B4 (v5k), default 0 = off (legacy)
  svp_ice_fix        = controldata["model_settings"].get("svp_ice_fix", 0).asInt();
  richards_infil_fix = controldata["model_settings"].get("richards_infil_fix", 0).asInt();
  stefan_fronts_fix  = controldata["model_settings"].get("stefan_fronts_fix", 0).asInt();
  snow_albedo_fix    = controldata["model_settings"].get("snow_albedo_fix", 0).asInt();
  // Fix C1 (v5k), default 0 = off (legacy)
  soil_nmin_tol_fix  = controldata["model_settings"].get("soil_nmin_tol_fix", 0).asInt();
  if (soil_nmin_tol_fix < 0 || soil_nmin_tol_fix > 1) {
    std::cerr << "FATAL: model_settings soil_nmin_tol_fix must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  // Fix L1 (v5m), default 0 = off (legacy). Refs (L1): see Integrator.h.
  n_safetynet_conserve = controldata["model_settings"].get("n_safetynet_conserve", 0).asInt();
  if (n_safetynet_conserve < 0 || n_safetynet_conserve > 1) {
    std::cerr << "FATAL: model_settings n_safetynet_conserve must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  // v5n (G1, G2, diagnostics), default 0 = off (legacy). Refs: see Soil_Bgc.h / Integrator.h.
  n_extract_limiter = controldata["model_settings"].get("n_extract_limiter", 0).asInt();
  n_extract_tau     = controldata["model_settings"].get("n_extract_tau", 0.2).asDouble();
  soil_giveup_retry = controldata["model_settings"].get("soil_giveup_retry", 0).asInt();
  soil_giveup_diag  = controldata["model_settings"].get("soil_giveup_diag", 0).asInt();
  n_netnmin_rescale_guard = controldata["model_settings"].get("n_netnmin_rescale_guard", 0).asInt(); // Fix G3
  if (n_netnmin_rescale_guard < 0 || n_netnmin_rescale_guard > 1) {
    std::cerr << "FATAL: model_settings n_netnmin_rescale_guard must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (n_extract_limiter < 0 || n_extract_limiter > 1 || soil_giveup_retry < 0 || soil_giveup_retry > 1
      || soil_giveup_diag < 0 || soil_giveup_diag > 1) {
    std::cerr << "FATAL: model_settings n_extract_limiter, soil_giveup_retry and soil_giveup_diag must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (!(n_extract_tau > 0.0 && n_extract_tau <= 1.0)) {
    std::cerr << "FATAL: model_settings n_extract_tau must be in (0, 1] (month fraction)" << std::endl;
    exit(EXIT_FAILURE);
  }
  // v5o: N1-N3 fixes ON by default (2/1/1); legacy = explicit 0 (bit-identical to v5n).
  // Refs (v5o defaults N1-N3 on): no published source for the fixes (code bugs, see Soil_Bgc.h / Vegetation_Bgc.h);
  // a closed TEM N cycle (no unbooked N sources/sinks): McGuire et al. 1992, doi:10.1029/92GB00219; Raich et al. 1991,
  // doi:10.2307/1941899. Legacy values (0/0/0) stay selectable and reproduce v5n bit-for-bit. See docs/daily-mode-REFERENCES.md.
  n_rootless_fix     = controldata["model_settings"].get("n_rootless_fix", 2).asInt(); // Fix N1 (v5k), default 2 since v5o
  n_som_move_fix     = controldata["model_settings"].get("n_som_move_fix", 1).asInt();     // Fix N2 (v5k), default 1 since v5o
  n_nonvasc_extract  = controldata["model_settings"].get("n_nonvasc_extract", 1).asInt();  // Fix N3 (v5k), default 1 since v5o
  if (n_som_move_fix < 0 || n_som_move_fix > 1 || n_nonvasc_extract < 0 || n_nonvasc_extract > 1) {
    std::cerr << "FATAL: model_settings n_som_move_fix and n_nonvasc_extract must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (n_rootless_fix < 0 || n_rootless_fix > 2) {
    std::cerr << "FATAL: model_settings.n_rootless_fix must be 0, 1 or 2" << std::endl;
    exit(EXIT_FAILURE);
  }
  // v5p: daily_c_split (output only), default 0; only the numbers 0 and 1 are valid.
  // Refs: see DailyCSplit.h / ModelData.h.
  {
    Json::Value dcs = controldata["model_settings"].get("daily_c_split", 0);
    if (!dcsplit::valid_switch(dcs)) {
      std::cerr << "FATAL: model_settings.daily_c_split must be 0 or 1" << std::endl;
      exit(EXIT_FAILURE);
    }
    daily_c_split = dcs.asInt();
  }
  // v5q: biological N fixation (bnf_scheme 0 = off, default). Refs: see Bnf.h.
  {
    std::string err = bnf::parse(controldata["model_settings"], n_rootless_fix, bnf);
    if (!err.empty()) {
      std::cerr << err << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // v5q-ndrain: n_leach_fix (0 = legacy, default; 1 = pools lose exactly NLOST). Refs: see Soil_Bgc.h.
  {
    Json::Value nl = controldata["model_settings"].get("n_leach_fix", 0);
    if (!bnf::valid_int(nl, 0, 1)) {
      std::cerr << "FATAL: model_settings.n_leach_fix must be 0 (legacy) or 1 (leaching removed = NLOST)" << std::endl;
      exit(EXIT_FAILURE);
    }
    n_leach_fix = nl.asInt();
  }
  // v5q-negpool: soil_negpool_fix (0 = legacy, default). Refs: see Integrator.h.
  {
    Json::Value np = controldata["model_settings"].get("soil_negpool_fix", 0);
    if (!bnf::valid_int(np, 0, 1)) {
      std::cerr << "FATAL: model_settings.soil_negpool_fix must be 0 (legacy) or 1 (no negative soil N pools)" << std::endl;
      exit(EXIT_FAILURE);
    }
    soil_negpool_fix = np.asInt();
  }
  // v5q-d2m-gate: d2m_transfer_gate (0 = legacy, default). Refs: see Soil_Bgc.h.
  {
    Json::Value dg = controldata["model_settings"].get("d2m_transfer_gate", 0);
    if (!bnf::valid_int(dg, 0, 2)) {
      std::cerr << "FATAL: model_settings.d2m_transfer_gate must be 0 (legacy), 1 (thaw-gated) or 2 (off)" << std::endl;
      exit(EXIT_FAILURE);
    }
    d2m_transfer_gate = dg.asInt();
  }
  // v5q-layern-fix: layer_n_remap_fix (0 = legacy, default). Refs: see Ground.h.
  {
    Json::Value lf = controldata["model_settings"].get("layer_n_remap_fix", 0);
    if (!bnf::valid_int(lf, 0, 1)) {
      std::cerr << "FATAL: model_settings.layer_n_remap_fix must be 0 (legacy) or 1 (layer N follows C in layer split/combine)" << std::endl;
      exit(EXIT_FAILURE);
    }
    layer_n_remap_fix = lf.asInt();
  }
  // v5q-layer-bookkeep: layer_bookkeep_fix (0 = legacy, default). Refs: see Ground.h.
  {
    Json::Value bk = controldata["model_settings"].get("layer_bookkeep_fix", 0);
    if (!bnf::valid_int(bk, 0, 1)) {
      std::cerr << "FATAL: model_settings.layer_bookkeep_fix must be 0 (legacy) or 1 (dynamic soil-layer C/N bookkeeping fixes)" << std::endl;
      exit(EXIT_FAILURE);
    }
    layer_bookkeep_fix = bk.asInt();
    if (layer_bookkeep_fix == 1 && layer_n_remap_fix != 1) {
      std::cerr << "FATAL: model_settings.layer_bookkeep_fix = 1 requires layer_n_remap_fix = 1 (layer N must follow C)" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // v5q-thawfrac-decomp: thawfrac_decomp (0 = legacy, default). Refs: see Soil_Bgc.h.
  {
    Json::Value tf = controldata["model_settings"].get("thawfrac_decomp", 0);
    if (!bnf::valid_int(tf, 0, 1)) {
      std::cerr << "FATAL: model_settings.thawfrac_decomp must be 0 (legacy) or 1 (decomposition modifier blended by the layer thawed fraction)" << std::endl;
      exit(EXIT_FAILURE);
    }
    thawfrac_decomp = tf.asInt();
  }
  // v5q-thawfrac-lmoist: thawfrac_moist_mode (0 = sws/f, default; 1 = layer-mean; 1 needs thawfrac_decomp = 1). Refs: see Soil_Bgc.h.
  {
    Json::Value tm = controldata["model_settings"].get("thawfrac_moist_mode", 0);
    if (!bnf::valid_int(tm, 0, 1)) {
      std::cerr << "FATAL: model_settings.thawfrac_moist_mode must be 0 (thawed-part moisture min(1, sws/f)) or 1 (layer-mean moisture)" << std::endl;
      exit(EXIT_FAILURE);
    }
    thawfrac_moist_mode = tm.asInt();
    if (thawfrac_moist_mode == 1 && thawfrac_decomp != 1) {
      std::cerr << "FATAL: model_settings.thawfrac_moist_mode = 1 requires thawfrac_decomp = 1 (it only acts in the thawed-fraction modifier)" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // v5q-leach-curn: n_leach_curn_fix (0 = legacy, default; 1 needs n_leach_fix = 1). Refs: see Soil_Bgc.h.
  {
    Json::Value lc = controldata["model_settings"].get("n_leach_curn_fix", 0);
    if (!bnf::valid_int(lc, 0, 1)) {
      std::cerr << "FATAL: model_settings.n_leach_curn_fix must be 0 (legacy) or 1 (N drainage leaching from the current drainage-zone AVLN)" << std::endl;
      exit(EXIT_FAILURE);
    }
    n_leach_curn_fix = lc.asInt();
    if (n_leach_curn_fix == 1 && n_leach_fix != 1) {
      std::cerr << "FATAL: model_settings.n_leach_curn_fix = 1 requires n_leach_fix = 1 (pools lose exactly NLOST)" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // ET-split fixes S1-S4 (v5l), default 0 = off (legacy); see ModelData.h / Vegetation_Env.h
  et_fpc_once_fix    = controldata["model_settings"].get("et_fpc_once_fix", 0).asInt();
  nv_snow_bury_fix   = controldata["model_settings"].get("nv_snow_bury_fix", 0).asInt();
  nv_snow_bury_depth = controldata["model_settings"].get("nv_snow_bury_depth", 0.0).asDouble();
  rain_int_lai_fix   = controldata["model_settings"].get("rain_int_lai_fix", 0).asInt();
  nv_rain_surface    = controldata["model_settings"].get("nv_rain_surface", 0).asInt();
  if (!ModelData::v5l_switches_valid(et_fpc_once_fix, nv_snow_bury_fix, nv_snow_bury_depth,
                                     rain_int_lai_fix, nv_rain_surface)) {
    std::cerr << "FATAL: model_settings et_fpc_once_fix, nv_snow_bury_fix and rain_int_lai_fix"
              << " must be 0 or 1, nv_rain_surface 0, 1 or 2, nv_snow_bury_depth >= 0 (m)" << std::endl;
    exit(EXIT_FAILURE);
  }
  // Option 4 (v5j): SEB top boundary, default seb_scheme 0 = off (see ModelData.h)
  seb_scheme    = controldata["model_settings"].get("seb_scheme", 0).asInt();
  seb_z0_snow   = controldata["model_settings"].get("seb_z0_snow", 0.001).asDouble();
  seb_z0_soil   = controldata["model_settings"].get("seb_z0_soil", 0.01).asDouble();
  seb_zwind     = controldata["model_settings"].get("seb_zwind", 10.0).asDouble();
  seb_ztemp     = controldata["model_settings"].get("seb_ztemp", 2.0).asDouble();
  seb_emis_snow = controldata["model_settings"].get("seb_emis_snow", 0.99).asDouble();
  seb_emis_soil = controldata["model_settings"].get("seb_emis_soil", 0.97).asDouble();
  seb_wind      = controldata["model_settings"].get("seb_wind", -1.0).asDouble();
  seb_minwind   = controldata["model_settings"].get("seb_minwind", 0.5).asDouble();
  seb_stab_b    = controldata["model_settings"].get("seb_stab_b", 10.0).asDouble();
  seb_tran_frac = controldata["model_settings"].get("seb_tran_frac", 1.0).asDouble();
  seb_pres      = controldata["model_settings"].get("seb_pres", 100000.0).asDouble();
  seb_lw_source = controldata["model_settings"].get("seb_lw_source", 0).asInt();
  seb_diag      = controldata["model_settings"].get("seb_diag", 0).asInt();
  seb_snow_tc   = controldata["model_settings"].get("seb_snow_tc", 0).asInt();
  if (climate_timestep.empty()) {
    climate_timestep = "monthly";
  }
  if (climate_timestep == "monthly") {
    daily_proj_climate = false;
  } else if (climate_timestep == "daily") {
    daily_proj_climate = true;
  } else {
    BOOST_LOG_SEV(glg, fatal) << "Invalid model_settings.climate_timestep: '"
                              << climate_timestep << "' (must be 'monthly' or 'daily')";
    std::cerr << "FATAL: Invalid model_settings.climate_timestep: '"
              << climate_timestep << "' (must be 'monthly' or 'daily')" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (seb_scheme > 0 && !daily_proj_climate) {
    std::cerr << "FATAL: model_settings.seb_scheme > 0 requires climate_timestep 'daily'" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (seb_scheme < 0 || seb_scheme > 2 || seb_snow_tc < 0 || seb_snow_tc > 1 || seb_lw_source < 0 || seb_lw_source > 1) {
    std::cerr << "FATAL: model_settings seb_scheme must be 0/1/2, seb_snow_tc and seb_lw_source 0/1" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (svp_ice_fix < 0 || svp_ice_fix > 1 || richards_infil_fix < 0 || richards_infil_fix > 1 ||
      stefan_fronts_fix < 0 || stefan_fronts_fix > 1 || snow_albedo_fix < 0 || snow_albedo_fix > 1) {
    std::cerr << "FATAL: model_settings svp_ice_fix, richards_infil_fix, stefan_fronts_fix and snow_albedo_fix must be 0 or 1" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (seb_scheme == 1) {
    std::cerr << "WARNING: model_settings.seb_scheme = 1 is EXPERIMENTAL (SEB on snow-free days too)" << std::endl;
  }

  // These must be directories. If user forgets to add the trailing slash, 
  // do it here. Might be better to do something with boost::filesystem for 
  // portability.
  output_dir.append("/");
  parameter_dir.append("/");

  // Unused (11/23/2015)
  //changeclimate = controldata["model_settings"]["dynamic_climate"].asInt();
  //changeco2     = controldata["model_settings"]["varied_co2"].asInt();
  //useseverity   = controldata["model_settings"]["fire_severity_as_input"].asInt();

}

/** Update all the appropriate fields in ModelData from an ArgHandler object.

    Pass const * so that access to ArgHandler is read-only.
*/
void ModelData::update(ArgHandler const * arghandler) {

  BOOST_LOG_SEV(glg, debug) << "Updating ModelData from an ArgHandler...";

  this->pr_yrs = arghandler->get_pr_yrs();
  this->eq_yrs = arghandler->get_eq_yrs();
  this->sp_yrs = arghandler->get_sp_yrs();
  this->tr_yrs = arghandler->get_tr_yrs();
  this->sc_yrs = arghandler->get_sc_yrs();
  this->pid_tag = arghandler->get_pid_tag();
  this->last_n_json_files = arghandler->get_last_n_json_files();
  this->archive_all_json = arghandler->get_archive_all_json();
  this->tar_caljson = arghandler->get_tar_caljson();

  // it it was set on the command line, then use that value, otherwise,
  // use the value
  if (arghandler->get_inter_stage_pause()) {
    this->inter_stage_pause = arghandler->get_inter_stage_pause();
  }

  // User wants to override the veg map
  if (arghandler->get_force_cmt() >= 0) {
    this->force_cmt = arghandler->get_force_cmt();
  }
}


ModelData::ModelData():force_cmt(-1) {
  et_fpc_once_fix = 0; nv_snow_bury_fix = 0; rain_int_lai_fix = 0; nv_rain_surface = 0; // v5l (off)
  n_safetynet_conserve = 0; // v5m L1 (off)
  n_extract_limiter = 0; soil_giveup_retry = 0; soil_giveup_diag = 0; n_extract_tau = 0.2; n_netnmin_rescale_guard = 0; // v5n (off)
  n_rootless_fix = 2; n_som_move_fix = 1; n_nonvasc_extract = 1; // v5o: N1-N3 fixes on by default (legacy = 0)
  daily_c_split = 0; // v5p (off)
  n_leach_fix = 0; // v5q-ndrain (off)
  soil_negpool_fix = 0; // v5q-negpool (off)
  d2m_transfer_gate = 0; // v5q-d2m-gate (legacy)
  layer_n_remap_fix = 0; // v5q-layern-fix (legacy)
  layer_bookkeep_fix = 0; // v5q-layer-bookkeep (legacy)
  thawfrac_decomp = 0; // v5q-thawfrac-decomp (legacy)
  thawfrac_moist_mode = 0; // v5q-thawfrac-lmoist (sws/f)
  n_leach_curn_fix = 0; // v5q-leach-curn (legacy)
  nv_snow_bury_depth = 0.0;
  set_envmodule(false);
  set_bgcmodule(false);
  set_dynamic_lai_module(false);
  set_dslmodule(false);
  set_dsbmodule(false);
}


bool ModelData::get_envmodule() {
  return this->envmodule;
}
void ModelData::set_envmodule(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting envmodule to " << s;
  this->envmodule = temutil::onoffstr2bool(s);
}
void ModelData::set_envmodule(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting envmodule to " << v;
  this->envmodule = v;
}

bool ModelData::get_bgcmodule() {
  return this->bgcmodule;
}
void ModelData::set_bgcmodule(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting bgcmodule to " << s;
  this->bgcmodule = temutil::onoffstr2bool(s);
}
void ModelData::set_bgcmodule(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting bgcmodule to " << v;
  this->bgcmodule = v;
}

bool ModelData::get_dynamic_lai_module() {
  return this->dynamic_lai_module;
}
void ModelData::set_dynamic_lai_module(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting dynamic_lai_module to " << s;
  this->dynamic_lai_module = temutil::onoffstr2bool(s);
}
void ModelData::set_dynamic_lai_module(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting dynamic_lai_module to " << v;
  this->dynamic_lai_module = v;
}

bool ModelData::get_dslmodule() {
  return this->dslmodule;
}
void ModelData::set_dslmodule(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting dslmodule to " << s;
  this->dslmodule = temutil::onoffstr2bool(s);
}
void ModelData::set_dslmodule(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting dslmodule to " << v;
  this->dslmodule = v;
}

bool ModelData::get_dsbmodule() {
  return this->dsbmodule;
}
void ModelData::set_dsbmodule(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting dsbmodule to " << s;
  this->dsbmodule = temutil::onoffstr2bool(s);
}
void ModelData::set_dsbmodule(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting dsbmodule to " << v;
  this->dsbmodule = v;
}

bool ModelData::get_nfeed() {
  return this->nfeed;
}
void ModelData::set_nfeed(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.nfeed to " << s;
  this->nfeed = temutil::onoffstr2bool(s);
}
void ModelData::set_nfeed(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.nfeed to " << v;
  this->nfeed = v;
}

bool ModelData::get_avlnflg() {
  return this->avlnflg;
}
void ModelData::set_avlnflg(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.avlnflg to " << s;
  this->avlnflg = temutil::onoffstr2bool(s);
}
void ModelData::set_avlnflg(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.avlnflg to " << v;
  this->avlnflg = v;
}

bool ModelData::get_baseline() {
  return this->baseline;
}
void ModelData::set_baseline(const std::string &s) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.baseline to " << s;
  this->baseline = temutil::onoffstr2bool(s);
}
void ModelData::set_baseline(const bool v) {
  BOOST_LOG_SEV(glg, info) << "Setting ModelData.baseline to " << v;
  this->baseline = v;
}

std::string ModelData::describe_module_settings() {
  std::stringstream s;
  s << table_row(22, "envmodule", this->get_envmodule());
  s << table_row(22, "bgcmodule", this->get_bgcmodule());
  s << table_row(22, "dynamic_lai_module", this->get_dynamic_lai_module());
  s << table_row(22, "dslmodule", this->get_dslmodule());
  s << table_row(22, "dsbmodule", this->get_dsbmodule());
  s << table_row(22, "baseline", this->get_baseline());
  s << table_row(22, "nfeed", this->get_nfeed());
  s << table_row(22, "avlnflg", this->get_avlnflg());
  return s.str();
}


/** Construct empty netCDF output files.
 *
 *  This function reads in output selections from a csv file specified
 *  in the config file. It creates an OutputSpec object and an empty
 *  NetCDF file for each line in the output_spec.csv file.
 *
 *  In addition, this function can copy the grid mapping from the input
 *  file(s)
*/
void ModelData::create_netCDF_output_files(int ysize, int xsize,
    const std::string& stage, int stage_year_count, bool copy_grip_mapping) {

  boost::filesystem::path output_base = output_dir;

  // Load output specification file
  BOOST_LOG_SEV(glg, debug) << "Loading output specification file "<<output_spec_file;
  std::ifstream output_csv(output_spec_file.c_str());

  std::string s;
  std::getline(output_csv, s); // Discard first line - header strings

  std::string token;    // Substrings between commas
  std::string name;     // CSV file variable name
  std::string timestep; // Yearly, monthly, or daily
  std::string invalid_option = "invalid"; // This marks an invalid selection

  // Handle for the NetCDF file
  int ncid;

  // NetCDF file dimensions, applicable to all output files
  int timeD;    // unlimited dimension
  int xD;
  int yD;

  // NetCDF file dimensions; different dims are applicable for different vars.
  int pftD;
  int pftpartD;
  int layerD;

  // NetCDF variable handle
  int Var;
  int tcVar; // time coordinate variable
  int gmVar; // the grid mapping variable

  // 1D Coordinate
  int vartype1D_dimids[1];

  // 3D Ecosystem
  int vartype3D_dimids[3];

  // 4D Soil
  int vartypeSoil4D_dimids[4];

  // 4D Veg - PFT or compartment but not both
  int vartypeVeg4D_dimids[4];

  // 5D Veg - PFT and PFT compartment
  int vartypeVeg5D_dimids[5];
 
  // Ingest output specification file, create OutputSpec for each entry.
  while(std::getline(output_csv, s)){ 

    std::istringstream ss(s);

    std::string units;
    std::string desc;

    OutputSpec new_spec;
    new_spec.pft = false;
    new_spec.compartment = false;
    new_spec.layer = false;
    new_spec.yearly = false;
    new_spec.monthly = false;
    new_spec.daily = false;
    new_spec.dim_count = 3; // All variables have time, y, x

    for(int ii=0; ii<10; ii++){
      std::getline(ss, token, ',');

      //If the token is the string that indicates an invalid option in
      //  the output_spec file, skip it.
      if(token.compare(invalid_option) == 0){
        continue;
      }
      //If the token is completely whitespace, we assume an erroneous
      //  entry in the output_spec file and skip it.
      if(temutil::all_whitespace(token)){
        continue;
      }

      if(ii==0){ // Variable name
        name = token; 
      }
      else if(ii==1){ // Short description
        desc = token;
      }
      else if(ii==2){ // Units
        units = token;
      }
      else if(ii==3){//Yearly
        if(token.length()>0){
          timestep = "yearly";
          new_spec.yearly = true;
        }
      }
      else if(ii==4){ // Monthly
        if(token.length()>0){
          timestep = "monthly";
          new_spec.monthly = true;
          new_spec.yearly = false;
        }
      }
      else if(ii==5){ // Daily
        if(token.length()>0){
          timestep = "daily";
          new_spec.daily = true;
          new_spec.monthly = false;
          new_spec.yearly = false;
        }
      }
      else if(ii==6){ // PFT
        if(token.length()>0){
          new_spec.pft = true;
          new_spec.dim_count++;
        }
      }
      else if(ii==7){ // Compartment
        if(token.length()>0){
          new_spec.compartment = true;
          new_spec.dim_count++;
        }
      }
      else if(ii==8){ // Layer
        //This will trigger if the layer option has actually been
        // specified by the user or if it is a variable only available
        // by layer with the "forced" marker in the csv file.
        if(token.length()>0){
          new_spec.layer = true;
          new_spec.dim_count++;
        }
      }
      else if(ii==9){ // Data type
        if(token.length()>0){
          if(token.find("int")!=std::string::npos){
            new_spec.data_type = NC_INT;
          }
          else if(token.find("float")!=std::string::npos){
            new_spec.data_type = NC_FLOAT;
          }
          else{
            new_spec.data_type = NC_DOUBLE;
          }
        }
      }
    } // end looping over tokens (aka columns) in a line

    // v5p: the *_DS outputs need model_settings.daily_c_split = 1 (FATAL otherwise).
    {
      std::string dcs_err = dcsplit::request_error(name,
          new_spec.yearly || new_spec.monthly || new_spec.daily, new_spec.daily, this->daily_c_split);
      if (!dcs_err.empty()) {
        std::cerr << dcs_err << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // Only create a file if a timestep is specified for the variable.
    // Otherwise, assume the user does not want that variable output.
    if(new_spec.yearly || new_spec.monthly || new_spec.daily){

      // File location information for reconstructing a complete path
      // and filename during output.
      new_spec.file_path = output_base.string();
      new_spec.filename_prefix = name + "_" + timestep;
      new_spec.var_name = name;

      // Temporary name for file creation.
      std::string creation_filename = name + "_" + timestep + "_" + stage + ".nc";

      BOOST_LOG_SEV(glg, debug)<<"Variable: "<<name<<". Timestep: "<<timestep;

      // filename with local path
      boost::filesystem::path output_filepath = output_base / creation_filename;
      // convert path to string for simplicity in the following function calls
      std::string creation_filestr = output_filepath.string();

#ifdef WITHMPI
      // Creating PARALLEL NetCDF file
      BOOST_LOG_SEV(glg, debug)<<"Creating a parallel output NetCDF file " << creation_filestr;
      temutil::nc( nc_create_par(creation_filestr.c_str(), NC_CLOBBER|NC_NETCDF4|NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid), creation_filestr );
#else
      // Creating NetCDF file
      BOOST_LOG_SEV(glg, debug) << "Creating an output NetCDF file " << creation_filestr;
      temutil::nc( nc_create(creation_filestr.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid), creation_filestr );
#endif

      BOOST_LOG_SEV(glg, debug) << "Adding file-level attributes";
      temutil::nc( nc_put_att_text(ncid, NC_GLOBAL, "Git_SHA", strlen(GIT_SHA), GIT_SHA ) );

      temutil::nc( nc_put_att_text(ncid, NC_GLOBAL, "run_name", this->run_name.length(), this->run_name.c_str() ) );
      temutil::nc( nc_put_att_text(ncid, NC_GLOBAL, "run_description", this->run_description.length(), this->run_description.c_str() ) );
      if (dcsplit::is_ds_name(name)) { // v5p: label the statistical daily split
        const std::string dsnote = "Statistical daily split of the model's monthly C flux "
          "(model_settings.daily_c_split = 1, output only): X_d = M*w_d/sum(w), w = the model's own "
          "environmental scalar with that day's drivers; the days of a month sum to the monthly value. "
          "Not a daily process model. See docs/daily-mode-REFERENCES.md.";
        temutil::nc( nc_put_att_text(ncid, NC_GLOBAL, "daily_c_split_note", dsnote.length(), dsnote.c_str()) );
      }

      //Calculating total timesteps
      int stage_timestep_count = 0;
      if(timestep == "yearly"){
        stage_timestep_count = stage_year_count;
      }
      else if(timestep == "monthly"){
        stage_timestep_count = stage_year_count*12;
      }
      else if(timestep == "daily"){
        stage_timestep_count = stage_year_count*365;
      }

      BOOST_LOG_SEV(glg, debug) << "Adding dimensions...";
      // All variables will have dimensions: time, y, x
      temutil::nc( nc_def_dim(ncid, "time", stage_timestep_count, &timeD) );
      temutil::nc( nc_def_dim(ncid, "y", ysize, &yD) );
      temutil::nc( nc_def_dim(ncid, "x", xsize, &xD) );

      // System-wide variables
      if(new_spec.dim_count==3){
        vartype3D_dimids[0] = timeD;
        vartype3D_dimids[1] = yD;
        vartype3D_dimids[2] = xD;

        temutil::nc( nc_def_var(ncid, name.c_str(), new_spec.data_type, 3, vartype3D_dimids, &Var) );
#ifdef WITHMPI
        //Instruct HDF5 to use independent parallel access for this variable
        temutil::nc( nc_var_par_access(ncid, Var, NC_INDEPENDENT) );
#endif
      }

      // PFT specific dimensions
      else if(new_spec.pft && !new_spec.compartment){
        temutil::nc( nc_def_dim(ncid, "pft", NUM_PFT, &pftD) );

        vartypeVeg4D_dimids[0] = timeD;
        vartypeVeg4D_dimids[1] = pftD;
        vartypeVeg4D_dimids[2] = yD;
        vartypeVeg4D_dimids[3] = xD;

        temutil::nc( nc_def_var(ncid, name.c_str(), new_spec.data_type, 4, vartypeVeg4D_dimids, &Var) );
#ifdef WITHMPI
        //Instruct HDF5 to use independent parallel access for this variable
        temutil::nc( nc_var_par_access(ncid, Var, NC_INDEPENDENT) );
#endif
      }

      // PFT compartment only
      else if(!new_spec.pft && new_spec.compartment){
        temutil::nc( nc_def_dim(ncid, "pftpart", NUM_PFT_PART, &pftpartD) );

        vartypeVeg4D_dimids[0] = timeD;
        vartypeVeg4D_dimids[1] = pftpartD;
        vartypeVeg4D_dimids[2] = yD;
        vartypeVeg4D_dimids[3] = xD;

        temutil::nc( nc_def_var(ncid, name.c_str(), new_spec.data_type, 4, vartypeVeg4D_dimids, &Var) );
#ifdef WITHMPI
        //Instruct HDF5 to use independent parallel access for this variable
        temutil::nc( nc_var_par_access(ncid, Var, NC_INDEPENDENT) );
#endif
      }

      // PFT and PFT compartments
      else if(new_spec.pft && new_spec.compartment){ 
        temutil::nc( nc_def_dim(ncid, "pft", NUM_PFT, &pftD) );
        temutil::nc( nc_def_dim(ncid, "pftpart", NUM_PFT_PART, &pftpartD) );

        vartypeVeg5D_dimids[0] = timeD;
        vartypeVeg5D_dimids[1] = pftpartD;
        vartypeVeg5D_dimids[2] = pftD;
        vartypeVeg5D_dimids[3] = yD;
        vartypeVeg5D_dimids[4] = xD;

        temutil::nc( nc_def_var(ncid, name.c_str(), new_spec.data_type, 5, vartypeVeg5D_dimids, &Var) );
#ifdef WITHMPI
        //Instruct HDF5 to use independent parallel access for this variable
        temutil::nc( nc_var_par_access(ncid, Var, NC_INDEPENDENT) );
#endif
      }

      // Soil specific dimensions
      else if(new_spec.layer){
        temutil::nc( nc_def_dim(ncid, "layer", MAX_SOI_LAY, &layerD) );

        vartypeSoil4D_dimids[0] = timeD;
        vartypeSoil4D_dimids[1] = layerD;
        vartypeSoil4D_dimids[2] = yD;
        vartypeSoil4D_dimids[3] = xD;

        temutil::nc( nc_def_var(ncid, name.c_str(), new_spec.data_type, 4, vartypeSoil4D_dimids, &Var) );
#ifdef WITHMPI
        //Instruct HDF5 to use independent parallel access for this variable
        temutil::nc( nc_var_par_access(ncid, Var, NC_INDEPENDENT) );
#endif
      }

      //  monthly tr :  days since 1901, copy attributes and data from hist climate
      //  yearly tr :  days since 1901, copy attributes from hist climate; generate data...
      //
      //  monthly sc : days since 2001, copy attributes and data from proj climate
      //  yearly sc : days since 2001, copy attributes from proj climate, generate data...

      if (stage == "tr" && (timestep == "yearly" || timestep == "monthly")) {
        vartype1D_dimids[0] = timeD;
        temutil::nc( nc_def_var(ncid, "time", NC_DOUBLE, 1, vartype1D_dimids, &tcVar) );

        int hist_climate_ncid;
        int hist_climate_tcV;

        BOOST_LOG_SEV(glg, debug) << "Opening historic climate file: "
                                  << this->hist_climate_file;
        temutil::nc( nc_open(this->hist_climate_file.c_str(), NC_NOWRITE, &hist_climate_ncid), this->hist_climate_file );
        temutil::nc( nc_inq_varid(hist_climate_ncid, "time", &hist_climate_tcV));

        // Copy attributes for time variable
        temutil::nc( nc_copy_att(hist_climate_ncid, hist_climate_tcV, "units", ncid, tcVar));
        temutil::nc( nc_copy_att(hist_climate_ncid, hist_climate_tcV, "calendar", ncid, tcVar));

        // perhaps write calendar attribute if it does not exist?
        //std::string s = "365_day";
        //temutil::nc( nc_put_att_text(ncid, tcVar, "calendar", s.length(), s.c_str()) );

      }

      if (stage == "sc" && (timestep == "yearly" || timestep == "monthly")) {
        vartype1D_dimids[0] = timeD;
        temutil::nc( nc_def_var(ncid, "time", NC_DOUBLE, 1, vartype1D_dimids, &tcVar) );

        int proj_climate_ncid;
        int proj_climate_tcV;

        BOOST_LOG_SEV(glg, debug) << "Opening projected climate file: "
                                  << this->proj_climate_file;
        temutil::nc( nc_open(this->proj_climate_file.c_str(), NC_NOWRITE, &proj_climate_ncid), this->proj_climate_file );
        temutil::nc( nc_inq_varid(proj_climate_ncid, "time", &proj_climate_tcV));

        temutil::nc( nc_copy_att(proj_climate_ncid, proj_climate_tcV, "units", ncid, tcVar));
        temutil::nc( nc_copy_att(proj_climate_ncid, proj_climate_tcV, "calendar", ncid, tcVar));

        // perhaps write calendar attribute if it does not exist?
        //std::string s = "365_day";
        //temutil::nc( nc_put_att_text(ncid, tcVar, "calendar", s.length(), s.c_str()) );

      }

      BOOST_LOG_SEV(glg, debug) << "Adding variable-level attributes from output spec.";

      // Swap out generic 'time' unit from output spec file with appropriate
      // timestep denomination.
      if (units.find("time") != std::string::npos) {
        BOOST_LOG_SEV(glg, debug) << "Updating units string! ";
        std::string real_timestep;
        if (timestep == "yearly") { real_timestep = "year"; }
        if (timestep == "monthly") { real_timestep = "month"; }
        if (timestep == "daily") { real_timestep = "day"; }
        units.replace(units.begin() + units.find("time"),
                      units.end(),
                      real_timestep.begin(),
                      real_timestep.end());

      }
      BOOST_LOG_SEV(glg, debug) << "Using units string: " << units.c_str();
      temutil::nc( nc_put_att_text(ncid, Var, "units", units.length(), units.c_str()) );
      temutil::nc( nc_put_att_text(ncid, Var, "long_name", desc.length(), desc.c_str()) );

      if(new_spec.data_type == NC_INT){
        temutil::nc( nc_put_att_int(ncid, Var, "_FillValue", NC_INT, 1, &MISSING_I) );
      }
      if(new_spec.data_type == NC_FLOAT){
        temutil::nc( nc_put_att_float(ncid, Var, "_FillValue", NC_FLOAT, 1, &MISSING_F) );
      }
      else if(new_spec.data_type == NC_DOUBLE){
        temutil::nc( nc_put_att_double(ncid, Var, "_FillValue", NC_DOUBLE, 1, &MISSING_D) );
      }

      if (copy_grip_mapping) {

        // open the vegetation input file
        int gmsrcid;
        BOOST_LOG_SEV(glg, debug) << "Opening vegetation file: "
                                  << this->veg_class_file;
        temutil::nc( nc_open(this->veg_class_file.c_str(), NC_NOWRITE, &gmsrcid), this->veg_class_file );

        // Figure out which id is for grid mapping variable
        int srcgmvid = -1;
        srcgmvid = temutil::get_gridmapping_vid(gmsrcid);

        if (srcgmvid >= 0) {

          temutil::copy_var(gmsrcid, srcgmvid, ncid);

        } else {
          BOOST_LOG_SEV(glg, warn) << "A grid_mapping variable is not"
                                   << " available for copying!! Your output"
                                   << " files will not be spatially referenced!";
        }

        // close the file we were reading the grid mapping from
        temutil::nc( nc_close(gmsrcid) );
      }

      /* End Define Mode (not strictly necessary for netcdf 4) */
      BOOST_LOG_SEV(glg, debug) << "Trying to leaving 'define mode'...";
      try {
        temutil::nc( nc_enddef(ncid) );
      } catch (const temutil::NetCDFDefineModeException& e) {
        BOOST_LOG_SEV(glg, info) << "Error ending define mode: " << e.what();
      }

      /* Fill out the time coordinate variable */
      if ((stage == "tr" || stage == "sc") && timestep == "yearly") {
        BOOST_LOG_SEV(glg, debug) << "Time coordinate variable, tr or sc, yearly.";

        int tcV;
        temutil::nc( nc_inq_varid(ncid, "time", &tcV));

        int runyrs;
        if (stage == "tr") { runyrs = this->tr_yrs; }
        if (stage == "sc") { runyrs = this->sc_yrs; }

        std::vector<int> time_coord_values(runyrs, 0);
        for (int i=0; i < runyrs; i++) {
          time_coord_values.at(i) = 365 * i; // assumes attribute calendar:'365_day'
        }
        size_t start[1];
        size_t count[1];

        start[0] = 0;
        count[0] = time_coord_values.size();

        temutil::nc( nc_put_vara_int(ncid, tcV, start, count, &time_coord_values[0]) );

      }

      if ((stage == "tr" || stage == "sc") && timestep == "monthly") {
        BOOST_LOG_SEV(glg, debug) << "Time coordinate variable, tr or sc, monthly.";

        int runyrs = 0;
        std::vector<int> full_time_coord;

        if (stage == "tr") {
          runyrs = this->tr_yrs;
          full_time_coord = temutil::get_timeseries2(this->hist_climate_file, "time", 0);
        }
        if (stage == "sc") {
          runyrs = this->sc_yrs;
          full_time_coord = temutil::get_timeseries2(this->proj_climate_file, "time", 0);
        }

        int tcV;
        temutil::nc( nc_inq_varid(ncid, "time", &tcV));

        size_t start[1];
        size_t count[1];

        start[0] = 0;
        count[0] = runyrs * 12;

        if (stage == "sc" && this->daily_proj_climate
            && full_time_coord.size() < size_t(runyrs) * 365) {
          std::cerr << "FATAL (daily climate mode): projected climate file has "
                    << full_time_coord.size() << " daily time values; " << runyrs
                    << " scenario years need " << runyrs * 365 << std::endl;
          exit(EXIT_FAILURE);
        }
        if (stage == "sc" && this->daily_proj_climate) {
          // Daily climate mode: the projected climate file has one time value
          // per day (365_day calendar), so label each monthly output step with
          // the file's time value for the first day of that month.
          std::vector<int> month_time_coord(runyrs * 12, 0);
          for (int iy = 0; iy < runyrs; iy++) {
            int doy = 0;
            for (int im = 0; im < 12; im++) {
              month_time_coord.at(iy * 12 + im) = full_time_coord.at(iy * 365 + doy);
              doy += DINM[im];
            }
          }
          temutil::nc( nc_put_vara_int(ncid, tcV, start, count, &month_time_coord[0]) );
        } else {
          temutil::nc( nc_put_vara_int(ncid, tcV, start, count, &full_time_coord[0]) );
        }

      }


      /* Close file. */
      BOOST_LOG_SEV(glg, debug) << "Closing new file...";
      temutil::nc( nc_close(ncid) );

      // Add OutputSpec objects to the map tracking the appropriate timestep
      if(new_spec.daily){
        daily_netcdf_outputs.insert(std::map<std::string, OutputSpec>::value_type(name, new_spec));
      }

      else if(new_spec.monthly){
        monthly_netcdf_outputs.insert(std::map<std::string, OutputSpec>::value_type(name, new_spec));
        //monthly_netcdf_outputs.insert({name, filename}); c++11
      }

      else if(new_spec.yearly){
        yearly_netcdf_outputs.insert(std::map<std::string, OutputSpec>::value_type(name, new_spec));
      }

    } // End file creation section (if timestep is specified)
  } // End looping over the output spec file.
}





