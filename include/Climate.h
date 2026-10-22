#ifndef _CLIMATE_H_
#define _CLIMATE_H_

#include <string>
#include <vector>

class Climate {
public:

  Climate();
  Climate(const std::string& fname, const std::string& co2fname, int y, int x);

  // Misc. climate variables
  // The following two values determine the span of historic climate
  //  years used to produce a baseline average climate for EQ
  //  TODO: SP repeated climate is not currently affected by these.
  int baseline_start;
  int baseline_end;
  // These values will change during a run when switching from
  //  historic to projected climate data
  int tseries_start_year;
  int tseries_end_year;

  // Opt-in daily climate mode (see ModelData::climate_timestep).
  // daily_proj: set from config; when true, load_proj_climate() expects a
  //   daily (365 steps/yr) file.
  // daily_series_loaded: true once a daily series has actually been loaded
  //   into tair/prec/nirr/vapo (and cld/par); only then does
  //   prepare_daily_driving_data() take the direct daily path.
  bool daily_proj = false;

  // Fix (a) (v5i), OFF by default: model_settings.vp_hpa_fix. The vapor_press
  // input is in hPa (docs model_overview.rst; CRU and the TVC/HPC MIP inputs,
  // whose files carry units = "hPa"), but svp is computed in Pa, so the legacy
  // vpd = svp(Pa) - vp(hPa) is ~svp and rh = vp/svp is ~1/100 of its true value.
  // 1 = convert vp to Pa (x100) for vpd (Penman-Monteith ET, stomatal VPD
  // limitation) and for EnvData d_atmd.vp (non-vascular rh in Cohort). The
  // vapo / vapo_d series themselves (and the VAPO output) stay in hPa.
  // 0 (default) = legacy units mix.
  static int vp_hpa_fix;
  // Refs: input vapour pressure in hPa (CRU TS vap: Harris et al. 2020, doi:10.1038/s41597-020-0453-3); VPD = es - ea in consistent units (FAO-56, Allen et al. 1998); Penman-Monteith (Monteith 1965, PMID 5321565).
  // Fix B1 (v5k), OFF by default: model_settings.svp_ice_fix. For t <= 0 C the
  // legacy saturation vapour pressure uses 17.27/(272.62+t), i.e. the Tetens
  // water coefficient with the WMO ice denominator (+51 % vs ice at -20 C).
  // 1 = WMO (2008) Magnus form over ice, 6.112 exp(22.46 t/(272.62+t)) hPa, for
  // t <= 0 C, in calculate_saturated_vapor_pressure() and in the SEB esat.
  static int svp_ice_fix;
  // Refs (B1): WMO-No. 8 (2008) CIMO Guide Annex 4.B, e_i = 6.112 exp(22.46 t/(272.62+t)) hPa; cf. Buck 1981, doi:10.1175/1520-0450(1981)020<1527:NEFCVP>2.0.CO;2; test values Murphy & Koop 2005, doi:10.1256/qj.04.94.
  bool daily_series_loaded = false;

  // Daily (365) extraterrestrial radiation, used in daily mode only.
  std::vector<float> girr_daily;

  // Option 4 (v5j) SEB: only filled when seb_extra is true.
  // lwdown (W/m2) and wind (m/s): full daily series (optional variables
  // 'lwdown' / 'wind' in the daily projected-climate file; LW parameterized
  // when absent or seb_lw_param == 1; wind = seb_wind_const when absent).
  bool seb_extra = false;
  int seb_lw_param = 0;
  float seb_wind_const = 3.0;
  bool seb_lw_from_file = false;
  bool seb_wind_from_file = false;
  std::vector<float> lwdown;
  std::vector<float> wind;
  std::vector<float> lwdown_d;
  std::vector<float> wind_d;
  void seb_load_extra(const std::string& fname, int y, int x);

  // driving variables
  std::vector<float> co2;
  std::vector<float> tair;
  std::vector<float> prec;
  std::vector<float> nirr;
  std::vector<float> vapo;

  // calculated driving variables
  std::vector<float> girr;
  std::vector<float> cld;
  std::vector<float> par;

  // simplified climate, averaged over some number of years
  std::vector<float> avgX_tair;
  std::vector<float> avgX_prec;
  std::vector<float> avgX_nirr;
  std::vector<float> avgX_vapo;

  // repetitive climate, repeated over a specified number of years
  //std::vector<float> rptX_tair;
  //std::vector<float> rptX_prec;
  //std::vector<float> rptX_nirr;
  //std::vector<float> rptX_vapo;

  // Daily containers

  // this is not really a daily value, but for scope/access
  // reasons it is easier to treat it as such. Every day of the year has the same co2 value...
  float co2_d;

  //  ->> should be interpolated from the monthly containers
  std::vector<float> tair_d;
  std::vector<float> prec_d;
  std::vector<float> nirr_d;
  std::vector<float> vapo_d;
  std::vector<float> rain_d;
  std::vector<float> snow_d;
  std::vector<float> girr_d;
  std::vector<float> cld_d;
  std::vector<float> par_d;

  // More daily containers, but unlike the others, these are not
  // ever stored at the monthly level.
  std::vector<float> rhoa_d;
  std::vector<float> svp_d;
  std::vector<float> vpd_d;
  std::vector<float> dersvp_d;
  std::vector<float> abshd_d;
  
  std::vector<float> monthly2daily(const std::vector<float>& mly_vals);

  void prepare_eq_daily_driving_data(int iy, const std::string& stage);
  void prepare_daily_driving_data(int, const std::string&);

  void monthlycontainers2log();
  void dailycontainers2log();

  void load_proj_climate(const std::string&, int, int);
  void load_proj_co2(const std::string&);
  // Scenario stage: index CO2 by calendar year (see Climate.cpp).
  void index_proj_co2_by_year(const std::string& proj_climate_fname,
                              const std::string& proj_co2_fname,
                              int sc_yrs);
  float co2_for_year_index(int iy) const;
  int co2_idx_offset = 0; // added to the model year index when reading co2


  void prep_avg_climate();

private:

  void load_from_file(const std::string& fname, int y, int x);
  void load_daily_from_file(const std::string& fname, int y, int x);
  void prepare_daily_driving_data_from_daily(int iy);
  void derive_daily_secondary_vars();

  void split_precip();

  std::vector<float> avg_over(const std::vector<float> & var, const int start_yr, const int end_yr);
  std::vector<float> rpt_years(const std::vector<float> &, const int);

  std::vector<float> interpolate_daily(const std::vector<float> & var);

  std::vector<float> eq_range(const std::vector<float>& data);
  std::vector<float> interpolation_range(const std::vector<float>& data, int year);

};

#endif /* _CLIMATE_H_ */


