#ifndef COHORT_H_
#define COHORT_H_

#include "Climate.h"

#include "Ground.h"
#include "Vegetation.h"

#include "Vegetation_Env.h"
#include "Vegetation_Bgc.h"

#include "Snow_Env.h"
#include "Soil_Env.h"
#include "SoilParent_Env.h"
#include "Soil_Bgc.h"

#include "WildFire.h"
#include "Thermokarst.h"
#include "ThermokarstData.h"

#include "CohortData.h"

#include "EnvData.h"
#include "BgcData.h"
#include "FireData.h"

#include "RestartData.h"

#include "CohortLookup.h"

#include "Integrator.h"

// headers for run
#include "ModelData.h"
#include "DailyCSplit.h"

class Cohort {
public :
  Cohort();
  Cohort(int y, int x, ModelData* modeldatapointer);
  ~Cohort();
  
  int y;
  int x;

  float lon;
  float lat;

  // model running status
  int errorid;
  bool failed;    // when an exception is caught, set failed to be true

  /*
    Note: FRI is a member of CohortData because it is checked in
    Soil_Bgc::prepareintegration(...), and at that point there is no access to
    the members/fields of a Cohort...
  */

  // old? can I deprecate these??
  //double pfsize[NUM_FSIZE];
  //double pfseason[NUM_FSEASON];
  
  //inputs
  CohortLookup chtlu;

  // domain
  Vegetation veg;
  Ground ground;
  
  // new domain
  Climate climate;
  bool snowcls_wind_warned = false; // snow_scheme = 1: warned once that the forcing has no wind

  // processes
  Vegetation_Env vegenv[NUM_PFT];
  Snow_Env snowenv;
  Soil_Env soilenv;
  SoilParent_Env solprntenv;

  Vegetation_Bgc vegbgc[NUM_PFT];
  Soil_Bgc soilbgc;

  WildFire fire;

  // data
  EnvData ed[NUM_PFT];
  BgcData bd[NUM_PFT];
  EnvData * edall;
  BgcData * bdall;

  FirData year_fd[12]; //Monthly fire data, for all PFTs and soil
  FirData * fd;   //Fire data for an individual month 

  // ---- RTS / thermokarst (model_settings.thermokarst; port of rts-regrow 7c9fd634)
  Thermokarst thermokarst;          // constructed from IO.hist_exp_thermokarst_file only when thermokarst = 1
  ThermokarstData year_tkdata[12];  // monthly RTS fluxes/diagnostics (outputs TK*)
  ThermokarstData tk_own;           // used when setProcessData() gets no ThermokarstData (unit tests)
  ThermokarstData * tkdata = &tk_own;
  double tk_prev_ald = UIN_D;       // previous year's ALD (m), updated in December only when thermokarst = 1

  ModelData * md;

  CohortData cd;
  RestartData restartdata;
  

//  void NEW_load_climate_from_file(int y, int x);
//  void NEW_load_veg_class_from_file(int y, int x);
//  void NEW_load_fire_from_file(int y, int x);

  void initialize_internal_pointers();

  void setModelData(ModelData* md);
  void setProcessData(EnvData * alledp, BgcData * allbdp, FirData *fdp);
  void setProcessData(EnvData * alledp, BgcData * allbdp, FirData *fdp, ThermokarstData * tdp);

  void initialize_state_parameters();
  //void prepareAllDrivingData();
  //void prepareDayDrivingData(const int & yrcnt, const int &usedatmyr);
  void updateMonthly(const int & yrcnt, const int & currmind,
                     const int & dinmcurr, std::string stage);
  
  void set_state_from_restartdata();
  void set_fire_state_from_restartdata();   // fire-balance fixes
  void set_restartdata_from_state();

  // Overwrites/fills climate containers with data from projected climate file
  void load_proj_climate(const std::string&);
  // Overwrites/fills co2 container with data from projected co2 input file
  void load_proj_co2(const std::string& proj_co2_file);
  // Overwrites/fills explicit fire data containers with data from the projected fire input file.
  void load_proj_explicit_fire(const std::string& proj_exp_fire_file);
  void load_proj_explicit_thermokarst(const std::string& proj_exp_thermokarst_file); // thermokarst = 1 only
  void set_tk_state_from_restartdata();   // RTS state (thermokarst = 1 only)

  // Integrator safety-net / give-up counters (read-only access for reporting)
  const Integrator& get_solintegrator() const { return solintegrator; }
  int get_veg_blackhole_count() const {
    int n = 0;
    for (int ip = 0; ip < NUM_PFT; ip++) { n += vegintegrator[ip].blackhole_count; }
    return n;
  }

  // v5p daily_c_split (output only): daily ecosystem totals of the current
  // month (g/m2/day, dcsplit::Var order, days 0..dinm-1) and the month's
  // fallback flags. Written by updateMonthly_DailyCSplit, read by Runner.
  double ds_daily[dcsplit::NVAR][31] = {};
  int ds_flags = 0;
  long ds_fallback_months = 0;  // months with any fallback flag (cumulative)
  long ds_remap_months = 0;     // months with the midpoint layer remap (cumulative)
  void load_ch4(const std::string& ch4_file);

  // RTS regrow / floor-RH restart state (layout version 1; see the private tk_rg_* members)
  void tk_rg_pack(double * s) const;
  void tk_rg_unpack(const double * s);

private:
  // v5p daily_c_split: soil layer structure of the Env loop and of the Bgc step
  int ds_env_numsl = 0, ds_bgc_numsl = 0;
  double ds_env_z[MAX_SOI_LAY] = {}, ds_env_dz[MAX_SOI_LAY] = {};
  double ds_bgc_z[MAX_SOI_LAY] = {}, ds_bgc_dz[MAX_SOI_LAY] = {};
  int ds_bgc_type[MAX_SOI_LAY] = {};
  void ds_begin_month();
  void ds_snapshot_soil(int& numsl, double* z, double* dz, int* type);
  void updateMonthly_DailyCSplit(const int & dinm);

  Integrator vegintegrator[NUM_PFT];
  Integrator solintegrator;


  void updateMonthly_DIMveg(const int & currmind, const bool & dynamic_lai_module);
  void updateMonthly_DIMgrd(const int & currmind, const bool & dslmodule);

  void updateMonthly_Env(const int & currmind, const int & dinmcurr);
  void updateMonthly_Bgc(const int & currmind, const int & dinmcurr);
  void updateMonthly_Dsb(const int & yrcnt, const int & currmind, std::string stage);

  // Fire is a type of disturbance
  void updateMonthly_Fir(const int & year, const int & midx, std::string stage);
  // ---- RTS (thermokarst = 1 only)
  void updateMonthly_Thermokarst(const int &year, const int &midx, std::string stage);
  void updateMonthly_RTSMethane(const int & midx);   // RTS CH4 (thermokarst_ch4)
  void updateMonthly_RTSFloorPrep();                  // per-layer RH modifiers before soil integration
  void updateMonthly_RTSFloorExport();                // monthly lateral export of floor mineral C/N
  void updateMonthly_RTSMethane_lateral();           // CO2/CH4 split of TKLATC
  int tk_cur_month = 0;               // calendar month of the current step (RTS floor rmin)
  double tk_rh_hist[10][12] = {{0}};  // pre-initiation monthly RHSOM ring buffer (years x months)
  int tk_rh_nyr[12] = {0};            // number of recorded years per calendar month (<= 10)
  int tk_rh_next[12] = {0};           // next ring slot per calendar month
  // RTS post-slump regrowth (RTSRegrow.h); this state and tk_rh_* go to the restart as
  // tk_regrow_state (thermokarst = 1 only; tk_rg_pack / tk_rg_unpack, RestartData::TK_RG_NSTATE)
  bool tk_rg_have = false;            // pre-slump state recorded
  double tk_rg_cpre[NUM_PFT][3] = {{0}}, tk_rg_npre[NUM_PFT][3] = {{0}};
  int tk_rg_type[NUM_PFT] = {0};
  double tk_rg_s = 1.0, tk_rg_gpp = 0.0, tk_rg_olpre = 0.0, tk_rg_olstab = 0.0, tk_rg_oladd = 0.0, tk_rg_ncarry = 0.0;
  double tk_rg_gpp12[12] = {0}; int tk_rg_gpp12_n = 0, tk_rg_gpp12_i = 0;  // trailing-12-month GPP ring (N taper)
  double tk_rg_ingpp = 0.0, tk_rg_ingpp12[12] = {0};  // pre-slump / trailing N-unlimited GPP (dynamic taper reference)
  double tk_rg_taper = 1.0;           // current N-supply taper factor in [0,1] (1 until GPP recovers)
  double tk_gpp_ring[60] = {0}, tk_ingpp_ring[60] = {0}; int tk_gpp_n = 0, tk_gpp_i = 0;
  double tk_rg_olc() const;           // organic (non-mineral) soil layer SOM C, g/m2
  void updateMonthly_RTSRegrowVeg();  // shrub boost only (>=1) before the veg integration; no growth ceiling in v2
  void updateMonthly_RTSRegrow();     // reseed + organic-layer rebuild (thermokarst step)
  void rebuildBdallVeg();
  soidiag_thermokarst tk_ch4_month; // this month's CH4 (kept across the slump step, which resets tkdata)
  void tk_set_month_flags();        // base fixes the RTS needs, on only after an RTS initiation

  // update root distribution
  void getSoilFineRootFrac_Monthly();
  double assignSoilLayerRootFrac(const double & topz, const double & botz,
                                 const double csumrootfrac[MAX_ROT_LAY],
                                 const double dzrotlay[MAX_ROT_LAY]);

  //
  void assignAtmEd2pfts_daily();
  void assignGroundEd2pfts_daily();
  void getSoilTransfactor4all_daily();
  void getEd4allveg_daily();
  void getEd4allgrnd_daily();
  void getEd4land_daily();

  void assignSoilBd2pfts_monthly();
  void getBd4allveg_monthly();

};
#endif /*COHORT_H_*/
