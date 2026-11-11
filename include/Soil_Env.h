#ifndef SOIL_ENV_H_
#define SOIL_ENV_H_

#include "Stefan.h"
#include "Richards.h"
#include "TemperatureUpdator.h"

#include "CohortData.h"
#include "EnvData.h"
#include "FireData.h"
#include "RestartData.h"

#include "errorcode.h"
#include "parameters.h"
#include "layerconst.h"
#include "CohortLookup.h"

#include "Ground.h"
#include "EnthalpySolver.h"

class Snow_Env;

// Option 4 (v5j) SEB: diagnostics of the last daily SEB solve.
// Energy-closure diagnostic (model_settings.energy_diag 1; output only): enthalpy of the
// snow + soil + rock column (J/m2, relative to liquid water at 0 C) before (h0) and after (h1)
// the daily thermal step updateDailyGroundT(), the solver's top/bottom boundary heat
// (gtop/gbot, J/m2, + = into the column) and the thermal-step residual
// res = (h1 - h0) - (gtop + gbot). nfront = number of freeze/thaw fronts after the step.
struct EnergyDiag {
  double h0 = 0.0, h1 = 0.0, gtop = 0.0, gbot = 0.0, res = 0.0, tsurf = 0.0;
  int nfront = 0, nsnow = 0, used = 0;
  // snowmk: enthalpy change of the snow layers inside Stefan::updateFronts (meltingSnowLayer marks
  // melt as liq += ice; the melt water is removed later by Snow_Env, whose degree-day melt supplies
  // that energy outside this step). Excluded from res (phase 2d): res = dH - (gtop + gbot) - snowmk.
  double snowmk = 0.0;
  // Option C (thermal_solver = 1 only; all J/m2 per day, + = into the column):
  //  src    radiation/SEB snow-melt energy applied as a heat source to the top snow cell
  //  clip   enthalpy of ice expelled by the pore-volume cap (and checkWaterValidity rounding)
  //  init   switch-on / restart projection onto the freezing curve at fixed T (first day only)
  //  equil  change of H by putting cells back on the curve (H-conserving; round-off only)
  //  bsnow, bwater, blayer, bfire: enthalpy changes booked by the non-thermal modules since the
  //         previous thermal step (Snow_Env::updateDailyM, Soil_Env::updateDailySM, monthly
  //         layer thickness/redivision, fire), measured around each call (Cohort.cpp)
  //  unattr change of H between two thermal steps not booked by any module (must be 0)
  //  resday full-day residual = res + equil + unattr (pass/fail <= 1e-6 W/m2 in the harness)
  double src = 0.0, clip = 0.0, init = 0.0, equil = 0.0, bsnow = 0.0, bwater = 0.0, blayer = 0.0,
         bfire = 0.0, unattr = 0.0, resday = 0.0;
  // Stage C2 (energy-conserving remap): bsnow and bwater are the PHYSICAL boundary enthalpy of the
  //  snow and soil-water windows (snowfall in as ice at min(Ts,0), sublimation out at the top snow T,
  //  melt water out at 0 C; infiltration in at max(Tair,0) for rain / 0 C for melt, drainage, runoff
  //  and ET out at the mixed enthalpy of the moved water). rsnow, rwater: correction the remap
  //  applied relative to the fixed-T module result (= the energy the legacy fixed-T moves would
  //  have created, + = the module result was too cold). wlayer: water mass change (kg m-2) of the
  //  monthly layer/BGC window (should be 0). remapfb: water windows that fell back to booking.
  double rsnow = 0.0, rwater = 0.0, wlayer = 0.0;
  long remapfb = 0;
  int solves = 0, fallbacks = 0, enth = 0;
};

struct SebDiag {
  double ta, ts, swin, alb, lwsky, lwin, lwout, h, le, g, qmelt, wind, fc, lam, t1;
  int snow, iters, used;
};

class Soil_Env {
public:

  Soil_Env();
  ~Soil_Env();

  soipar_env envpar;

  double root_water_up[MAX_SOI_LAY];

  Richards richards;
  Stefan stefan;
  TemperatureUpdator tempupdator;

  // Supply-limited ET bookkeeping (per day, mm): water the minimum-liquid floor had to add in
  // checkSoilLiquidWaterValidity(), and soil evaporation that was reported but never removed
  // from the soil (frozen-top branch). Cohort subtracts these from the reported ET.
  double floor_added_mm;
  // DRAFT (daily-structural): see ModelData::frozen_infil_eice (<0 = off)
  double frozen_infil_eice;
  // Refs: Swenson et al. 2012, doi:10.1029/2012MS000165; CLM5 Tech Note eq. 2.7.48 (Omega = 6); see docs/daily-mode-REFERENCES.md.
  double frozen_infil_mm; // water infiltrated into frozen soil today (mm)
  int soil_evap_beta;     // DRAFT: see ModelData::soil_evap_beta (0 = off)
  // Refs: Lee & Pielke 1992, doi:10.1175/1520-0450(1992)031<0480:ETSSSH>2.0.CO;2; theta_fc from Clapp & Hornberger 1978, doi:10.1029/WR014i004p00601.
  double evap_unapplied_mm;
  // Fix (c) (v5j): per-layer flag (index as root_water_up) = 1 if today's root
  // uptake from that layer was actually applied by updateDailySM, else 0.
  // Filled every call; only used when model_settings.nonvasc_trans_withdraw = 1.
  double uptake_applied[MAX_SOI_LAY];
  // Option 2 (v5h): lateral drainage of thawed near-surface
  // organic layers towards the peat retention capacity; see ModelData.h.
  double org_drain_rate, org_drain_thre_f, org_drain_thre_h, org_drain_depth;
  // Refs: Quinton & Marsh 1999, doi:10.1002/(SICI)1099-1085(199911)13:16<2563::AID-HYP942>3.0.CO;2-D; Letts et al. 2000, doi:10.1080/07055900.2000.9649643; cf. CLM5 Tech Note sec. 2.7.4.
  double orgLateralDrain(); // returns mm drained today

  // Option 4 (v5j) SEB (see ModelData::seb_scheme). 0 = off.
  int seb_scheme;
  // Refs (SEB): bulk transfer, Brutsaert 1982, doi:10.1007/978-94-017-1497-6; simplified stability functions after Louis 1979, doi:10.1007/BF00117978
  //   (stable) and Dyer 1974, doi:10.1007/BF00240838 (unstable); LW: Prata 1996, doi:10.1002/qj.49712253306; Crawford & Duchon 1999,
  //   doi:10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2; snow emissivity Warren 1982, doi:10.1029/RG020i001p00067; esat WMO-No. 8 (2008);
  //   safeguarded Newton, Press et al. 2007 (Numerical Recipes 3rd ed., sec. 9.4); see docs/daily-mode-REFERENCES.md.
  // Fix B2 (v5k), OFF by default: model_settings.richards_infil_fix. When
  // Richards::update returns early, the legacy code adds the day's infiltration
  // to the puddle/runoff of Richards' private EnvData copy, which nothing reads,
  // so the water is lost (and soil evaporation / uptake are reported although
  // not applied). 1 = the infiltration goes to this cohort's puddle (up to
  // ponding_max_mm) and the rest to surface runoff (qover), as in the frozen-top
  // branch; evaporation is booked as unapplied and the uptake as not applied.
  int richards_infil_fix = 0;
  // Refs (B2): no published source; code bug (Richards::update early return, upstream 419a4eb3, fix commit 07f6c8d0); restores water mass conservation (continuity in Richards 1931, doi:10.1063/1.1745010).
  void rebookUnappliedRichards(const double infil, const double evap);
  double seb_z0_snow, seb_z0_soil, seb_zwind, seb_ztemp, seb_emis_snow,
         seb_emis_soil, seb_wind, seb_minwind, seb_stab_b, seb_tran_frac,
         seb_pres;
  double seb_fc_last; // last canopy LW-obscured fraction (from SW transmission)
  Snow_Env* snowenvp; // for the snow albedo
  SebDiag sebdiag;
  // Energy-closure diagnostic, see EnergyDiag. model_settings.energy_diag: 0 = off (default),
  // 1 = fill energydiag every day (Cohort writes output_dir/energy_daily.csv). Never changes
  // the model state. Refs: enthalpy form with latent heat of fusion, e.g. Bonan 2019
  // (Climate Change and Terrestrial Ecosystem Modeling, ch. 5), Jafarov et al. 2012,
  // doi:10.5194/tc-6-613-2012 (dvmdostem/GIPL soil thermal context); no published source for
  // the dvmdostem-specific residual.
  int energy_diag = 0;
  EnergyDiag energydiag;
  // Option C: enthalpy thermal solver (see ModelData.h; Soil_Env_enthalpy.cpp). thermal_solver
  // 0 = legacy (default; none of the members below is used), 1 = enthalpy solver.
  int thermal_solver = 0, thermal_sfcc = 1, thermal_face_mode = 0, thermal_nl_solver = 0,
      thermal_substeps = 1;
  double geothermal_flux_Wm2 = 0.0;   // bottom heat flux into the column (W m-2), this cohort
  // between-step enthalpy bookings (J m-2), added by Cohort around the non-thermal modules
  double enth_bk_snow = 0.0, enth_bk_water = 0.0, enth_bk_layer = 0.0, enth_bk_fire = 0.0;
  double enth_h_last = 0.0;          // column enthalpy at the end of the last thermal step
  bool enth_have_last = false;
  bool enth_init_pending = true;     // project onto the curve at fixed T on the next step
  long enth_fallbacks = 0;           // steps that needed dt halving (non-convergence)
  long enth_failures = 0;            // steps still not converged after halving
  double enth_init_maxdliq = 0.0;    // max |liq change| of the switch-on projection (kg m-2)
  // Stage C2: energy-conserving remap of the non-thermal windows (thermal_solver 1 only)
  struct EnthSnap {
    std::vector<Layer*> L; std::vector<int> snow; std::vector<double> T, liq, ice, H;
    double hcol = 0.0, ttop = 0.0; bool topsnow = false;
  };
  void enthSnap(EnthSnap& s) const;
  // soil water window: mass moved between cells carries its enthalpy (pool mixing), infiltration
  // enters at h_in (J/kg), the excess leaves at the pool enthalpy; T set from H at fixed liq/ice
  void enthRemapWaterCore(const EnthSnap& s, double infil_kg, double h_in);
  void enthRemapWater(const EnthSnap& s);                  // reads ed (qinfl, rain/melt, Tair)
  // snow window: physical boundary terms, correction spread over the snow cells by heat capacity
  void enthRemapSnowCore(const EnthSnap& s, double tsurface, double sublim_kg, double melt_kg);
  void enthRemapSnow(const EnthSnap& s, const double tdrv); // reads ed (sublim, melt, n-factor)
  static void enthSetT(Layer* l, double h);                // T from H at fixed liq/ice (exact)
  double columnWaterE() const;                             // liq + ice of all non-rock cells (kg m-2)
  double enth_rm_snow = 0.0, enth_rm_water = 0.0, enth_wk_layer = 0.0;
  long enth_remap_fallbacks = 0;
  enthalpy::ColumnSolver enth_solver;
  // enthalpy of one layer as used by the enthalpy solver (J m-2): Cs*T + liq*c_w*T + ice*(c_i*T - L_f),
  // Cs = vhcsolid*(1-poro)*dz (soil), vhcsolid*dz (rock), 0 (snow)
  static double enthLayerCs(Layer* l);
  static double enthLayerH(Layer* l);
  double columnEnthalpyE() const;
  // solver cell of a layer (texture/organic fraction for SFCC 1 as Soil_Bgc::rhFrozenEfr)
  void enthalpyCell(Layer* l, enthalpy::Cell& c) const;
  // radiation / SEB snow-melt energy (W m-2) applied to the top snow cell (0 without snow)
  double enthMeltSource(const double tsurface) const;
  // Core of one daily enthalpy step on the ground column (no ed/cd access; unit-testable):
  // equilibrate, solve, scatter, melt-excess hand-off, ice cap, fronts and layer states.
  void enthalpyThermalStep(const double tsurface, const double qmelt_Wm2, EnergyDiag& q);
  void updateDailyGroundT_enthalpy(const double tsurface);
  void enthCheckFinite(const char* where) const;
  // Enthalpy of one layer, J/m2: (solid heat capacity)*T + liq*c_liq*T + ice*(c_ice*T - L_f),
  // solid part = getMixVolHeatCapa()*dz - liq*c_liq - ice*c_ice (0 for snow).
  static double layerEnthalpy(Layer* l);
  double columnEnthalpy() const;
  double snowEnthalpy() const;   // snow layers only (energy_diag)
  double solveSEB(const double & tdrv);

  void setGround(Ground* grndp);
  void setCohortData(CohortData* cdp);
  void setEnvData(EnvData* edp);
  void setCohortLookup(CohortLookup * chtlup);

  void resetDiagnostic();   /*! reset diagnostic variables to initial values */

  void initializeParameter();
  void initializeState();
  void set_state_from_restartdata(const RestartData & rdata);

  void updateDailyGroundT(const double & tdrv, const double & dayl);
  void updateDailySM(double weighted_veg_tran);

  double getSoilTransFactor(double r_e_ij[MAX_SOI_LAY], Layer* fstsoill,
                            const double vrootfr[MAX_SOI_LAY]);

  void retrieveDailyTM(Layer* toplayer, Layer* lstsoill);
  void checkSoilLiquidWaterValidity(Layer *topsoill, int topind);

  double getWaterTable(Layer* fstsoil);

private:

  Ground * ground;
  CohortData * cd;
  EnvData * ed;
  CohortLookup* chtlu;

  double ponding_max_mm;
  // wetland (phase 3) diagnostics, mm/day: run-on inflow added to infiltration today (0 unless
  // cd->runon_status); groundwater recharge is Richards::qgw_mm.
  double runon_in_mm = 0.0;
  double frozenInfiltration(Layer* fstsoill, double supply_mm);//max ponding (surface water storage) (mm)

  void updateDailySurfFlux(Layer* frontl, const double & dayl);
  void updateDailySoilThermal4Growth(Layer* fstsoill, const double &tsurface);
  void updateLayerStateAfterThermal(Layer* fstsoill, Layer *lstsoill,
                                    Layer* botlayer);

  void retrieveDailyFronts();

  double getEvaporation(const double & dayl, const double &rad);
  double getPenMonET(const double & ta, const double& vpd, const double &irad,
                     const double &rv, const double & rh);
  double getRunoff(Layer* fstsoill, Layer* drainl,
                   const double & rnth,
                   const double & melt);

  // the following codes not used anymore
  double getInflFrozen(Layer *fstminl, const double &  rnth,
                       const double & melt);
  double updateLayerTemp5Lat(Layer* currl, const double & infil);

};

#endif /*SOIL_ENV_H_*/
