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
