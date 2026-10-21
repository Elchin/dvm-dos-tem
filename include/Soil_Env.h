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
  double frozen_infil_mm; // water infiltrated into frozen soil today (mm)
  int soil_evap_beta;     // DRAFT: see ModelData::soil_evap_beta (0 = off)
  double evap_unapplied_mm;
  // Fix (c) (v5j): per-layer flag (index as root_water_up) = 1 if today's root
  // uptake from that layer was actually applied by updateDailySM, else 0.
  // Filled every call; only used when model_settings.nonvasc_trans_withdraw = 1.
  double uptake_applied[MAX_SOI_LAY];
  // Option 2 (v5h): lateral drainage of thawed near-surface
  // organic layers towards the peat retention capacity; see ModelData.h.
  double org_drain_rate, org_drain_thre_f, org_drain_thre_h, org_drain_depth;
  double orgLateralDrain(); // returns mm drained today

  // Option 4 (v5j) SEB (see ModelData::seb_scheme). 0 = off.
  int seb_scheme;
  double seb_z0_snow, seb_z0_soil, seb_zwind, seb_ztemp, seb_emis_snow,
         seb_emis_soil, seb_wind, seb_minwind, seb_stab_b, seb_tran_frac,
         seb_pres;
  double seb_fc_last; // last canopy LW-obscured fraction (from SW transmission)
  Snow_Env* snowenvp; // for the snow albedo
  SebDiag sebdiag;
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
