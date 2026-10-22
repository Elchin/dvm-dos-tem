#ifndef VEGETATION_ENV_H_
#define VEGETATION_ENV_H_
#include <cmath>

#include "CohortLookup.h"
#include "states.h"
#include "fluxes.h"
#include "diagnostics.h"
#include "parameters.h"

#include "CohortData.h"
#include "EnvData.h"
#include "FireData.h"
#include "RestartData.h"

class Vegetation_Env {
public:

  Vegetation_Env();
  ~Vegetation_Env();

  int ipft;

  vegpar_env envpar;

  void setCohortLookup(CohortLookup * chtlup);
  void setEnvData(EnvData* edatap);
  void setFirData(FirData* fdp);

  void initializeParameter();
  void initializeState();
  void set_state_from_restartdata(const RestartData & rdata);

  void updateRadiation(double leaf_area_index, double foliar_projected_cover);
  void updateWaterBalance(const double &daylhr, double leaf_area_index, double foliar_projected_cover);

  // ---- ET-split fixes (v5l). All default 0 = legacy (bit-identical). Set from
  // ModelData (model_settings) in Cohort::initialize_internal_pointers().
  //
  // S1 et_fpc_once_fix = 1: transpiration applies the foliar projected cover
  //   (FPC) once. Legacy: updateRadiation() scales SW and PAR by FPC (ground-area
  //   basis), the Penman-Monteith transpiration uses that FPC-scaled radiation
  //   with a single-leaf conductance, and Cohort then multiplies tran by FPC again
  //   ("ALT fix 1", commit 28167031), so the radiation term scales with FPC^2 and
  //   the stomatal light factor m_ppfd sees FPC-diluted PAR. With the switch the
  //   transpiration radiation term and m_ppfd use the radiation per unit covered
  //   area (swinter/FPC, pardown/FPC); Cohort's single tran *= FPC then puts
  //   transpiration on the ground-area basis. Wet-canopy evaporation and the
  //   *_pet diagnostics are unchanged.
  int fpc_once_fix;
  // S2 nv_snow_bury_fix = 1: a non-vascular PFT (moss/lichen) intercepts no snow
  //   when the snowpack (previous day's snow depth, m) is deeper than
  //   nv_snow_bury_depth (m); snow it still stores falls to the snowpack (sdrip).
  //   Default depth 0.0: in dvm-dos-tem the moss layer is the top of the soil
  //   column and snow layers are built on top of it, so any snowpack lies on
  //   the moss/lichen surface.
  int nv_snow_bury_fix;
  double nv_snow_bury_depth;
  // S3 rain_int_lai_fix = 1: rain interception capacity = ircoef [mm per unit
  //   LAI] * LAI (mm, ground area), filled by today's rain on the PFT (rnfl*FPC)
  //   up to the capacity left by yesterday's stored water; stored + new water
  //   evaporates first (wet canopy), transpiration uses the remaining daylight,
  //   and what is left stays in storage for the next day (excess over capacity
  //   drips). Legacy: max_int = 0.21*P + 0.5 mm/day per PFT (no LAI, ircoef unused).
  int rain_int_lai_fix;
  // S4 nv_rain_surface: rain "intercepted" by a non-vascular PFT is moss-surface
  //   water, not canopy water. 1 = same interception/evaporation physics, but the
  //   evaporation is booked as d_v2a.evap_moss (output EVAP_MOSS, included in EET)
  //   instead of canopy evaporation d_v2a.evap; 2 = no non-vascular interception:
  //   the rain reaches the ground surface (soil/moss layer infiltration) and
  //   evaporates through soil evaporation (stored water drips to the ground).
  int nv_rain_surface;
  // per-day context set by Cohort before updateRadiation()/updateWaterBalance()
  int nonvasc;        // cd.d_veg.nonvascular[ipft]
  double snow_depth;  // edall->d_snws.snowthick (m), start of day

private:

  EnvData * ed;
  FirData * fd;
  CohortLookup * chtlu;

  //function
  double getRainInterception(const double & rain, const double & lai);
  double getSnowInterception(const double & snow, const double & lai);
  double getLeafStomaCond(const double & ta, const double &  parin,
                          const double & vpdin, const double& btran,
                          double & m_ppfd, double & m_vpd );
  double getCanopySubl(const double & rac, const double & sinter,
                       const double & lai );

  double getPenMonET(const double & ta, const double& vpd, const double &irad,
                     const double &rv, const double & rh);

};
#endif /*VEGETATION_ENV_H_*/
