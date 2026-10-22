#ifndef SNOW_ENV_H_
#define SNOW_ENV_H_

#include "errorcode.h"

#include "CohortData.h"
#include "EnvData.h"
#include "RestartData.h"

#include "Ground.h"

#include "CohortLookup.h"

class Snow_Env {
public:
  Snow_Env();
  ~Snow_Env();

  double wind;     /*! the wind speed class (0, 1)? of a vegetation type*/

  // parameters for snow physics
  snwpar_env snowenvpar;
  // Fix B4 (v5k), OFF by default: model_settings.snow_albedo_fix. Legacy
  // getAlbedoVis() returns a hard-coded 0.9 below -10 C (above snwalbmax, e.g.
  // 0.8, with a jump at -10 C); 1 = snwalbmax below -10 C.
  static int albedo_fix;

  void setGround(Ground* grndp);

  void setCohortLookup(CohortLookup* chtlup);
  void setCohortData(CohortData* cdp);
  void setEnvData(EnvData* edp);

  void initializeParameter();
  void initializeNewSnowState();
  void set_state_from_restartdata(const RestartData & rdata);

  void updateDailyM(const double & tdrv);
  void checkSnowLayersT(Layer* frontl);
  friend class Soil_Env; // Option 4 (v5j) SEB: snow albedo in Soil_Env::solveSEB

private:

  CohortLookup * chtlu;
  CohortData * cd;
  EnvData * ed;

  Ground * ground;

  void updateSnowEd(Layer * frontl);
  double meltSnowLayersAfterT(Layer * frontl);

  void updateDailySurfFlux( Layer* frontl, const double & tdrv);

  double getSublimation(double const & rn, double const & swe,
                        double const & ta, double & radmelt);

  // DRAFT (daily-structural): radiation-driven melt of the top snow layers on
  // days with Tair >= -1 C (mm/day). Previously this was booked as sublimation
  // (lost to the atmosphere and counted in EET); now it is routed to melt water.
  double rad_melt;

  double getAlbedoVis(const double & tem); //get albedo of visible
                                           //  radiation of snow
  double getAlbedoNir(const double & tem); //get albedo of Nir radiation of snow

};

#endif /*SNOW_ENV_H_*/
