
#ifndef GROUND_H_
#define GROUND_H_

#include <sstream>
#include <cmath>
#include <iostream>
#include <deque>
using namespace std;

#include "DoubleLinkedList.h"
#include "Snow.h"
#include "Moss.h"
#include "Organic.h"
#include "MineralInfo.h"
#include "SoilParent.h"

#include "Layer.h"
#include "MossLayer.h"
#include "SnowLayer.h"
#include "MineralLayer.h"
#include "OrganicLayer.h"
#include "ParentLayer.h"

#include "CohortLookup.h"

#include "BgcData.h"
#include "CohortData.h"
#include "FireData.h"
#include "RestartData.h"

#include "parameters.h"
#include "TEMUtilityFunctions.h"

class Ground: public DoubleLinkedList {

public :
  Ground();
  Ground(MineralInfo mi);

  // Option 4 (v5j) SEB: today's SEB result, set by
  // Soil_Env::updateDailyGroundT and used by Snow_Env::updateDailyM.
  bool seb_on = false;     // SEB set today's top boundary
  double seb_ts = 0.0;     // skin temperature (C)
  double seb_melt = 0.0;   // snow melt from the SEB energy surplus (mm/day)
  double seb_sublim = 0.0; // snow sublimation (mm/day)

  ~Ground();

  std::string layer_report_string();
  std::string layer_report_string(const std::string& colunm_groups);

  // A ground (snow-soil<moss-peat-mineral>-soilparent column has the following
  //   5 types of horizons, each of which has a number of layers defined below
  Snow snow;
  Moss moss;
  Organic organic;
  MineralInfo mineralinfo;
  SoilParent soilparent; //soil parent materials, including rock

  // dimension parameters
  snwpar_dim snowdimpar;
  soipar_dim soildimpar;

  // layer boundary for each type of horizon

  Layer* fstsoill; // first layer of soil column: snow/rock horizons excluded
  Layer* lstsoill; // last layer of soil column: snow/rock horizons excluded

  Layer* fstmossl; // first moss layer
  Layer* lstmossl; // last moss layer
  Layer* fstshlwl; // first fibric os layer
  Layer* lstshlwl; // last fibric os layer
  Layer* fstdeepl; // first amorphous os layer
  Layer* lstdeepl; // last amorphous os layer
  Layer* fstminel; // first mineral layer
  Layer* lstminel; // last mineral layer

  // freezing/thawing fronts
  double frntz[MAX_NUM_FNT];
  int frnttype[MAX_NUM_FNT];
  deque<double> frontsz; //fronts depth in order from top to bottom:
                         //  distance from ground soil surface
  deque<int> frontstype; //fronts type in order as above:
                         //  1 = freezing front - front with upper
                         //    frozen/lower unfrozen,
                         //  -1 = thawing front - front with upper
                         //    unfrozen/lower frozen.
  // SO freezing/thawing fronts are alternatively in order,
  // i.e., a freezing front must be followed by a thawing front, or vice versa

  Layer* fstfntl; /*! first snow/soil layer containing phase change front */
  Layer* lstfntl; /*! last snow/soil layer containing phase change front */

  // this is the last layer in which Richards Equation applied for
  Layer* drainl; /*! layer in which the subsurface drainage occurs */
  // Water-leak fix (v5f): soil water that the layer
  // bookkeeping has to remove (freeze expansion beyond pore space, organic layer
  // shrinkage). Booked as subsurface drainage on the next Soil_Env::updateDailySM
  // instead of being deleted (mm).
  double expelled_water_mm = 0.0;

  double draindepth; //the subsurface drainage depth (m), which much
                     //  be in somewhere of 'drainl'
  //this depth is one of following situation: watertable, barrier
  //  (e.g., frozen), or soil bottom)

  int ststate; // thermal state of whole column,
  // 0: partially frozen, 1: totally frozen,, -1: totally unfrozen

  void setBgcData(BgcData *bdp);
  void setCohortLookup(CohortLookup *chtlu);

  void initParameter();
  void initDimension();
  void initLayerStructure(snwstate_dim * snowdim, soistate_dim * soildim);
  void set_state_from_restartdata(snwstate_dim *snowdim, soistate_dim *soildim,
                             const RestartData & rdata);
  void resortGroundLayers();

  // snow layers
  bool constructSnowLayers(const double & dsmass, const double & tdrv);
  bool divideSnowLayers();
  bool combineSnowLayers();

  void updateSnowLayerPropertiesDaily();
  void updateSnowHorizon();
  void checkSnowLayer();

  // soil layers
  void redivideSoilLayers();

  // v5q-layern-fix: model_settings.layer_n_remap_fix (0 = legacy, default).
  //  1 = (a) splitOneSoilLayer, organic layers: the new lower layer gets orgn/avln in the same
  //          fraction as its C (legacy: all N of the split layer goes to the lower part, the
  //          upper part keeps C with N = 0); (b) new deep layer formed from the last fibric
  //          layer: N moves in the fraction of the C moved (legacy: all N moves, fibric N = 0);
  //      (c) combineTwoSoilLayersL2U: avln is summed (legacy typo 'usl->avln =+lsl->avln'
  //          assigns the lower layer's value, so the upper layer's avln is lost).
  //      Soil_Bgc::assignNitrogen{Bd2Layer,Layer2Bd}Monthly carry layer N with the C pools
  //      around redivideSoilLayers (Cohort::updateMonthly_DIMgrd and the post-burn path).
  // Refs: layer structure and the organic-layer C<->thickness relation of the dynamic organic
  //   soil (DOS) module: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010,
  //   doi:10.1029/2010JG001302. The fix itself is a mass/stoichiometry bookkeeping
  //   (conservation) bug fix: no published source for the remap rule (N follows C).
  int layer_n_remap_fix = 0;
  // fraction of N that goes with the C part c_part of a layer holding c_whole (clamped to
  // [0,1]); 'fallback' (thickness fraction) if c_whole <= 0. Unit-tested.
  static double layerNSplitFrac(const double c_part, const double c_whole, const double fallback);

  // thawing/freezing fronts
  void setFstLstFrontLayers();
  void updateWholeFrozenStatus(); //update if whole soil frozen(1)/
                                  //  unfrozen(-1)/partially-frozen(0)

  // update water drainage layer and depth
  void setDrainL();

  // soil burning caused soil structure change
  double adjustSoilAfterburn();

  // organic layer thickness dynamics with C content
  void updateOslThickness5Carbon(Layer* fstsoill);

  //
  void retrieveSnowDimension(snwstate_dim * snowdim);
  void retrieveSoilDimension(soistate_dim * soildim);//This is required if
                                                     //  anything changed in
                                                     //  the dimension

  //
  void checkWaterValidity();

  double getCarbonForDepthRange(double topz, double bottomz);

  double getVWCAtDepth(double vwcdepth);
  double getTempAtDepth(double temperaturez);

private :

  bool rocklayercreated;

  BgcData * bd;
  CohortLookup * chtlu;

  void initRockLayers();
  void initSnowSoilLayers();

  void setFstLstSoilLayer();
  void setFstLstMossLayers();
  void setFstLstShlwLayers();
  void setFstLstDeepLayers();
  void setFstLstMineLayers();
  void updateLayerIndex();
  void updateLayerZ();

  void updateSoilHorizons();

  void redivideMossLayers(const int &mosstype);
  void redivideShlwLayers();
  void redivideDeepLayers();

  void splitOneSoilLayer(SoilLayer*usl, SoilLayer* lsl,
                         const double & updeptop, const double &lsldz);
  void combineTwoSoilLayersU2L(SoilLayer* usl, SoilLayer* lsl);
  void combineTwoSoilLayersL2U(SoilLayer* lsl, SoilLayer* usl);

  void adjustFrontsAfterThickchange(const double & depth,
                                    const double & thickchange);
  void getLayerFrozenstatusByFronts(Layer * soill);

  double thicknessFromCarbon(const double carbon, const double coefA, const double coefB);
  double carbonFromThickness(const double thickness, const double coefA, const double coefB);

  void getOslCarbon5Thickness(SoilLayer* sl, const double &plctop,
                              const double &plcbot);
  void getOslThickness5Carbon(SoilLayer* sl, const double &plztop,
                              const double &plzbot);

  void checkFrontsValidity();

  void cleanSnowSoilLayers();
  void cleanAllLayers();
  void cleanRockLayers();

};

#endif /*GROUND_H_*/
