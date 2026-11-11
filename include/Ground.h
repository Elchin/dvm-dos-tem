
#ifndef GROUND_H_
#define GROUND_H_

#include <sstream>
#include <cmath>
#include <iostream>
#include "SnowClass.h"
#include <deque>
#include <vector>
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
  // model_settings.snow_scheme (default 0 = legacy dvmdostem snow, bit-identical). 1 = Sturm et al.
  // (1995) snow classes as in SiBCASA (Schaefer et al. 2009, doi:10.1029/2008JF001125): new-snow
  // density from T and wind (eqs 4-7), Anderson (1976)/CLM4 compaction instead of the 4-day relaxation
  // to snwdenmax, two-layer density stratification (eqs 1-3) superposed on the snow layers (SWE and
  // sum(ice*T) conserved), fixed depth-hoar k in the bottom layer. See SnowClass.h.
  int snow_scheme = 0;
  SnowClassPar snowcls;
  double snowcls_tair = 0.;   // today's air temperature (C), set by Cohort before Snow_Env::updateDailyM
  double snowcls_wind = 0.;   // today's wind speed (m/s), idem
  double newSnowDensity() const;          // snowdimpar.newden (scheme 0) or rho_new (scheme 1)
  void compactSnowLayersClass();          // scheme 1: Anderson/CLM compaction of each snow layer (1 day)
  void stratifySnowLayers(bool redistribute); // scheme 1: two-layer profile; false = only set fdh/kdh
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
  // snow_empty_cleanup (model_settings, default 0 = legacy): remove snow layers that hold no water
  // (|ice| < SNOW_EMPTY_MASS_EPS and |liq| < SNOW_EMPTY_MASS_EPS kg m-2, the threshold below which
  // checkWaterValidity() zeroes ice/liq anyway) and have no thickness (|dz| <= SNOW_EMPTY_DZ_EPS m).
  // Such layers are left when the thermal step melts a lower snow layer completely while the
  // layer above keeps ice: constructSnowLayers() only drops empty layers from the top, and
  // combineSnowLayers() tests the pre-compaction dz (dz -> 0 only in updateThick()), so the layer
  // reaches the next day's thermal step with zero heat capacity. Any residual mass (<= eps) goes
  // to the snow layer above (else the one below, else extramass): mass is conserved exactly; it
  // carries no heat capacity worth booking (< 2e-9 kg m-2). Returns true if a layer was removed
  // (caller must resortGroundLayers()). Counters are cumulative over the run (diagnostics only).
  static constexpr double SNOW_EMPTY_MASS_EPS = 1.0e-9;   // kg m-2 (as checkWaterValidity)
  static constexpr double SNOW_EMPTY_DZ_EPS   = 1.0e-9;   // m
  bool removeEmptySnowLayers();
  int snow_empty_cleanup = 0;
  long snow_empty_removed = 0;      // layers removed
  long snow_empty_days = 0;         // days with >= 1 removal
  double snow_empty_mass = 0.0;     // kg m-2 moved to a neighbour / extramass

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

  // v5q-layer-bookkeep: model_settings.layer_bookkeep_fix (0 = legacy, default, bit-identical;
  //  1 = fix; 1 requires layer_n_remap_fix 1, checked in ModelData). Conservation bug fixes of
  //  the dynamic soil-layer bookkeeping, no physics/parameter change:
  //  (1) redivideShlwLayers, no fibric layer left: the new fibric layer's C is moved out of the
  //      layer below it (all four C pools, orgn and avln in the same fraction; at most all of
  //      that layer's C; no new layer if it holds no C). Legacy: rawc = aboveground litterfall
  //      C of the month, not subtracted from any pool (C created; orgn = avln = 0).
  //  (2) redivideDeepLayers, deepnum == 0: C and N of the removed humic layer(s) are added to
  //      the next soil layer below (the layer above if there is none). Legacy: lost.
  //  (3) BgcData::prvltrfcnque (12-month litterfall C:N per layer index, sets kd in
  //      Soil_Bgc::updateKdyrly4all): carried to the new layers after the January re-division
  //      in Cohort::updateMonthly_DIMgrd (same horizon type, relative depth in the horizon).
  //      Legacy: stays at the old layer index after layers are split/combined.
  //  (4) splitOneSoilLayer, organic: the new lower layer gets at most the C of the layer being
  //      split, so the upper part keeps C >= 0. Legacy: no cap (upper C can go negative when
  //      the layer thickness is not on the C-thickness curve).
  // Refs: layer structure and the organic-layer C<->thickness relation of the dynamic organic
  //   soil (DOS) module: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010,
  //   doi:10.1029/2010JG001302. The fixes themselves are mass-conservation bookkeeping
  //   bug fixes: no published source.
  int layer_bookkeep_fix = 0;
  // Option C (thermal_solver = 1): enthalpy-conserving soil layer split/merge (split: both parts
  // keep the parent's T and specific liq/ice; merge: heat-capacity-weighted T with summed
  // liq/ice). 0 = legacy (default; gradient/thickness-weighted T, front-based water split).
  int enthalpy_remap = 0;
  static double enthHeatCap(const Layer* l); // Cs + liq c_w + ice c_i (J m-2 K-1), soil/rock
  // (1): fraction of the layer below moved into a new fibric layer that needs c_want; c_have =
  //      C of the layer below. min(c_want, c_have)/c_have, 0 if c_have <= 0 or c_want <= 0.
  static double shlwFromBelowFrac(const double c_want, const double c_have);
  // (2): add all C/N pools of 'from' to 'to'
  static void addLayerCN(Layer* to, const Layer* from);
  // (3): litterfall C:N history snapshot before the re-division, remap after it
  struct LtrCNSnap { int cls; double rtop, rbot, amid; std::deque<double> q; };
  std::vector<LtrCNSnap> ltrcn_snap;
  static int layerHorizonClass(const Layer* l);   // 0 moss, 1 fibric, 2 humic, 3 mineral, 4 other
  static int ltrcnSourceIndex(const std::vector<LtrCNSnap>& s, const int cls, const double rmid,
                              const double amid);
  void snapshotLitterCNQueues();
  void remapLitterCNQueues();

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

  std::array<double, MAX_SOI_LAY> get_ch4_bylayer();
  void set_ch4_bylayer(const double* newch4);

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
