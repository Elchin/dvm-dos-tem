
/*!\file
 * Implementation of Richard's law for soil water dynamics*/

#ifndef RICHARDS_H_
#define RICHARDS_H_

#include <cmath>
#include <limits>

#include "CohortData.h"
#include "EnvData.h"

#include "Layer.h"
#include "SoilLayer.h"
#include "CrankNicholson.h"

// Richards work arrays are indexed by Layer::indl, which counts snow layers above the
// soil column (up to MAX_SNW_LAY), and prepareSoilColumn() reads index drainind+1.
const int MAX_RICH_IND = MAX_SNW_LAY + MAX_SOI_LAY + 2;

class Richards {
public :
  Richards();
  ~Richards();

  int num_al;   // active layers' layer numbers
  double qdrain; // mm/day
  double qinfil; // mm/day
  // Fix B2 (v5k): true when the last update() returned early (drainage depth at
  // or above the top of its domain: no infiltration, evaporation or uptake applied).
  bool early_return = false;
  double qevap; // mm/day
  double e_ice; //ice impedance factor (CLM5)

  // model_settings.richards_dt_fix, OFF by default (0 = legacy, bit-identical).
  // Bug R1: when a trial substep is rejected (max_tridiag_error > toler_upper) and dtsub is
  // halved, legacy re-solves with the dt_dz = dtsub/dz of the rejected (longer) substep, but
  // only advances dtdone by the halved dtsub. The accepted deltathetaliq therefore carries the
  // fluxes of a longer interval, so the day's net flux is applied for more than SEC_IN_DAY
  // (water is created or destroyed; 3.3 mm/d at 20 mm/d infiltration in the unit test).
  // 1 = recompute dt_dz from the halved dtsub before every retry, so the tridiagonal system
  // (computeLHS/computeRHS), the error estimate (fluxNet0) and dtdone use the same dtsub.
  // Refs: CLM5 tech note (Lawrence et al. 2018) sect. 7.3.3 / CTSM SoilWaterMovementMod
  // (adaptive time stepping recomputes dt-dependent coefficients each trial).
  int dt_fix = 0;
  // Diagnostics only (never affect results): cumulative counts over the object's life.
  long n_calls = 0;        // update() calls that entered the substep loop
  long n_rejections = 0;   // rejected trial substeps (dtsub halvings)
  long n_reject_calls = 0; // update() calls with >= 1 rejection

  // regional_wetland_deepc (wetland phase 3): run-on / groundwater; set per cohort in
  // Soil_Env::initializeParameter (wetland off: max_ponding 100 mm = the legacy literal,
  // inflow_factor 0, both flags false).
  double max_ponding = 100.0;
  double inflow_factor = 0.0;
  bool groundwater_status = false;
  bool runon_status = false;

  //Array to hold lateral drainage values for output to file.
  double layer_drain[MAX_RICH_IND];
  double percolation[MAX_RICH_IND];
  double layer_input[MAX_RICH_IND];
  double qgw_mm = 0.0;   // wetland (phase 3) diagnostic: groundwater recharge added in the last update (mm)   // groundwater input by Richards index (mm/day)

  void update(Layer *fstsoill, Layer* bdrainl, const double & bdraindepth,
              const double & fbaseflow, const double & watertab,
              double root_water_up[MAX_SOI_LAY],
              const double & evap, const double & infil,
              const double & cell_slope);

  void setCohortData(CohortData* cdp);
  void setEnvData(EnvData* edp);

  CrankNicholson cn;

private:

  EnvData ed;

  void prepareSoilColumn(Layer *currsoill, int drainind);
  void clearRichardsArrays();
  void computeHydraulicProperties(Layer *fstsoill, int drainind);
  void computeMoistureFluxesAndDerivs(Layer *fstsoill, int topind, int drainind);
  void computeLHS(Layer *fstsoill, int topind, int drainind);
  void computeRHS(Layer *fstsoill, int topind, int drainind);

  Layer * drainl; //bottom active layer
  double z_watertab; //water table depth in mm

  //Note that the Richards code uses arrays whose indices follow the 1-based indexing
  //of soil layers. That is, element zero of these arrays will not be used. An extra
  //element is added to the length of the arrays (MAX_SOI_LAY+1) to allow for this empty
  //placeholder. Element 1 will generally represent a moss layer, which are not currently
  //used in percolation and lateral drainage. Therefore, the first element of the arrays that
  //is filled will usually be element 2, corresponding to the index of the first fibric organic layer.

  //Hydraulic properties and water fluxes for calculation of percolation and lateral drainage
  double qtrans[MAX_RICH_IND]; //root water uptake (incoming) (mm/s)
  double fluxNet0[MAX_RICH_IND]; //column water flux based on solver solution; from CLM fortran code
  double fluxNet1[MAX_RICH_IND]; //estimated column water flux based on water qin/qout/trans/evap/infil; from CLM fortran code
  double qin[MAX_RICH_IND]; //water flux into soil layer (mm/s); see CLM 5 eq 7.74
  double qout[MAX_RICH_IND]; //water flux out of soil layer (mm/s); see CLM 5 eq 7.75
  double s1[MAX_RICH_IND]; //liquid saturation at the layer interface (unitless); from CLM fortran code
  double s2[MAX_RICH_IND]; //liquid saturation at the layer node (unitless); from CLM fortran code
  double imped[MAX_RICH_IND]; // ice impedance (unitless); see CLM 5 eq 7.48
  double hk[MAX_RICH_IND]; //hydraulic conductivity (mm/s); see CLM 5 eq 7.47
  double dhkdw[MAX_RICH_IND]; //d(hk)/d(vol_liq); see CLM 5 eqs 7.83-7.84
  double dsmpdw[MAX_RICH_IND]; //d(smp)/d(vol_liq); see CLM 5 eqs 7.80-7.82
  double smp[MAX_RICH_IND]; //soil matric potential (mm); see CLM 5 eq 7.53
  double dqidw0[MAX_RICH_IND]; //d(qin)/d(vol_liq(i-1)); see CLM 5 eq 7.76
  double dqidw1[MAX_RICH_IND]; //d(qin)/d(vol_liq(i)); see CLM 5 eq 7.77
  double dqodw1[MAX_RICH_IND]; //d(qout)/d(vol_liq(i)); see CLM 5 eq 7.78
  double dqodw2[MAX_RICH_IND]; //d(qout)/d(vol_liq(i+1)); see CLM 5 eq 7.79

  //Incoming layer property values copied to arrays for convenience.
  double Bsw[MAX_RICH_IND]; //bsw hornberger constant (unitless) (by horizon type)
  double ksat[MAX_RICH_IND]; //saturated hydraulic conductivity (mm/s) (by horizon type)
  double psisat[MAX_RICH_IND]; //saturated soil matric potential (mm/s) (by horizon type)

  //Physical values calculated and/or collected to arrays for convenience
  double vol_liq[MAX_RICH_IND];//volumetric liquid water (mm3 liquid water/mm3 soil)
  double vol_ice[MAX_RICH_IND];//volumetric ice (mm3 ice/mm3 soil)
  double vol_h2o[MAX_RICH_IND]; //total volumetric water (mm3 liq + mm3 ice)/mm3 soil)
  double vol_poro[MAX_RICH_IND]; //total layer porosity
  double liq_poro[MAX_RICH_IND]; //porosity available to be filled by liquid (vol_poro - vol_ice)
  double soi_liq[MAX_RICH_IND]; //layer liquid water (kg/m2 or mm)
  double ice_frac[MAX_RICH_IND]; //fraction of porosity filled by ice (vol_ice/vol_poro)
  double effminliq[MAX_RICH_IND]; //minimum layer liquid adjusted for frozen layer fraction (kg m-2 or mm)
  double effmaxliq[MAX_RICH_IND]; //maximum layer liquid adjusted for frozen layer fraction (kg m-2 or mm)
  double z_h[MAX_RICH_IND]; //depth of layer bottom interface in mm; see CLM 5 Figure 7.2
  double dzmm[MAX_RICH_IND]; // layer thickness in mm
  double nodemm[MAX_RICH_IND]; //depth of center of layer in mm; CLM 5 Figure 7.2 'z(i)'

  //Arrays related to iteration and the tridiagonal matrix solution for richards percolation
  double deltathetaliq[MAX_RICH_IND]; //Results from tridiagonal calculation
  double amx[MAX_RICH_IND]; //left off-diagonal of tridiagonal matrix for richards solution
  double bmx[MAX_RICH_IND]; //diagonal of tridiagonal matrix input for richards solution
  double cmx[MAX_RICH_IND]; //right off-diagonal of tridiagonal matrix input for richards solution
  double rmx[MAX_RICH_IND]; //forcing term of input for tridiagonal richards solution
  double dt_dz[MAX_RICH_IND]; //temporary time/thickness variable used in iteration
  double tridiag_error[MAX_RICH_IND]; //layer-wise error in tridiagonal solver solution per iteration; from CLM 5 eq 7.98

  //Values related to iteration and the tridiagonal matrix solution for richards percolation
  double max_tridiag_error; //maximum layer-wise error in tridiagonal solver solution per iteration
  double delta_t; //(seconds) total time to be completed during iteration (one day)
  double dtmin;  // minimum timestep length (seconds). Shorter values may be more accurate but take longer for total iteration.
  double toler_upper;  // Tolerance to halve length of substep
  double toler_lower;  // Tolerance to double length of substep

};

#endif /*RICHARDS_H_*/
