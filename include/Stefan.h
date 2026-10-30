/*! \file
 * this class updates the position of freezing and thawing fronts using Stefan Algorithm (ATWoo42004a)
 */
#ifndef STEFAN_H_
#define STEFAN_H_

#include <cmath>
#include <vector>

#include "Ground.h"
#include "errorcode.h"
#include "layerconst.h"

#include "CrankNicholson.h"

class Stefan {
public:
  Stefan();
  ~Stefan();

  int itsumall;

  void setGround(Ground* grndp);
  void initpce();
  void updateFronts(const double & tdrv, const double &timestep);
  // Fix B3 (v5k), OFF by default: model_settings.stefan_fronts_fix. 1 = reset
  // the minimum front distance on every pass of combineExtraFronts() (legacy:
  // set once, which never terminates with more than one surplus pair).
  static int fronts_fix;
  // v5q-humic-stefan: model_settings.stefan_front_resist_fix (0 = legacy, default, bit-identical;
  // 1 = for a THAWING front, the soil layer that contains the front contributes only its thawed
  // part (1 - frozenfrac) * dz / k_unf to the thermal resistance above the front (legacy: the full
  // dz, i.e. the still-frozen part below the front is counted as resistance above it); if the front
  // passes through the whole layer in the same step, the remaining part is added before the next
  // layer, so deeper layers see the full layer as in legacy. The freezing direction is unchanged
  // (possible extension). Integers 0/1 only, anything else FATAL.
  // Refs: modified Stefan front tracking: Woo et al. 2004, doi:10.1029/2004GL019475; its use in
  //   DOS-TEM: Yi et al. 2009, doi:10.1029/2008JG000841. The fix (resistance only above the front)
  //   is a discretization bug fix: no published source.
  static int front_resist_fix;
  // Pure helper (unit-tested): resistance added now (res_now) and, if the front passes the whole
  // layer this step, the extra resistance to add afterwards (res_after).
  static void frontLayerResistance(int fix, double tdrv, bool isSoil, int frozen, int freezing,
                                   double frozenfrac, double dz, double tkres,
                                   double & res_now, double & res_after);
  // Refs (B3): no published source; loop-termination bug in combineExtraFronts (upstream 419a4eb3, fix commit 07f6c8d0); algorithm: Woo et al. 2004, doi:10.1029/2004GL019475.

private:
  bool all_fronts_swept = false; // set in updateFronts(): today's front swept out every existing front


  Ground *ground;
  Layer *botdrvl;

  void meltingSnowLayer(double const & tkfront, double & dse,
                        double & sumresabv, const double & tdrv,
                        Layer* currl, const double & timestep);

  //20180820 This function requires the timestep in order to
  // modify the energy needed to thaw soil layers by adding the
  // energy needed to raise the layer to 0 degrees.
  void processNewFrontSoilLayerDown(const int &frozenstate,
                                    double const & sumrescum,
                                    double const & tkfront, double & dse,
                                    double & newfntdz, Layer* currl,
                                    const double & timestep);

  void frontsDequeDown(const double & newfntz, const int & newfnttype);

  //20180820 Unlike the downward process, the upward process does
  // not need the timestep because layer temperature is already
  // allowed for.
  void processNewFrontSoilLayerUp(const int &frozenstate,
                                  double const & sumrescum,
                                  double const & tkfront ,
                                  double & dse, double & newfntdz,
                                  Layer* currl);

  void frontsDequeUp(const double & newfntz, const int & newfnttype);

  double prepareBottomDriving();

  //get the degree seconds needed to fully freeze/thaw  one or part of one layer
  double getDegSecNeeded(const double & dz, const double & volwat,
                         const double & tk , const double & sumresabv);
  //calculate partial depth based on extra degree seconds
  double getPartialDepth(const double & volwat,const double & tk,
                         const double & sumresabv, const double & dse);

  void combineExtraFronts();
  void rebuildFrontsFromLayerStates();
  void updateLayerFrozenState(Layer* toplayer, const int freezing1);
  void updateWaterAfterFront(Layer* toplayer);


};

#endif /*STEFAN_H_*/
