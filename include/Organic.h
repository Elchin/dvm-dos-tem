/*! \file
 */
#ifndef ORGANIC_H_
#define ORGANIC_H_

#include "iostream"
using namespace std;

#include "OrganicLayer.h"

class Organic {
public:
  Organic();
  ~Organic();

  int shlwnum;
  double shlwthick;
  double shlwc; //used as a tracker to determine if a shallow Organic
                //  horizon exists.

  int deepnum;
  double deepthick;
  double deepc; //used as a tracker to determine if a deep Organic
                //  horizon exists.

  double shlwdz[MAX_SLW_LAY];
  double deepdz[MAX_DEP_LAY];

  bool shlwchanged;
  double lstshlwdz;

  void ShlwThickScheme(const double & totthickness);
  void DeepThickScheme(const double & totthickness);

  // v5q-humic-stefan: model_settings.humic_max_thick (m; 0 = off/legacy, default, bit-identical).
  // > 0: after the legacy split, if any deep (humic) sublayer is thicker than the cap, the deep
  // horizon is re-split into n = min(MAX_DEP_LAY, max(deepnum, ceil(deepthick/cap))) equal layers
  // (restart/output layout unchanged, MAX_DEP_LAY = 3). Known limit: the cap cannot be met once
  // deepthick > MAX_DEP_LAY * cap (then 3 equal layers of deepthick/3 > cap).
  // Refs: humic layering and the organic C<->thickness relation of the dynamic organic soil module:
  //   Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010, doi:10.1029/2010JG001302.
  //   The cap itself is a vertical-resolution (discretization) choice: no published source.
  double humic_max_thick;
  // Pure helper (unit-tested); returns true if the split was changed.
  static bool applyHumicMaxThick(double cap, double deepthick, double deepdz[], int & deepnum);

  void assignShlwThicknesses(int type[], double dz[],const int & maxnum);
  void assignDeepThicknesses(int type[], double dz[],const int & maxnum);

private:

};
#endif /*ORGANIC_H_*/
