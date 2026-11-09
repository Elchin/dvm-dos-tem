/*! \file SJCurve.cpp  (see SJCurve.h) */
#include "../include/SJCurve.h"
#include <cmath>

namespace sjcurve {

double phiPure(const double t, const double b) {
  if (t > 0.0) { return 1.0; }
  return pow((0.1 - t) / 0.01, b);   // SJ2016 Eq. 4, Tref = 0.1 degC, T* = 0.01 degC
}

double layerPhi(const double t, const double forg, const double fc, const double fsi,
                const double fsa) {
  const double phio = phiPure(t, -1.0);   // Table 1: organic b = -1.0
  double fo = fmin(1.0, fmax(0.0, forg));
  const double c = fmax(0.0, fc), si = fmax(0.0, fsi), sa = fmax(0.0, fsa);
  const double sum = c + si + sa;
  if (!(sum > 0.0)) { return phio; }   // no mineral texture: organic
  const double phim = (c * phiPure(t, -0.3) + si * phiPure(t, -0.5) + sa * phiPure(t, -0.9)) / sum;
  return (1.0 - fo) * phim + fo * phio;   // SJ2016 Eq. 5
}

double orgFraction(const double somc_gm2, const double dz) {
  if (!(dz > 0.0)) { return 1.0; }
  const double mo = 2.0 * fmax(0.0, somc_gm2) / 1000.0;   // kg organic matter m-2 (50 % C)
  return fmin(1.0, mo / (dz * 140.0));                     // SJ2016 Eq. 6, rho_o = 140 kg m-3
}

}
