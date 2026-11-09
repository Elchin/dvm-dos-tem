/*! \file SJCurve.h
 * Schaefer & Jafarov (2016) unfrozen-water curve, the ONE shared implementation used by
 * soil respiration (Soil_Bgc::sjPhiPure / sjLayerPhi / sjOrgFraction, rh_frozen_mode) and by
 * the enthalpy thermal solver (EnthalpySolver, thermal_solver = 1, thermal_sfcc = 1).
 * The arithmetic is exactly the former Soil_Bgc code (moved, not changed), so the legacy
 * respiration path stays bit-identical.
 * Refs: Schaefer & Jafarov 2016, doi:10.5194/bg-13-1991-2016 (Eqs. 4-6, Table 1).
 */
#ifndef SJCURVE_H_
#define SJCURVE_H_

namespace sjcurve {
// 1 for t > 0, else ((0.1 - t)/0.01)^b  (SJ2016 Eq. 4, Tref = 0.1 degC, T* = 0.01 degC)
double phiPure(const double t, const double b);
// texture/organic blend (SJ2016 Eq. 5): clay/silt/sand b = -0.3/-0.5/-0.9, organic b = -1.0
double layerPhi(const double t, const double forg, const double fc, const double fsi,
                const double fsa);
// organic fraction (SJ2016 Eq. 6, rho_o = 140 kg m-3, 50 % C), capped at 1
double orgFraction(const double somc_gm2, const double dz);
}

#endif
