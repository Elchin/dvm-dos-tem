#ifndef RTSFLOOR_H_
#define RTSFLOOR_H_
/** RTS-light v0.2: post-slump floor carbon treatment (three independent switches, all off by default
 *  -> bitwise v0.1). m = months since RTS initiation (<0: not initiated, nothing is applied).
 *  (1) rh_cap:   labile-limited mineralization of the exposed mineral (deep) C. Per mineral layer and month
 *                RH_l' = min(RH_l, k(m)/12 * C_l),  k(m) = cap_inf + (cap0 - cap_inf) * exp(-m / (12 * cap_tau))
 *                (1/yr): initial labile-limited rate cap0 (Knoblauch et al. 2021: ~0.8% of thawed C in yr 1)
 *                moving to the stabilised-floor rate cap_inf over ~cap_tau yr (labile depletion, floor stabilisation)
 *  (2) export:   lateral export of exposed mineral C/N by floor mass wasting / mudflows, monthly fraction
 *                e(m) = lat0/12 * exp(-m / (12 * lat_tau)) of every mineral-layer SOM pool and orgN -> TKLATC/TKLATN
 *  (3) anox:     oxygen limitation of RH on the saturated floor: RH_l *= 1 - (1 - anox_f) * sat_l,
 *                sat_l = saturated fraction of layer l below the water table, raised by the post-slump ponding
 *                of RTSMethane (sat -> sat + (1-sat) pond).
 *  Defaults/sources: docs_src/rts_ch4.md. */
namespace Json { class Value; }

struct RTSFloorParams {
  bool rh_cap = false;   double cap0 = 0.045; double cap_inf = 0.002; double cap_tau = 1.5; // 1/yr, 1/yr, yr (CALIBRATED, docs)
  bool exp_on = false;   double lat0 = 0.02;  double lat_tau = 5.0;   // 1/yr, yr
  bool anox = false;     double anox_f = 0.3;                        // anaerobic/aerobic RH ratio
  bool any() const { return rh_cap || exp_on || anox; }
  void load(const Json::Value & ms);
};

class RTSFloor {
public:
  static double capRate(const RTSFloorParams & p, const int m);   // k(m), 1/yr
  // returns the scaled RH of one layer (rh >= 0, c = layer SOM C, gC/m2)
  static double capLayerRH(const RTSFloorParams & p, const int m, const double rh, const double c);
  static double exportFraction(const RTSFloorParams & p, const int m);  // monthly fraction
  static double anoxFactor(const RTSFloorParams & p, const double sat);
};
#endif
