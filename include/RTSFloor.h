#ifndef RTSFLOOR_H_
#define RTSFLOOR_H_
/** RTS-light v0.2: post-slump floor carbon treatment (three independent switches, all ON by default
 *  in v0.2; all three false -> bitwise v0.1). m = months since RTS initiation (<0: not initiated, nothing is applied).
 *  (1) rh_cap:   labile-limited mineralization of the exposed mineral (deep) C. Per mineral layer and month
 *                RH_l' = min(RH_l, k(m)/12 * C_l),  k(m) = cap_inf + (cap0 - cap_inf) * exp(-m / (12 * cap_tau))
 *                (1/yr): initial labile-limited rate cap0 (Knoblauch et al. 2021: ~0.8% of thawed C in yr 1)
 *                moving to the stabilised-floor rate cap_inf over ~cap_tau yr (labile depletion, floor stabilisation)
 *  (2) export:   lateral export of exposed mineral C/N by floor mass wasting / mudflows, monthly fraction
 *                e(m) = lat0/12 * exp(-m / (12 * lat_tau)) of every mineral-layer SOM pool and orgN -> TKLATC/TKLATN
 *  (3) anox:     oxygen limitation of RH on the saturated floor: RH_l *= 1 - (1 - anox_f) * sat_l,
 *                sat_l = saturated fraction of layer l below the water table, raised by the post-slump ponding
 *                of RTSMethane (sat -> sat + (1-sat) pond).
 *  (4) rmin:     minimum floor respiration (column): RHSOM(month) >= r(m) * RH_pre(calendar month), with
 *                r(m) = 1 - (1 - rmin0) * exp(-m / (12 * rmin_tau)) and RH_pre the mean column RHSOM of the same
 *                calendar month over the last rmin_years years before initiation (running, per calendar month,
 *                so the seasonal cycle is kept). The deficit is taken from the existing soil SOM pools (layer
 *                weight = layer SOM C x the layer's moisture/temperature scalers, pools in proportion to their C,
 *                at most rmin_maxfrac of a layer's C per month): no C is created, closure is unaffected.
 *                Floor/undisturbed chamber ratios 0.28-0.55 for dry floors, stabilised slumps ~1 (Jensen et al.
 *                2014, doi:10.1088/1748-9326/9/2/025001; Cassidy et al. 2016; russia-inputs/ch4_compare/
 *                chamber_data_dryupland.md). Keys: thermokarst_floor_rmin, rts_floor_rmin0, rts_floor_rmin_tau.
 *  Defaults/sources: docs_src/rts_ch4.md. */
namespace Json { class Value; }

struct RTSFloorParams {
  bool rh_cap = true;    double cap0 = 0.045; double cap_inf = 0.002; double cap_tau = 1.5; // 1/yr, 1/yr, yr (CALIBRATED, docs)
  bool exp_on = true;    double lat0 = 0.02;  double lat_tau = 5.0;   // 1/yr, yr
  bool anox = true;      double anox_f = 0.3;                        // anaerobic/aerobic RH ratio
  bool rmin_on = true;   double rmin0 = 0.3;  double rmin_tau = 10.0;  // -, yr (see (4))
  int rmin_years = 5;                                                 // pre-slump climatology window, yr
  bool any() const { return rh_cap || exp_on || anox || rmin_on; }
  void load(const Json::Value & ms);
};

class RTSFloor {
public:
  static double capRate(const RTSFloorParams & p, const int m);   // k(m), 1/yr
  // returns the scaled RH of one layer (rh >= 0, c = layer SOM C, gC/m2)
  static double capLayerRH(const RTSFloorParams & p, const int m, const double rh, const double c);
  static double exportFraction(const RTSFloorParams & p, const int m);  // monthly fraction
  static double anoxFactor(const RTSFloorParams & p, const double sat);
  static double rminFactor(const RTSFloorParams & p, const int m);   // r(m), 1 before initiation
  static constexpr double rmin_maxfrac = 0.05;  // max fraction of a layer's SOM C added per month
  /** Raises the column RH (per-layer pool RH arrays, n layers) to at least target (g/month) by adding RH
   *  from the pools c*[] in proportion to w_l * C_l (C_l = layer SOM C), each layer limited to
   *  rmin_maxfrac * C_l. Returns the RH added (<= target - current RH; 0 if already above). */
  static double applyMinimum(const double target, const int n, const double * w,
                             const double * crawc, const double * csoma, const double * csompr, const double * csomcr,
                             double * rrawc, double * rsoma, double * rsompr, double * rsomcr);
};
#endif
