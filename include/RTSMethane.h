#ifndef RTSMETHANE_H_
#define RTSMETHANE_H_
/** RTS-light v0.2 CH4, scheme (b): monthly port of the production / oxidation /
 * plant-transport terms of Soil_Bgc::CH4Flux (upstream origin/combined_ch4_hydro,
 * Fan et al. 2013 lineage) with ONE CH4 pool per soil layer (restart state) and
 * NO diffusion solver, NO hydrology/interception changes (vegetation physics and
 * soil water/heat unchanged).  Per thawed layer l, month of `days` days:
 *  prod_l = kfrac * sat_l * prodq10^((T-tref)/10) * sum_k kd_k C_k       (k = rawc, soma, sompr, somcr;
 *           kd_k = the model's own monthly decomposition constants; each pool gives <= avail_frac C_k)
 *           (removed from the SOM pools in proportion to kd_k C_k: an extra,
 *            anaerobic decomposition pathway, as upstream)
 *  P_l   += prod_l
 *  ox_l   = min(P_l, (1-sat_l) * Vmax dz days * P_l/(km dz + P_l) * oxidq10^((T-tref)/10))   (Michaelis-Menten)
 *  pt_l   = P_l * (1 - exp(-plant_k * fLAI * min(1, froot_l/froot_ref) * days)); rhizo_ox of it oxidized
 *  esc_l  = P_l * (1 - exp(-days/tau_l)),  tau_l = tau_sat (sat_l>=0.5) or tau_unsat
 *           escaping CH4 from below the WT passes the oxic zone: oxidized share f_ox(WTD_eff)*sat_l
 * Frozen layers (T <= 0): no production, transport or oxidation -> CH4 is stored
 * and released after thaw (shoulder-season pulse).  At RTS initiation the whole
 * column CH4 pool degasses (counted as emission in the slump month).
 * C balance:  dC_soil+veg = ... - prod ;  dP = prod - oxid - emis  ->
 *  dC_total(incl. P) = GPP - RA - RH_CO2 - RHDWD - CH4emis - lateral - RTS respired,
 *  RH_CO2 = RHSOM + oxid.
 * RTS additions (ponding, lateral CH4 split) as in scheme (a).
 * Units: gC/m2/month (pools gC/m2), m, degC, days.  Defaults/sources: docs_src/rts_ch4.md
 */
#include <vector>
namespace Json { class Value; }

struct RTSMethaneParams {
  // scheme (b) production (calibrated scale factor kfrac, see docs)
  double kfrac = 0.29;      // CH4 production per unit potential decomposition; CALIBRATED: 5 gC/m2/yr wet sedge tundra (docs)
  double prodq10 = 3.0;     // Q10 of methanogenesis
  double tref = 10.0;       // degC
  double avail_frac = 0.05; // production limited to this fraction of the layer SOC per month
  // oxidation (Michaelis-Menten, upstream values converted: Vmax 20 umol/L/h, km 5 umol/L)
  double oxid_vmax = 5.76;  // gC/m3/day
  double oxid_km = 0.06;    // gC/m3
  double oxidq10 = 1.9;
  // plant transport
  double plant_k = 0.1;     // 1/day at LAI >= lai_ref and froot >= froot_ref
  double lai_ref = 1.0;
  double froot_ref = 0.1;
  double rhizo_ox = 0.5;    // share of plant-transported CH4 oxidized in the rhizosphere (upstream 50%)
  // escape (diffusion + ebullition, reduced form)
  double tau_sat = 20.0;    // days, layers below the water table
  double tau_unsat = 1.0;   // days, layers above the water table
  // transit oxidation vs water-table depth (shared with scheme a)
  double oxid_min = 0.2, oxid_max = 0.95, oxid_zscale = 0.1;
  // RTS additions
  double pond_sat0 = 0.5, pond_tau = 10.0, lat_ch4_frac = 0.02;
  // background anoxic/inundated fraction of the cell (0 = model water table only;
  // 1 = wet tundra with WT at the surface; e.g. from a wetland map such as BAWLD)
  double bg_sat = 0.0;
  void load(const Json::Value & ms);
};

struct RTSMethaneLayer {
  double ztop = 0.0, dz = 0.1, tem = 0.0;
  double c[4] = {0, 0, 0, 0};   // rawc, soma, sompr, somcr (gC/m2)
  double kd[4] = {0, 0, 0, 0};  // monthly decomposition constants
  double froot = 0.0;           // fine-root fraction in layer
  double pool = 0.0;            // CH4 pool (gC/m2), updated in place
  double dc[4] = {0, 0, 0, 0};  // out: C removed from each SOM pool (gC/m2)
};

struct RTSMethaneFluxes {
  double prod = 0.0, oxid = 0.0, emis = 0.0, plant = 0.0, pool = 0.0;
};

class RTSMethane {
public:
  static double saturatedFraction(const double ztop, const double dz, const double wtd);
  static double oxidFraction(const RTSMethaneParams & p, const double wtd);
  static double pondSaturation(const RTSMethaneParams & p, const int months_since_init);
  static RTSMethaneFluxes compute(std::vector<RTSMethaneLayer> & col, const RTSMethaneParams & p,
                                  const double wtd, const double pond, const double lai,
                                  const double days);
};
#endif
