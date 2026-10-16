#ifndef RTSMETHANE_H_
#define RTSMETHANE_H_
/** RTS-light v0.2 CH4, scheme (a): reduced-form monthly CH4 as a partition of RH.
 *
 * O(layers) per month, no column solver, no CH4 state.
 *  prod_l  = min(1, ch4_frac * rel_q10^((T_l - tref)/10)) * sat_l * RH_l      (T_l > 0)
 *            sat_l = saturated fraction of layer l below the (effective) water table
 *  plant   = prod * plant_frac * min(1, LAI/lai_ref)          (bypasses oxidation)
 *  oxid    = (prod - plant) * f_ox(WTD_eff)
 *            f_ox = oxid_min + (oxid_max - oxid_min) * (1 - exp(-WTD_eff/oxid_zscale))
 *  emis    = prod - oxid
 * CH4 C is taken out of RH (it is a partition, no extra C leaves the soil pools):
 *  RH_CO2  = RH - prod + oxid ,  RH_CO2 + emis = RH  (so the C balance is
 *  dC = GPP - RA - RH_CO2 - RHDWD - CH4emis - lateral - RTS respired).
 * RTS additions (both schemes):
 *  pond(m) = pond_sat0 * exp(-m / (12 * pond_tau)), m = months since RTS initiation
 *            (0 before initiation). sat_l -> sat_l + (1 - sat_l) * pond and
 *            WTD_eff = WTD * (1 - pond): a saturated / ponded slump floor that
 *            drains and revegetates over ~pond_tau years.
 *  mg(m)   = 1 - (1 - mg0) * exp(-m / (12 * mg_tau)) (1 before initiation): CH4
 *            production on the fresh slump floor is suppressed (no methanogen
 *            community in exposed permafrost) and recovers over ~mg_tau years.
 *  lat_ch4 = lat_ch4_frac * TKLATC: the part of the laterally exported C assumed
 *            to be mineralized as CH4 downstream (diagnostic split of TKLATC; the
 *            remainder is CO2/DOC/POC). It does not change the cell C balance.
 * Units: gC/m2/month, m, degC.  Default values and sources: docs_src/rts_ch4.md.
 */
#include <vector>
namespace Json { class Value; }

struct RTSMethaneParams {
  // scheme (a)
  double ch4_frac = 0.135;  // CH4 share of anaerobic RH at tref; CALIBRATED: 5 gC/m2/yr wet sedge tundra (docs)
  double rel_q10 = 1.5;     // extra Q10 of methanogenesis relative to RH
  double tref = 10.0;       // degC
  double plant_frac = 0.3;  // share of production emitted via aerenchyma at LAI >= lai_ref
  double lai_ref = 1.0;     // m2/m2
  // oxidation vs water-table depth (both schemes; transit oxidation in b)
  double oxid_min = 0.2;    // WT at surface
  double oxid_max = 0.95;   // deep WT
  double oxid_zscale = 0.1; // m
  // RTS additions (both schemes)
  double pond_sat0 = 0.2;   // ponding saturation factor right after initiation (0-1); chamber-constrained, see docs
  double pond_tau = 5.0;    // e-folding time of ponding (yr)
  // post-slump methanogen recovery (exposed Pleistocene/mineral floor lacks active
  // methanogens, Knoblauch et al. 2021): prod is scaled by
  // mg(m) = 1 - (1 - mg0) * exp(-m / (12 * mg_tau)) for m >= 0 months since initiation
  double mg0 = 0.05;        // methanogen activity factor right after initiation (0-1); CALIBRATED (docs)
  double mg_tau = 10.0;     // e-folding recovery time (yr)
  double lat_ch4_frac = 0.02; // part of laterally exported C mineralized as CH4 downstream
  // background anoxic/inundated fraction of the cell (0 = model water table only;
  // 1 = wet tundra with WT at the surface; e.g. from a wetland map such as BAWLD)
  double bg_sat = 0.0;
  void load(const Json::Value & ms); // thermokarst_ch4_<name> keys in model_settings
};

struct RTSMethaneLayer {
  double ztop = 0.0, dz = 0.1, tem = 0.0;
  double rh = 0.0;          // heterotrophic respiration of this layer (gC/m2/month)
};

struct RTSMethaneFluxes {
  double prod = 0.0, oxid = 0.0, emis = 0.0, plant = 0.0, pool = 0.0;
};

class RTSMethane {
public:
  static double saturatedFraction(const double ztop, const double dz, const double wtd);
  static double oxidFraction(const RTSMethaneParams & p, const double wtd);
  static double pondSaturation(const RTSMethaneParams & p, const int months_since_init);
  static double methanogenFactor(const RTSMethaneParams & p, const int months_since_init);
  static RTSMethaneFluxes compute(const std::vector<RTSMethaneLayer> & col,
                                  const RTSMethaneParams & p, const double wtd,
                                  const double pond, const double lai);
};
#endif
