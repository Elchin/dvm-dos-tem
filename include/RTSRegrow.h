#ifndef RTSREGROW_H_
#define RTSREGROW_H_
/** RTS post-slump vegetation recovery (reduced form; switch thermokarst_regrow, default OFF -> bitwise v0.2).
 *  Spec: russia-inputs/rts_recovery_equations.md; docs: docs_src/rts_ch4.md ("Post-slump regrowth").
 *  t = years since stabilization = (months since initiation - 12 t_active) / 12.
 *  (0) active slump (t < 0): the 99% vegetation kill is kept; growth limited to (1 - veg_remove_frac) of pre-slump C.
 *  (1) reseed at t = 0: every PFT is raised to f_reseed x its pre-slump C (tissue split and C:N of the pre-slump
 *      plant). External propagule import, reported in TKREGROWC/TKREGROWN (C/N balances close).
 *  (2) growth limit (nudge, not an override): per PFT p, L_p(t) = max(f_reseed, 1/(1+exp(-k_p (t - s t50_p)))),
 *      s = tau_site / tau_ref, tau_site = 1.35 GPP_pre^-1.68 (Xia et al. 2026: yr, regional mean annual GPP in
 *      kgC m-2 yr-1; GPP_pre = model GPP over the <= 5 yr before initiation), tau_ref = 7 yr. If the PFT's C exceeds
 *      B_p(t) L_p(t) C_pre,p, its cmax is scaled by that ratio (GPP limited, not removed); released for L_p >= release.
 *  (3) shrub boost on deciduous shrubs: cmax and nmax x B(t) = 1 + (boost0 - 1) exp(-t / boost_T).
 *  (4) organic-layer rebuild: up to dC/dt = ol_rate1 (t <= ol_T) / ol_rate2 gC m-2 yr-1 of litter C added to the
 *      top organic (else top) soil layer while the organic C is below the cumulative target (and below pre-slump
 *      organic C), with N at C:N ol_cn taken from the floor N supply (external input, TKREGROWC/TKREGROWN).
 *  (5) floor N supply: n_rate1 (t <= n_T) / n_rate2 gN m-2 yr-1 (Pizano et al. 2014), minus the nominal OL N of (4),
 *      enters soil available N through the atmospheric-input hook m_a2soi.avlninput (like deposition; NINPUT output,
 *      'avlninput' column of the NBAL_DIAG line).
 *  (6) net N mineralization uses the RH before the floor RH cap (cap applies to CO2 only) while regrow is on. */
namespace Json { class Value; }

enum RTSRegrowType { RG_GRAM = 0, RG_SHRUB = 1, RG_MOSS = 2, RG_EVERGREEN = 3, RG_LICHEN = 4, RG_NTYPE = 5 };

struct RTSRegrowParams {
  bool on = false;
  double t_active = 1.0;     // yr from initiation to stabilization
  double f_reseed = 0.07;    // fraction of pre-slump PFT C at stabilization
  double tau_ref = 7.0, tau_min = 2.0, tau_max = 150.0;  // yr
  double t50[RG_NTYPE] = { 4.0, 30.0, 30.0, 100.0, 100.0 };   // yr (low Arctic, before site scaling)
  double k[RG_NTYPE]   = { 0.8, 0.15, 0.10, 0.03, 0.03 };     // 1/yr
  double boost0 = 1.3, boost_T = 80.0;                                  // -, yr
  double ol_rate1 = 32.0, ol_rate2 = 6.5, ol_T = 60.0, ol_cn = 30.0;  // gC/m2/yr, yr, g/g
  double ol_frac_raw = 0.0;  // part of the OL input added as litter (rawc); rest as slow SOM (sompr): Pizano rate is net accumulation
  double ol_maxmult = 1.0;  // monthly input cap = ol_maxmult x rate/12 (1: gross input at the Pizano net rate)
  double n_rate1 = 2.2, n_rate2 = 0.4, n_T = 60.0;                    // gN/m2/yr, yr
  bool limit_on = true;  double release = 0.95;
  bool navln_exempt = true;
  void load(const Json::Value & ms);
};

class RTSRegrow {
public:
  static int pftType(int ifwoody, int ifdeciwoody, int nonvascular);
  static double tau(const RTSRegrowParams & p, double gpp_kg);       // yr
  static double siteScale(const RTSRegrowParams & p, double gpp_kg); // tau / tau_ref
  static double tYears(const RTSRegrowParams & p, int months);       // years since stabilization
  static int stabMonth(const RTSRegrowParams & p);
  static double targetFrac(const RTSRegrowParams & p, int type, double t, double s, double hold);
  static double boost(const RTSRegrowParams & p, int type, double t);
  static double olRate(const RTSRegrowParams & p, double t);   // gC/m2/yr (0 for t<0)
  static double nRate(const RTSRegrowParams & p, double t);    // gN/m2/yr (0 for t<0)
  static double avlnSupplyMonth(const RTSRegrowParams & p, double t); // gN/m2/month to avln
  static double limitMult(double c_now, double c_limit);       // cmax multiplier in (0,1]
};
#endif
