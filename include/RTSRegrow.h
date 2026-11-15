#ifndef RTSREGROW_H_
#define RTSREGROW_H_
/** RTS post-slump vegetation recovery v2 (reduced form; switch thermokarst_regrow, default ON since v0.3;
 *  thermokarst_regrow=false -> bitwise v0.2).
 *  Spec: russia-inputs/rts_recovery_equations.md; docs: docs_src/rts_ch4.md ("Post-slump regrowth").
 *  t = years since stabilization = (months since initiation - 12 t_active) / 12.
 *  v1's per-PFT logistic growth ceiling (and its active-slump hold) is REMOVED in v2: it held GPP below the
 *  model's natural recovery for 12-36 yr and harmed cells that recover unaided (CMT05/07). No cmax/nmax factor
 *  below 1 is ever applied, so vegetation equations recover at their natural rate. Legacy keys
 *  rts_regrow_limit / rts_regrow_t50_* / rts_regrow_k_* / rts_regrow_release are ignored if present.
 *  (1) reseed at t = 0: every PFT is raised to f_reseed x its pre-slump C (tissue split and C:N of the pre-slump
 *      plant). External propagule import, reported in TKREGROWC/TKREGROWN (C/N balances close).
 *  (2) shrub boost on deciduous shrubs only: cmax and nmax x B(t) = 1 + (boost0 - 1) exp(-t / boost_T) >= 1,
 *      so it can only raise, never cap, shrub growth (kept from v1 for that reason).
 *  (3) organic-layer rebuild: up to dC/dt = ol_rate1 (t <= ol_T) / ol_rate2 gC m-2 yr-1 of litter C added to the
 *      top organic (else top) soil layer while the organic C is below the cumulative target (and below pre-slump
 *      organic C), with N at C:N ol_cn (external input, TKREGROWC/TKREGROWN; capped with the C at removed C).
 *  (4) floor N supply with GPP-recovery taper: base n_rate1 (t <= n_T) / n_rate2 gN m-2 yr-1 (Pizano et al.
 *      2014), minus the nominal OL N of (3), enters soil available N through m_a2soi.avlninput (like
 *      deposition; NINPUT output, 'avlninput' column of NBAL_DIAG). The available-N part is scaled by
 *      taperFactor = 1 - r^n_taper_pow clamped to [0,1], with r = trailing-12-month GPP / GPP_ref(t) and
 *      GPP_ref(t) = GPP_pre * INGPP(t) / INGPP_pre: the pre-slump GPP (model GPP over <= 5 yr before
 *      initiation) scaled by the change in N-unlimited GPP (INGPP, same trailing window), so the
 *      reference tracks the undisturbed climate/CO2 trend instead of freezing at the initiation year.
 *      It is a smooth plateau taper -- essentially the full supply while GPP is well below the reference
 *      (1 - 0.5^4 = 0.94 at half recovery), then falling steeply to zero once current GPP reaches the
 *      (contemporary) pre-slump level, so the supply cannot sustain the v1 no-limit over-fertilisation
 *      (end GPP 1.07-1.25 x RTS-off). (Static-reference tapers, linear 1 - r and 1 - r^4 against GPP_pre
 *      alone, were tested first in v2 and rejected: they starved slow cells, whose GPP plateaued at the
 *      absolute pre-slump level while undisturbed GPP kept rising; see regrowth_diag/regrow/README.md.)
 *      Set rts_regrow_n_taper=false to restore the untapered v1 supply.
 *  (5) net N mineralization uses the RH before the floor RH cap (cap applies to CO2 only) while regrow is on.
 *  tau() (Xia et al. 2026 site timescale) is retained only as a diagnostic logged at initiation. */
namespace Json { class Value; }

enum RTSRegrowType { RG_GRAM = 0, RG_SHRUB = 1, RG_MOSS = 2, RG_EVERGREEN = 3, RG_LICHEN = 4, RG_NTYPE = 5 };

struct RTSRegrowParams {
  bool on = true;   // v0.3 default; false -> bitwise v0.2 (tested)
  double t_active = 1.0;     // yr from initiation to stabilization
  double f_reseed = 0.07;    // fraction of pre-slump PFT C at stabilization
  double tau_ref = 7.0, tau_min = 2.0, tau_max = 150.0;  // yr (diagnostic only in v2)
  double boost0 = 1.3, boost_T = 80.0;                                  // -, yr
  double ol_rate1 = 32.0, ol_rate2 = 6.5, ol_T = 60.0, ol_cn = 30.0;  // gC/m2/yr, yr, g/g
  double ol_frac_raw = 0.0;  // part of the OL input added as litter (rawc); rest as slow SOM (sompr): Pizano rate is net accumulation
  double ol_maxmult = 1.0;  // monthly input cap = ol_maxmult x rate/12 (1: gross input at the Pizano net rate)
  double n_rate1 = 2.2, n_rate2 = 0.4, n_T = 60.0;                    // gN/m2/yr, yr
  bool n_taper = true;  double n_taper_pow = 4.0;  // GPP-recovery taper of the available-N supply (v2)
  bool navln_exempt = true;
  void load(const Json::Value & ms);
};

class RTSRegrow {
public:
  static int pftType(int ifwoody, int ifdeciwoody, int nonvascular);
  static double tau(const RTSRegrowParams & p, double gpp_kg);       // yr (diagnostic)
  static double tYears(const RTSRegrowParams & p, int months);       // years since stabilization
  static int stabMonth(const RTSRegrowParams & p);
  static double boost(const RTSRegrowParams & p, int type, double t); // >= 1, shrubs only
  static double olRate(const RTSRegrowParams & p, double t);   // gC/m2/yr (0 for t<0)
  static double nRate(const RTSRegrowParams & p, double t);    // gN/m2/yr (0 for t<0)
  static double taperFactor(const RTSRegrowParams & p, double gpp_ratio); // [0,1]; 1 low GPP, 0 recovered
  static double avlnSupplyMonth(const RTSRegrowParams & p, double t, double taper = 1.0); // gN/m2/month to avln
};
#endif
