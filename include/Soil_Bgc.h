#ifndef SOIL_BGC_H_
#define SOIL_BGC_H_

#include "Bnf.h"
#include "errorcode.h"

#include "CohortData.h"
#include "EnvData.h"
#include "FireData.h"
#include "BgcData.h"
#include "RestartData.h"

#include "CohortLookup.h"
#include "ModelData.h"

#include "parameters.h"

#include "Ground.h"

#include <cmath>

using namespace std;

class Soil_Bgc {
  friend class Integrator; // soil safety net reads the monthly inputs
public:
  Soil_Bgc();
  ~Soil_Bgc();

  soistate_bgc tmp_sois;   // the previous soistate_bgc
  soistate_bgc del_sois;   // the change of soistate_bgc

  soi2soi_bgc del_soi2soi;
  soi2atm_bgc del_soi2a;
  soi2lnd_bgc del_soi2l;

  atm2soi_bgc del_a2soi;

  soipar_bgc bgcpar;
  soipar_cal calpar;

  void setCohortData(CohortData* cdp);
  void setEnvData(EnvData* edp);
  void setBgcData(BgcData* bdp);
  void setFirData(FirData* fdp);
  void setCohortLookup(CohortLookup* chtlup);

  void setGround(Ground * ground);

  void initializeParameter();

  void initializeState();
  void set_state_from_restartdata(const RestartData & resdata);

  void assignCarbonBd2LayerMonthly();
  void assignCarbonLayer2BdMonthly();
  // v5q-layern-fix (layer_n_remap_fix 1 only): carry layer orgn/avln with the C pools around
  // Ground::redivideSoilLayers (legacy copies only C, so after a layer split/combine the bd N
  // stays at the old layer index while C moves). Refs: see Ground.h.
  void assignNitrogenBd2LayerMonthly();
  void assignNitrogenLayer2BdMonthly();

  void deltac();
  void deltan();
  void deltastate();

  void prepareIntegration(const bool &mdnfeedback, const bool &mdavlnflg,
                          const bool &mdbaseline);
  void afterIntegration();
  void clear_del_structs();

  int get_nfeed();
  void set_nfeed(int);

  int get_avlnflg();
  void set_avlnflg(int);

  int get_baseline();
  void set_baseline(int);

  // Fix N1 (v5k), ON by default since v5o (default 2; 0 = legacy, still selectable):
  // model_settings.n_rootless_fix. In soil layers
  // without roots (rtlfalfrac == 0) deltastate() sets the inorganic N change and
  // net mineralisation to 0, but only after the organic N change was computed with
  // the net mineralisation removed, so that N leaves ORGN and is never added to
  // AVLN (lost from the column; unreported N sink). 0 = legacy;
  // 1 = no net mineralisation in rootless layers (the N stays in ORGN, consistent
  //     with the zeroed NETNMIN and AVLN changes);
  // 2 = the mineralised N goes to AVLN of that layer (no root extraction there;
  //     leaching as for other layers).
  int n_rootless_fix = 2;   // v5o default (was 0)
  // Refs (N1): no published source; code bug in Soil_Bgc::deltastate (upstream 419a4eb3, fix commit 29d72af8); restores N mass conservation.
  // Fix N2 (v5k), ON by default since v5o (default 1; 0 = legacy, still selectable):
  // model_settings.n_som_move_fix. deltastate()
  // moves organic N with the SOM C that is moved down the profile, but (a) the N
  // added to the top humic layer(s) is overwritten ('=' instead of '+=') when the
  // mobile part of the humic layer's own decomposition products is exported, so
  // the shallow->deep N transfer disappears every step (the dominant unreported
  // N sink at HPC nfeed), (b) a humic receiving layer gets the full shallow N
  // transfer instead of its thickness share (dcaddfrac, as for C) and (c) each
  // fibric donor gives N in proportion to the CUMULATIVE C moved from all fibric
  // layers above it, not its own. 1 = N moves with C: per-layer donor share,
  // dcaddfrac share to receivers, export accumulated (+=).
  int n_som_move_fix = 1;   // v5o default (was 0)
  // v5q-d2m-gate: deep(humic)->mineral SOM transfer of deltastate() step (II)3)-4),
  // model_settings.d2m_transfer_gate (FATAL unless integer 0/1/2):
  //  0 = legacy (default; code path unchanged): the mobile part rh*fsoma*som2co2 of the
  //      decomposition of the type-2 layers in the top xtopdlthick (0.10 m) of the humic horizon
  //      moves to the top xtopmlthick (0.20 m) of the mineral horizon, frozen or not; receiver N
  //      is assigned (=); C/N is lost where the mineral column is thinner than 0.20 m.
  //  1 = gated: receiver thickness counts only where it lies above zgate = max(this year's
  //      running max ALD, last year's ALD) (d2mGateDepth; no transfer if zgate <= 0).
  //      f = eligible receiver thickness / 0.20 (d2mGateReceivers) scales the donor removal
  //      (C and N, N at the donor C:N); a donor layer either gives C and N together or nothing;
  //      receivers get it in proportion to their eligible thickness (weights sum to 1, so
  //      nothing is lost), N with +=. f = 0 (no thawed receiver) = no transfer at all.
  //  2 = transfer off (physics-only / no-transfer spin-up phase).
  // Refs: Koven et al. 2009, doi:10.1029/2009GL040150 (no advective SOM transport in
  //   permafrost); Koven et al. 2013, doi:10.5194/bg-10-7109-2013 and Lawrence et al. 2019 (CLM5),
  //   doi:10.1029/2018MS001583: vertical SOM transport only within the active layer, gated on
  //   max(altmax, altmax_lastyear), none where ALT = 0 (CTSM SoilBiogeochemLittVertTransp code).
  //   The f-scaling of the dvmdostem transfer and the restart variable: no published source.
  int d2m_transfer_gate = 0;
  double d2m_moved_c = 0.0, d2m_moved_n = 0.0, d2m_f = 0.0;   // gate 1: last deltastate() (diagnostic)
  // pure helpers (unit-tested). Gate depth: max(running max ALD, previous-year ALD), MISSING -> 0.
  static double d2mGateDepth(const double ald_runmax, const double ald_prevyr);
  // receiver weights w[0..n-1] (type-3 layers inside the top xtopmlthick, eligible part above
  // zgate; sum 1 if any eligible, else all 0); returns f = eligible thickness / xtopmlthick.
  static double d2mGateReceivers(const int n, const int type[], const double z[], const double dz[],
                                 const double xtopmlthick, const double zgate, double w[]);
  // Refs (N2): no published source; code bug in Soil_Bgc::deltastate (upstream 419a4eb3, fix commit 20b2723e); restores N mass conservation (N moves with C).
  // v5q-thawfrac-decomp: model_settings.thawfrac_decomp (FATAL unless integer 0/1).
  //  0 = legacy (default; code path unchanged): deltac() treats a layer as wholly frozen when its
  //      monthly mean ts < 0 (rhmoist = rhmoistfrozen) and wholly thawed otherwise, so a layer that
  //      is half thawed but whose mean ts is pinned just below 0 C by latent heat has no RH where
  //      rhmoistfrozen = 0 (CMT51); deltan() mineralizes N only where ts > 0.
  //  1 = for each soil layer, with f = 1 - m_sois.frozenfrac (the frozen thickness fraction set from
  //      the Stefan freeze/thaw fronts, Ground::getLayerFrozenstatusByFronts, monthly mean; clamped
  //      to [0,1]; legacy ts sign if not finite/missing):
  //        E = f * rhmoist(sws_th) * Q10^(max(ts,0)/10) + (1-f) * rhmoistfrozen * Q10^(min(ts,0)/10)
  //      thawed part: temperature max(ts, 0) (>= 0 C, the fronts are the 0 C isotherm; the layer-mean
  //      ts of a partly frozen layer sits at ~0 C) and liquid saturation sws_th = min(1, sws/f) for f < 1 (the
  //      liquid of a partly frozen layer is held in its unfrozen part, as in Soil_Env / Richards
  //      liq = dz*(1-frozenfrac)*lwc); frozen part: the legacy frozen path at min(ts, 0).
  //      E multiplies kd*C of every pool (rawc, soma, sompr, somcr; wood debris inherits the first
  //      non-moss layer's E as in legacy); bd rhmoist[il] = E/rhq10[il] (rhq10 diagnostic legacy).
  //      N: net mineralization/immobilization (getNetmin) uses the thawed-part RH share
  //      f*E_th/E of the layer RH with the thaw gate open, the frozen part keeps the legacy closed
  //      gate, so at f = 1 (ts > 0) and f = 0 (ts < 0) deltan() equals legacy. f = 0 / 1 give
  //      the legacy frozen / thawed rate; linear in f in between. No calibrated parameter changed.
  // Refs: vertically resolved, per-layer decomposition with frozen-soil limitation: Koven et al. 2013,
  //   doi:10.5194/bg-10-7109-2013; Lawrence et al. 2019 (CLM5), doi:10.1029/2018MS001583. Layer
  //   frozen fraction from the two-directional Stefan freeze/thaw fronts: Yi et al. 2009,
  //   doi:10.1029/2008JG000841. The linear thawed-fraction blend of the dvmdostem modifier
  //   (sub-layer mixing of the thawed and frozen parts): no published source.
  int thawfrac_decomp = 0;
  // v5q-thawfrac-lmoist: model_settings.thawfrac_moist_mode (FATAL unless integer 0/1; 1 requires
  //   thawfrac_decomp = 1). Thawed-part liquid saturation sws_th used by thawfracModifier:
  //   0 = min(1, sws/f) for f < 1, sws for f = 1 (default; v5q-thawfrac-decomp unchanged: the layer
  //       liquid is assumed to sit in the unfrozen part);
  //   1 = layer-mean sws for every f (no sub-layer redistribution of the layer liquid).
  //   Only the moisture argument of rhmoist() in the thawed part changes; f, the temperatures, the
  //   frozen part, the N thawed share and all parameters are as thawfrac_decomp 1. Equal at f = 1.
  // Refs: layer-resolved decomposition with a per-layer moisture scalar from the layer liquid-water
  //   state and frozen-soil limitation: Koven et al. 2013, doi:10.5194/bg-10-7109-2013; Lawrence et
  //   al. 2019 (CLM5), doi:10.1029/2018MS001583. Layer frozen fraction from the Stefan freeze/thaw
  //   fronts: Yi et al. 2009, doi:10.1029/2008JG000841. The choice between thawed-part (sws/f) and
  //   layer-mean moisture for a partly frozen layer: no published source.
  int thawfrac_moist_mode = 0;
  // pure helper (unit-tested): thawed-part saturation for fthaw in (0,1] and mode 0/1 (see above).
  static double thawedPartSws(const double sws, const double fthaw, const int mode);
  double tfd_thawshare[MAX_SOI_LAY] = {0};   // thawfrac_decomp 1: thawed share of the layer RH (last deltac)
  // pure helper (unit-tested): thawed fraction f from a frozen fraction (clamped to [0,1]; not finite
  // or outside [-0.5, 1.5] (MISSING) -> legacy indicator ts < 0 ? 0 : 1).
  static double thawFraction(const double frozenfrac, const double ts);
  // thawfrac_decomp 1: E (returned) and its thawed share (f*E_th/E, 0 if E == 0) for one layer.
  double thawfracModifier(const double fthaw, const double ts, const double sws, double & thawshare);
  // v5q-rh-frozen: as above, with the frozen-part rate blended per layer:
  //   efr = (1 - w_al)*rhmoistfrozen*Q10^(min(ts,0)/10) + w_al*efr_al  (w_al = 0: the 4-argument form).
  double thawfracModifier(const double fthaw, const double ts, const double sws, double & thawshare,
                          const double w_al, const double efr_al);
  // v5q-rh-frozen: model_settings.rh_frozen_mode (FATAL unless integer 0/1; 1 requires thawfrac_decomp
  //   = 1). Frozen-part (1 - f) decomposition rate of the thawed-fraction modifier:
  //   0 = rhmoistfrozen * Q10^(min(ts,0)/10) in every frozen layer part (default, legacy; TVC CMT51 has
  //       rhmoistfrozen = 0, so no frozen-soil RH at all; > 0 respires the whole permafrost column).
  //   1 = for the fraction w_al of the layer above the permafrost table zt (w_al = clamp((zt - z)/dz, 0, 1),
  //       zt = max(this month's ALD, this year's running-max ALD, previous year's ALD), the same
  //       previous-year ALD as d2m_transfer_gate): efr_al = phi(T) * rhmoist(s_tot) * Q10^(T/10),
  //       T = min(ts, 0); below zt the legacy rate. Schaefer & Jafarov (2016): only the carbon in the
  //       thin unfrozen water films ("thawed" carbon, fraction phi = liquid water fraction) decays below
  //       0 degC, with the normal kinetic temperature and moisture factors (this replaces an extra
  //       frozen Q10, which made permafrost C decay). phi(T): Eq. 4 power law phi_i = ((Tref - T)/T*)^b_i,
  //       Tref = 0.1 degC, T* = 0.01 degC, b = -0.3 clay, -0.5 silt, -0.9 sand, -1.0 organic (Table 1);
  //       Eq. 5 phi = (1 - forg)(fc phi_c + fsi phi_si + fsa phi_sa) + forg phi_o; Eq. 6 forg = Mo/(dz rho_o),
  //       Mo = 2 x layer SOM C (50 % C), rho_o = 140 kg m-3 (capped at 1). dvmdostem frozen layer parts
  //       hold no liquid (the Stefan scheme puts the layer liquid in the thawed part), so phi comes from
  //       the power law, not from the layer liquid. s_tot = sws + iwc/poro (liquid + ice pore saturation,
  //       the SiBCASA moisture argument in SJ2016), evaluated with the dvmdostem rhmoist() curve and the
  //       CMT moistmin/opt/max. Deviations from SJ2016 (no published source): dvmdostem moss/peat layers
  //       (type 0-2) use forg = 1; mineral layers without texture use forg = 1; the film-frozen ->
  //       thawed pool transfers of SJ2016 are not modelled (phi scales the rate of the single layer
  //       pool, the dvmdostem pools are not split); the frozen-part RH keeps the legacy closed
  //       N-mineralisation gate (tfd_thawshare = thawed share), as for rhmoistfrozen > 0.
  // Refs: Schaefer & Jafarov 2016, doi:10.5194/bg-13-1991-2016 (Eqs. 3-6, Table 1);
//   Romanovsky & Osterkamp 2000, doi:10.1002/1099-1530(200007/09)11:3<219::AID-PPP352>3.0.CO;2-7;
//   Nicolsky et al. 2009, doi:10.1016/j.coldregions.2008.03.003; Kurylyk & Watanabe 2013,
//   doi:10.1016/j.advwatres.2013.07.016 (unfrozen-water power law); active-layer gate as
//   d2m_transfer_gate: Koven et al. 2013, doi:10.5194/bg-10-7109-2013, Lawrence et al. 2019,
//   doi:10.1029/2018MS001583; frozen-soil respiration observed: Mikan et al. 2002,
//   doi:10.1016/S0038-0717(02)00168-2, Natali et al. 2019, doi:10.1038/s41558-019-0592-8.
  //   The permafrost-table restriction and the forg = 1 choices: no published source.
  int rh_frozen_mode = 0;
  // rh_frozen_mode 1: per soil layer, set from the Ground layers at prepareIntegration (if ground set):
  double rhf_fclay[MAX_SOI_LAY] = {0}, rhf_fsilt[MAX_SOI_LAY] = {0}, rhf_fsand[MAX_SOI_LAY] = {0};
  double rhf_por[MAX_SOI_LAY] = {0};
  int rhf_texok[MAX_SOI_LAY] = {0};   // 1 = mineral layer with texture (fractions sum > 0)
  void rhFrozenSetLayers();
  // rh_frozen_mode 1: active-layer frozen-part rate efr_al of soil layer il at layer temperature ts.
  double rhFrozenEfr(const int il, const double ts);
  // pure helpers (unit-tested), SJ2016 Eqs. 4-6, Table 1:
  static double sjPhiPure(const double t, const double b);         // 1 for t > 0, else ((0.1 - t)/0.01)^b
  static double sjLayerPhi(const double t, const double forg, const double fc, const double fsi,
                           const double fsa);                        // Eq. 5 (texture renormalised)
  static double sjOrgFraction(const double somc_gm2, const double dz);  // Eq. 6, capped at 1
  static double rhFrozenTableDepth(const double ald_mon, const double ald_yr, const double ald_prevyr);
  static double fracAboveDepth(const double z, const double dz, const double zt);
  static double totalSaturation(const double sws, const double iwc, const double por);
  // Fix G1 (v5n), OFF by default: model_settings.n_extract_limiter.
  // The root N extraction rtnextract[il] is fixed for the whole month (it equals
  // the vegetation N uptake, computed from the month-start AVLN and integrated
  // before the soil), so inside the soil integrator it is a zero-order sink on
  // AVLN. When AVLN runs out part-way through the month (immobilisation,
  // leaching, warm-month uptake ~ the whole pool) AVLN crosses zero in finite
  // time; RKF45 rejects every step that makes a pool negative, halves dt down to
  // 0.5^20 and finally gives up (blackhole -> safety net). 1 = flux limiter:
  // the AVLN part of the extraction is r*min(1, AVLN/(tau*r)) (tau =
  // n_extract_tau, month fraction), i.e. unchanged while AVLN covers tau months
  // of demand, and the rest r - that is taken from ORGN of the same layer
  // (organic N uptake), so the total extraction (= vegetation uptake) and N
  // mass are unchanged. Applied only where deltastate() applies extraction.
  static int n_extract_limiter;
  static double n_extract_tau;   // month fraction, default 0.2 (v5n tests: 0.02, 0.05, 0.2) (0 < tau <= 1)
  // G2 (v5n, Integrator::soil_giveup_retry): limiter on for one retry month
  bool nx_limiter_active = false;
  bool nxLimiterOn() const { return n_extract_limiter == 1 || nx_limiter_active; }
  // G1 helper (pure, unit-tested): AVLN part of an extraction r (g/m2/month)
  // from a pool A (g/m2): r for r <= 0 or A >= tau*r, else r*max(A,0)/(tau*r).
  static double limitedExtraction(const double r, const double A, const double tau);
  // Fix G3 (v5n), OFF by default: model_settings.n_netnmin_rescale_guard.
  // deltan() closes the inorganic N budget, when supply (month-start AVLN of the
  // drainage zone + input + net mineralisation) cannot cover the extraction, by
  // multiplying every layer's net mineralisation by nscale = nminadj/totnetnmin.
  // When the profile is net immobilising (totnetnmin < 0) and nminadj > 0,
  // nscale < 0: the sign of every layer flips, so mineralising (rooted) layers
  // are turned into immobilising ones and their AVLN is drained at a fixed
  // rate. 1 = skip the rescale when nscale < 0 (net mineralisation unchanged);
  // positive nscale is applied as in legacy. Meant to be used with G1.
  // Refs (G3): no published source for the fix (code: Soil_Bgc::deltan, the N-budget rescale of net N mineralisation inherited from TEM / DOS-TEM); principle: a scale factor applied to a signed flux must be >= 0, otherwise mineralising layers become immobilising and vice versa; net N mineralisation / immobilisation as the balance of gross mineralisation and microbial uptake (Schimel & Bennett 2004, doi:10.1890/03-8002); immobilisation limited by available mineral N in CLM-CN (Thornton et al. 2007, doi:10.1029/2006GB002868).
  static int n_netnmin_rescale_guard;
  // G3 helper (pure, unit-tested): the legacy rescale; returns the raw nscale
  // (1.0 if not computed); with guard == 1 a negative nscale is not applied.
  static double rescaleNetNmin(double netnmin[], const int n, const double nminadj,
                               const double totnetnmin, const int guard);
  // v5n diagnostics (written by every deltan(); read only by the give-up diag)
  double dg_nscale = 1.0, dg_totnetnmin_raw = 0.0, dg_totnextract = 0.0, dg_totdzavln = 0.0;
  // Refs (G1/G2): plant N uptake is a saturating (Michaelis-Menten) function of soil available N in TEM (Raich et al. 1991, doi:10.2307/1941899; McGuire et al. 1992, doi:10.1029/92GB00219); demand-based uptake is down-regulated when soil mineral N cannot supply it (Thornton et al. 2007, doi:10.1029/2006GB002868); arctic/boreal plants also take up organic N (Nasholm et al. 1998, doi:10.1038/31921; Schimel & Bennett 2004, doi:10.1890/03-8002); a zero-order sink drives a non-negative ODE state below zero in finite time, which no step-size reduction can avoid, positivity needs a state-dependent sink (Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011; Burchard et al. 2003, doi:10.1016/S0168-9274(03)00101-6).
  // v5p daily_c_split (output only, see DailyCSplit.h): daily RH weight of a
  // soil layer, rhmoist*rhq10 as in deltac (rhmoistfrozen when ts < 0).
  // Refs: see DailyCSplit.h; no published source for the exact method.
  double ds_rh_weight(double ts, double sws);
  // v5q-ndrain: N leaching bookkeeping fix, model_settings.n_leach_fix (0 = legacy, default; 1 = fix).
  // deltan() books the inorganic N leaching avlnlost = (totdzavln/totdzliq)*fnloss*qdrain
  // + (avln[0]/liq[0])*fnloss*qover (NLOST, integrated as I_AVLNLOSS), but legacy deltastate()
  // removes ndrain = avlnlost*(liq[il]/totdzliq*qdrain)*fdz from the drainage-zone layers, i.e.
  // about avlnlost*qdrain in total (qdrain in mm/month), takes nothing for the qover part and
  // nothing when totdzliq <= 0.01. 1 = the pools lose exactly avlnlost: the qover part from
  // layer 0, the drainage part from the drainage-zone layers in proportion to their current
  // (integrator-state) AVLN
  // (avln[il]*fdz/totdzavln = each layer drains a share liq[il]*fdz/totdzliq of the water at
  // its own concentration, so no layer is charged for N it does not hold); fallbacks: by
  // liq*fdz/totdzliq if totdzavln <= 0, layer 0 if totdzliq <= 0 too. Same for avlnflg 0.
  // Refs: no published source; code bug in Soil_Bgc::deltastate (also in upstream
  //   dvm-dos-tem c0fa75cb, src/Soil_Bgc.cpp l.1004-1018); restores N mass conservation
  //   (NLOST = AVLN removed). Observed arctic DIN export context: McNamara et al. 2008,
  //   doi:10.1002/hyp.6920; Khosh et al. 2017, doi:10.1002/2016jg003377.
  int n_leach_fix = 0;
  double leach_drain = 0.0, leach_over = 0.0;   // n_leach_fix 1: parts of avlnlost (set in deltan)
  // v5q-leach-curn: model_settings.n_leach_curn_fix (0 = legacy, default; 1 = fix; needs n_leach_fix 1).
  // Legacy deltan(): drainage part of avlnlost = totdzavln/totdzliq*fnloss*qdrain with totdzavln
  // built once per month in prepareIntegration() from the month-start AVLN (bd->m_sois.avln). It
  // is a constant-rate (zero-order) sink during the month: with n_leach_fix 1 it is charged to
  // the drainage-zone layers even after their AVLN is used up (e.g. a frozen moss layer 0 that
  // holds the whole drainage zone, emptied by the overland part under a cap raised by netnmin
  // of deeper layers), AVLN goes negative and adapt() gives up (V5Q_STEP8 doc section 13).
  // 1 = the drainage part uses the current (integrator-state, tmp_sois) AVLN of the same zone
  // (drainZoneAvln, same fdz as prepareIntegration): rate = k*AVLN_zone(t), which goes smoothly
  // to 0 as the zone empties (the overland part already uses tmp_sois.avln[0]). The cap, its
  // netnmin rescale, the drain/overland split and leachLayers are unchanged (cap still uses the
  // month-start totdzavln). totdzliq, qdrain, fnloss unchanged.
  // Refs: positivity of a non-negative ODE state needs a sink that vanishes with the state
  //   (Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011; Burchard et al. 2003,
  //   doi:10.1016/S0168-9274(03)00101-6); leaching context as n_leach_fix (McNamara et al. 2008,
  //   doi:10.1002/hyp.6920; Khosh et al. 2017, doi:10.1002/2016jg003377). The fix itself:
  //   no published source (code defect, see Soil_Bgc.h).
  int n_leach_curn_fix = 0;
  // pure helper (unit-tested): sum over the drainage zone of max(0, avln[i])*fdz[i], zone and fdz
  // exactly as prepareIntegration() builds totdzavln.
  static double drainZoneAvln(const int n, const double z[], const double dz[], const double avln[],
                              const double draindepth);
  // pure helper (unit-tested): per-layer leaching removal ndrain[0..n-1] for n_leach_fix 1;
  // returns sum(ndrain) (= drain + over up to round-off).
  static double leachLayers(const int n, const double z[], const double dz[], const double liq[],
                            const double avln[], const double draindepth, const double drain,
                            const double over, double ndrain[]);
  // v5q: biological N fixation (model_settings.bnf_scheme, 0 = off): see Bnf.h.
  // bnfBeginMonth() (Cohort, before the vegetation) computes this month's BNF
  // from the previous model year; prepareIntegration() puts it into
  // m_a2soi.avlninput (AVLN of the top layer, needs n_rootless_fix = 2);
  // bnfRecordMonth() (Cohort, after the soil) stores this month's NPP / ET.
  // Nothing of this runs with bnf_scheme = 0.
  // Refs: Meyerholt et al. 2016, doi:10.5194/bg-13-1491-2016; Cleveland et al. 1999,
  //   doi:10.1029/1999GB900014; Thornton et al. 2007, doi:10.1029/2006GB002868 (BNF to soil mineral N).
  bnf::Params bnfp;
  bnf::YearMemory bnfmem;
  double bnf_soil_month = 0.0;   // g/m2/month into AVLN
  double bnf_veg_month = 0.0;    // Phase 2 (bnf_to_plant 1): g/m2/month to plant labile N (Cohort shares it among PFTs)
  double bnf_veg_delivered = 0.0;// Phase 2: amount the vegetation integrators delivered (veg give-up: fraction bh_t0)
  void bnfRecordPft(const int ip, const double npp);
  void bnfBeginMonth(const int m, const bool nfeedon);
  void bnfRecordMonth(const int m, const double npp, const double et);

private:

  // NOTE: seems like these 3 variables are supposed to shadow the
  // corresponding ModelData members. Not sure if this should be
  // strictly enforced somehow?

  int nfeed;  // soil-plant-air N module switch

  int avlnflg; // open-N cycle switch; otherwise, N budget method used to
               // balance the N I/O from the ecosystem

  int baseline;  // open-N cycle switch; otherwise, N budget method used to
                 // balance the N I/O from the ecosystem

  double d2wdebrisc;
  double d2wdebrisn;
  double mossdeathc;
  double mossdeathn;
  double ltrflc[MAX_SOI_LAY];     //litterfall C into each soil layer
  double ltrfln[MAX_SOI_LAY];     //litterfall N into each soil layer

  double rtnextract[MAX_SOI_LAY];  // root N extraction from each soil layer

  double kdshlw;
  double kddeep;
  double decay;
  double nup;

  double totdzliq;     // total drainage zone liq water
  double totdzavln;    // total drainage zone avln

  double totnetnmin;   // total net N mineralization
  double totnextract;  // total root N extract (uptake) from soil

  CohortData * cd;
  EnvData * ed;
  BgcData * bd;
  FirData * fd;

  CohortLookup * chtlu;
  Ground * ground;

  void initSoilCarbon(double & initshlwc, double & initdeepc,
                      double & initminec);
  void initOslayerCarbon(double & shlwc, double & deepc);
  void initMslayerCarbon(double & minec);

  double getRhmoist(const double &vsm, const double & moistmin,
                    const double & moistmax, const double & moistopt);
  double getRhq10(const double & tsoil);

  double getNimmob(const double & soilh2o, const double & soilorgc,
                   const double & soilorgn, const double & availn,
                   const double & ksoil, const double kn2);

  double getNetmin(const double & nimmob, const double & soilorgc,
                   const double & soilorgn,
                   const double & rh, const double & tcnsoil,
                   const double & decay, const double & nup, const double & soilts );


  double getKnsoilmoist(const double & vsm);   //

  void updateKdyrly4all();

  double getKdyrly(double& yrltrcn, const double lcclnc, const double & kdc);

};
#endif /*SOIL_BGC_H_*/
