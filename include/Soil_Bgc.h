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
  // Refs (N2): no published source; code bug in Soil_Bgc::deltastate (upstream 419a4eb3, fix commit 20b2723e); restores N mass conservation (N moves with C).
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
