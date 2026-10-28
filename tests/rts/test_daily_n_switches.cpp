// Unit tests for the three N/solver switches ported from hpc-tvc daily-climate-mode onto
// rts-regrow (same lightweight CHECK framework as test_rts_units.cpp; no external framework).
// Build/run: make rts-unit-test   (run from the repository root)
//
//   switch (daily commits)                       tests
//   N3 n_nonvasc_extract                         test_n3_*
//      (5812ff19 / 20b2723e / d3d3fbba)          refs: Raich et al. 1991 doi:10.2307/1941899;
//                                                McGuire et al. 1992 doi:10.1029/92GB00219;
//                                                Turetsky 2003 doi:10.1639/05
//   G1 n_extract_limiter / n_extract_tau         test_g1_*
//      (904a4a7c / 5c597aec)                     refs: Thornton et al. 2007 doi:10.1029/2006GB002868;
//                                                Nasholm et al. 1998 doi:10.1038/31921;
//                                                Schimel & Bennett 2004 doi:10.1890/03-8002
//   soil give-up counters + soil_giveup_diag     test_giveup_counter_*
//      (diagnostic, cf. hpc-tvc daily solver     refs: Shampine 2005 doi:10.1016/j.amc.2004.12.011;
//       give-up reporting)                       Burchard et al. 2003 doi:10.1016/S0168-9274(03)00101-6
//   config keys / defaults / validation          test_switch_config_*
//   end-to-end (pytest, test_rts_regression.py): test_soil_giveup_diag_synthetic and
//   test_soil_giveup_diag_cmt70_forced_giveups (diag on == off bit-identically; counters ==
//   per-event / per-year / summary reports on a real forced give-up case).
//
// Every test checks the switch-off path too: with the switch off the result must equal the
// legacy code path exactly (==, not a tolerance).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <functional>
#include <stdexcept>
#include <unistd.h>
#include <sys/wait.h>
#include <json/json.h>
#include <boost/filesystem.hpp>
#include "TEMLogger.h"
#include "TEMUtilityFunctions.h"
#define private public
#define protected public
#include "ModelData.h"
#include "Soil_Bgc.h"
#include "Vegetation_Bgc.h"
#include "Integrator.h"
#undef private
#undef protected
#include "CohortData.h"
#include "EnvData.h"
#include "BgcData.h"

static int g_fail = 0, g_pass = 0;
#define CHECK(cond) do { if (cond) { g_pass++; } else { g_fail++; \
  std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)
#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

static int run_isolated(const std::function<int()>& fn, unsigned timeout_s = 30) {
  std::fflush(stdout); std::fflush(stderr);
  pid_t pid = fork();
  if (pid == 0) {
    if (!std::freopen("/dev/null", "w", stderr)) {}
    if (!std::freopen("/dev/null", "w", stdout)) {}
    alarm(timeout_s);
    int rc = 3;
    try { rc = fn(); } catch (const std::invalid_argument&) { rc = 5; } catch (...) { rc = 4; }
    _exit(rc);
  }
  int st = 0;
  waitpid(pid, &st, 0);
  if (WIFEXITED(st)) return WEXITSTATUS(st);
  if (WIFSIGNALED(st)) return 128 + WTERMSIG(st);
  return -1;
}

// ---------------------------------------------------------------- config keys
static void test_switch_config_defaults() {
  ModelData md;
  CHECK(md.n_nonvasc_extract == false);
  CHECK(md.n_extract_limiter == false);
  CHECK(md.n_extract_tau == 0.2);
  CHECK(md.soil_giveup_diag == false);
  Vegetation_Bgc vb; Soil_Bgc sb; Integrator in;
  CHECK(vb.n_nonvasc_extract == false);
  CHECK(sb.n_extract_limiter == false);
  CHECK(sb.n_extract_tau == 0.2);
  CHECK(in.sbgc_months == 0 && in.sbgc_giveups == 0 && in.sbgc_retry_fails == 0);
  CHECK(in.sbgc_rkf_calls == 0 && in.sbgc_forced_steps == 0);
}

static void test_switch_config_parse() {
  if (!boost::filesystem::exists("config/config.js")) {
    std::printf("  SKIP config parse test (run from repo root)\n"); return;
  }
  Json::Value cfg = temutil::parse_control_file("config/config.js");
  const Json::Value& ms = cfg["model_settings"];
  // shipped config: keys present and off
  CHECK(ms.isMember("n_nonvasc_extract") && ms["n_nonvasc_extract"].asBool() == false);
  CHECK(ms.isMember("n_extract_limiter") && ms["n_extract_limiter"].asBool() == false);
  CHECK(ms.isMember("n_extract_tau") && ms["n_extract_tau"].asDouble() == 0.2);
  CHECK(ms.isMember("soil_giveup_diag") && ms["soil_giveup_diag"].asBool() == false);
  int rc = run_isolated([&]() {
    ModelData md(cfg);
    return (!md.n_nonvasc_extract && !md.n_extract_limiter && md.n_extract_tau == 0.2 && !md.soil_giveup_diag) ? 0 : 1;
  });
  CHECK(rc == 0);
  Json::Value on = cfg;
  on["model_settings"]["n_nonvasc_extract"] = true;
  on["model_settings"]["n_extract_limiter"] = 1;   // 0/1 accepted as well
  on["model_settings"]["n_extract_tau"] = 0.5;
  on["model_settings"]["soil_giveup_diag"] = true;
  rc = run_isolated([&]() {
    ModelData md(on);
    return (md.n_nonvasc_extract && md.n_extract_limiter && md.n_extract_tau == 0.5 && md.soil_giveup_diag) ? 0 : 1;
  });
  CHECK(rc == 0);
  // keys absent (older configs): defaults
  Json::Value old = cfg;
  old["model_settings"].removeMember("n_nonvasc_extract");
  old["model_settings"].removeMember("n_extract_limiter");
  old["model_settings"].removeMember("n_extract_tau");
  old["model_settings"].removeMember("soil_giveup_diag");
  rc = run_isolated([&]() {
    ModelData md(old);
    return (!md.n_nonvasc_extract && !md.n_extract_limiter && md.n_extract_tau == 0.2 && !md.soil_giveup_diag) ? 0 : 1;
  });
  CHECK(rc == 0);
  // tau must be a month fraction in (0, 1]
  const double bad[] = {0.0, -0.1, 1.5};
  for (double t : bad) {
    Json::Value b = cfg; b["model_settings"]["n_extract_tau"] = t;
    rc = run_isolated([&]() { ModelData md(b); return 0; });
    CHECK(rc == 5);   // std::invalid_argument
  }
  Json::Value one = cfg; one["model_settings"]["n_extract_tau"] = 1.0;
  rc = run_isolated([&]() { ModelData md(one); return md.n_extract_tau == 1.0 ? 0 : 1; });
  CHECK(rc == 0);
}

// ---------------------------------------------------------------- N3 fixtures
// Column of the hpc-tvc N3 test: moss | fibric | humic | mineral, AVLN .1/.2/.3/.4, liq 10.
struct VegFix {
  CohortData cd; EnvData ed; BgcData bd; Vegetation_Bgc vb;
  VegFix(bool n3, int nonvasc = 1) {
    vb.setCohortData(&cd); vb.setEnvData(&ed); vb.setBgcData(&bd);
    vb.ipft = 0; vb.set_nfeed(true); vb.n_nonvasc_extract = n3;
    cd.m_veg.nonvascular[0] = nonvasc;
    cd.m_soil.numsl = 4;
    const int type[4] = {0, 1, 2, 3}; const double avln[4] = {0.1, 0.2, 0.3, 0.4};
    for (int il = 0; il < MAX_SOI_LAY; ++il) {
      cd.m_soil.type[il] = il < 4 ? type[il] : 3;
      bd.m_sois.avln[il] = il < 4 ? avln[il] : 0.0;
      ed.m_sois.liq[il] = il < 4 ? 10.0 : 0.0;
      vb.fracnuptake[il] = 0.0; bd.m_soi2v.nextract[il] = -99.0;
      cd.m_soil.frootfrac[il][0] = 0.0;
    }
    for (int i = 0; i < NUM_PFT_PART; ++i) {
      bd.m_vegs.c[i] = 0.; bd.m_vegs.strn[i] = 0.; bd.m_a2v.ingpp[i] = 0.; bd.m_a2v.innpp[i] = 0.;
      bd.m_a2v.gpp[i] = 0.; bd.m_a2v.npp[i] = 0.; bd.m_v2a.rg[i] = 0.; bd.m_v2a.rm[i] = 0.;
      bd.m_v2soi.ltrfalc[i] = 0.; bd.m_v2soi.ltrfaln[i] = 0.; bd.m_v2v.nmobil[i] = 0.;
      bd.m_v2v.nresorb[i] = 0.; bd.m_soi2v.snuptake[i] = 0.;
      vb.calpar.nfall[i] = 0.; vb.tmp_vegs.c[i] = 0.; vb.tmp_vegs.strn[i] = 0.;
      vb.del_v2soi.ltrfalc[i] = 0.; vb.del_a2v.innpp[i] = 0.; vb.del_a2v.ingpp[i] = 0.;
      vb.del_v2a.rm[i] = 0.; vb.bgcpar.c2neven[i] = 30.; vb.bgcpar.cpart[i] = 0.;
    }
    bd.m_vegs.labn = 0.; bd.m_soi2v.lnuptake = 0.3;   // total actual uptake 0.3 g N/m2/month
    vb.tmp_vegs.labn = 0.;
    vb.calpar.nmax = 10.0; cd.m_vegd.ffoliage[0] = 1.0;  // potential non-vascular uptake 10
  }
  double sum_nextract() const { double s = 0; for (int il = 0; il < MAX_SOI_LAY; ++il) s += bd.m_soi2v.nextract[il]; return s; }
};

static void test_n3_uptake_cap() {
  VegFix off(false), on(true);
  off.vb.deltanfeed(); on.vb.deltanfeed();
  CHECK(off.vb.del_soi2v.innuptake == 10.0);               // legacy: soil AVLN ignored
  CHECK_NEAR(on.vb.del_soi2v.innuptake, 0.95 * (0.1 + 0.2), 1e-12);  // moss+fibric AVLN only
  // no AVLN/water in the top layers -> second pass over the whole column
  VegFix dry(true);
  dry.ed.m_sois.liq[0] = 0.; dry.bd.m_sois.avln[1] = 0.;
  dry.vb.deltanfeed();
  CHECK_NEAR(dry.vb.del_soi2v.innuptake, 0.95 * (0.2 * 0 + 0.1 * 0 + 0.3 + 0.4), 1e-12);
  // vascular PFT unaffected by the switch
  VegFix voff(false, 0), von(true, 0);
  for (VegFix* f : {&voff, &von}) { f->cd.m_soil.frootfrac[2][0] = 1.0; f->bd.m_vegd.raq10 = 1.0;
    f->vb.bgcpar.knuptake = 1.0; f->vb.deltanfeed(); }
  CHECK(voff.vb.del_soi2v.innuptake == von.vb.del_soi2v.innuptake);
}

static void test_n3_extraction() {
  VegFix off(false), on(true);
  off.vb.afterIntegration(); on.vb.afterIntegration();
  // off: non-vascular uptake is not taken from any layer (legacy, N created from nothing)
  for (int il = 0; il < MAX_SOI_LAY; ++il) CHECK(off.bd.m_soi2v.nextract[il] == 0.0 || il >= 4);
  CHECK(off.sum_nextract() == 0.0 || off.bd.m_soi2v.nextract[4] == -99.0);
  double s_off = 0; for (int il = 0; il < 4; ++il) s_off += off.bd.m_soi2v.nextract[il];
  CHECK(s_off == 0.0);
  // on: whole uptake extracted, from moss+fibric, weighted by AVLN x liq (0.1:0.2)
  double s_on = 0; for (int il = 0; il < 4; ++il) s_on += on.bd.m_soi2v.nextract[il];
  CHECK_NEAR(s_on, 0.3, 1e-12);
  CHECK_NEAR(on.bd.m_soi2v.nextract[0], 0.1, 1e-12);
  CHECK_NEAR(on.bd.m_soi2v.nextract[1], 0.2, 1e-12);
  CHECK(on.bd.m_soi2v.nextract[2] == 0.0 && on.bd.m_soi2v.nextract[3] == 0.0);
  // fallback: top layers dry -> weights over the whole column (0.3*10 : 0.4*10)
  VegFix dry(true); dry.ed.m_sois.liq[0] = 0.; dry.ed.m_sois.liq[1] = 0.;
  dry.vb.afterIntegration();
  CHECK(dry.bd.m_soi2v.nextract[0] == 0.0 && dry.bd.m_soi2v.nextract[1] == 0.0);
  CHECK_NEAR(dry.bd.m_soi2v.nextract[2], 0.3 * 3.0 / 7.0, 1e-12);
  CHECK_NEAR(dry.bd.m_soi2v.nextract[3], 0.3 * 4.0 / 7.0, 1e-12);
  // vascular PFT: switch has no effect
  VegFix voff(false, 0), von(true, 0);
  for (VegFix* f : {&voff, &von}) { f->vb.fracnuptake[2] = 1.0; f->vb.afterIntegration(); }
  bool same = true;
  for (int il = 0; il < 4; ++il) same = same && voff.bd.m_soi2v.nextract[il] == von.bd.m_soi2v.nextract[il];
  CHECK(same); CHECK_NEAR(von.bd.m_soi2v.nextract[2], 0.3, 1e-12);
}

// ---------------------------------------------------------------- Soil_Bgc fixture
// Column as in test_rts_units.cpp burial test (fibric | humic | mineral | mineral), no
// decomposition, no litter, no inputs -> the only N flux is the root/moss extraction.
struct SoilFix {
  CohortData cd; EnvData ed; BgcData bd; Soil_Bgc sb;
  int n;
  SoilFix(int nl = 4) : n(nl) {
    sb.setCohortData(&cd); sb.setEnvData(&ed); sb.setBgcData(&bd);
    cd.m_soil.numsl = n;
    const int type[4] = {1, 2, 3, 3}; const double dz[4] = {0.05, 0.10, 0.3, 0.5};
    double z = 0.0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      cd.m_soil.type[i] = i < 4 ? type[i] : 3; cd.m_soil.dz[i] = i < 4 ? dz[i] : 0.5;
      cd.m_soil.z[i] = z; z += cd.m_soil.dz[i];
      cd.m_soil.frootfrac[i][0] = 0.0;
      ed.m_sois.ts[i] = 5.0; ed.m_sois.liq[i] = 0.0; ed.m_soid.sws[i] = 0.5;
      bd.m_soid.knmoist[i] = 0.0;
      sb.bgcpar.kdrawc[i] = 0.; sb.bgcpar.kdsoma[i] = 0.; sb.bgcpar.kdsompr[i] = 0.; sb.bgcpar.kdsomcr[i] = 0.;
      sb.del_soi2a.rhrawc[i] = 0.; sb.del_soi2a.rhsoma[i] = 0.; sb.del_soi2a.rhsompr[i] = 0.; sb.del_soi2a.rhsomcr[i] = 0.;
      sb.ltrflc[i] = 0.; sb.ltrfln[i] = 0.; sb.rtnextract[i] = 0.;
      sb.del_soi2soi.netnmin[i] = 0.; sb.del_soi2soi.nimmob[i] = 0.;
      bd.m_v2soi.rtlfalfrac[i] = 0.; bd.m_soi2v.nextract[i] = 0.;
      bd.m_sois.rawc[i] = sb.tmp_sois.rawc[i] = 100.0; bd.m_sois.soma[i] = sb.tmp_sois.soma[i] = 50.0;
      bd.m_sois.sompr[i] = sb.tmp_sois.sompr[i] = 80.0; bd.m_sois.somcr[i] = sb.tmp_sois.somcr[i] = 120.0;
      bd.m_sois.orgn[i] = sb.tmp_sois.orgn[i] = 10.0; bd.m_sois.avln[i] = sb.tmp_sois.avln[i] = 1.0;
    }
    cd.m_soil.deepthick = 0.10;
    sb.bgcpar.som2co2 = 0.5; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.bgcpar.moistmin = 0.; sb.bgcpar.moistmax = 1.; sb.bgcpar.moistopt = 0.5; sb.bgcpar.rhq10 = 2.;
    sb.bgcpar.rhmoistfrozen = 0.; sb.bgcpar.kn2 = 0.; sb.bgcpar.nmincnsoil = 20.; sb.calpar.micbnup = 0.;
    sb.bgcpar.fnloss = 0.;
    sb.set_nfeed(1); sb.set_avlnflg(0); sb.set_baseline(0);
    sb.d2wdebrisc = 0.; sb.d2wdebrisn = 0.; sb.totdzliq = 0.; sb.totdzavln = 0.; sb.del_soi2l.avlnlost = 0.;
    sb.tmp_sois.wdebrisc = 0.; sb.tmp_sois.wdebrisn = 0.; sb.del_soi2a.rhwdeb = 0.;
    bd.m_sois.wdebrisc = 0.; bd.m_sois.wdebrisn = 0.;
    bd.m_a2soi.avlninput = 0.; bd.m_a2soi.orgninput = 0.;
  }
  void extract(int il, double r) { sb.rtnextract[il] = r; bd.m_soi2v.nextract[il] = r; }
  double dN() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.avln[i] + sb.del_sois.orgn[i]; return s; }
  double dA() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.avln[i]; return s; }
  double N() const { double s = 0; for (int i = 0; i < n; ++i) s += bd.m_sois.avln[i] + bd.m_sois.orgn[i]; return s; }
};

static void test_n3_soil_n_decreases_by_uptake() {
  // feed the N3 extraction of the veg fixture into the soil budget
  for (int sw = 0; sw < 2; ++sw) {
    VegFix v(sw == 1); v.vb.afterIntegration();
    SoilFix s;
    for (int il = 0; il < 4; ++il) s.extract(il, v.bd.m_soi2v.nextract[il]);
    s.sb.deltastate();
    if (sw == 0) {
      CHECK(s.dN() == 0.0);   // off: soil N untouched by the moss uptake
      CHECK(s.dA() == 0.0);
    } else {
      CHECK_NEAR(s.dN(), -0.3, 1e-12);   // on: soil loses exactly the non-vascular uptake
      CHECK_NEAR(s.sb.del_sois.avln[0], -0.1, 1e-12);
      CHECK_NEAR(s.sb.del_sois.avln[1], -0.2, 1e-12);
    }
  }
}

// ---------------------------------------------------------------- G1
static void test_g1_limited_extraction_shape() {
  const double tau = 0.2;
  CHECK(Soil_Bgc::limitedExtraction(1.0, 0.5, tau) == 1.0);    // A >= tau r: unlimited (== r exactly)
  CHECK(Soil_Bgc::limitedExtraction(1.0, 0.2, tau) == 1.0);    // boundary A = tau r
  CHECK_NEAR(Soil_Bgc::limitedExtraction(1.0, 0.1, tau), 0.5, 1e-12);  // A/tau
  CHECK(Soil_Bgc::limitedExtraction(1.0, 0.0, tau) == 0.0);
  CHECK(Soil_Bgc::limitedExtraction(1.0, -0.1, tau) == 0.0);   // negative AVLN: no AVLN draw
  CHECK(Soil_Bgc::limitedExtraction(0.0, 0.3, tau) == 0.0);
  CHECK(Soil_Bgc::limitedExtraction(0.7, 1e3, tau) == 0.7);
  // continuous and monotone in A
  double prev = -1; bool mono = true;
  for (int k = 0; k <= 100; ++k) { double v = Soil_Bgc::limitedExtraction(1.0, 0.003 * k, tau); mono = mono && v >= prev; prev = v; }
  CHECK(mono);
  // the point of G1: dA/dt = -rx(A) keeps A >= 0 over a month where the zero-order sink
  // (rx = r) would overdraw it (explicit Euler, 1000 sub-steps)
  double a_lim = 0.1, a_raw = 0.1; const double r = 0.5, h = 1e-3;
  for (int k = 0; k < 1000; ++k) { a_lim -= h * Soil_Bgc::limitedExtraction(r, a_lim, tau); a_raw -= h * r; }
  CHECK(a_lim > 0.0); CHECK(a_raw < 0.0);
}

static void test_g1_deltastate_split() {
  // r = 1 from layer 1 (humic), AVLN 0.1 < tau r: 0.5 from AVLN, 0.5 from ORGN
  SoilFix off, on; on.sb.n_extract_limiter = true;
  for (SoilFix* s : {&off, &on}) { s->sb.tmp_sois.avln[1] = 0.1; s->extract(1, 1.0); s->sb.deltan(); s->sb.deltastate(); }
  CHECK(off.sb.del_sois.avln[1] == -1.0); CHECK(off.sb.del_sois.orgn[1] == 0.0);   // legacy
  CHECK_NEAR(on.sb.del_sois.avln[1], -0.5, 1e-12); CHECK_NEAR(on.sb.del_sois.orgn[1], -0.5, 1e-12);
  CHECK(off.dN() == -1.0); CHECK_NEAR(on.dN(), -1.0, 1e-12);   // total N removed unchanged
  CHECK(off.sb.totnextract == 1.0); CHECK_NEAR(on.sb.totnextract, 0.5, 1e-12);  // AVLN budget
  // plenty of AVLN: limiter is an exact no-op
  SoilFix off2, on2; on2.sb.n_extract_limiter = true;
  for (SoilFix* s : {&off2, &on2}) { s->extract(0, 0.1); s->extract(2, 0.15); s->sb.deltan(); s->sb.deltastate(); }
  bool same = off2.sb.totnextract == on2.sb.totnextract;
  for (int i = 0; i < 4; ++i) same = same && off2.sb.del_sois.avln[i] == on2.sb.del_sois.avln[i]
                                            && off2.sb.del_sois.orgn[i] == on2.sb.del_sois.orgn[i];
  CHECK(same);
  // tau = 1 (whole month): stronger limit
  SoilFix t1; t1.sb.n_extract_limiter = true; t1.sb.n_extract_tau = 1.0;
  t1.sb.tmp_sois.avln[1] = 0.1; t1.extract(1, 1.0); t1.sb.deltan(); t1.sb.deltastate();
  CHECK_NEAR(t1.sb.del_sois.avln[1], -0.1, 1e-12); CHECK_NEAR(t1.sb.del_sois.orgn[1], -0.9, 1e-12);
}

// ---------------------------------------------------------------- give-up counters
// Forced give-up: one rooted layer whose AVLN (0.1) is overdrawn by a zero-order root
// extraction of 0.5 g/month (the Gydan CMT70 failure mode). The first integration hits the
// negative-AVLN pool check from t ~ 0.2 on, forces maxitmon steps and blackholes; the retry
// (negative AVLN tolerated) completes.
struct IntFix {
  SoilFix s; Integrator in;
  IntFix(bool g1) : s(2) {
    s.sb.n_extract_limiter = g1;
    s.bd.m_sois.avln[0] = 0.1; s.extract(0, 0.5);
    in.setBgcData(&s.bd); in.setSoil_Bgc(&s.sb);
  }
};

static void test_giveup_counter_forced_case() {
  IntFix f(false);
  const double n0 = f.s.N();
  f.in.updateMonthlySbgc(2);
  CHECK(f.in.sbgc_months == 1);
  CHECK(f.in.sbgc_giveups == 1);
  CHECK(f.in.sbgc_retry_fails == 0);
  CHECK(f.in.n_avln_retries == 1);
  CHECK(f.in.sbgc_last == 1);
  CHECK(f.in.sbgc_forced_steps >= f.in.maxitmon);     // blackhole = maxitmon forced steps
  CHECK(f.in.sbgc_rkf_calls > f.in.sbgc_forced_steps);
  CHECK(f.in.sbgc_last_bh_time > 0.15f && f.in.sbgc_last_bh_time < 0.25f);  // AVLN gone at t = 0.2
  CHECK(f.in.sbgc_last_bh_time_retry == 1.0f);
  CHECK_NEAR(f.s.N(), n0 - 0.5, 1e-4);                // retry: whole extraction, N conserved
  // a second (identical) month accumulates
  f.s.bd.m_sois.avln[0] = 0.1; f.s.bd.m_sois.orgn[0] = 10.0;
  f.in.updateMonthlySbgc(2);
  CHECK(f.in.sbgc_months == 2 && f.in.sbgc_giveups == 2 && f.in.sbgc_retry_fails == 0);
  // G1 on: the same month integrates without give-up, same N removed
  IntFix g(true);
  g.in.updateMonthlySbgc(2);
  CHECK(g.in.sbgc_months == 1 && g.in.sbgc_giveups == 0 && g.in.sbgc_retry_fails == 0);
  CHECK(g.in.sbgc_last == 0 && g.in.sbgc_forced_steps == 0);
  CHECK(g.s.bd.m_sois.avln[0] > 0.0);
  CHECK_NEAR(g.s.N(), n0 - 0.5, 1e-4);
}

static void test_giveup_counter_retry_fails() {
  // negative litter C into layer 1 drives RAWC < 0 (a pool the retry does not tolerate):
  // first integration and retry both blackhole -> remainder path
  IntFix f(false);
  f.s.bd.m_sois.rawc[1] = 1.0; f.s.sb.ltrflc[1] = -5.0;
  f.in.updateMonthlySbgc(2);
  CHECK(f.in.sbgc_months == 1);
  CHECK(f.in.sbgc_giveups == 1);
  CHECK(f.in.sbgc_retry_fails == 1);
  CHECK(f.in.sbgc_last == 2);
  CHECK(f.in.sbgc_last_bh_time < 1.0f && f.in.sbgc_last_bh_time_retry < 1.0f);
  CHECK(f.in.sbgc_forced_steps >= 2 * f.in.maxitmon);
}

static void test_giveup_counter_clean_month() {
  IntFix f(false); f.s.extract(0, 0.01);   // AVLN 0.1 covers the month
  f.in.updateMonthlySbgc(2);
  CHECK(f.in.sbgc_months == 1 && f.in.sbgc_giveups == 0 && f.in.sbgc_retry_fails == 0);
  CHECK(f.in.sbgc_last == 0 && f.in.n_avln_retries == 0 && f.in.sbgc_forced_steps == 0);
  CHECK(f.in.sbgc_rkf_calls >= 1);
}

int main() {
  test_switch_config_defaults();
  test_switch_config_parse();
  test_n3_uptake_cap();
  test_n3_extraction();
  test_n3_soil_n_decreases_by_uptake();
  test_g1_limited_extraction_shape();
  test_g1_deltastate_split();
  test_giveup_counter_forced_case();
  test_giveup_counter_retry_fails();
  test_giveup_counter_clean_month();
  std::printf("Daily-port N switch unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail ? 1 : 0;
}
