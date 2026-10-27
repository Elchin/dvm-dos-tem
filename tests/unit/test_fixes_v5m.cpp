// v5m: switchable fix L1 (default OFF = legacy).
//   L1 model_settings.n_safetynet_conserve -> Integrator::n_safetynet_conserve
// When the soil integrator gives up (adapt() resets the soil to month fraction
// t0) the safety net delivers the remaining inputs but legacy does not apply
// the remaining (1 - t0) of the root N extraction that the vegetation has
// already taken up for the whole month -> N is created. With L1 on, that N is
// debited from soil ORGN (rooted layers by extraction share, then any layer),
// then AVLN, and any rest is carried as debt to the next soil month.
// Refs (L1): McGuire et al. 1992, doi:10.1029/92GB00219; Nasholm et al. 1998,
//   doi:10.1038/31921; Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>

namespace {
struct L1Guard {              // static switch: always restore legacy after a test
  explicit L1Guard(int v) { Integrator::n_safetynet_conserve = v; }
  ~L1Guard() { Integrator::n_safetynet_conserve = 0; }
};
struct SnFix {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb; Integrator in;
  explicit SnFix(int numsl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    in.setBgcData(&bd); in.setSoil_Bgc(&sb);
    cd.m_soil.numsl = numsl;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      cd.m_soil.type[i] = (i == 0) ? 0 : 2;
      sb.ltrflc[i] = 0.; sb.ltrfln[i] = 0.; sb.rtnextract[i] = 0.;
      bd.m_v2soi.rtlfalfrac[i] = 0.;
    }
    sb.d2wdebrisc = 0.; sb.d2wdebrisn = 0.;
    bd.m_a2soi.orgninput = 0.; bd.m_a2soi.avlninput = 0.;
    sb.nfeed = 1; sb.n_rootless_fix = 0;
    for (int i = 0; i < NUMEQ; ++i) in.y[i] = 0.f;
  }
  double soilN() const {
    double s = in.y[Integrator::I_WDEBRISN];
    for (int il = 0; il < MAX_SOI_LAY; ++il) s += in.y[Integrator::I_L_ORGN + il] + in.y[Integrator::I_L_AVLN + il];
    return s;
  }
  // typical event column: 4 layers, roots in 1-3, extraction > AVLN
  void event_column() {
    double ex[4] = {0.0, 0.30, 0.15, 0.05}, rf[4] = {0.0, 0.5, 0.3, 0.2};
    double on[4] = {5.0, 40.0, 120.0, 300.0}, av[4] = {0.0, 0.02, 0.01, 0.005};
    for (int i = 0; i < 4; ++i) {
      sb.rtnextract[i] = ex[i]; bd.m_v2soi.rtlfalfrac[i] = rf[i];
      in.y[Integrator::I_L_ORGN + i] = on[i]; in.y[Integrator::I_L_AVLN + i] = av[i];
      sb.ltrfln[i] = 0.01 * i; sb.ltrflc[i] = 0.5 * i;
    }
    sb.d2wdebrisn = 0.02; sb.d2wdebrisc = 1.0; bd.m_a2soi.orgninput = 0.003; bd.m_a2soi.avlninput = 0.004;
  }
};
Json::Value ms(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
}

TEST_SUITE("fixes_v5m") {

TEST_CASE("L1 default 0 (legacy) in Integrator and ModelData(); ModelData(Json) parses, invalid exits") {
  CHECK(Integrator::n_safetynet_conserve == 0);
  ModelData md; CHECK(md.n_safetynet_conserve == 0);
  ModelData m0(ms("svp_ice_fix", 0)); CHECK(m0.n_safetynet_conserve == 0);
  ModelData m1(ms("n_safetynet_conserve", 1)); CHECK(m1.n_safetynet_conserve == 1);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_safetynet_conserve", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_safetynet_conserve", -1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_safetynet_conserve", 1)); }) == 0);
  Integrator in;
  CHECK(in.safetynet_ndebit == 0.0); CHECK(in.safetynet_nneed == 0.0);
  CHECK(in.n_debt == 0.0); CHECK(in.n_debit_total == 0.0);
}

TEST_CASE("L1 off: legacy creates (1 - t0) * root extraction; state identical to legacy formula") {
  L1Guard g(0);
  SnFix f(4); f.event_column(); f.in.bh_t0 = 0.3f;
  float y0[NUMEQ]; for (int i = 0; i < NUMEQ; ++i) y0[i] = f.in.y[i];
  double n0 = f.soilN();
  f.in.applySoilSafetyNet();
  // legacy: only inputs are added (no debit), soil N rises by exactly the inputs
  CHECK(f.soilN() - n0 == doctest::Approx(f.in.safetynet_n).epsilon(1e-5));
  const double frac = 0.7;
  double exp_in = frac * (0.01 + 0.02 + 0.03 + 0.02 + 0.003 + 0.0); // avlninput not delivered: rtlfalfrac[0]=0
  CHECK(f.in.safetynet_n == doctest::Approx(exp_in).epsilon(1e-6));
  for (int il = 0; il < 4; ++il) {   // ORGN: + litter N only; AVLN unchanged; nothing debited
    double add = (il > 0 ? frac * 0.01 * il : 0.0) + (il == 0 ? frac * 0.003 : 0.0);
    CHECK(f.in.y[Integrator::I_L_ORGN + il] == doctest::Approx(y0[Integrator::I_L_ORGN + il] + add).epsilon(1e-6));
    CHECK(f.in.y[Integrator::I_L_AVLN + il] == y0[Integrator::I_L_AVLN + il]);
  }
  CHECK(f.in.safetynet_ndebit == 0.0); CHECK(f.in.safetynet_nneed == 0.0); CHECK(f.in.n_debt == 0.0);
  // the legacy imbalance (veg took the full month, soil debited only t0)
  double leak = frac * (0.30 + 0.15 + 0.05);
  CHECK(leak == doctest::Approx(0.35));
}

TEST_CASE("L1 on: soil N change = inputs - (1 - t0) * root extraction (veg + soil N conserved)") {
  L1Guard g(1);
  SnFix f(4); f.event_column(); f.in.bh_t0 = 0.3f;
  float y0[NUMEQ]; for (int i = 0; i < NUMEQ; ++i) y0[i] = f.in.y[i];
  double n0 = f.soilN();
  f.in.applySoilSafetyNet();
  const double frac = 0.7, E = 0.30 + 0.15 + 0.05;
  CHECK(f.in.safetynet_nneed == doctest::Approx(frac * E).epsilon(1e-6));
  // donor pools are float (ORGN ~300 g/m2 -> ulp ~3e-5): debit exact to float round-off
  CHECK(f.in.safetynet_ndebit == doctest::Approx(frac * E).epsilon(1e-4));
  CHECK(f.in.n_debt == 0.0);
  // vegetation already holds E (whole month); the soil lost t0*E during the
  // integrated part and must lose frac*E now: total N change = inputs only.
  CHECK(f.soilN() - n0 == doctest::Approx(f.in.safetynet_n - frac * E).epsilon(1e-5));
  // donors: ORGN of the rooted layers in proportion to extraction; AVLN untouched
  for (int il = 1; il < 4; ++il) {
    double add = frac * 0.01 * il;
    double deb = frac * f.sb.rtnextract[il];
    CHECK(f.in.y[Integrator::I_L_ORGN + il] == doctest::Approx(y0[Integrator::I_L_ORGN + il] + add - deb).epsilon(1e-6));
    CHECK(f.in.y[Integrator::I_L_AVLN + il] == y0[Integrator::I_L_AVLN + il]);
  }
  CHECK(f.in.y[Integrator::I_L_ORGN + 0] == doctest::Approx(y0[Integrator::I_L_ORGN] + frac * 0.003).epsilon(1e-6));
  // C delivery identical to legacy (C has no leak)
  L1Guard g0(0);
  SnFix h(4); h.event_column(); h.in.bh_t0 = 0.3f; h.in.applySoilSafetyNet();
  CHECK(h.in.safetynet_c == f.in.safetynet_c);
  for (int il = 0; il < MAX_SOI_LAY; ++il) CHECK(h.in.y[Integrator::I_L_RAWC + il] == f.in.y[Integrator::I_L_RAWC + il]);
  CHECK(h.in.y[Integrator::I_WDEBRISC] == f.in.y[Integrator::I_WDEBRISC]);
}

TEST_CASE("L1 on: rootless layers count only with n_rootless_fix 2 (as in deltastate); deposition too") {
  L1Guard g(1);
  SnFix f(3);
  f.sb.rtnextract[0] = 0.2; f.sb.rtnextract[1] = 0.1; f.bd.m_v2soi.rtlfalfrac[1] = 1.0;  // layer 0 rootless
  for (int il = 0; il < 3; ++il) f.in.y[Integrator::I_L_ORGN + il] = 10.f;
  f.bd.m_a2soi.avlninput = 0.05; f.in.bh_t0 = 0.5f;
  f.in.applySoilSafetyNet();
  CHECK(f.in.safetynet_nneed == doctest::Approx(0.05));        // 0.5 * 0.1, layer 0 excluded
  CHECK(f.in.safetynet_n == doctest::Approx(0.0));             // no deposit into rootless top layer
  SnFix h(3); h.sb.n_rootless_fix = 2;
  h.sb.rtnextract[0] = 0.2; h.sb.rtnextract[1] = 0.1; h.bd.m_v2soi.rtlfalfrac[1] = 1.0;
  for (int il = 0; il < 3; ++il) h.in.y[Integrator::I_L_ORGN + il] = 10.f;
  h.bd.m_a2soi.avlninput = 0.05; h.in.bh_t0 = 0.5f;
  double n0 = h.soilN();
  h.in.applySoilSafetyNet();
  CHECK(h.in.safetynet_nneed == doctest::Approx(0.15));
  CHECK(h.in.safetynet_n == doctest::Approx(0.025));           // deposit delivered (deltastate does)
  CHECK(h.soilN() - n0 == doctest::Approx(0.025 - 0.15).epsilon(1e-5));
  // legacy (off) with n_rootless_fix 2: no deposit (unchanged legacy behaviour)
  L1Guard g0(0);
  SnFix k(3); k.sb.n_rootless_fix = 2; k.bd.m_a2soi.avlninput = 0.05; k.in.bh_t0 = 0.5f;
  k.in.applySoilSafetyNet();
  CHECK(k.in.safetynet_n == 0.0);
}

TEST_CASE("L1 on, edge: all soil N pools zero -> nothing debited now, debt carried and paid next month") {
  L1Guard g(1);
  SnFix f(3);
  f.sb.rtnextract[1] = 0.4; f.bd.m_v2soi.rtlfalfrac[1] = 1.0; f.in.bh_t0 = 0.75f;
  f.in.applySoilSafetyNet();
  CHECK(f.in.safetynet_nneed == doctest::Approx(0.1));
  CHECK(f.in.safetynet_ndebit == 0.0);
  CHECK(f.in.n_debt == doctest::Approx(0.1));
  CHECK(f.soilN() == 0.0);
  // next month: mineralisation/litter refilled the pools; debt paid ORGN first
  f.in.y[Integrator::I_L_ORGN + 1] = 0.06f; f.in.y[Integrator::I_L_ORGN + 2] = 0.02f;
  f.in.y[Integrator::I_L_AVLN + 1] = 0.05f;
  double n0 = f.soilN();
  f.in.payNDebt();
  CHECK(n0 - f.soilN() == doctest::Approx(0.1).epsilon(1e-5));
  CHECK(f.in.y[Integrator::I_L_ORGN + 1] == 0.0f); CHECK(f.in.y[Integrator::I_L_ORGN + 2] == 0.0f);
  CHECK(f.in.y[Integrator::I_L_AVLN + 1] == doctest::Approx(0.03).epsilon(1e-5));
  CHECK(f.in.n_debt == 0.0);
  CHECK(f.in.n_debit_total == doctest::Approx(0.1).epsilon(1e-5));
}

TEST_CASE("L1 on, edge: all soil N pools negative -> pools untouched (never made more negative), debt carried") {
  L1Guard g(1);
  SnFix f(3);
  for (int il = 0; il < 3; ++il) { f.in.y[Integrator::I_L_ORGN + il] = -0.5f; f.in.y[Integrator::I_L_AVLN + il] = -0.1f; }
  f.sb.rtnextract[1] = 0.2; f.sb.rtnextract[2] = 0.2;
  f.bd.m_v2soi.rtlfalfrac[1] = 0.5; f.bd.m_v2soi.rtlfalfrac[2] = 0.5; f.in.bh_t0 = 0.0f;
  f.in.applySoilSafetyNet();
  for (int il = 0; il < 3; ++il) {
    CHECK(f.in.y[Integrator::I_L_ORGN + il] == -0.5f); CHECK(f.in.y[Integrator::I_L_AVLN + il] == -0.1f);
  }
  CHECK(f.in.safetynet_ndebit == 0.0);
  CHECK(f.in.n_debt == doctest::Approx(0.4));
  f.in.payNDebt();                       // still nothing to take
  CHECK(f.in.n_debt == doctest::Approx(0.4));
}

TEST_CASE("L1 on, edge: ORGN short -> AVLN donates the rest; t0 = 1 or no extraction -> no debit") {
  L1Guard g(1);
  SnFix f(2);
  f.sb.rtnextract[0] = 1.0; f.bd.m_v2soi.rtlfalfrac[0] = 1.0; f.cd.m_soil.type[0] = 2;
  f.in.y[Integrator::I_L_ORGN] = 0.3f; f.in.y[Integrator::I_L_ORGN + 1] = 0.2f; f.in.y[Integrator::I_L_AVLN + 1] = 1.0f;
  f.in.bh_t0 = 0.0f;
  double n0 = f.soilN();
  f.in.applySoilSafetyNet();
  CHECK(f.in.safetynet_ndebit == doctest::Approx(1.0).epsilon(1e-6));
  CHECK(n0 - f.soilN() == doctest::Approx(1.0).epsilon(1e-6));
  CHECK(f.in.y[Integrator::I_L_ORGN] == 0.0f); CHECK(f.in.y[Integrator::I_L_ORGN + 1] == 0.0f);
  CHECK(f.in.y[Integrator::I_L_AVLN + 1] == doctest::Approx(0.5).epsilon(1e-6));
  SnFix h(2); h.sb.rtnextract[0] = 1.0; h.bd.m_v2soi.rtlfalfrac[0] = 1.0; h.in.y[Integrator::I_L_ORGN] = 5.f;
  h.in.bh_t0 = 1.0f; h.in.applySoilSafetyNet();
  CHECK_FALSE(h.in.safetynet_fired); CHECK(h.in.y[Integrator::I_L_ORGN] == 5.f);
  SnFix k(2); k.in.y[Integrator::I_L_ORGN] = 5.f; k.in.bh_t0 = 0.5f; k.in.applySoilSafetyNet();
  CHECK(k.in.safetynet_nneed == 0.0); CHECK(k.in.y[Integrator::I_L_ORGN] == 5.f);
  SnFix m(2); m.sb.nfeed = 0; m.sb.rtnextract[0] = 1.0; m.bd.m_v2soi.rtlfalfrac[0] = 1.0;   // nfeed off: no N bookkeeping
  m.in.y[Integrator::I_L_ORGN] = 5.f; m.in.bh_t0 = 0.5f; m.in.applySoilSafetyNet();
  CHECK(m.in.y[Integrator::I_L_ORGN] == 5.f); CHECK(m.in.safetynet_nneed == 0.0);
}

TEST_CASE("L1 debitProportional: weights, caps, remainder, round-off, degenerate inputs") {
  float p[4] = {1.0f, 2.0f, 0.0f, -1.0f};
  double w[4] = {1.0, 1.0, 5.0, 5.0};             // layers 2, 3 have no positive pool -> ignored
  double d = Integrator::debitProportional(1.0, p, 4, w);
  CHECK(d == doctest::Approx(1.0).epsilon(1e-7));
  CHECK(p[0] == doctest::Approx(0.5f)); CHECK(p[1] == doctest::Approx(1.5f));
  CHECK(p[2] == 0.0f); CHECK(p[3] == -1.0f);
  float q[2] = {0.1f, 10.0f}; double wq[2] = {1.0, 0.0};
  d = Integrator::debitProportional(0.5, q, 2, wq);  // cap at 0.1, rest from remaining pool
  CHECK(d == doctest::Approx(0.5).epsilon(1e-6));
  CHECK(q[0] == 0.0f); CHECK(q[1] == doctest::Approx(9.6f));
  float r[2] = {0.2f, 0.3f};
  d = Integrator::debitProportional(2.0, r, 2, NULL);   // more than available -> empty, returns available
  CHECK(d == doctest::Approx(0.5).epsilon(1e-6)); CHECK(r[0] == 0.0f); CHECK(r[1] == 0.0f);
  float s[1] = {3.0f};
  CHECK(Integrator::debitProportional(0.0, s, 1, NULL) == 0.0);
  CHECK(Integrator::debitProportional(-1.0, s, 1, NULL) == 0.0);
  CHECK(Integrator::debitProportional(1.0, s, 0, NULL) == 0.0);
  CHECK(s[0] == 3.0f);
  // returned amount equals the actual float change (bookkeeping exact to round-off)
  float t[3] = {123.456f, 0.001f, 77.7f}; double tw[3] = {0.3, 0.2, 0.1};
  double before = (double)t[0] + t[1] + t[2];
  d = Integrator::debitProportional(0.0123, t, 3, tw);
  CHECK(before - ((double)t[0] + t[1] + t[2]) == doctest::Approx(d).epsilon(1e-12));
  CHECK(d == doctest::Approx(0.0123).epsilon(1e-4));
}

} // TEST_SUITE fixes_v5m
