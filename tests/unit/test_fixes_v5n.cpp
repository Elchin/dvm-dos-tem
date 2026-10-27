// v5n draft: soil-integrator give-up fixes (all default OFF = legacy).
//   G1 model_settings.n_extract_limiter (+ n_extract_tau) -> Soil_Bgc::n_extract_limiter / n_extract_tau
//   G2 model_settings.soil_giveup_retry -> Integrator::soil_giveup_retry
//   diagnostics model_settings.soil_giveup_diag -> Integrator::soil_giveup_diag (stderr/csv only)
// Mechanism: the root N extraction rtnextract is constant within the soil month
// (zero-order sink on AVLN). Immobilisation and leaching vanish as AVLN -> 0 but
// extraction does not, so AVLN reaches zero in finite time; RKF45 rejects every
// step that makes AVLN negative and adapt() gives up after 100 steps at dt_min.
// Refs (G1/G2): Raich et al. 1991, doi:10.2307/1941899; McGuire et al. 1992,
//   doi:10.1029/92GB00219; Thornton et al. 2007, doi:10.1029/2006GB002868;
//   Nasholm et al. 1998, doi:10.1038/31921; Schimel & Bennett 2004, doi:10.1890/03-8002;
//   Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011; Burchard et al. 2003,
//   doi:10.1016/S0168-9274(03)00101-6.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
struct G1Guard {   // static switches: always restore legacy after a test
  G1Guard(int on, double tau) { Soil_Bgc::n_extract_limiter = on; Soil_Bgc::n_extract_tau = tau; }
  ~G1Guard() { Soil_Bgc::n_extract_limiter = 0; Soil_Bgc::n_extract_tau = 0.05; }
};
// same set-up as test_fixes_v5k.cpp SoilNFix (deltastate bookkeeping, no leaching)
struct SoilNFix5n {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  explicit SoilNFix5n(double avln) {
    const int types[6] = {1, 1, 2, 2, 3, 3}; const double dz[6] = {0.05, 0.05, 0.05, 0.05, 0.1, 0.1};
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.n_som_move_fix = 1; sb.n_rootless_fix = 0;
    cd.m_soil.numsl = 6; cd.m_soil.deepthick = 0.1; sb.totdzliq = 0.0;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < 6;
      cd.m_soil.type[i] = on ? types[i] : 3; cd.m_soil.dz[i] = on ? dz[i] : 0.1; cd.m_soil.z[i] = z; z += cd.m_soil.dz[i];
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = avln;
      sb.del_soi2a.rhrawc[i] = 3.0; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = on ? 2.0 : 0.0; sb.ltrfln[i] = on ? 0.1 : 0.0; sb.rtnextract[i] = on ? 0.2 : 0.0;
      sb.del_soi2soi.netnmin[i] = on ? 0.05 : 0.0; bd.m_v2soi.rtlfalfrac[i] = on ? 0.2 : 0.0;
      ed.m_sois.liq[i] = 10.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0; sb.del_soi2l.avlnlost = 0; ed.m_sois.draindepth = 0;
  }
  double dN() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double expected() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.ltrfln[i] - sb.rtnextract[i]; return s; }
};
Json::Value ms(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
}

TEST_SUITE("fixes_v5n") {

TEST_CASE("v5n defaults are legacy (0) in the statics and ModelData(); ModelData(Json) parses, invalid exits") {
  CHECK(Soil_Bgc::n_extract_limiter == 0); CHECK(Soil_Bgc::n_extract_tau == 0.05);
  CHECK(Integrator::soil_giveup_retry == 0); CHECK(Integrator::soil_giveup_diag == 0);
  ModelData md;
  CHECK(md.n_extract_limiter == 0); CHECK(md.soil_giveup_retry == 0); CHECK(md.soil_giveup_diag == 0); CHECK(md.n_extract_tau == 0.05);
  ModelData m0(ms("svp_ice_fix", 0));
  CHECK(m0.n_extract_limiter == 0); CHECK(m0.soil_giveup_retry == 0); CHECK(m0.soil_giveup_diag == 0); CHECK(m0.n_extract_tau == 0.05);
  CHECK(ModelData(ms("n_extract_limiter", 1)).n_extract_limiter == 1);
  CHECK(ModelData(ms("soil_giveup_retry", 1)).soil_giveup_retry == 1);
  CHECK(ModelData(ms("soil_giveup_diag", 1)).soil_giveup_diag == 1);
  CHECK(ModelData(ms("n_extract_tau", 0.2)).n_extract_tau == 0.2);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_extract_limiter", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("soil_giveup_retry", -1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("soil_giveup_diag", 3)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_extract_tau", 0.0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_extract_tau", 1.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("n_extract_tau", 1.0)); }) == 0);
  Integrator in;
  CHECK(in.retry_count == 0); CHECK(in.retry_fail_count == 0); CHECK(!in.retry_fired);
  CHECK(in.rkf45_calls == 0); CHECK(in.forced_months == 0);
  Soil_Bgc sb; CHECK(!sb.nx_limiter_active); CHECK(!sb.nxLimiterOn());
}

TEST_CASE("G1 limitedExtraction: r unchanged while AVLN >= tau*r, linear in AVLN below, 0 at AVLN <= 0") {
  const double r = 0.4, tau = 0.05;                       // tau*r = 0.02 g/m2
  CHECK(Soil_Bgc::limitedExtraction(r, 1.0, tau) == r);    // exactly r (bit-identical branch)
  CHECK(Soil_Bgc::limitedExtraction(r, 0.02, tau) == doctest::Approx(r).epsilon(1e-14));   // continuous at the knee
  CHECK(Soil_Bgc::limitedExtraction(r, 0.01, tau) == doctest::Approx(0.2));
  CHECK(Soil_Bgc::limitedExtraction(r, 0.0, tau) == 0.0);
  CHECK(Soil_Bgc::limitedExtraction(r, -1e-9, tau) == 0.0);
  CHECK(Soil_Bgc::limitedExtraction(0.0, 0.0, tau) == 0.0);    // no extraction
  CHECK(Soil_Bgc::limitedExtraction(-0.1, 0.0, tau) == -0.1);  // non-positive r passed through
  double prev = -1;
  for (int k = 0; k <= 100; ++k) {                          // monotone, bounded by r and by A/tau
    double A = 0.0005 * k, x = Soil_Bgc::limitedExtraction(r, A, tau);
    CHECK(x >= prev); CHECK(x <= r); CHECK(x <= A / tau + 1e-15); prev = x;
  }
}

TEST_CASE("G1 mechanism: zero-order extraction drives AVLN negative within the month, the limiter keeps it >= 0") {
  // dA/dt = -r (legacy, zero order) vs dA/dt = -limitedExtraction(r, A) (G1), A0 < r*1 month
  const double r = 0.4, A0 = 0.1, tau = 0.05; const int n = 20000; const double h = 1.0 / n;
  double Al = A0, Ag = A0, tzero = -1, taken = 0;
  for (int k = 0; k < n; ++k) {
    Al -= r * h; if (Al < 0 && tzero < 0) tzero = (k + 1) * h;
    double x = Soil_Bgc::limitedExtraction(r, Ag, tau); Ag -= x * h; taken += x * h;
    CHECK(Ag >= 0.0);
  }
  CHECK(tzero == doctest::Approx(A0 / r).epsilon(1e-3));   // legacy: AVLN = 0 at t = A0/r = 0.25 month
  CHECK(Ag > 0.0); CHECK(Ag < 1e-3);                       // G1: exponential approach, time scale tau
  CHECK(taken == doctest::Approx(A0 - Ag).epsilon(1e-9));  // AVLN part; ORGN supplies r - taken
}

TEST_CASE("G1 off / on with ample AVLN: deltastate bit-identical to legacy") {
  SoilNFix5n a(5.0); a.sb.deltastate();
  G1Guard g(1, 0.05);
  SoilNFix5n b(5.0); b.sb.deltastate();                    // AVLN 5 >= tau*r = 0.01 in every layer
  for (int i = 0; i < 6; ++i) {
    CHECK(a.sb.del_sois.avln[i] == b.sb.del_sois.avln[i]);
    CHECK(a.sb.del_sois.orgn[i] == b.sb.del_sois.orgn[i]);
  }
}

TEST_CASE("G1 on, AVLN short: AVLN supplies r*AVLN/(tau*r), ORGN the rest; soil N change unchanged") {
  const double avl = 0.004, tau = 0.05;                    // tau*r = 0.01 > AVLN -> rx = 0.08
  SoilNFix5n a(avl); a.sb.deltastate();                    // legacy
  G1Guard g(1, tau);
  SoilNFix5n b(avl); b.sb.deltastate();
  const double rx = 0.2 * avl / (tau * 0.2);
  for (int i = 0; i < 6; ++i) {
    CAPTURE(i);
    CHECK(b.sb.del_sois.avln[i] == doctest::Approx(a.sb.del_sois.avln[i] + (0.2 - rx)).epsilon(1e-12));
    CHECK(b.sb.del_sois.orgn[i] == doctest::Approx(a.sb.del_sois.orgn[i] - (0.2 - rx)).epsilon(1e-12));
  }
  CHECK(b.dN() == doctest::Approx(a.dN()).epsilon(1e-12));
  CHECK(b.dN() == doctest::Approx(b.expected()).epsilon(1e-9));   // N conserved (= litter N - extraction)
}

TEST_CASE("G1 applies only where deltastate applies extraction (rooted, or rootless with n_rootless_fix 2)") {
  G1Guard g(1, 0.05);
  for (int fx = 0; fx <= 2; fx += 2) {
    CAPTURE(fx);
    SoilNFix5n a(0.001), b(0.001);
    a.sb.n_rootless_fix = fx; b.sb.n_rootless_fix = fx;
    a.bd.m_v2soi.rtlfalfrac[5] = 0.0; b.bd.m_v2soi.rtlfalfrac[5] = 0.0;   // layer 5 rootless
    Soil_Bgc::n_extract_limiter = 0; a.sb.deltastate();
    Soil_Bgc::n_extract_limiter = 1; b.sb.deltastate();
    if (fx == 0) {   // legacy rootless branch: AVLN change 0, extraction ignored -> limiter not applied
      CHECK(b.sb.del_sois.avln[5] == a.sb.del_sois.avln[5]); CHECK(b.sb.del_sois.orgn[5] == a.sb.del_sois.orgn[5]);
    } else {
      CHECK(b.sb.del_sois.orgn[5] < a.sb.del_sois.orgn[5]);
    }
    CHECK(b.dN() == doctest::Approx(a.dN()).epsilon(1e-12));
  }
}

TEST_CASE("G2 per-month flag: nx_limiter_active turns the limiter on with the switch off, and only then") {
  G1Guard g(0, 0.05);
  SoilNFix5n a(0.004); a.sb.deltastate();
  SoilNFix5n b(0.004); b.sb.nx_limiter_active = true; CHECK(b.sb.nxLimiterOn()); b.sb.deltastate();
  Soil_Bgc::n_extract_limiter = 1;
  SoilNFix5n c(0.004); c.sb.deltastate();
  for (int i = 0; i < 6; ++i) {
    CHECK(b.sb.del_sois.avln[i] == c.sb.del_sois.avln[i]); CHECK(b.sb.del_sois.orgn[i] == c.sb.del_sois.orgn[i]);
    CHECK(b.sb.del_sois.avln[i] != a.sb.del_sois.avln[i]);
  }
}

} // TEST_SUITE
