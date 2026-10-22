// Area 11: soil integrator safety net and N-conserving root N extraction.
// Area 12b (NLOST) is checked at run level by the regression harness (N closure).
#include "doctest.h"
#include "test_helpers.h"

namespace {
struct BgcFixture {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb; Integrator in;
  explicit BgcFixture(int numsl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    cd.yrsdist = 1000; cd.fri = 100;     // no fire re-deposition branch in prepareIntegration
    fd.fire_a2soi.orgn = 0.0;
    in.setBgcData(&bd); in.setSoil_Bgc(&sb);
    cd.m_soil.numsl = numsl;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      cd.m_soil.type[i] = 0; cd.m_soil.z[i] = 0.1 * i; cd.m_soil.dz[i] = 0.1;
      sb.ltrflc[i] = 0.; sb.ltrfln[i] = 0.; bd.m_sois.avln[i] = 0.; bd.m_soi2v.nextract[i] = 0.;
      bd.m_v2soi.rtlfalfrac[i] = 0.; ed.m_soid.sws[i] = 0.5; ed.m_sois.liq[i] = 10.;
    }
    for (int i = 0; i < NUMEQ; ++i) in.y[i] = 0.f;
    ed.m_sois.draindepth = 0.3;
  }
};
}

TEST_SUITE("integrator") {

TEST_CASE("soil safety net delivers the remaining (1 - t0) of the month's inputs; C and N conserved") {
  BgcFixture f(4);
  // moss (type 0), fibric, humic, mineral
  f.cd.m_soil.type[0] = 0; f.cd.m_soil.type[1] = 1; f.cd.m_soil.type[2] = 2; f.cd.m_soil.type[3] = 3;
  double lc[4] = {5.0, 20.0, 6.0, 2.0}, ln[4] = {0.1, 0.4, 0.12, 0.04};
  for (int i = 0; i < 4; ++i) { f.sb.ltrflc[i] = lc[i]; f.sb.ltrfln[i] = ln[i]; }
  f.sb.d2wdebrisc = 3.0; f.sb.d2wdebrisn = 0.03;
  f.bd.m_a2soi.orgninput = 0.02; f.bd.m_a2soi.avlninput = 0.05; f.bd.m_v2soi.rtlfalfrac[0] = 0.1;
  f.sb.nfeed = 1;
  float y0[NUMEQ]; for (int i = 0; i < NUMEQ; ++i) { f.in.y[i] = 1.0f; y0[i] = 1.0f; }
  f.in.bh_t0 = 0.4f;
  f.in.applySoilSafetyNet();
  const double frac = 0.6;
  double expC = frac * (20.0 + 6.0 + 2.0 + 3.0);                 // moss layer excluded
  double expN = frac * (0.4 + 0.12 + 0.04 + 0.03 + 0.02 + 0.05);
  CHECK(f.in.safetynet_fired);
  CHECK(f.in.safetynet_count == 1);
  CHECK(f.in.safetynet_t0 == doctest::Approx(0.4));
  CHECK(f.in.safetynet_c == doctest::Approx(expC).epsilon(1e-6));
  CHECK(f.in.safetynet_n == doctest::Approx(expN).epsilon(1e-6));
  // state increments equal the reported deliveries
  double dC = 0, dN = 0;
  for (int il = 0; il < MAX_SOI_LAY; ++il) {
    dC += (f.in.y[Integrator::I_L_RAWC + il] - y0[Integrator::I_L_RAWC + il]);
    dN += (f.in.y[Integrator::I_L_ORGN + il] - y0[Integrator::I_L_ORGN + il]) + (f.in.y[Integrator::I_L_AVLN + il] - y0[Integrator::I_L_AVLN + il]);
  }
  dC += f.in.y[Integrator::I_WDEBRISC] - y0[Integrator::I_WDEBRISC];
  dN += f.in.y[Integrator::I_WDEBRISN] - y0[Integrator::I_WDEBRISN];
  CHECK(dC == doctest::Approx(f.in.safetynet_c).epsilon(1e-5));
  CHECK(dN == doctest::Approx(f.in.safetynet_n).epsilon(1e-5));
  CHECK(f.in.y[Integrator::I_L_RAWC + 0] == doctest::Approx(1.0));   // moss untouched
  // SOM pools other than RAWC untouched
  for (int il = 0; il < MAX_SOI_LAY; ++il) { CHECK(f.in.y[Integrator::I_L_SOMA + il] == 1.0f); CHECK(f.in.y[Integrator::I_L_SOMCR + il] == 1.0f); }
}

TEST_CASE("soil safety net: nfeed off delivers C only; t0 = 1 delivers nothing") {
  BgcFixture f(2);
  f.cd.m_soil.type[0] = 1; f.cd.m_soil.type[1] = 3;
  f.sb.ltrflc[0] = 10; f.sb.ltrfln[0] = 1; f.sb.d2wdebrisc = 2; f.sb.d2wdebrisn = 0.5;
  f.sb.nfeed = 0; f.in.bh_t0 = 0.0f;
  f.in.applySoilSafetyNet();
  CHECK(f.in.safetynet_c == doctest::Approx(12.0));
  CHECK(f.in.safetynet_n == 0.0);
  BgcFixture g(2); g.cd.m_soil.type[0] = 1; g.sb.ltrflc[0] = 10; g.sb.nfeed = 1; g.in.bh_t0 = 1.0f;
  g.in.applySoilSafetyNet();
  CHECK_FALSE(g.in.safetynet_fired);
  CHECK(g.in.safetynet_count == 0);
}

TEST_CASE("Integrator counters start at zero") {
  Integrator in;
  CHECK(in.blackhole_count == 0); CHECK(in.safetynet_count == 0); CHECK_FALSE(in.safetynet_fired);
}

} // TEST_SUITE integrator

TEST_SUITE("root_n") {

static void run_extraction(BgcFixture& f) { f.sb.prepareIntegration(true, true, false); }

TEST_CASE("root N extraction: per-layer excess moved to layers with headroom; total = vegetation uptake") {
  BgcFixture f(4);
  for (int i = 0; i < 4; ++i) f.cd.m_soil.type[i] = 3;
  double ex[4] = {0.5, 0.1, 0.1, 0.0}, av[4] = {0.2, 0.3, 0.5, 0.2};
  for (int i = 0; i < 4; ++i) { f.bd.m_soi2v.nextract[i] = ex[i]; f.bd.m_sois.avln[i] = av[i]; }
  run_extraction(f);
  double tot = 0;
  for (int i = 0; i < 4; ++i) { tot += f.sb.rtnextract[i]; CHECK(f.sb.rtnextract[i] <= av[i] + 1e-12); CHECK(f.sb.rtnextract[i] >= 0.0); }
  CHECK(tot == doctest::Approx(0.7));
  // hand calculation: excess 0.3 shared by headroom (0, 0.2, 0.4, 0.2)
  CHECK(f.sb.rtnextract[0] == doctest::Approx(0.2));
  CHECK(f.sb.rtnextract[1] == doctest::Approx(0.1 + 0.3 * 0.2 / 0.8));
  CHECK(f.sb.rtnextract[2] == doctest::Approx(0.1 + 0.3 * 0.4 / 0.8));
  CHECK(f.sb.rtnextract[3] == doctest::Approx(0.0 + 0.3 * 0.2 / 0.8));
}

TEST_CASE("root N extraction: column demand > column avln -> vegetation extraction unchanged (no cap)") {
  BgcFixture f(3);
  for (int i = 0; i < 3; ++i) f.cd.m_soil.type[i] = 3;
  double ex[3] = {0.5, 0.4, 0.3}, av[3] = {0.1, 0.2, 0.1};
  for (int i = 0; i < 3; ++i) { f.bd.m_soi2v.nextract[i] = ex[i]; f.bd.m_sois.avln[i] = av[i]; }
  run_extraction(f);
  for (int i = 0; i < 3; ++i) CHECK(f.sb.rtnextract[i] == doctest::Approx(ex[i]));
}

TEST_CASE("root N extraction: no excess -> unchanged; negative avln treated as 0") {
  BgcFixture f(3);
  for (int i = 0; i < 3; ++i) f.cd.m_soil.type[i] = 3;
  double ex[3] = {0.05, 0.1, 0.0}, av[3] = {0.2, 0.3, -0.1};
  for (int i = 0; i < 3; ++i) { f.bd.m_soi2v.nextract[i] = ex[i]; f.bd.m_sois.avln[i] = av[i]; }
  run_extraction(f);
  for (int i = 0; i < 3; ++i) CHECK(f.sb.rtnextract[i] == doctest::Approx(ex[i]));
  // property test: random cases conserve the total whenever demand <= supply
  unsigned s = 12345u;
  auto rnd = [&]() { s = s * 1103515245u + 12345u; return ((s >> 8) & 0xffff) / 65536.0; };
  for (int k = 0; k < 200; ++k) {
    std::unique_ptr<BgcFixture> gp(new BgcFixture(6)); BgcFixture& g = *gp;   // large: keep off the stack
    double te = 0, ta = 0;
    for (int i = 0; i < 6; ++i) { g.cd.m_soil.type[i] = 3; g.bd.m_soi2v.nextract[i] = rnd() * 0.3; g.bd.m_sois.avln[i] = rnd() * 0.5; te += g.bd.m_soi2v.nextract[i]; ta += g.bd.m_sois.avln[i]; }
    run_extraction(g);
    double tot = 0; for (int i = 0; i < 6; ++i) tot += g.sb.rtnextract[i];
    CHECK(tot == doctest::Approx(te).epsilon(1e-12));
    if (te <= ta) for (int i = 0; i < 6; ++i) CHECK(g.sb.rtnextract[i] <= g.bd.m_sois.avln[i] + 1e-12);
  }
}

} // TEST_SUITE root_n
