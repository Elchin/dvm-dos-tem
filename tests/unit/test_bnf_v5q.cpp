// v5q: biological N fixation (BNF), model_settings.bnf_scheme (0 = off, default).
//   1 PRO   B = c*(1 - exp(d*NPP)),  c = 1.8 gN/m2/yr, d = -0.003 m2 yr/gC
//   2 AET   B = max(0, a*ET + b),    a = 0.00234 gN/mm/m2, b = -0.0172 gN/m2/yr
//   3 CONST B = bnf_const;  + bnf_asym (free-living, constant);  taper bnf_avln_cap
// from the previous model year, monthly weight = previous-year NPP shape, into
// m_a2soi.avlninput (AVLN of the top layer; needs n_rootless_fix = 2).
// Refs: Meyerholt et al. 2016, doi:10.5194/bg-13-1491-2016 (Eqs. A1, A2, Table A1);
//   Cleveland et al. 1999, doi:10.1029/1999GB900014; Thornton et al. 2007,
//   doi:10.1029/2006GB002868. See include/Bnf.h, docs/daily-mode-REFERENCES.md.
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/Bnf.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value msb(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
Json::Value msb2(const char* k, Json::Value v, const char* k2, Json::Value v2) {
  Json::Value c; c["model_settings"][k] = v; c["model_settings"][k2] = v2; return c;
}
struct BnfSoil {   // soil column for one Soil_Bgc::deltastate() call (cf. test_fixes_v5k SoilNFix)
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  explicit BnfSoil(int n) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    cd.m_soil.numsl = n; cd.m_soil.deepthick = 0.1; sb.totdzliq = 0.0;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < n;
      cd.m_soil.type[i] = 3; cd.m_soil.dz[i] = 0.1; cd.m_soil.z[i] = z; z += 0.1;
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = 1; bd.m_sois.avln[i] = on ? 0.5 : 0.0;
      sb.del_soi2a.rhrawc[i] = 3.0; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = on ? 2.0 : 0.0; sb.ltrfln[i] = on ? 0.1 : 0.0; sb.rtnextract[i] = on ? 0.02 : 0.0;
      sb.del_soi2soi.netnmin[i] = on ? 0.05 : 0.0; bd.m_v2soi.rtlfalfrac[i] = on ? 0.5 : 0.0;
      ed.m_sois.liq[i] = 10.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0; sb.del_soi2l.avlnlost = 0; ed.m_sois.draindepth = 0;
  }
  double dN() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double expected() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.ltrfln[i] - sb.rtnextract[i]; return s; }
};
}

TEST_SUITE("bnf_v5q") {

TEST_CASE("BNF default off: ModelData() and ModelData(Json) without keys; paper defaults (Table A1)") {
  ModelData md; CHECK(md.bnf.scheme == 0);
  ModelData mj(msb("svp_ice_fix", 0)); CHECK(mj.bnf.scheme == 0);
  CHECK(mj.bnf.pro_c == 1.8); CHECK(mj.bnf.pro_d == -0.003);
  CHECK(mj.bnf.aet_a == 0.00234); CHECK(mj.bnf.aet_b == -0.0172);
  CHECK(mj.bnf.cnst == 0.0); CHECK(mj.bnf.asym == 0.0); CHECK(mj.bnf.avln_cap == 0.0);
  Soil_Bgc sb; CHECK(sb.bnfp.scheme == 0); CHECK(sb.bnf_soil_month == 0.0); CHECK(!sb.bnfmem.have_prev);
  atm2soi_bgc a; CHECK(a.bnf == 0.0);
}

TEST_CASE("BNF switch validation: 0..3 integers only; parameters; n_rootless_fix 2 required") {
  for (int s = 0; s <= 3; ++s) { CAPTURE(s); CHECK(ModelData(msb("bnf_scheme", s)).bnf.scheme == s); }
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", 4)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", -1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", 1.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", "1")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_scheme", 1)); }) == 0);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 1, "n_rootless_fix", 0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 2, "n_rootless_fix", 1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 0, "n_rootless_fix", 0)); }) == 0);
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_asym", 0.1)); }) == EXIT_FAILURE);      // needs scheme > 0
  CHECK(ut::exit_status_of([] { ModelData m(msb("bnf_const", 0.1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 1, "bnf_asym", -0.1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 1, "bnf_pro_d", 0.01)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msb2("bnf_scheme", 1, "bnf_pro_c", "x")); }) == EXIT_FAILURE);
  ModelData m3(msb2("bnf_scheme", 3, "bnf_const", 0.25)); CHECK(m3.bnf.cnst == 0.25);
  ModelData ma(msb2("bnf_scheme", 1, "bnf_asym", 0.05)); CHECK(ma.bnf.asym == 0.05);
}

TEST_CASE("BNF formulas vs Meyerholt et al. 2016 Eqs. A1/A2 (Table A1) and Cleveland et al. 1999 units") {
  bnf::Params p;
  CHECK(bnf::pro(0.0, p.pro_c, p.pro_d) == 0.0);
  CHECK(bnf::pro(-50.0, p.pro_c, p.pro_d) == 0.0);
  CHECK(bnf::pro(1000.0, p.pro_c, p.pro_d) == doctest::Approx(1.8 * (1.0 - std::exp(-3.0))));
  CHECK(bnf::pro(1.0e6, p.pro_c, p.pro_d) == doctest::Approx(1.8));               // saturates at c
  CHECK(bnf::pro(100.0, p.pro_c, p.pro_d) == doctest::Approx(0.46654).epsilon(1e-4));
  // AET: 0.00234 gN/mm/m2 * ET - 0.0172 == Cleveland's 0.0234 kgN/ha/mm * ET - 0.172 kgN/ha (1 kg/ha = 0.1 g/m2)
  for (double et : {50.0, 139.5, 500.0, 1500.0}) {
    CAPTURE(et);
    CHECK(bnf::aet(et, p.aet_a, p.aet_b) == doctest::Approx(0.1 * (0.0234 * et - 0.172)));
  }
  CHECK(bnf::aet(500.0, p.aet_a, p.aet_b) == doctest::Approx(1.1528));
  CHECK(bnf::aet(5.0, p.aet_a, p.aet_b) == 0.0);                                  // below 7.35 mm: clamped at 0
  p.scheme = 1; CHECK(bnf::scheme_annual(p, 200.0, 0.0) == doctest::Approx(1.8 * (1 - std::exp(-0.6))));
  p.scheme = 2; CHECK(bnf::scheme_annual(p, 200.0, 300.0) == doctest::Approx(0.6848));
  p.scheme = 3; p.cnst = 0.3; CHECK(bnf::scheme_annual(p, 200.0, 300.0) == 0.3);
  p.scheme = 0; CHECK(bnf::scheme_annual(p, 200.0, 300.0) == 0.0);
}

TEST_CASE("BNF month weights, taper, previous-year memory (first year 0, months sum to the annual rate)") {
  double npp[12] = {-1, -0.5, 0, 2, 10, 30, 40, 30, 8, 0, -1, -2};
  double s = 0; for (int m = 0; m < 12; ++m) s += bnf::month_weight(npp, m);
  CHECK(s == doctest::Approx(1.0));
  CHECK(bnf::month_weight(npp, 0) == 0.0);
  CHECK(bnf::month_weight(npp, 6) == doctest::Approx(40.0 / 120.0));
  double zero[12] = {0}; CHECK(bnf::month_weight(zero, 3) == doctest::Approx(1.0 / 12.0));
  CHECK(bnf::cap_factor(5.0, 0.0) == 1.0);
  CHECK(bnf::cap_factor(0.5, 2.0) == doctest::Approx(0.75));
  CHECK(bnf::cap_factor(3.0, 2.0) == 0.0);
  bnf::Params p; p.scheme = 1; p.asym = 0.02;
  bnf::YearMemory y;
  for (int m = 0; m < 12; ++m) {           // first year: nothing in memory -> 0
    CHECK(bnf::month_bnf(p, y, m, 0.0) == 0.0);
    y.record(m, npp[m], 10.0);
  }
  CHECK(y.have_prev); CHECK(y.prev_npp_yr() == doctest::Approx(115.5)); CHECK(y.prev_et_yr() == doctest::Approx(120.0));
  double tot = 0;
  for (int m = 0; m < 12; ++m) { tot += bnf::month_bnf(p, y, m, 0.0); y.record(m, 0.0, 0.0); }
  CHECK(tot == doctest::Approx(bnf::pro(115.5, 1.8, -0.003) + 0.02));
  // the year just recorded (all zero) is now the previous year: uniform weight, PRO(0) = 0, asym only
  CHECK(bnf::month_bnf(p, y, 4, 0.0) == doctest::Approx(0.02 / 12.0));
  bnf::Params off; CHECK(bnf::month_bnf(off, y, 4, 0.0) == 0.0);
}

TEST_CASE("BNF Soil_Bgc: off computes nothing; nfeed 0 gives 0; on uses the column AVLN taper") {
  BnfSoil f(3);
  f.sb.bnfBeginMonth(5, true); CHECK(f.sb.bnf_soil_month == 0.0);
  f.sb.bnfRecordMonth(11, 50.0, 100.0); CHECK(!f.sb.bnfmem.have_prev);        // off: memory untouched
  f.sb.bnfp.scheme = 3; f.sb.bnfp.cnst = 1.2;
  for (int m = 0; m < 12; ++m) f.sb.bnfRecordMonth(m, 1.0, 1.0);
  f.sb.bnfBeginMonth(0, false); CHECK(f.sb.bnf_soil_month == 0.0);
  f.sb.bnfBeginMonth(0, true);  CHECK(f.sb.bnf_soil_month == doctest::Approx(0.1));
  f.sb.bnfp.avln_cap = 3.0;   // column AVLN 1.5 -> factor 0.5
  f.sb.bnfBeginMonth(0, true);  CHECK(f.sb.bnf_soil_month == doctest::Approx(0.05));
}

TEST_CASE("BNF conservation: deltastate adds the input to AVLN of a rootless top layer (n_rootless_fix 2)") {
  BnfSoil f(3); f.sb.n_rootless_fix = 2;
  f.bd.m_v2soi.rtlfalfrac[0] = 0.0; f.sb.rtnextract[0] = 0.0;   // rootless top layer
  f.bd.m_a2soi.avlninput = 0.07;
  f.sb.deltastate();
  CHECK(f.dN() == doctest::Approx(f.expected() + 0.07).epsilon(1e-12));
  BnfSoil g(3); g.sb.n_rootless_fix = 2; g.bd.m_v2soi.rtlfalfrac[0] = 0.0; g.sb.rtnextract[0] = 0.0;
  g.sb.deltastate();
  CHECK(f.sb.del_sois.avln[0] - g.sb.del_sois.avln[0] == doctest::Approx(0.07).epsilon(1e-12));
  // rooted top layer: same
  BnfSoil r(3); r.sb.n_rootless_fix = 2; r.bd.m_a2soi.avlninput = 0.07; r.sb.deltastate();
  CHECK(r.dN() == doctest::Approx(r.expected() + 0.07).epsilon(1e-12));
}

TEST_CASE("BNF safety net: on, the remaining (1 - t0) input reaches a rootless top layer; off = legacy") {
  for (int scheme = 0; scheme <= 1; ++scheme) {
    CAPTURE(scheme);
    CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb; Integrator in;
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    in.setBgcData(&bd); in.setSoil_Bgc(&sb);
    cd.m_soil.numsl = 3;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      cd.m_soil.type[i] = 3; sb.ltrflc[i] = 0.; sb.ltrfln[i] = 0.; sb.rtnextract[i] = 0.; bd.m_v2soi.rtlfalfrac[i] = 0.;
    }
    bd.m_v2soi.rtlfalfrac[1] = 1.0;   // top layer rootless
    sb.d2wdebrisc = 0.; sb.d2wdebrisn = 0.; bd.m_a2soi.orgninput = 0.; bd.m_a2soi.avlninput = 0.08;
    sb.nfeed = 1; sb.n_rootless_fix = 2; sb.bnfp.scheme = scheme;
    for (int i = 0; i < NUMEQ; ++i) in.y[i] = 0.f;
    in.bh_t0 = 0.25f;
    in.applySoilSafetyNet();
    if (scheme == 0) { CHECK(in.safetynet_n == 0.0); CHECK(in.y[Integrator::I_L_AVLN] == 0.f); }
    else { CHECK(in.safetynet_n == doctest::Approx(0.06)); CHECK(in.y[Integrator::I_L_AVLN] == doctest::Approx(0.06)); }
  }
}

TEST_CASE("BNF output spec row: BNF exists, g/m2/time, monthly and yearly only") {
  std::ifstream f("../../config/output_spec.csv");
  REQUIRE(f.good());
  std::string line; bool found = false;
  while (std::getline(f, line)) {
    if (line.rfind("BNF,", 0) == 0) {
      found = true;
      CHECK(line.find(",g/m2/time,") != std::string::npos);
      CHECK(line.find(",,,invalid,invalid,invalid,invalid,") != std::string::npos);
    }
  }
  CHECK(found);
}

} // TEST_SUITE bnf_v5q
