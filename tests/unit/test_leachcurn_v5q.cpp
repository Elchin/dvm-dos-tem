// v5q-leach-curn: N drainage leaching from the current drainage-zone AVLN,
// model_settings.n_leach_curn_fix (0 = legacy, default; 1 = fix; requires n_leach_fix = 1).
// Legacy: the drainage part of avlnlost uses totdzavln built once per month from the month-start
// AVLN, i.e. a constant-rate sink during the month that keeps drawing on an emptied layer (soil
// AVLN < 0 -> adapt() give-up -> safety net; V5Q_STEP8_D2M_GATE.md section 13). 1 = the rate uses
// the current (tmp_sois) AVLN of the same zone and goes to 0 with it. See include/Soil_Bgc.h.
// Refs: a sink that does not vanish with a non-negative state drives it below zero in finite time;
//   positivity needs a state-dependent sink: Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011;
//   Burchard et al. 2003, doi:10.1016/S0168-9274(03)00101-6. Leaching context (as n_leach_fix):
//   McNamara et al. 2008, doi:10.1002/hyp.6920; Khosh et al. 2017, doi:10.1002/2016jg003377.
//   The fix itself: no published source (code defect).
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>
#include <algorithm>

namespace {
Json::Value lcm(Json::Value v, int nl = 1) {
  Json::Value c; c["model_settings"]["n_leach_fix"] = nl; c["model_settings"]["n_leach_curn_fix"] = v; return c;
}

// Column as test_leach_v5q LeachSoil. deep_ts > 0 opens net mineralization in layers >= 1.
struct CurSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  CurSoil(int n_, int sw, bool set_switch = true) : n(n_) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2;
    if (set_switch) { sb.n_leach_curn_fix = sw; }
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.bgcpar.fnloss = 0.1; sb.bgcpar.kn2 = 1.0; sb.bgcpar.nmincnsoil = 1.0; sb.calpar.micbnup = 0.5;
    sb.decay = 0.0;
    cd.m_soil.numsl = n; cd.m_soil.deepthick = 0.1;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < n;
      cd.m_soil.type[i] = 3; cd.m_soil.dz[i] = 0.1; cd.m_soil.z[i] = z; z += 0.1;
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = on ? 0.4 : 0.0; bd.m_sois.avln[i] = on ? 0.4 : 0.0;
      sb.del_soi2a.rhrawc[i] = 3.0; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = 0; sb.ltrfln[i] = 0; sb.rtnextract[i] = on ? 0.002 : 0.0; bd.m_soi2v.nextract[i] = sb.rtnextract[i];
      bd.m_v2soi.rtlfalfrac[i] = on ? 0.5 : 0.0; bd.m_soid.knmoist[i] = 0.5;
      ed.m_sois.liq[i] = on ? 20.0 : 0.0; ed.m_sois.ts[i] = -1.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0.01; sb.del_soi2l.avlnlost = 0;
    ed.m_sois.draindepth = 0.25;
    ed.m_soi2l.qdrain = 30.0; ed.m_soi2l.qover = 5.0;
    totals();
  }
  // Section-13 geometry: 1-cm frozen moss layer 0 (no SOM) holds the whole drainage zone, tiny
  // liquid, big snowmelt runoff; deeper layers mineralize N (raises the column-wide cap).
  void moss0() {
    cd.m_soil.dz[0] = 0.01; double z = 0.01;
    for (int i = 1; i < MAX_SOI_LAY; ++i) { cd.m_soil.z[i] = z; z += cd.m_soil.dz[i]; }
    sb.tmp_sois.rawc[0] = sb.tmp_sois.soma[0] = sb.tmp_sois.sompr[0] = sb.tmp_sois.somcr[0] = 0; sb.tmp_sois.orgn[0] = 0;
    sb.del_soi2a.rhrawc[0] = sb.del_soi2a.rhsoma[0] = sb.del_soi2a.rhsompr[0] = sb.del_soi2a.rhsomcr[0] = 0;
    sb.rtnextract[0] = 0; bd.m_soi2v.nextract[0] = 0; bd.m_v2soi.rtlfalfrac[0] = 0;
    sb.tmp_sois.avln[0] = 0.0244; bd.m_sois.avln[0] = 0.0244;
    ed.m_sois.liq[0] = 0.66; ed.m_sois.ts[0] = -5.0;
    for (int i = 1; i < n; ++i) { ed.m_sois.ts[i] = 1.0; sb.del_soi2a.rhrawc[i] = 30.0; }
    sb.bgcpar.nmincnsoil = 10.0;
    bd.m_a2soi.avlninput = 0.0;
    ed.m_sois.draindepth = 0.0097; ed.m_soi2l.qover = 121.0; ed.m_soi2l.qdrain = 3.41;
    totals();
  }
  void totals() {   // as Soil_Bgc::prepareIntegration() (month-start AVLN)
    sb.totdzliq = 0; sb.totdzavln = 0;
    for (int i = 0; i < cd.m_soil.numsl; i++) {
      double f = 0;
      if (cd.m_soil.z[i] + cd.m_soil.dz[i] <= ed.m_sois.draindepth) f = 1;
      else if (cd.m_soil.z[i] < ed.m_sois.draindepth) f = (ed.m_sois.draindepth - cd.m_soil.z[i]) / cd.m_soil.dz[i];
      else break;
      sb.totdzliq += fmax(0., ed.m_sois.liq[i]) * f; sb.totdzavln += fmax(0., bd.m_sois.avln[i]) * f;
    }
  }
  void step() { sb.deltan(); sb.deltastate(); }
  double removed() {   // AVLN removed by leaching in this derivative evaluation
    double s = 0;
    for (int i = 0; i < n; ++i)
      s += (i == 0 ? bd.m_a2soi.avlninput : 0.0) + sb.del_soi2soi.netnmin[i] - sb.rtnextract[i] - sb.del_sois.avln[i];
    return s;
  }
  double zone_now() { return Soil_Bgc::drainZoneAvln(n, cd.m_soil.z, cd.m_soil.dz, sb.tmp_sois.avln, ed.m_sois.draindepth); }
  // explicit-Euler month (nsteps), returns {min AVLN[0], integrated NLOST, pool loss, budget residual}
  struct M { double min0, nlost, pool_loss, resid; };
  M month(int nsteps) {
    const double dt = 1.0 / nsteps;
    double av0 = 0, av1 = 0, inp = 0, ext = 0, nmin = 0, nlost = 0, mn = 1e300;
    for (int i = 0; i < n; ++i) av0 += sb.tmp_sois.avln[i];
    for (int s = 0; s < nsteps; ++s) {
      step();
      nlost += sb.del_soi2l.avlnlost * dt; inp += bd.m_a2soi.avlninput * dt;
      for (int i = 0; i < n; ++i) { ext += sb.rtnextract[i] * dt; nmin += sb.del_soi2soi.netnmin[i] * dt;
                                    sb.tmp_sois.avln[i] += sb.del_sois.avln[i] * dt; }
      mn = std::min(mn, sb.tmp_sois.avln[0]);
    }
    for (int i = 0; i < n; ++i) av1 += sb.tmp_sois.avln[i];
    return M{mn, nlost, av0 + inp + nmin - ext - av1, (av1 - av0) - (inp + nmin - ext - nlost)};
  }
};
}

TEST_SUITE("leachcurn_v5q") {

TEST_CASE("n_leach_curn_fix default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc") {
  ModelData md; CHECK(md.n_leach_curn_fix == 0);
  Json::Value c; c["model_settings"]["n_leach_fix"] = 1; ModelData mj(c); CHECK(mj.n_leach_curn_fix == 0);
  Soil_Bgc sb; CHECK(sb.n_leach_curn_fix == 0);
}

TEST_CASE("n_leach_curn_fix validation: 0/1 integers only, 1 needs n_leach_fix 1, anything else FATAL") {
  CHECK(ModelData(lcm(0)).n_leach_curn_fix == 0);
  CHECK(ModelData(lcm(1)).n_leach_curn_fix == 1);
  CHECK(ModelData(lcm(Json::Value(0), 0)).n_leach_curn_fix == 0);
  CHECK(ut::exit_status_of([] { ModelData m(lcm(Json::Value(1), 0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lcm(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lcm(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lcm(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lcm(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lcm("1")); }) == EXIT_FAILURE);
}

TEST_CASE("drainZoneAvln == prepareIntegration totdzavln (full + partial layers, negatives as 0, zone inside layer 0)") {
  CurSoil f(5, 1);
  CHECK(Soil_Bgc::drainZoneAvln(5, f.cd.m_soil.z, f.cd.m_soil.dz, f.bd.m_sois.avln, 0.25) == f.sb.totdzavln);
  double z[4] = {0, .1, .2, .3}, dz[4] = {.1, .1, .1, .1}, av[4] = {-1.0, 0.3, 0.2, 9.0};
  CHECK(Soil_Bgc::drainZoneAvln(4, z, dz, av, 0.25) == doctest::Approx(0.0 + 0.3 + 0.2 * 0.5).epsilon(1e-15));
  double av2[4] = {0.0244, 1, 1, 1}, z2[4] = {0, .01, .05, .1}, dz2[4] = {.01, .04, .05, .1};
  CHECK(Soil_Bgc::drainZoneAvln(4, z2, dz2, av2, 0.0097) == doctest::Approx(0.0244 * 0.97).epsilon(1e-14));
  CHECK(Soil_Bgc::drainZoneAvln(4, z2, dz2, av2, 0.0) == 0.0);
}

TEST_CASE("switch 0 is the legacy code path, bit for bit (unset == 0; drain part from month-start totdzavln)") {
  CurSoil a(5, 0, false), b(5, 0, true);
  for (CurSoil *f : {&a, &b}) { for (int i = 0; i < 5; ++i) f->sb.tmp_sois.avln[i] = 0.1 + 0.05 * i; }   // state != month start
  a.step(); b.step();
  CHECK(a.sb.del_soi2l.avlnlost == b.sb.del_soi2l.avlnlost);
  CHECK(a.sb.leach_drain == b.sb.leach_drain); CHECK(a.sb.leach_over == b.sb.leach_over);
  for (int i = 0; i < 5; ++i) { CHECK(a.sb.del_sois.avln[i] == b.sb.del_sois.avln[i]); }
  const double over = a.sb.tmp_sois.avln[0] / a.ed.m_sois.liq[0] * 0.1 * 5.0;
  CHECK(a.sb.del_soi2l.avlnlost - over == doctest::Approx(a.sb.totdzavln / a.sb.totdzliq * 0.1 * 30.0).epsilon(1e-13));
  // switch 1 on the same state: drain part from the current zone AVLN
  CurSoil c(5, 1); for (int i = 0; i < 5; ++i) c.sb.tmp_sois.avln[i] = 0.1 + 0.05 * i; c.step();
  CHECK(c.sb.del_soi2l.avlnlost - over == doctest::Approx(c.zone_now() / c.sb.totdzliq * 0.1 * 30.0).epsilon(1e-13));
  CHECK(c.sb.del_soi2l.avlnlost != doctest::Approx(a.sb.del_soi2l.avlnlost));
  // state == month start: switch 1 == switch 0
  CurSoil d0(5, 0), d1(5, 1); d0.step(); d1.step();
  CHECK(d1.sb.del_soi2l.avlnlost == doctest::Approx(d0.sb.del_soi2l.avlnlost).epsilon(1e-15));
}

TEST_CASE("switch 1: leaching goes smoothly to 0 as the zone empties (linear in the zone AVLN, 0 at 0); legacy does not") {
  for (int sw = 0; sw <= 1; ++sw) {
    CAPTURE(sw);
    double L1 = 0;
    for (double lam : {1.0, 0.5, 0.1, 1e-3, 1e-6, 0.0}) {
      CAPTURE(lam);
      CurSoil f(5, sw); f.bd.m_a2soi.avlninput = 0;
      for (int i = 0; i < 5; ++i) { f.sb.rtnextract[i] = 0; f.bd.m_soi2v.nextract[i] = 0; }
      for (int i = 0; i < 3; ++i) f.sb.tmp_sois.avln[i] = 0.4 * lam;      // zone layers 0-2 (layer 2 half)
      f.step();
      const double L = f.removed();
      CHECK(L == doctest::Approx(f.sb.del_soi2l.avlnlost).epsilon(1e-12));
      if (lam == 1.0) { L1 = L; REQUIRE(L1 > 0); continue; }
      if (sw == 1) { CHECK(L == doctest::Approx(lam * L1).epsilon(1e-12)); }
      else         { CHECK(L > 0.9 * f.sb.totdzavln / f.sb.totdzliq * 0.1 * 30.0); }   // month-start drain stays
      if (sw == 1 && lam == 0.0) { CHECK(L == 0.0); for (int i = 0; i < 3; ++i) CHECK(f.sb.del_sois.avln[i] == 0.0); }
    }
  }
}

TEST_CASE("section-13 case (frozen moss layer 0 = drainage zone, cap raised by deep netnmin): legacy goes negative, switch 1 never exceeds the available N") {
  for (int sw = 0; sw <= 1; ++sw) {
    CAPTURE(sw);
    CurSoil f(5, sw); f.moss0();
    f.sb.deltan(); REQUIRE(f.sb.totnetnmin > f.sb.tmp_sois.avln[0]);   // cap > layer-0 stock
    CurSoil g(5, sw); g.moss0();
    const double a0 = g.sb.tmp_sois.avln[0];
    CurSoil::M m = g.month(20000);
    if (sw == 0) { CHECK(m.min0 < 0.0); }                              // legacy: the safety-net trigger
    else {
      CHECK(m.min0 >= 0.0);
      CHECK(m.nlost <= a0 * (1.0 + 1e-12));                            // all leaching is from layer 0
      CHECK(m.nlost > 0.5 * a0);                                       // and it does leach
    }
  }
}

TEST_CASE("switch 1 N closure over a month: AVLN lost by the pools == integrated NLOST; column AVLN budget closes") {
  for (int geo = 0; geo <= 1; ++geo) {
    CAPTURE(geo);
    CurSoil f(5, 1); if (geo == 1) f.moss0();
    CurSoil::M m = f.month(3000);
    CHECK(m.nlost > 0.0);
    CHECK(m.pool_loss == doctest::Approx(m.nlost).epsilon(1e-11));
    CHECK(std::fabs(m.resid) < 1e-12);
  }
}

} // TEST_SUITE leachcurn_v5q
