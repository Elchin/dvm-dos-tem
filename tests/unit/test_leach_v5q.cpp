// v5q-ndrain: N leaching bookkeeping fix, model_settings.n_leach_fix (0 = legacy, default; 1 = fix).
// Legacy: deltan() books avlnlost (-> I_AVLNLOSS -> NLOST) = drain part + overland part, but
// deltastate() removes sum(ndrain) = avlnlost*qdrain (qdrain mm/month) from the drainage zone,
// nothing for the overland part and nothing when totdzliq <= 0.01. n_leach_fix 1: the pools lose
// exactly avlnlost (overland part from layer 0, drainage part by current AVLN in the zone).
// Refs: no published source (code bug, N mass conservation; same code in upstream dvm-dos-tem
//   c0fa75cb src/Soil_Bgc.cpp l.1004-1018). See include/Soil_Bgc.h, docs/daily-mode-REFERENCES.md.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value msl(Json::Value v) { Json::Value c; c["model_settings"]["n_leach_fix"] = v; return c; }

// Soil column for deltan() + deltastate() (cf. test_bnf_v5q BnfSoil). ts < 0 -> netnmin = 0, all
// layers rooted, no extraction limiter: d(AVLN_il)/dt = ninput_il - ndrain_il - rtnextract_il.
struct LeachSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  explicit LeachSoil(int n, int fix, double fnloss = 0.1) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.avlnflg = 1; sb.n_leach_fix = fix; sb.n_rootless_fix = 2;
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.bgcpar.fnloss = fnloss; sb.bgcpar.kn2 = 1.0; sb.bgcpar.nmincnsoil = 1.0; sb.calpar.micbnup = 0.5;
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
    ed.m_sois.draindepth = 0.25;                     // 2 full layers + half of layer 2
    ed.m_soi2l.qdrain = 30.0; ed.m_soi2l.qover = 5.0;  // mm/month
    totals();
  }
  void totals() {   // as Soil_Bgc::prepareIntegration()
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
  // AVLN removed by leaching in this derivative evaluation
  double removed() {
    double s = 0;
    for (int i = 0; i < cd.m_soil.numsl; ++i)
      s += (i == 0 ? bd.m_a2soi.avlninput : 0.0) + sb.del_soi2soi.netnmin[i] - sb.rtnextract[i] - sb.del_sois.avln[i];
    return s;
  }
};
}

TEST_SUITE("leach_v5q") {

TEST_CASE("n_leach_fix default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc") {
  ModelData md; CHECK(md.n_leach_fix == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.n_leach_fix == 0);
  Soil_Bgc sb; CHECK(sb.n_leach_fix == 0);
}

TEST_CASE("n_leach_fix validation: 0/1 integers only, anything else FATAL") {
  CHECK(ModelData(msl(0)).n_leach_fix == 0);
  CHECK(ModelData(msl(1)).n_leach_fix == 1);
  CHECK(ut::exit_status_of([] { ModelData m(msl(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msl(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msl(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msl(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msl("1")); }) == EXIT_FAILURE);
}

TEST_CASE("leachLayers: AVLN-weighted in the drainage zone (partial layer by fdz), overland on layer 0") {
  double z[4] = {0, .1, .2, .3}, dz[4] = {.1, .1, .1, .1}, liq[4] = {5, 20, 20, 20}, av[4] = {0.0, 0.3, 0.2, 9.0}, nd[4];
  double t = Soil_Bgc::leachLayers(4, z, dz, liq, av, 0.25, 0.04, 0.01, nd);
  const double sa = 0.0 + 0.3 + 0.2 * 0.5;
  CHECK(t == doctest::Approx(0.05).epsilon(1e-14));
  CHECK(nd[0] == doctest::Approx(0.01));                 // N-free top layer: only the overland part
  CHECK(nd[1] == doctest::Approx(0.04 * 0.3 / sa));
  CHECK(nd[2] == doctest::Approx(0.04 * 0.1 / sa));
  CHECK(nd[3] == 0.0);                                    // below the drainage depth
  // fallbacks: no AVLN in the zone -> by liq*fdz; no water either -> layer 0; negative pools count 0
  double av0[4] = {-1, 0, 0, 5};
  Soil_Bgc::leachLayers(4, z, dz, liq, av0, 0.25, 0.04, 0.0, nd);
  CHECK(nd[0] == doctest::Approx(0.04 * 5 / 35.0)); CHECK(nd[2] == doctest::Approx(0.04 * 10 / 35.0)); CHECK(nd[3] == 0.0);
  double l0[4] = {0, 0, 0, 0};
  t = Soil_Bgc::leachLayers(4, z, dz, l0, av0, 0.25, 0.04, 0.01, nd);
  CHECK(nd[0] == doctest::Approx(0.05)); CHECK(t == doctest::Approx(0.05));
}

TEST_CASE("fix 1: pools lose exactly the booked avlnlost (drain + overland), also when uncapped / capped / totdzliq <= 0.01") {
  for (int c = 0; c < 4; ++c) {
    CAPTURE(c);
    LeachSoil f(5, 1);
    if (c == 1) { f.ed.m_soi2l.qover = 0; }
    if (c == 2) { f.sb.bgcpar.fnloss = 50.0; }                  // cap binds: avlnlost = totdzavln - nextract + input
    if (c == 3) { for (int i = 0; i < 5; ++i) f.ed.m_sois.liq[i] = 0.002; f.totals(); }   // totdzliq = 0.005
    f.step();
    REQUIRE(f.sb.del_soi2l.avlnlost > 0.0);
    CHECK(f.removed() == doctest::Approx(f.sb.del_soi2l.avlnlost).epsilon(1e-12));
    CHECK(f.sb.leach_drain + f.sb.leach_over == doctest::Approx(f.sb.del_soi2l.avlnlost).epsilon(1e-14));
    if (c == 2) {
      double tx = 0; for (int i = 0; i < 5; ++i) tx += f.sb.rtnextract[i];
      CHECK(f.sb.del_soi2l.avlnlost == doctest::Approx(f.sb.totdzavln - tx + 0.01));
    }
  }
}

TEST_CASE("legacy (fix 0) bookkeeping documented: removal = drain part * qdrain, overland part never removed") {
  LeachSoil f(5, 0); f.step();
  const double booked = f.sb.del_soi2l.avlnlost;
  const double over = f.sb.tmp_sois.avln[0] / f.ed.m_sois.liq[0] * 0.1 * 5.0;
  const double drain = booked - over;
  CHECK(f.removed() == doctest::Approx(booked * 30.0).epsilon(1e-12));   // sum(liq*fdz)/totdzliq = 1
  CHECK(f.removed() != doctest::Approx(booked));
  CHECK(drain == doctest::Approx(f.sb.totdzavln / f.sb.totdzliq * 0.1 * 30.0));
  LeachSoil g(5, 0); for (int i = 0; i < 5; ++i) g.ed.m_sois.liq[i] = 0.002; g.totals(); g.step();
  CHECK(g.sb.del_soi2l.avlnlost > 0.0);
  CHECK(g.removed() == doctest::Approx(0.0));                            // nothing removed, NLOST booked
}

TEST_CASE("fix 1, avlnflg 0 (budget residual): pools conserve AVLN, removal = avlnlost") {
  LeachSoil f(5, 1); f.sb.avlnflg = 0; f.step();
  CHECK(f.removed() == doctest::Approx(f.sb.del_soi2l.avlnlost).epsilon(1e-12));
  double s = 0; for (int i = 0; i < 5; ++i) s += f.sb.del_sois.avln[i];
  CHECK(s == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("monthly (30 explicit daily steps): sum of AVLN lost by the pools == integrated NLOST (fix 1), not (fix 0)") {
  for (int fix = 0; fix <= 1; ++fix) {
    CAPTURE(fix);
    LeachSoil f(5, fix);
    const double dt = 1.0 / 30.0;
    double nlost = 0, av0 = 0, av1 = 0, inp = 0, ext = 0;
    for (int i = 0; i < 5; ++i) av0 += f.sb.tmp_sois.avln[i];
    for (int d = 0; d < 30; ++d) {
      f.step();
      nlost += f.sb.del_soi2l.avlnlost * dt;                 // what I_AVLNLOSS integrates
      inp += f.bd.m_a2soi.avlninput * dt;
      for (int i = 0; i < 5; ++i) { ext += f.sb.rtnextract[i] * dt; f.sb.tmp_sois.avln[i] += f.sb.del_sois.avln[i] * dt; }
    }
    for (int i = 0; i < 5; ++i) av1 += f.sb.tmp_sois.avln[i];
    const double lost_by_pools = av0 + inp - ext - av1;
    CHECK(nlost > 0.0);
    if (fix == 1) { CHECK(lost_by_pools == doctest::Approx(nlost).epsilon(1e-12)); }
    else          { CHECK(lost_by_pools > 5.0 * nlost); }
  }
}

} // TEST_SUITE leach_v5q
