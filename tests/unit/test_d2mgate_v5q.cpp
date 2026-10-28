// v5q-d2m-gate: deep(humic)->mineral SOM transfer gate, model_settings.d2m_transfer_gate
// (0 = legacy, default; 1 = only into receiver thickness above max(running max ALD, previous-year
// ALD), C and N together at the donor C:N, receiver N +=; 2 = transfer off).
// Refs: Koven et al. 2009, doi:10.1029/2009GL040150; Koven et al. 2013, doi:10.5194/bg-10-7109-2013;
//   Lawrence et al. 2019 (CLM5), doi:10.1029/2018MS001583 (transport only within the active layer,
//   gated on max(altmax, altmax_lastyear)). f-scaling of the dvmdostem transfer: no published source.
//   See include/Soil_Bgc.h, docs/daily-mode-REFERENCES.md.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value msg(Json::Value v) { Json::Value c; c["model_settings"]["d2m_transfer_gate"] = v; return c; }

// Column: 0 fibric (type 1, 0.05 m), 1-2 humic (type 2, 0.05 m each = the top 0.10 m donors),
// 3.. mineral (type 3, 0.10 m each; top at 0.15 m). Receivers (legacy): layers 3 and 4.
struct GateSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  GateSoil(int gate, int nl = 9, double mdz = 0.1, int nfeed = 1) : n(nl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = nfeed; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2; sb.d2m_transfer_gate = gate;
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.bgcpar.fnloss = 0.1; sb.bgcpar.kn2 = 1.0; sb.bgcpar.nmincnsoil = 1.0; sb.calpar.micbnup = 0.5;
    cd.m_soil.numsl = nl; cd.m_soil.deepthick = 0.1;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < nl;
      int t = (i == 0) ? 1 : (i <= 2 ? 2 : 3);
      double dz = (i <= 2) ? 0.05 : mdz;
      cd.m_soil.type[i] = t; cd.m_soil.dz[i] = dz; cd.m_soil.z[i] = z; z += dz;
      sb.tmp_sois.rawc[i] = 100 + 3 * i; sb.tmp_sois.soma[i] = 50 + i; sb.tmp_sois.sompr[i] = 200 - 5 * i; sb.tmp_sois.somcr[i] = 300 + 7 * i;
      sb.tmp_sois.orgn[i] = 20 + 0.5 * i; sb.tmp_sois.avln[i] = on ? 0.4 : 0.0; bd.m_sois.avln[i] = on ? 0.4 : 0.0;
      sb.del_soi2a.rhrawc[i] = 3.0 + 0.1 * i; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1 + 0.01 * i; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = 0; sb.ltrfln[i] = 0; sb.rtnextract[i] = on ? 0.002 : 0.0; bd.m_soi2v.nextract[i] = sb.rtnextract[i];
      bd.m_v2soi.rtlfalfrac[i] = on ? 0.5 : 0.0; bd.m_soid.knmoist[i] = 0.5;
      ed.m_sois.liq[i] = on ? 20.0 : 0.0; ed.m_sois.ts[i] = -1.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0.01; sb.del_soi2l.avlnlost = 0;
    ed.m_sois.draindepth = 0.25; ed.m_soi2l.qdrain = 30.0; ed.m_soi2l.qover = 5.0;
    ed.y_soid.ald = 1.0; cd.d2m_ald_prevyr = 0.0;
    sb.totdzliq = 0; sb.totdzavln = 0;
    for (int i = 0; i < nl; i++) {
      double f = 0;
      if (cd.m_soil.z[i] + cd.m_soil.dz[i] <= ed.m_sois.draindepth) f = 1;
      else if (cd.m_soil.z[i] < ed.m_sois.draindepth) f = (ed.m_sois.draindepth - cd.m_soil.z[i]) / cd.m_soil.dz[i];
      else break;
      sb.totdzliq += ed.m_sois.liq[i] * f; sb.totdzavln += bd.m_sois.avln[i] * f;
    }
  }
  GateSoil& ald(double runmax, double prev) { ed.y_soid.ald = runmax; cd.d2m_ald_prevyr = prev; return *this; }
  GateSoil& step() { sb.deltan(); sb.deltastate(); return *this; }
  double c(int i) const { return sb.del_sois.rawc[i] + sb.del_sois.soma[i] + sb.del_sois.sompr[i] + sb.del_sois.somcr[i]; }
  double sumc() const { double s = 0; for (int i = 0; i < n; ++i) s += c(i); return s; }
  double sumn() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.orgn[i]; return s; }
  double legacy_mobile() const {   // what legacy moves: donors 1, 2
    double s = 0; for (int i = 1; i <= 2; ++i) s += (sb.del_soi2a.rhrawc[i] + sb.del_soi2a.rhsoma[i] + sb.del_soi2a.rhsompr[i] + sb.del_soi2a.rhsomcr[i]) * 0.3 * 0.2;
    return s;
  }
};
}

TEST_SUITE("d2mgate_v5q") {

TEST_CASE("d2m_transfer_gate default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc, CohortData") {
  ModelData md; CHECK(md.d2m_transfer_gate == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.d2m_transfer_gate == 0);
  Soil_Bgc sb; CHECK(sb.d2m_transfer_gate == 0);
  CohortData cd; CHECK(cd.d2m_ald_prevyr == 0.0);
}

TEST_CASE("d2m_transfer_gate validation: integers 0/1/2 only, anything else FATAL") {
  CHECK(ModelData(msg(0)).d2m_transfer_gate == 0);
  CHECK(ModelData(msg(1)).d2m_transfer_gate == 1);
  CHECK(ModelData(msg(2)).d2m_transfer_gate == 2);
  CHECK(ut::exit_status_of([] { ModelData m(msg(3)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("d2mGateDepth: max(running max ALD, previous-year ALD), MISSING/negative -> 0") {
  CHECK(Soil_Bgc::d2mGateDepth(0.4, 0.6) == 0.6);
  CHECK(Soil_Bgc::d2mGateDepth(0.7, 0.6) == 0.7);
  CHECK(Soil_Bgc::d2mGateDepth(MISSING_D, 0.5) == 0.5);
  CHECK(Soil_Bgc::d2mGateDepth(MISSING_D, MISSING_D) == 0.0);
  CHECK(Soil_Bgc::d2mGateDepth(0.0, 0.0) == 0.0);
}

TEST_CASE("d2mGateReceivers: open / partial / at the top / ALD 0 / thin mineral column") {
  int ty[6] = {1, 2, 2, 3, 3, 3};
  double z[6] = {0, .05, .10, .15, .25, .35}, dz[6] = {.05, .05, .05, .1, .1, .1}, w[6];
  double f = Soil_Bgc::d2mGateReceivers(6, ty, z, dz, 0.2, 2.0, w);           // open
  CHECK(f == doctest::Approx(1.0)); CHECK(w[3] == doctest::Approx(0.5)); CHECK(w[4] == doctest::Approx(0.5)); CHECK(w[5] == 0.0);
  CHECK(w[0] == 0.0); CHECK(w[1] == 0.0); CHECK(w[2] == 0.0);
  f = Soil_Bgc::d2mGateReceivers(6, ty, z, dz, 0.2, 0.32, w);                 // partial: 0.10 + 0.07
  CHECK(f == doctest::Approx(0.85)); CHECK(w[3] == doctest::Approx(0.1 / 0.17)); CHECK(w[4] == doctest::Approx(0.07 / 0.17));
  CHECK(w[3] + w[4] == doctest::Approx(1.0));
  f = Soil_Bgc::d2mGateReceivers(6, ty, z, dz, 0.2, 0.15, w);                 // exactly at the mineral top
  CHECK(f == 0.0); for (int i = 0; i < 6; ++i) CHECK(w[i] == 0.0);
  f = Soil_Bgc::d2mGateReceivers(6, ty, z, dz, 0.2, 0.0, w);                  // ALD 0: no transfer
  CHECK(f == 0.0); for (int i = 0; i < 6; ++i) CHECK(w[i] == 0.0);
  f = Soil_Bgc::d2mGateReceivers(4, ty, z, dz, 0.2, 2.0, w);                  // only 0.10 m of mineral
  CHECK(f == doctest::Approx(0.5)); CHECK(w[3] == doctest::Approx(1.0));
}

TEST_CASE("gate 1 fully open (>= 0.20 m mineral) == legacy per layer (C and N), to round-off") {
  GateSoil a(0); a.step();
  GateSoil b(1); b.ald(2.0, 0.0).step();
  CHECK(b.sb.d2m_f == doctest::Approx(1.0));
  CHECK(b.sb.d2m_moved_c == doctest::Approx(a.legacy_mobile()).epsilon(1e-13));
  for (int i = 0; i < a.n; ++i) {
    CAPTURE(i);
    CHECK(b.c(i) == doctest::Approx(a.c(i)).epsilon(1e-12));
    CHECK(b.sb.del_sois.orgn[i] == doctest::Approx(a.sb.del_sois.orgn[i]).epsilon(1e-12));
  }
}

TEST_CASE("gate 1 closed (zgate above the mineral top, or ALD 0) == gate 2 exactly, nothing moved") {
  for (int k = 0; k < 3; ++k) {
    CAPTURE(k);
    GateSoil off(2); off.step();
    GateSoil g(1);
    if (k == 0) g.ald(0.10, 0.12);
    if (k == 1) g.ald(MISSING_D, 0.0);
    if (k == 2) g.ald(0.0, 0.0);
    g.step();
    CHECK(g.sb.d2m_f == 0.0); CHECK(g.sb.d2m_moved_c == 0.0); CHECK(g.sb.d2m_moved_n == 0.0);
    for (int i = 0; i < g.n; ++i) {
      CAPTURE(i);
      CHECK(g.sb.del_sois.rawc[i] == off.sb.del_sois.rawc[i]); CHECK(g.sb.del_sois.soma[i] == off.sb.del_sois.soma[i]);
      CHECK(g.sb.del_sois.sompr[i] == off.sb.del_sois.sompr[i]); CHECK(g.sb.del_sois.somcr[i] == off.sb.del_sois.somcr[i]);
      CHECK(g.sb.del_sois.orgn[i] == off.sb.del_sois.orgn[i]);
    }
  }
}

TEST_CASE("previous-year ALD opens the gate when this year's running max is still shallow") {
  GateSoil g(1); g.ald(0.05, 2.0).step();
  CHECK(g.sb.d2m_f == doctest::Approx(1.0));
}

TEST_CASE("gate 1 partial: moved C = f * legacy, N at the donor C:N, column C and N conserved vs gate 2") {
  GateSoil off(2); off.step();
  GateSoil g(1); g.ald(0.32, 0.0).step();
  CHECK(g.sb.d2m_f == doctest::Approx(0.85));
  CHECK(g.sb.d2m_moved_c == doctest::Approx(0.85 * g.legacy_mobile()).epsilon(1e-13));
  double en = 0;
  for (int i = 1; i <= 2; ++i) {
    double tm = 0.85 * (g.sb.del_soi2a.rhrawc[i] + g.sb.del_soi2a.rhsoma[i] + g.sb.del_soi2a.rhsompr[i] + g.sb.del_soi2a.rhsomcr[i]) * 0.06;
    double ts = g.sb.tmp_sois.rawc[i] + g.sb.tmp_sois.soma[i] + g.sb.tmp_sois.sompr[i] + g.sb.tmp_sois.somcr[i];
    en += tm / ts * g.sb.tmp_sois.orgn[i];
    CHECK(off.c(i) - g.c(i) == doctest::Approx(tm).epsilon(1e-10));
    CHECK(off.sb.del_sois.orgn[i] - g.sb.del_sois.orgn[i] == doctest::Approx(tm / ts * g.sb.tmp_sois.orgn[i]).epsilon(1e-10));
  }
  CHECK(g.sb.d2m_moved_n == doctest::Approx(en).epsilon(1e-13));
  CHECK(g.c(3) - off.c(3) == doctest::Approx(g.sb.d2m_moved_c * 0.1 / 0.17).epsilon(1e-10));
  CHECK(g.c(4) - off.c(4) == doctest::Approx(g.sb.d2m_moved_c * 0.07 / 0.17).epsilon(1e-10));
  CHECK(g.sb.del_sois.orgn[4] - off.sb.del_sois.orgn[4] == doctest::Approx(g.sb.d2m_moved_n * 0.07 / 0.17).epsilon(1e-10));
  CHECK(g.c(5) == off.c(5)); CHECK(g.sb.del_sois.orgn[5] == off.sb.del_sois.orgn[5]);
  CHECK(g.sumc() == doctest::Approx(off.sumc()).epsilon(1e-13));
  CHECK(g.sumn() == doctest::Approx(off.sumn()).epsilon(1e-13));
}

TEST_CASE("thin mineral column (0.10 m): legacy loses half of the transfer, gate 1 conserves C and N") {
  GateSoil off(2, 4); off.step();
  GateSoil leg(0, 4); leg.step();
  GateSoil g(1, 4); g.ald(2.0, 0.0).step();
  CHECK(g.sb.d2m_f == doctest::Approx(0.5));
  CHECK(off.sumc() - leg.sumc() == doctest::Approx(0.5 * leg.legacy_mobile()).epsilon(1e-10));
  CHECK(g.sumc() == doctest::Approx(off.sumc()).epsilon(1e-13));
  CHECK(g.sumn() == doctest::Approx(off.sumn()).epsilon(1e-13));
}

TEST_CASE("donor with totsomc <= mobile C: neither C nor N leaves it (C and N together or not at all)") {
  GateSoil off(2); GateSoil g(1);
  for (GateSoil* s : {&off, &g}) { s->sb.tmp_sois.rawc[2] = 0.01; s->sb.tmp_sois.soma[2] = 0; s->sb.tmp_sois.sompr[2] = 0; s->sb.tmp_sois.somcr[2] = 0; }
  off.step(); g.ald(2.0, 0).step();
  CHECK(g.c(2) == off.c(2)); CHECK(g.sb.del_sois.orgn[2] == off.sb.del_sois.orgn[2]);
  CHECK(g.c(1) < off.c(1));
  CHECK(g.sumc() == doctest::Approx(off.sumc()).epsilon(1e-13));
  CHECK(g.sumn() == doctest::Approx(off.sumn()).epsilon(1e-13));
}

TEST_CASE("nfeed 0: gate 1 moves C only, column C conserved") {
  GateSoil off(2, 9, 0.1, 0); off.step();
  GateSoil g(1, 9, 0.1, 0); g.ald(0.32, 0).step();
  CHECK(g.sb.d2m_moved_n == 0.0);
  CHECK(g.sb.d2m_moved_c > 0.0);
  CHECK(g.sumc() == doctest::Approx(off.sumc()).epsilon(1e-13));
}

}
