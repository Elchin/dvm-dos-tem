// wetland_deep_C phase 2: model_settings.deep_carbon (1 = deep soil C burial scheme of upstream
// deep_carbon_v2 7f607479 with the fork's conservation fix a59ce891, DEFAULT; 0 = v5q burial path).
// Refs: deep_carbon_v2 commit messages (fractions: no published source); Yi et al. 2010,
//   doi:10.1029/2010JG001302; Hugelius et al. 2014, doi:10.5194/bg-11-6573-2014;
//   Koven et al. 2013, doi:10.5194/bg-10-7109-2013. See include/Soil_Bgc.h.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value msg(Json::Value v) { Json::Value c; c["model_settings"]["deep_carbon"] = v; return c; }

// Column as test_d2mgate_v5q: 0 fibric (0.05 m), 1-2 humic (0.05 m each, deepthick 0.10),
// 3.. mineral (0.10 m each; 6 layers = 0.6 m < deepc_xtopmlthick 5.2 -> all are receivers).
struct DcSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  DcSoil(int dc, int gate, int nl = 9, double mdz = 0.1, int nfeed = 1) : n(nl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = nfeed; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2; sb.d2m_transfer_gate = gate; sb.deep_carbon = dc; sb.n_som_move_fix = 1;
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
  DcSoil& ald(double runmax, double prev) { ed.y_soid.ald = runmax; cd.d2m_ald_prevyr = prev; return *this; }
  DcSoil& step() { sb.deltan(); sb.deltastate(); return *this; }
  double c(int i) const { return sb.del_sois.rawc[i] + sb.del_sois.soma[i] + sb.del_sois.sompr[i] + sb.del_sois.somcr[i]; }
  double sumc() const { double s = 0; for (int i = 0; i < n; ++i) s += c(i); return s; }
  double sumn() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.orgn[i]; return s; }
  DcSoil& noop() { for (int k = 0; k < 4; ++k) { sb.deepc_s2df[k] = 0.0; sb.deepc_d2mf[k] = 0.0; } return *this; }
  double p(int k, int i) const { return k == 0 ? sb.del_sois.rawc[i] : k == 1 ? sb.del_sois.soma[i] : k == 2 ? sb.del_sois.sompr[i] : sb.del_sois.somcr[i]; }
  double tc(int i) const { return sb.tmp_sois.rawc[i] + sb.tmp_sois.soma[i] + sb.tmp_sois.sompr[i] + sb.tmp_sois.somcr[i]; }
};
}

TEST_SUITE("deepc") {

TEST_CASE("deep_carbon: ModelData default 1 (with and without the key), Soil_Bgc member default 0") {
  ModelData md; CHECK(md.deep_carbon == 1);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.deep_carbon == 1);
  Soil_Bgc sb; CHECK(sb.deep_carbon == 0);
  CHECK(sb.deepc_s2df[0] == 0.5); CHECK(sb.deepc_s2df[1] == 1.0); CHECK(sb.deepc_s2df[2] == 1.0); CHECK(sb.deepc_s2df[3] == 1.0);
  CHECK(sb.deepc_d2mf[0] == 0.5); CHECK(sb.deepc_d2mf[3] == 1.0);
  CHECK(sb.deepc_xtopmlthick == 5.2); CHECK(sb.deepc_minec_maxdepth == 5.0);
}

TEST_CASE("deep_carbon validation: integers 0/1 only, anything else FATAL") {
  CHECK(ModelData(msg(0)).deep_carbon == 0);
  CHECK(ModelData(msg(1)).deep_carbon == 1);
  CHECK(ut::exit_status_of([] { ModelData m(msg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("deep_carbon 0 ignores the deep-C fractions (legacy path untouched)") {
  DcSoil a(0, 0); a.step();
  DcSoil b(0, 0); b.noop(); b.sb.deepc_xtopmlthick = 0.3; b.step();
  for (int i = 0; i < a.n; ++i) { CAPTURE(i); for (int k = 0; k < 4; ++k) CHECK(a.p(k, i) == b.p(k, i)); CHECK(a.sb.del_sois.orgn[i] == b.sb.del_sois.orgn[i]); }
  DcSoil d(1, 0); d.step();
  bool differs = false; for (int i = 0; i < a.n; ++i) differs = differs || (d.c(i) != a.c(i));
  CHECK(differs);
}

TEST_CASE("deep_carbon 1, gate 0: per-pool fractions of the positive increments, thickness shares, column C and N conserved") {
  DcSoil z(1, 0); z.noop().step();           // pre-transfer increments (no-op fractions)
  DcSoil g(1, 0); g.step();
  double s2d[4] = {0, 0, 0, 0}, d2m[4] = {0, 0, 0, 0};
  for (int k = 0; k < 4; ++k) {
    const double fr = g.sb.deepc_s2df[k];
    const double pre = z.p(k, 0);
    s2d[k] = pre > 0 ? pre * fr : 0.0;
    CHECK(g.p(k, 0) == doctest::Approx(pre - s2d[k]).epsilon(1e-12));
    for (int i = 1; i <= 2; ++i) {           // humic: receive half each, then donate d2mf of positive part
      const double rec = z.p(k, i) + 0.5 * s2d[k];
      const double out = rec > 0 ? rec * g.sb.deepc_d2mf[k] : 0.0;
      d2m[k] += out;
      CHECK(g.p(k, i) == doctest::Approx(rec - out).epsilon(1e-12));
    }
  }
  for (int i = 3; i < g.n; ++i) { CAPTURE(i); for (int k = 0; k < 4; ++k) CHECK(g.p(k, i) == doctest::Approx(z.p(k, i) + d2m[k] / 6.0).epsilon(1e-12)); }
  CHECK(g.sb.d2m_moved_c == doctest::Approx(d2m[0] + d2m[1] + d2m[2] + d2m[3]).epsilon(1e-13));
  // N: the fibric layer loses its own buried C at its C:N (a59ce891)
  const double s2dc = s2d[0] + s2d[1] + s2d[2] + s2d[3];
  CHECK(z.sb.del_sois.orgn[0] - g.sb.del_sois.orgn[0] == doctest::Approx(s2dc / g.tc(0) * g.sb.tmp_sois.orgn[0]).epsilon(1e-10));
  CHECK(g.sumc() == doctest::Approx(z.sumc()).epsilon(1e-13));
  CHECK(g.sumn() == doctest::Approx(z.sumn()).epsilon(1e-13));
}

TEST_CASE("deep_carbon 1: C and N conserved for gates 0/1(open, partial, closed)/2, nfeed 1/0, n_som_move_fix 1/0, no mineral layer") {
  for (int nf = 0; nf <= 1; ++nf) for (int sm = 0; sm <= 1; ++sm) for (int gate = 0; gate <= 2; ++gate) for (int nl = 3; nl <= 9; nl += 3) {
    double alds[3] = {2.0, 0.32, 0.0};
    for (int a = 0; a < (gate == 1 ? 3 : 1); ++a) {
      CAPTURE(nf); CAPTURE(sm); CAPTURE(gate); CAPTURE(nl); CAPTURE(a);
      DcSoil z(1, gate, nl, 0.1, nf); z.sb.n_som_move_fix = sm; z.noop(); z.ald(alds[a], 0).step();
      DcSoil g(1, gate, nl, 0.1, nf); g.sb.n_som_move_fix = sm; g.ald(alds[a], 0).step();
      CHECK(g.sumc() == doctest::Approx(z.sumc()).epsilon(1e-13));
      if (nf == 1 && sm == 1) CHECK(g.sumn() == doctest::Approx(z.sumn()).epsilon(1e-13));
      if (gate == 2 || nl == 3 || (gate == 1 && a == 2)) { CHECK(g.sb.d2m_moved_c == 0.0); }
      else { CHECK(g.sb.d2m_moved_c > 0.0); }
    }
  }
}

TEST_CASE("deep_carbon 1, gate 1 partial: f = eligible thickness / 5.2, moved C = f * gate-0 donor C") {
  DcSoil o(1, 0); o.step();
  DcSoil g(1, 1); g.ald(0.32, 0).step();
  CHECK(g.sb.d2m_f == doctest::Approx(0.17 / 5.2));
  CHECK(g.sb.d2m_moved_c == doctest::Approx(0.17 / 5.2 * o.sb.d2m_moved_c).epsilon(1e-12));
  CHECK(g.c(5) - DcSoil(1, 2).step().c(5) == 0.0);    // below the gate: receives nothing
}

}
