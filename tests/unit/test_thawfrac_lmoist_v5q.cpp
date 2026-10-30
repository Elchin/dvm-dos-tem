// v5q-thawfrac-lmoist: thawed-part moisture of the thawed-fraction decomposition modifier,
// model_settings.thawfrac_moist_mode (0 = min(1, sws/f), default, v5q-thawfrac-decomp unchanged;
// 1 = layer-mean sws; 1 requires thawfrac_decomp = 1). See include/Soil_Bgc.h.
// Refs: Koven et al. 2013, doi:10.5194/bg-10-7109-2013; Lawrence et al. 2019 (CLM5),
//   doi:10.1029/2018MS001583; Yi et al. 2009, doi:10.1029/2008JG000841. The mode choice itself:
//   no published source.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value lmsg(Json::Value mode, Json::Value tfd = 1) {
  Json::Value c; c["model_settings"]["thawfrac_decomp"] = tfd; c["model_settings"]["thawfrac_moist_mode"] = mode; return c;
}

// Column as test_thawfrac_v5q TfSoil: 0 fibric, 1-2 humic, 3.. mineral.
struct LmSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  LmSoil(int tfd, int mode, double rmf, int nl = 6) : n(nl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2; sb.d2m_transfer_gate = 2;
    sb.thawfrac_decomp = tfd; sb.thawfrac_moist_mode = mode;
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    sb.bgcpar.fnloss = 0.1; sb.bgcpar.kn2 = 1.0; sb.bgcpar.nmincnsoil = 1.0; sb.calpar.micbnup = 0.5;
    sb.bgcpar.rhq10 = 2.0; sb.bgcpar.rhmoistfrozen = rmf;
    sb.bgcpar.moistmin = 0.0; sb.bgcpar.moistopt = 0.5; sb.bgcpar.moistmax = 1.0;
    cd.m_soil.numsl = nl; cd.m_soil.deepthick = 0.1;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < nl;
      int t = (i == 0) ? 1 : (i <= 2 ? 2 : 3);
      double dz = (i <= 2) ? 0.05 : 0.1;
      cd.m_soil.type[i] = t; cd.m_soil.dz[i] = dz; cd.m_soil.z[i] = z; z += dz;
      sb.bgcpar.kdrawc[i] = 0.05; sb.bgcpar.kdsoma[i] = 0.02; sb.bgcpar.kdsompr[i] = 0.005; sb.bgcpar.kdsomcr[i] = 0.0005;
      sb.tmp_sois.rawc[i] = 100 + 3 * i; sb.tmp_sois.soma[i] = 50 + i; sb.tmp_sois.sompr[i] = 200 - 5 * i; sb.tmp_sois.somcr[i] = 300 + 7 * i;
      sb.tmp_sois.orgn[i] = 20 + 0.5 * i; sb.tmp_sois.avln[i] = on ? 0.4 : 0.0; bd.m_sois.avln[i] = on ? 0.4 : 0.0;
      sb.ltrflc[i] = 0; sb.ltrfln[i] = 0; sb.rtnextract[i] = on ? 0.002 : 0.0; bd.m_soi2v.nextract[i] = sb.rtnextract[i];
      bd.m_v2soi.rtlfalfrac[i] = on ? 0.5 : 0.0; bd.m_soid.knmoist[i] = 0.5;
      ed.m_sois.liq[i] = on ? 20.0 : 0.0; ed.m_sois.ts[i] = 2.0; ed.m_sois.frozenfrac[i] = 0.0; ed.m_soid.sws[i] = 0.4;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0;
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
  LmSoil& lay(int i, double ts, double ff, double sws) { ed.m_sois.ts[i] = ts; ed.m_sois.frozenfrac[i] = ff; ed.m_soid.sws[i] = sws; return *this; }
  LmSoil& step() { sb.deltac(); sb.deltan(); sb.deltastate(); return *this; }
  double rh(int i) const { return sb.del_soi2a.rhrawc[i] + sb.del_soi2a.rhsoma[i] + sb.del_soi2a.rhsompr[i] + sb.del_soi2a.rhsomcr[i]; }
  double kc(int i) const { return 0.05 * sb.tmp_sois.rawc[i] + 0.02 * sb.tmp_sois.soma[i] + 0.005 * sb.tmp_sois.sompr[i] + 0.0005 * sb.tmp_sois.somcr[i]; }
  double dc(int i) const { return sb.del_sois.rawc[i] + sb.del_sois.soma[i] + sb.del_sois.sompr[i] + sb.del_sois.somcr[i]; }
  double sum_dc() const { double s = 0; for (int i = 0; i < n; ++i) s += dc(i); return s; }
  double sum_rh() const { double s = 0; for (int i = 0; i < n; ++i) s += rh(i); return s; }
  double sum_dn() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double n_ext() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.rtnextract[i]; return s; }
};
// rhmoist with moistmin 0, moistopt 0.5, moistmax 1 (as Soil_Bgc::getRhmoist), clamped at 0
double M(double sws) { double m = sws * (sws - 1.0) / (sws * (sws - 1.0) - (sws - 0.5) * (sws - 0.5)); return m < 0 ? 0 : m; }
double Q(double t) { return std::pow(2.0, t / 10.0); }
}

TEST_SUITE("thawfrac_lmoist_v5q") {

TEST_CASE("thawfrac_moist_mode default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc") {
  ModelData md; CHECK(md.thawfrac_moist_mode == 0);
  Json::Value c; c["model_settings"]["thawfrac_decomp"] = 1; ModelData mj(c); CHECK(mj.thawfrac_moist_mode == 0);
  Soil_Bgc sb; CHECK(sb.thawfrac_moist_mode == 0);
}

TEST_CASE("thawfrac_moist_mode validation: integers 0/1 only; 1 requires thawfrac_decomp 1; anything else FATAL") {
  CHECK(ModelData(lmsg(0)).thawfrac_moist_mode == 0);
  CHECK(ModelData(lmsg(0, 0)).thawfrac_moist_mode == 0);
  CHECK(ModelData(lmsg(1)).thawfrac_moist_mode == 1);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(1, 0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("thawedPartSws: mode 0 = min(1, sws/f) (sws at f = 1), mode 1 = sws") {
  CHECK(Soil_Bgc::thawedPartSws(0.3, 0.5, 0) == doctest::Approx(0.6));
  CHECK(Soil_Bgc::thawedPartSws(0.8, 0.5, 0) == 1.0);
  CHECK(Soil_Bgc::thawedPartSws(0.3, 1.0, 0) == 0.3);
  CHECK(Soil_Bgc::thawedPartSws(0.3, 0.1, 0) == 1.0);
  for (double f : {0.01, 0.1, 0.5, 0.99, 1.0}) { CAPTURE(f); CHECK(Soil_Bgc::thawedPartSws(0.3, f, 1) == 0.3); }
  // mode 0 expression unchanged (bit-identical to v5q-thawfrac-decomp)
  for (double f : {0.07, 0.33, 0.6}) for (double s : {0.01, 0.2, 0.45}) {
    CAPTURE(f); CAPTURE(s);
    CHECK(Soil_Bgc::thawedPartSws(s, f, 0) == fmin(1.0, s / f));
  }
}

TEST_CASE("thawfracModifier: mode 0 == thawed-part formula with sws/f; mode 1 == layer-mean sws; equal at f = 1 and f = 0") {
  for (double rmf : {0.0, 0.012}) {
    CAPTURE(rmf);
    for (double ts : {-0.004, 0.0, 1.5}) for (double f : {0.0, 0.2, 0.5, 0.8, 1.0}) for (double sws : {0.05, 0.25, 0.45}) {
      CAPTURE(ts); CAPTURE(f); CAPTURE(sws);
      LmSoil a(1, 0, rmf), b(1, 1, rmf);
      double sa = -1, sbh = -1;
      const double Ea = a.sb.thawfracModifier(f, ts, sws, sa);
      const double Eb = b.sb.thawfracModifier(f, ts, sws, sbh);
      const double fr = rmf * Q(fmin(ts, 0.0));
      const double s0 = (f < 1.0 && f > 0.0) ? fmin(1.0, sws / f) : sws;
      CHECK(Ea == doctest::Approx(f * M(s0) * Q(fmax(ts, 0.0)) + (1 - f) * fr).epsilon(1e-13));
      CHECK(Eb == doctest::Approx(f * M(sws) * Q(fmax(ts, 0.0)) + (1 - f) * fr).epsilon(1e-13));
      if (f == 1.0 || f == 0.0) { CHECK(Ea == Eb); CHECK(sa == sbh); }
    }
  }
  // the TVC-like case: thin thawed part of a wet layer; mode 0 caps at saturation (rhmoist 0 at moistmax 1)
  LmSoil a(1, 0, 0.0), b(1, 1, 0.0);
  double sh;
  CHECK(a.sb.thawfracModifier(0.3, -0.003, 0.4, sh) == 0.0);
  CHECK(b.sb.thawfracModifier(0.3, -0.003, 0.4, sh) == doctest::Approx(0.3 * M(0.4)).epsilon(1e-13));
}

TEST_CASE("deltac/deltan with mode 1: per-layer RH and netnmin from the layer-mean moisture; mode 0 unchanged") {
  LmSoil a(1, 0, 0.0), b(1, 1, 0.0);
  for (LmSoil* s : {&a, &b}) s->lay(2, -0.003, 0.5, 0.2).lay(3, -0.003, 0.75, 0.1).lay(4, 1.0, 0.0, 0.3);
  a.step(); b.step();
  CHECK(a.rh(2) == doctest::Approx(a.kc(2) * 0.5 * M(0.4)).epsilon(1e-12));
  CHECK(b.rh(2) == doctest::Approx(b.kc(2) * 0.5 * M(0.2)).epsilon(1e-12));
  CHECK(a.rh(3) == doctest::Approx(a.kc(3) * 0.25 * M(0.4)).epsilon(1e-12));
  CHECK(b.rh(3) == doctest::Approx(b.kc(3) * 0.25 * M(0.1)).epsilon(1e-12));
  CHECK(b.rh(4) == a.rh(4));                                        // f = 1: identical
  CHECK(b.sb.del_soi2soi.netnmin[4] == a.sb.del_soi2soi.netnmin[4]);
  CHECK(b.sb.tfd_thawshare[2] == 1.0);                              // rmf 0: all RH from the thawed part
  CHECK(b.sb.del_soi2soi.netnmin[2] != a.sb.del_soi2soi.netnmin[2]);
}

TEST_CASE("C and N conservation with mode 1 (mixed column): column dC = -RH, column dN = inputs - losses") {
  for (int mode : {0, 1}) for (double rmf : {0.0, 0.012}) {
    CAPTURE(mode); CAPTURE(rmf);
    LmSoil s(1, mode, rmf);
    s.lay(0, 3.0, 0.0, 0.5).lay(1, 1.0, 0.2, 0.4).lay(2, -0.002, 0.5, 0.3).lay(3, -0.004, 0.7, 0.2).lay(4, -1.0, 1.0, 0.05).lay(5, -2.0, 1.0, 0.02);
    s.step();
    CHECK(s.sum_rh() > 0.0);
    CHECK(s.sum_dc() == doctest::Approx(-s.sum_rh()).epsilon(1e-12));
    const double nin = s.bd.m_a2soi.orgninput + s.bd.m_a2soi.avlninput;
    const double nout = s.sb.del_soi2l.avlnlost + s.sb.del_soi2l.orgnlost + s.n_ext();
    CHECK(s.sum_dn() == doctest::Approx(nin - nout).epsilon(1e-12));
  }
}

TEST_CASE("mode ignored with thawfrac_decomp 0 (legacy path)") {
  LmSoil a(0, 0, 0.0), b(0, 1, 0.0);
  for (LmSoil* s : {&a, &b}) s->lay(1, 2.0, 0.4, 0.3).lay(2, -0.003, 0.5, 0.2);
  a.step(); b.step();
  for (int i = 0; i < a.n; ++i) {
    CAPTURE(i);
    CHECK(a.rh(i) == b.rh(i));
    CHECK(a.bd.m_soid.rhmoist[i] == b.bd.m_soid.rhmoist[i]);
    CHECK(a.sb.del_soi2soi.netnmin[i] == b.sb.del_soi2soi.netnmin[i]);
  }
}

} // TEST_SUITE thawfrac_lmoist_v5q
