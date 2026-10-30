// v5q-thawfrac-decomp: soil decomposition modifier blended by the layer thawed fraction,
// model_settings.thawfrac_decomp (0 = legacy, default; 1 = E = f*E_thawed + (1-f)*E_frozen,
// f = 1 - frozenfrac from the freeze/thaw fronts; thawed part at max(ts,0) C and liquid
// saturation min(1, sws/f); frozen part = legacy rhmoistfrozen path at min(ts,0); net N
// mineralization from the thawed RH share). See include/Soil_Bgc.h.
// Refs: Koven et al. 2013, doi:10.5194/bg-10-7109-2013; Lawrence et al. 2019 (CLM5),
//   doi:10.1029/2018MS001583; Yi et al. 2009, doi:10.1029/2008JG000841. The linear blend itself:
//   no published source.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value tmsg(Json::Value v) { Json::Value c; c["model_settings"]["thawfrac_decomp"] = v; return c; }

// Column of nl layers: 0 fibric, 1-2 humic, 3.. mineral; every layer gets ts, frozenfrac, sws.
struct TfSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  TfSoil(int sw, double rmf, int nl = 6) : n(nl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2; sb.d2m_transfer_gate = 2;
    sb.thawfrac_decomp = sw;
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
  TfSoil& lay(int i, double ts, double ff, double sws) { ed.m_sois.ts[i] = ts; ed.m_sois.frozenfrac[i] = ff; ed.m_soid.sws[i] = sws; return *this; }
  TfSoil& step() { sb.deltac(); sb.deltan(); sb.deltastate(); return *this; }
  double rh(int i) const { return sb.del_soi2a.rhrawc[i] + sb.del_soi2a.rhsoma[i] + sb.del_soi2a.rhsompr[i] + sb.del_soi2a.rhsomcr[i]; }
  double totc0(int i) const { return sb.tmp_sois.rawc[i] + sb.tmp_sois.soma[i] + sb.tmp_sois.sompr[i] + sb.tmp_sois.somcr[i]; }
  double dc(int i) const { return sb.del_sois.rawc[i] + sb.del_sois.soma[i] + sb.del_sois.sompr[i] + sb.del_sois.somcr[i]; }
  double sum_dc() const { double s = 0; for (int i = 0; i < n; ++i) s += dc(i); return s; }
  double sum_rh() const { double s = 0; for (int i = 0; i < n; ++i) s += rh(i); return s; }
  double sum_dn() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double n_ext() const { double s = 0; for (int i = 0; i < n; ++i) s += sb.rtnextract[i]; return s; }
};
// legacy modifier with the same parameters
double Eleg(double ts, double sws, double rmf) {
  double m = (ts < 0.0) ? rmf : (sws - 0.0) * (sws - 1.0) / ((sws - 0.0) * (sws - 1.0) - (sws - 0.5) * (sws - 0.5));
  if (m < 0) m = 0;
  return m * std::pow(2.0, ts / 10.0);
}
}

TEST_SUITE("thawfrac_v5q") {

TEST_CASE("thawfrac_decomp default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc") {
  ModelData md; CHECK(md.thawfrac_decomp == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.thawfrac_decomp == 0);
  Soil_Bgc sb; CHECK(sb.thawfrac_decomp == 0);
}

TEST_CASE("thawfrac_decomp validation: integers 0/1 only, anything else FATAL") {
  CHECK(ModelData(tmsg(0)).thawfrac_decomp == 0);
  CHECK(ModelData(tmsg(1)).thawfrac_decomp == 1);
  CHECK(ut::exit_status_of([] { ModelData m(tmsg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tmsg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tmsg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("thawFraction: 1 - frozenfrac clamped to [0,1]; missing/non-finite -> legacy ts indicator") {
  CHECK(Soil_Bgc::thawFraction(0.0, -1.0) == 1.0);
  CHECK(Soil_Bgc::thawFraction(1.0, 3.0) == 0.0);
  CHECK(Soil_Bgc::thawFraction(0.25, -0.003) == doctest::Approx(0.75));
  CHECK(Soil_Bgc::thawFraction(1.0 + 1e-12, 0.0) == 0.0);
  CHECK(Soil_Bgc::thawFraction(-1e-12, 0.0) == 1.0);
  CHECK(Soil_Bgc::thawFraction(MISSING_D, -0.1) == 0.0);
  CHECK(Soil_Bgc::thawFraction(MISSING_D, 0.1) == 1.0);
  CHECK(Soil_Bgc::thawFraction(UIN_D, 0.1) == 1.0);
  CHECK(Soil_Bgc::thawFraction(std::nan(""), -0.1) == 0.0);
}

TEST_CASE("f_thaw = 0 (frozenfrac 1, ts < 0) == legacy frozen rate, every pool and layer; N == legacy") {
  for (double rmf : {0.0, 0.012}) {
    CAPTURE(rmf);
    TfSoil a(0, rmf), b(1, rmf);
    for (TfSoil* s : {&a, &b}) for (int i = 0; i < s->n; ++i) s->lay(i, -0.5 - 0.3 * i, 1.0, 0.1);
    a.step(); b.step();
    for (int i = 0; i < a.n; ++i) {
      CAPTURE(i);
      CHECK(b.sb.del_soi2a.rhrawc[i] == doctest::Approx(a.sb.del_soi2a.rhrawc[i]).epsilon(1e-14));
      CHECK(b.sb.del_soi2a.rhsoma[i] == doctest::Approx(a.sb.del_soi2a.rhsoma[i]).epsilon(1e-14));
      CHECK(b.sb.del_soi2a.rhsompr[i] == doctest::Approx(a.sb.del_soi2a.rhsompr[i]).epsilon(1e-14));
      CHECK(b.sb.del_soi2a.rhsomcr[i] == doctest::Approx(a.sb.del_soi2a.rhsomcr[i]).epsilon(1e-14));
      CHECK(b.rh(i) == doctest::Approx(0.05 * b.sb.tmp_sois.rawc[i] * Eleg(-0.5 - 0.3 * i, 0.1, rmf)
                                       + 0.02 * b.sb.tmp_sois.soma[i] * Eleg(-0.5 - 0.3 * i, 0.1, rmf)
                                       + 0.005 * b.sb.tmp_sois.sompr[i] * Eleg(-0.5 - 0.3 * i, 0.1, rmf)
                                       + 0.0005 * b.sb.tmp_sois.somcr[i] * Eleg(-0.5 - 0.3 * i, 0.1, rmf)).epsilon(1e-12));
      CHECK(b.sb.tfd_thawshare[i] == 0.0);
      CHECK(b.sb.del_soi2soi.netnmin[i] == a.sb.del_soi2soi.netnmin[i]);
      CHECK(b.sb.del_sois.orgn[i] == doctest::Approx(a.sb.del_sois.orgn[i]).epsilon(1e-13));
      CHECK(b.sb.del_sois.avln[i] == doctest::Approx(a.sb.del_sois.avln[i]).epsilon(1e-13));
    }
    if (rmf == 0.0) CHECK(b.sum_rh() == 0.0);
  }
}

TEST_CASE("f_thaw = 1 (frozenfrac 0, ts > 0) == legacy thawed rate, every pool and layer; N == legacy") {
  TfSoil a(0, 0.0), b(1, 0.0);
  for (TfSoil* s : {&a, &b}) for (int i = 0; i < s->n; ++i) s->lay(i, 0.5 + 1.1 * i, 0.0, 0.15 + 0.12 * i);
  a.step(); b.step();
  for (int i = 0; i < a.n; ++i) {
    CAPTURE(i);
    CHECK(b.sb.del_soi2a.rhrawc[i] == doctest::Approx(a.sb.del_soi2a.rhrawc[i]).epsilon(1e-14));
    CHECK(b.sb.del_soi2a.rhsoma[i] == doctest::Approx(a.sb.del_soi2a.rhsoma[i]).epsilon(1e-14));
    CHECK(b.sb.del_soi2a.rhsompr[i] == doctest::Approx(a.sb.del_soi2a.rhsompr[i]).epsilon(1e-14));
    CHECK(b.sb.del_soi2a.rhsomcr[i] == doctest::Approx(a.sb.del_soi2a.rhsomcr[i]).epsilon(1e-14));
    CHECK(b.sb.tfd_thawshare[i] == 1.0);
    CHECK(a.rh(i) > 0.0);
    CHECK(b.sb.del_soi2soi.netnmin[i] == doctest::Approx(a.sb.del_soi2soi.netnmin[i]).epsilon(1e-13));
    CHECK(b.sb.del_sois.orgn[i] == doctest::Approx(a.sb.del_sois.orgn[i]).epsilon(1e-13));
    CHECK(b.sb.del_sois.avln[i] == doctest::Approx(a.sb.del_sois.avln[i]).epsilon(1e-13));
  }
}

TEST_CASE("rate linear in f_thaw between the frozen and thawed rates (thawed-part state fixed), netnmin linear") {
  // thawed part: saturation 0.6, ts = -0.003 (layer mean pinned just below 0 C) -> thawed part at 0 C
  for (double rmf : {0.0, 0.012}) {
    CAPTURE(rmf);
    const double sth = 0.6, ts = -0.003;
    TfSoil t0(1, rmf), t1(1, rmf);
    t0.lay(4, ts, 1.0, 0.0).step();          // f = 0
    t1.lay(4, ts, 0.0, sth).step();          // f = 1
    const double r0 = t0.rh(4), r1 = t1.rh(4);
    CHECK(r1 == doctest::Approx(t1.totc0(4) * 0 + (0.05 * t1.sb.tmp_sois.rawc[4] + 0.02 * t1.sb.tmp_sois.soma[4]
                                + 0.005 * t1.sb.tmp_sois.sompr[4] + 0.0005 * t1.sb.tmp_sois.somcr[4]) * Eleg(0.0, sth, rmf)).epsilon(1e-12));
    const double n1 = t1.sb.del_soi2soi.netnmin[4];
    CHECK(t0.sb.del_soi2soi.netnmin[4] == 0.0);
    CHECK(n1 != 0.0);
    for (double f : {0.1, 0.25, 0.5, 0.75, 0.9}) {
      CAPTURE(f);
      TfSoil m(1, rmf); m.lay(4, ts, 1.0 - f, sth * f).step();   // layer sws = f * thawed-part saturation
      CHECK(m.rh(4) == doctest::Approx(f * r1 + (1.0 - f) * r0).epsilon(1e-12));
      CHECK(m.sb.del_soi2soi.netnmin[4] == doctest::Approx(f * n1).epsilon(1e-12));   // N from the thawed part only
      CHECK(m.sb.tfd_thawshare[4] == doctest::Approx(f * r1 / (f * r1 + (1.0 - f) * r0)).epsilon(1e-12));
      // legacy (ts < 0): the whole layer frozen
      TfSoil l(0, rmf); l.lay(4, ts, 1.0 - f, sth * f).step();
      CHECK(l.rh(4) == doctest::Approx(r0).epsilon(1e-12));
      if (rmf == 0.0) { CHECK(l.rh(4) == 0.0); CHECK(m.rh(4) > 0.0); }
    }
  }
}

TEST_CASE("C and N conservation with switch 1 (mixed column): column dC = -RH, column dN = inputs - losses") {
  for (int sw : {0, 1}) {
    for (double rmf : {0.0, 0.012}) {
      CAPTURE(sw); CAPTURE(rmf);
      TfSoil s(sw, rmf);
      s.lay(0, 3.0, 0.0, 0.5).lay(1, 1.0, 0.2, 0.4).lay(2, -0.002, 0.5, 0.3).lay(3, -0.004, 0.7, 0.2).lay(4, -1.0, 1.0, 0.05).lay(5, -2.0, 1.0, 0.02);
      s.step();
      CHECK(s.sum_rh() > 0.0);
      CHECK(s.sum_dc() == doctest::Approx(-s.sum_rh()).epsilon(1e-12));
      const double nin = s.bd.m_a2soi.orgninput + s.bd.m_a2soi.avlninput;
      const double nout = s.sb.del_soi2l.avlnlost + s.sb.del_soi2l.orgnlost + s.n_ext();
      CHECK(s.sum_dn() == doctest::Approx(nin - nout).epsilon(1e-12));
    }
  }
}

TEST_CASE("switch 0 unchanged: frozenfrac has no effect, partly thawed layer with ts < 0 has the legacy frozen rate") {
  TfSoil a(0, 0.0), b(0, 0.0);
  for (int i = 0; i < a.n; ++i) { a.lay(i, -0.003, 0.0, 0.3); b.lay(i, -0.003, 0.5, 0.3); }
  a.lay(0, 2.0, 0.0, 0.3); b.lay(0, 2.0, 0.9, 0.3);
  a.step(); b.step();
  for (int i = 0; i < a.n; ++i) {
    CAPTURE(i);
    CHECK(a.rh(i) == b.rh(i));
    CHECK(a.sb.del_soi2a.rhsomcr[i] == b.sb.del_soi2a.rhsomcr[i]);
    CHECK(a.sb.del_soi2soi.netnmin[i] == b.sb.del_soi2soi.netnmin[i]);
    CHECK(a.sb.del_sois.orgn[i] == b.sb.del_sois.orgn[i]);
    CHECK(a.bd.m_soid.rhmoist[i] == b.bd.m_soid.rhmoist[i]);
    CHECK(a.bd.m_soid.rhmoist[i] == (i == 0 ? a.bd.m_soid.rhmoist[0] : 0.0));
  }
  CHECK(a.rh(0) > 0.0);
}

TEST_CASE("wood debris inherits the first non-moss layer modifier (switch 1: blended)") {
  TfSoil a(1, 0.0);
  a.sb.tmp_sois.wdebrisc = 50.0;
  a.lay(0, -0.003, 0.4, 0.3);
  a.sb.deltac();
  const double E0 = a.bd.m_soid.rhmoist[0] * a.bd.m_soid.rhq10[0];
  CHECK(E0 > 0.0);
  CHECK(a.sb.del_soi2a.rhwdeb == doctest::Approx(0.05 * 50.0 * E0).epsilon(1e-13));
  CHECK(E0 == doctest::Approx(0.6 * Eleg(0.0, 0.5, 0.0)).epsilon(1e-12));
}

} // TEST_SUITE thawfrac_v5q
