// v5q-rh-frozen: model_settings.rh_frozen_mode (0 = legacy frozen-part rate rhmoistfrozen*Q10^(min(ts,0)/10),
// default, bit-identical; 1 = above the permafrost table the frozen-part rate is the unfrozen-water fraction
// phi(T) x rhmoist(liquid+ice saturation) x Q10^(min(ts,0)/10); 1 requires thawfrac_decomp = 1).
// See include/Soil_Bgc.h.
// Refs: Schaefer & Jafarov 2016, doi:10.5194/bg-13-1991-2016 (Eqs. 4-6, Table 1); Romanovsky &
//   Osterkamp 2000, doi:10.1002/1099-1530(200007/09)11:3<219::AID-PPP352>3.0.CO;2-7; Nicolsky et al.
//   2009, doi:10.1016/j.coldregions.2008.03.003; Kurylyk & Watanabe 2013,
//   doi:10.1016/j.advwatres.2013.07.016; Koven et al. 2013, doi:10.5194/bg-10-7109-2013; Lawrence et
//   al. 2019, doi:10.1029/2018MS001583. The permafrost-table restriction: no published source.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value rfmsg(Json::Value mode, Json::Value tfd = 1) {
  Json::Value c; c["model_settings"]["thawfrac_decomp"] = tfd; c["model_settings"]["rh_frozen_mode"] = mode; return c;
}

// Column: 0 fibric (type 1), 1-2 humic (type 2), 3.. mineral (type 3, silt loam 20/60/20 clay/silt/sand).
struct RfSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  int n;
  RfSoil(int mode, double rmf = 0.0, int nl = 6) : n(nl) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.avlnflg = 1; sb.n_leach_fix = 1; sb.n_rootless_fix = 2; sb.d2m_transfer_gate = 1;
    sb.thawfrac_decomp = 1; sb.thawfrac_moist_mode = 0; sb.rh_frozen_mode = mode;
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
      ed.m_sois.liq[i] = on ? 20.0 : 0.0; ed.m_sois.ts[i] = -2.0; ed.m_sois.frozenfrac[i] = 1.0;
      ed.m_soid.sws[i] = 0.0; ed.m_soid.iwc[i] = (t == 3) ? 0.2 : 0.3;
      sb.rhf_por[i] = (t == 3) ? 0.45 : 0.9;
      if (t == 3) { sb.rhf_fclay[i] = 20; sb.rhf_fsilt[i] = 60; sb.rhf_fsand[i] = 20; sb.rhf_texok[i] = 1; }
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0.01; sb.del_soi2l.avlnlost = 0;
    ed.m_sois.draindepth = 0.0; ed.m_soi2l.qdrain = 0.0; ed.m_soi2l.qover = 0.0;
    // winter: no thaw this month or year; previous-year ALD 0.25 m (mid layer 3: 0.15..0.25 -> fully above)
    ed.m_soid.ald = 0.0; ed.y_soid.ald = 0.0; cd.d2m_ald_prevyr = 0.25;
    sb.totdzliq = 0; sb.totdzavln = 0;
  }
  RfSoil& step() { sb.deltac(); sb.deltan(); sb.deltastate(); return *this; }
  double rh(int i) const { return sb.del_soi2a.rhrawc[i] + sb.del_soi2a.rhsoma[i] + sb.del_soi2a.rhsompr[i] + sb.del_soi2a.rhsomcr[i]; }
  double kc(int i) const { return 0.05 * sb.tmp_sois.rawc[i] + 0.02 * sb.tmp_sois.soma[i] + 0.005 * sb.tmp_sois.sompr[i] + 0.0005 * sb.tmp_sois.somcr[i]; }
  double somc(int i) const { return sb.tmp_sois.rawc[i] + sb.tmp_sois.soma[i] + sb.tmp_sois.sompr[i] + sb.tmp_sois.somcr[i]; }
};
double M(double s) { double m = s * (s - 1.0) / (s * (s - 1.0) - (s - 0.5) * (s - 0.5)); return m < 0 ? 0 : m; }
double Q(double t) { return std::pow(2.0, t / 10.0); }
double P(double t, double b) { return t > 0 ? 1.0 : std::pow((0.1 - t) / 0.01, b); }
}

TEST_SUITE("rhfrozen_v5q") {

TEST_CASE("rh_frozen_mode default 0: ModelData(), ModelData(Json) without the key, Soil_Bgc") {
  ModelData md; CHECK(md.rh_frozen_mode == 0);
  Json::Value c; c["model_settings"]["thawfrac_decomp"] = 1; ModelData mj(c); CHECK(mj.rh_frozen_mode == 0);
  Soil_Bgc sb; CHECK(sb.rh_frozen_mode == 0);
}

TEST_CASE("rh_frozen_mode validation: integers 0/1 only; 1 requires thawfrac_decomp 1; anything else FATAL") {
  CHECK(ModelData(rfmsg(0)).rh_frozen_mode == 0);
  CHECK(ModelData(rfmsg(0, 0)).rh_frozen_mode == 0);
  CHECK(ModelData(rfmsg(1)).rh_frozen_mode == 1);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg(1, 0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(rfmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("sjPhiPure: SJ2016 Eq. 4 / Table 1 (phi = 1 above 0 degC, phi_crit at 0 degC, power law below)") {
  CHECK(Soil_Bgc::sjPhiPure(0.5, -0.3) == 1.0);
  CHECK(Soil_Bgc::sjPhiPure(1e-9, -1.0) == 1.0);
  // phi_crit (Table 1): clay 0.50, silt 0.32, sand 0.13, organic 0.10 (= 10^b)
  CHECK(Soil_Bgc::sjPhiPure(0.0, -0.3) == doctest::Approx(0.50).epsilon(0.005));
  CHECK(Soil_Bgc::sjPhiPure(0.0, -0.5) == doctest::Approx(0.32).epsilon(0.02));
  CHECK(Soil_Bgc::sjPhiPure(0.0, -0.9) == doctest::Approx(0.13).epsilon(0.02));
  CHECK(Soil_Bgc::sjPhiPure(0.0, -1.0) == doctest::Approx(0.10));
  CHECK(Soil_Bgc::sjPhiPure(-1.0, -1.0) == doctest::Approx(1.0 / 110.0));
  CHECK(Soil_Bgc::sjPhiPure(-5.0, -0.5) == doctest::Approx(std::pow(510.0, -0.5)));
  // monotone: colder -> less liquid; finer texture -> more liquid
  double prev = 2.0;
  for (double t : {0.0, -0.1, -0.5, -1.0, -3.0, -10.0, -30.0}) {
    CAPTURE(t);
    double p = Soil_Bgc::sjPhiPure(t, -0.5); CHECK(p < prev); CHECK(p > 0.0); prev = p;
    CHECK(Soil_Bgc::sjPhiPure(t, -0.3) > Soil_Bgc::sjPhiPure(t, -0.5));
    CHECK(Soil_Bgc::sjPhiPure(t, -0.5) > Soil_Bgc::sjPhiPure(t, -0.9));
    CHECK(Soil_Bgc::sjPhiPure(t, -0.9) > Soil_Bgc::sjPhiPure(t, -1.0));
  }
}

TEST_CASE("sjLayerPhi: SJ2016 Eq. 5 mix, texture renormalised, no texture = organic") {
  for (double t : {0.0, -0.7, -4.0}) {
    CAPTURE(t);
    CHECK(Soil_Bgc::sjLayerPhi(t, 1.0, 20, 60, 20) == doctest::Approx(P(t, -1.0)));
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.0, 1, 0, 0) == doctest::Approx(P(t, -0.3)));
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.0, 0, 0, 1) == doctest::Approx(P(t, -0.9)));
    double mix = 0.2 * P(t, -0.3) + 0.6 * P(t, -0.5) + 0.2 * P(t, -0.9);
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.0, 20, 60, 20) == doctest::Approx(mix));   // percent in, renormalised
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.0, 0.2, 0.6, 0.2) == doctest::Approx(mix));
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.4, 20, 60, 20) == doctest::Approx(0.6 * mix + 0.4 * P(t, -1.0)));
    CHECK(Soil_Bgc::sjLayerPhi(t, 0.3, 0, 0, 0) == doctest::Approx(P(t, -1.0)));
    CHECK(Soil_Bgc::sjLayerPhi(t, 1.7, 20, 60, 20) == doctest::Approx(P(t, -1.0)));   // forg clamped
  }
  CHECK(Soil_Bgc::sjLayerPhi(0.3, 0.2, 20, 60, 20) == 1.0);
}

TEST_CASE("sjOrgFraction: SJ2016 Eq. 6, Mo = 2 x C, rho_o = 140 kg m-3, capped at 1") {
  CHECK(Soil_Bgc::sjOrgFraction(7000.0, 0.1) == doctest::Approx(1.0));
  CHECK(Soil_Bgc::sjOrgFraction(3500.0, 0.1) == doctest::Approx(0.5));
  CHECK(Soil_Bgc::sjOrgFraction(700.0, 0.2) == doctest::Approx(0.05));
  CHECK(Soil_Bgc::sjOrgFraction(1e6, 0.1) == 1.0);
  CHECK(Soil_Bgc::sjOrgFraction(-5.0, 0.1) == 0.0);
  CHECK(Soil_Bgc::sjOrgFraction(100.0, 0.0) == 1.0);
}

TEST_CASE("rhFrozenTableDepth, fracAboveDepth, totalSaturation") {
  CHECK(Soil_Bgc::rhFrozenTableDepth(0.0, 0.0, 0.45) == 0.45);
  CHECK(Soil_Bgc::rhFrozenTableDepth(0.5, 0.3, 0.45) == 0.5);
  CHECK(Soil_Bgc::rhFrozenTableDepth(0.2, 0.6, 0.45) == 0.6);
  CHECK(Soil_Bgc::rhFrozenTableDepth(MISSING_D, MISSING_D, MISSING_D) == 0.0);
  CHECK(Soil_Bgc::fracAboveDepth(0.0, 0.1, 0.45) == 1.0);
  CHECK(Soil_Bgc::fracAboveDepth(0.4, 0.1, 0.45) == doctest::Approx(0.5));
  CHECK(Soil_Bgc::fracAboveDepth(0.5, 0.1, 0.45) == 0.0);
  CHECK(Soil_Bgc::fracAboveDepth(0.45, 0.1, 0.45) == 0.0);
  CHECK(Soil_Bgc::fracAboveDepth(0.3, 0.0, 0.45) == 1.0);
  CHECK(Soil_Bgc::totalSaturation(0.2, 0.27, 0.9) == doctest::Approx(0.5));
  CHECK(Soil_Bgc::totalSaturation(0.0, 0.5, 0.45) == 1.0);
  CHECK(Soil_Bgc::totalSaturation(0.3, 0.2, 0.0) == doctest::Approx(0.3));   // porosity unknown: liquid only
  CHECK(Soil_Bgc::totalSaturation(MISSING_D, 0.2, 0.4) == doctest::Approx(0.5));
}

TEST_CASE("thawfracModifier: w_al = 0 is the legacy expression bit for bit; w_al blends the frozen part") {
  Soil_Bgc sb; sb.thawfrac_moist_mode = 0; sb.bgcpar.rhq10 = 2.0; sb.bgcpar.moistmin = 0; sb.bgcpar.moistopt = 0.5; sb.bgcpar.moistmax = 1.0;
  for (double rmf : {0.0, 0.012, 0.2}) for (double f : {0.0, 0.3, 1.0}) for (double ts : {-3.0, -0.2, 0.0, 1.5}) {
    CAPTURE(rmf); CAPTURE(f); CAPTURE(ts);
    sb.bgcpar.rhmoistfrozen = rmf;
    double sh4 = -1, sh6 = -1;
    double e4 = sb.thawfracModifier(f, ts, 0.4, sh4);
    double e6 = sb.thawfracModifier(f, ts, 0.4, sh6, 0.0, 123.0);
    CHECK(e4 == e6); CHECK(sh4 == sh6);
    double swth = (f < 1.0) ? std::fmin(1.0, 0.4 / f) : 0.4;
    double eth = (f > 0) ? M(swth) * Q(ts > 0 ? ts : 0) : 0;
    double efr = rmf * Q(ts < 0 ? ts : 0);
    CHECK(e4 == doctest::Approx(f * eth + (1 - f) * efr));
    double sh = -1;
    double eb = sb.thawfracModifier(f, ts, 0.4, sh, 0.25, 0.03);
    CHECK(eb == doctest::Approx(f * eth + (1 - f) * (0.75 * efr + 0.25 * 0.03)));
    double e1 = sb.thawfracModifier(f, ts, 0.4, sh, 1.0, 0.03);
    CHECK(e1 == doctest::Approx(f * eth + (1 - f) * 0.03));
    CHECK(sh == doctest::Approx(e1 > 0 ? f * eth / e1 : 0.0));
  }
}

TEST_CASE("deltac mode 1: frozen layers above the table respire phi*rhmoist(s_tot)*Q10, below legacy (0)") {
  RfSoil s(1); s.sb.deltac();
  // organic layers 0-2 (z 0..0.15) and mineral layer 3 (0.15..0.25) are above zt = 0.25; layers 4-5 below
  double t = -2.0;
  double eo = P(t, -1.0) * M(0.3 / 0.9) * Q(t);
  for (int i = 0; i <= 2; ++i) { CAPTURE(i); CHECK(s.rh(i) == doctest::Approx(s.kc(i) * eo)); CHECK(s.sb.tfd_thawshare[i] == 0.0); }
  double forg3 = std::fmin(1.0, 2.0 * s.somc(3) / 1000.0 / (0.1 * 140.0));
  double mix = 0.2 * P(t, -0.3) + 0.6 * P(t, -0.5) + 0.2 * P(t, -0.9);
  double em = ((1 - forg3) * mix + forg3 * P(t, -1.0)) * M(0.2 / 0.45) * Q(t);
  CHECK(s.rh(3) == doctest::Approx(s.kc(3) * em));
  CHECK(s.rh(4) == 0.0); CHECK(s.rh(5) == 0.0);
  // partly above: zt = 0.30 -> layer 4 (0.25..0.35) half
  RfSoil h(1); h.cd.d2m_ald_prevyr = 0.30; h.sb.deltac();
  double forg4 = std::fmin(1.0, 2.0 * h.somc(4) / 1000.0 / (0.1 * 140.0));
  double em4 = ((1 - forg4) * mix + forg4 * P(t, -1.0)) * M(0.2 / 0.45) * Q(t);
  CHECK(h.rh(4) == doctest::Approx(0.5 * h.kc(4) * em4));
  CHECK(h.rh(5) == 0.0);
  // legacy rhmoistfrozen below the table is kept
  RfSoil l(1, 0.012); l.sb.deltac();
  CHECK(l.rh(5) == doctest::Approx(l.kc(5) * 0.012 * Q(t)));
  CHECK(l.rh(0) == doctest::Approx(l.kc(0) * eo));   // replaced above the table
}

TEST_CASE("deltac mode 1: saturated frozen part gives 0 (rhmoist(1) = 0), colder gives less, thawed part unchanged") {
  RfSoil a(1); a.ed.m_soid.iwc[0] = 0.9; a.sb.deltac(); CHECK(a.rh(0) == 0.0);
  RfSoil w(1); w.ed.m_sois.ts[0] = -0.5; w.sb.deltac();
  RfSoil c(1); c.ed.m_sois.ts[0] = -10.0; c.sb.deltac();
  CHECK(w.rh(0) > c.rh(0)); CHECK(c.rh(0) > 0.0);
  // partly frozen layer (f = 0.4, ts = 0.5): thawed part as mode 0, frozen part at T = 0 (phi_crit)
  RfSoil p0(0), p1(1);
  for (RfSoil* r : {&p0, &p1}) { r->ed.m_sois.ts[1] = 0.5; r->ed.m_sois.frozenfrac[1] = 0.6; r->ed.m_soid.sws[1] = 0.2; r->sb.deltac(); }
  double eth = M(0.2 / 0.4) * Q(0.5);
  CHECK(p0.rh(1) == doctest::Approx(p0.kc(1) * 0.4 * eth));
  double efr = P(0.0, -1.0) * M(0.2 + 0.3 / 0.9) * Q(0.0);
  CHECK(p1.rh(1) == doctest::Approx(p1.kc(1) * (0.4 * eth + 0.6 * efr)));
  CHECK(p1.sb.tfd_thawshare[1] == doctest::Approx(0.4 * eth / (0.4 * eth + 0.6 * efr)));
}

TEST_CASE("mode 1 with no active layer (zt = 0) and mode 0 are bit-identical (deltac, deltan, deltastate)") {
  for (double rmf : {0.0, 0.012}) {
    CAPTURE(rmf);
    RfSoil a(0, rmf), b(1, rmf);
    for (RfSoil* r : {&a, &b}) { r->cd.d2m_ald_prevyr = 0.0; r->ed.m_sois.ts[2] = 1.0; r->ed.m_sois.frozenfrac[2] = 0.5; r->ed.m_soid.sws[2] = 0.3; r->step(); }
    for (int i = 0; i < a.n; ++i) {
      CAPTURE(i);
      CHECK(a.rh(i) == b.rh(i));
      CHECK(a.sb.del_sois.orgn[i] == b.sb.del_sois.orgn[i]);
      CHECK(a.sb.del_sois.avln[i] == b.sb.del_sois.avln[i]);
      CHECK(a.sb.del_sois.rawc[i] == b.sb.del_sois.rawc[i]);
      CHECK(a.sb.tfd_thawshare[i] == b.sb.tfd_thawshare[i]);
    }
  }
  // mode 0 ignores the table entirely
  RfSoil c(0); c.cd.d2m_ald_prevyr = 5.0; c.sb.deltac();
  for (int i = 0; i < c.n; ++i) { CHECK(c.rh(i) == 0.0); }
}

TEST_CASE("mode 1: frozen-part RH keeps the legacy closed N-mineralisation gate (thawshare 0 -> no net mineralisation)") {
  RfSoil s(1); s.step();
  for (int i = 0; i <= 3; ++i) { CAPTURE(i); CHECK(s.rh(i) > 0.0); CHECK(s.sb.tfd_thawshare[i] == 0.0); CHECK(s.sb.del_soi2soi.netnmin[i] == 0.0); }
}

TEST_CASE("rhFrozenSetLayers without a Ground: zero texture/porosity (organic phi, liquid-only saturation)") {
  RfSoil s(1); s.sb.ground = NULL; s.sb.rhFrozenSetLayers();
  for (int i = 0; i < MAX_SOI_LAY; ++i) { CHECK(s.sb.rhf_por[i] == 0.0); CHECK(s.sb.rhf_texok[i] == 0); }
  s.ed.m_soid.sws[0] = 0.3; s.sb.deltac();
  CHECK(s.rh(0) == doctest::Approx(s.kc(0) * P(-2.0, -1.0) * M(0.3) * Q(-2.0)));
  CHECK(s.rh(3) == doctest::Approx(s.kc(3) * P(-2.0, -1.0) * M(0.0) * Q(-2.0)));
}

}  // TEST_SUITE
