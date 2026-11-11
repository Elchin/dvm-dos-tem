// v5q-humic-stefan: three independent switches, all default off (bit-identical to v5q-thawfrac-lmoist):
//  model_settings.humic_max_thick (m; 0 = off): max thickness of a deep (humic) sublayer, equal
//    re-split within MAX_DEP_LAY = 3 (Organic::applyHumicMaxThick); cannot be met if humic > 3 x cap.
//  model_settings.stefan_front_resist_fix (0/1): thawing front, only the thawed part of the front
//    layer counts as resistance above the front (Stefan::frontLayerResistance).
//  output ALD_FRONT (output only): yearly max of the daily deepest thawing front, not
//    permafrost-gated (EnvData m_soid/y_soid.aldfront).
// Refs: layering / organic C<->thickness: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010,
//   doi:10.1029/2010JG001302. Modified Stefan front tracking: Woo et al. 2004,
//   doi:10.1029/2004GL019475. The cap and the resistance fix are discretization fixes: no published source.
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/OrganicLayer.h"
#include "../../include/CohortLookup.h"
#include <json/json.h>
#include <cmath>
#include <cstring>

using ut::make_soil_layer;

namespace {
Json::Value msg(const char* key, Json::Value v) { Json::Value c; c["model_settings"][key] = v; return c; }

const double CA = 0.0230, CB = 1.2450;    // CMT51 humic (Yi et al. 2009 form)
const double FA = 0.0130, FB = 1.3750;    // CMT51 fibric
CohortLookup* lookup() {
  static CohortLookup ch;
  ch.poro_f = 0.95; ch.bulkden_f = 60000.; ch.tcsolid_f = 0.25; ch.hksat_f = 0.1;
  ch.poro_h = 0.85; ch.bulkden_h = 150000.; ch.tcsolid_h = 0.25; ch.hksat_h = 0.01;
  for (int pp = 0; pp < NUM_PFT_PART; ++pp) for (int ip = 0; ip < NUM_PFT; ++ip) ch.cfall[pp][ip] = 0.;
  for (int ip = 0; ip < NUM_PFT; ++ip) ch.nonvascular[ip] = 0;
  return &ch;
}
SoilLayer* org(int type, double dz, double tem, double liqv, double icev) {
  OrganicLayer* l = new OrganicLayer(dz, type, lookup());
  l->tem = tem; l->liq = liqv * dz * DENLIQ; l->ice = icev * dz * DENICE;
  l->frozen = (tem > 0) ? -1 : 1; l->frozenfrac = (tem > 0) ? 0. : 1.;
  return l;
}
void setc(Layer* l, double c, double cn, double avln) {
  l->rawc = 0.1 * c; l->soma = 0.2 * c; l->sompr = 0.3 * c; l->somcr = 0.4 * c;
  l->orgn = c / cn; l->avln = avln;
}
double C(const Layer* l) { return l->rawc + l->soma + l->sompr + l->somcr; }
bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

// legacy reference split (copy of the 6630518f DeepThickScheme thresholds, for the "off" check)
void legacy_split(double lst, double th, double dz[3], int& n) {
  if (th <= 0.) { dz[0] = dz[1] = dz[2] = MISSING_D; n = 0; return; }
  double t1 = (lst > 0.) ? 3 * lst : 0.02, t2 = (lst > 0.) ? 6 * lst : 0.06;
  bool one = (lst > 0.) ? (th < t1) : (th <= t1), two = (lst > 0.) ? (th < t2) : (th <= t2);
  if (one) { dz[0] = th; dz[1] = dz[2] = MISSING_D; n = 1; }
  else if (two) { dz[0] = 1./3. * th; dz[1] = 2./3. * th; dz[2] = MISSING_D; n = 2; }
  else { dz[0] = 1./6. * th; dz[1] = 2./6. * th; dz[2] = 3./6. * th; n = 3; }
}
} // namespace

TEST_SUITE("humic_stefan_v5q") {

TEST_CASE("defaults: humic_max_thick 0 and stefan_front_resist_fix 0 (ModelData(), ModelData(Json) without keys, Organic, Stefan)") {
  ModelData md; CHECK(md.humic_max_thick == 0.0); CHECK(md.stefan_front_resist_fix == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c);
  CHECK(mj.humic_max_thick == 0.0); CHECK(mj.stefan_front_resist_fix == 0);
  Organic o; CHECK(o.humic_max_thick == 0.0);
  Ground g; CHECK(g.organic.humic_max_thick == 0.0);
  CHECK(Stefan::front_resist_fix == 0);
}

TEST_CASE("humic_max_thick validation: 0 or a finite number in [0.02, 5.0] m, anything else FATAL") {
  CHECK(ModelData(msg("humic_max_thick", 0)).humic_max_thick == 0.0);
  CHECK(ModelData(msg("humic_max_thick", 0.0)).humic_max_thick == 0.0);
  CHECK(ModelData(msg("humic_max_thick", 0.15)).humic_max_thick == 0.15);
  CHECK(ModelData(msg("humic_max_thick", 0.02)).humic_max_thick == 0.02);
  CHECK(ModelData(msg("humic_max_thick", 5)).humic_max_thick == 5.0);
  for (Json::Value bad : {Json::Value(0.019), Json::Value(-0.15), Json::Value(5.01), Json::Value(1e300),
                          Json::Value(true), Json::Value("0.15"), Json::Value(Json::nullValue)}) {
    CAPTURE(bad.toStyledString());
    CHECK(ut::exit_status_of([bad] { ModelData m(msg("humic_max_thick", bad)); }) == EXIT_FAILURE);
  }
}

TEST_CASE("stefan_front_resist_fix validation: integers 0/1 only, anything else FATAL") {
  CHECK(ModelData(msg("stefan_front_resist_fix", 0)).stefan_front_resist_fix == 0);
  CHECK(ModelData(msg("stefan_front_resist_fix", 1)).stefan_front_resist_fix == 1);
  CHECK(ModelData(msg("stefan_front_resist_fix", 1.0)).stefan_front_resist_fix == 1);
  for (Json::Value bad : {Json::Value(2), Json::Value(-1), Json::Value(0.5), Json::Value(true), Json::Value("1")}) {
    CAPTURE(bad.toStyledString());
    CHECK(ut::exit_status_of([bad] { ModelData m(msg("stefan_front_resist_fix", bad)); }) == EXIT_FAILURE);
  }
}

// ---------------- humic_max_thick ----------------
TEST_CASE("DeepThickScheme, cap off (0) and cap >= every legacy sublayer: bit-identical to the legacy split") {
  for (double lst : {-9999., 0., 0.02, 0.05, 0.124}) for (int i = 0; i <= 600; ++i) {
    double th = 0.0025 * i; CAPTURE(lst); CAPTURE(th);
    double rdz[3]; int rn; legacy_split(lst, th, rdz, rn);
    for (double cap : {0.0, 5.0}) {
      Organic o; o.lstshlwdz = lst; o.humic_max_thick = cap; o.DeepThickScheme(th);
      if (cap == 5.0 && th > 5.0 * 3) continue;
      REQUIRE(o.deepnum == rn);
      for (int k = 0; k < 3; ++k) CHECK(same_bits(o.deepdz[k], rdz[k]));
    }
  }
}

TEST_CASE("DeepThickScheme, cap on: n <= 3, n >= legacy n, sum dz == deepthick, max dz <= cap when humic <= 3 cap, equal split when changed") {
  int changed = 0, limit = 0;
  for (double cap : {0.05, 0.10, 0.15, 0.20}) for (double lst : {-9999., 0.02, 0.05, 0.124}) for (int i = 1; i <= 600; ++i) {
    double th = 0.0025 * i; CAPTURE(cap); CAPTURE(lst); CAPTURE(th);
    double rdz[3]; int rn; legacy_split(lst, th, rdz, rn);
    Organic o; o.lstshlwdz = lst; o.humic_max_thick = cap; o.DeepThickScheme(th);
    REQUIRE(o.deepnum >= rn); REQUIRE(o.deepnum <= MAX_DEP_LAY);
    double s = 0., mx = 0.;
    for (int k = 0; k < o.deepnum; ++k) { CHECK(o.deepdz[k] > 0.); s += o.deepdz[k]; mx = std::max(mx, o.deepdz[k]); }
    for (int k = o.deepnum; k < 3; ++k) CHECK(o.deepdz[k] == MISSING_D);
    CHECK(s == doctest::Approx(th).epsilon(1e-14));
    double rmx = 0.; for (int k = 0; k < rn; ++k) rmx = std::max(rmx, rdz[k]);
    if (rmx <= cap) {                       // already satisfied: unchanged
      CHECK(o.deepnum == rn); for (int k = 0; k < 3; ++k) CHECK(same_bits(o.deepdz[k], rdz[k]));
    } else {
      ++changed;
      for (int k = 0; k < o.deepnum; ++k) CHECK(o.deepdz[k] == doctest::Approx(th / o.deepnum).epsilon(1e-15));
      CHECK(mx <= rmx);
      if (th <= 3 * cap * (1 + 1e-12)) { CHECK(mx <= cap * (1 + 1e-12)); CHECK(o.deepnum == std::max(rn, (int)std::ceil(th / cap - 1e-12))); }
      else { ++limit; CHECK(o.deepnum == 3); CHECK(mx > cap); }   // documented limit: humic > 3 x cap
    }
  }
  CHECK(changed > 100); CHECK(limit > 100);
}

TEST_CASE("TVC-like humic: 0.34 m with last fibric 0.124 m -> legacy 1 layer; cap 0.15 -> 3 x 0.1133; cap 0.20 -> 2 x 0.17") {
  Organic a; a.lstshlwdz = 0.124; a.DeepThickScheme(0.34);
  CHECK(a.deepnum == 1); CHECK(a.deepdz[0] == 0.34);
  Organic b; b.lstshlwdz = 0.124; b.humic_max_thick = 0.15; b.DeepThickScheme(0.34);
  CHECK(b.deepnum == 3); for (int k = 0; k < 3; ++k) CHECK(b.deepdz[k] == doctest::Approx(0.34 / 3));
  Organic c; c.lstshlwdz = 0.124; c.humic_max_thick = 0.20; c.DeepThickScheme(0.34);
  CHECK(c.deepnum == 2); CHECK(c.deepdz[0] == doctest::Approx(0.17)); CHECK(c.deepdz[1] == doctest::Approx(0.17));
  CHECK(c.deepdz[2] == MISSING_D);
  Organic d; d.lstshlwdz = 0.124; d.humic_max_thick = 0.10; d.DeepThickScheme(0.60);   // limit: 0.60 > 3 x 0.10
  CHECK(d.deepnum == 3); for (int k = 0; k < 3; ++k) CHECK(d.deepdz[k] == doctest::Approx(0.20));
}

TEST_CASE("applyHumicMaxThick no-op guards: cap <= 0 or NaN, deepthick <= 0, deepnum 0") {
  double dz[3] = {0.3, MISSING_D, MISSING_D}; int n = 1;
  CHECK(!Organic::applyHumicMaxThick(0.0, 0.3, dz, n)); CHECK(!Organic::applyHumicMaxThick(-1., 0.3, dz, n));
  CHECK(!Organic::applyHumicMaxThick(std::nan(""), 0.3, dz, n)); CHECK(!Organic::applyHumicMaxThick(0.1, 0.0, dz, n));
  int n0 = 0; CHECK(!Organic::applyHumicMaxThick(0.1, 0.3, dz, n0));
  CHECK(n == 1); CHECK(dz[0] == 0.3);
  CHECK(Organic::applyHumicMaxThick(0.1, 0.3, dz, n)); CHECK(n == 3);
}

namespace {
// fibric 0.124 / humic 0.34 / mineral 0.464 column; humic C on the Yi et al. (2009) C-thickness curve
// (as in the model, where the humic thickness is derived from its C). tf, tm = fibric / mineral T.
struct DeepCol {
  Ground g; SoilLayer *fib, *hum, *mn;
  DeepCol(double cap, double tem, double tf, double tm) {
    g.layer_n_remap_fix = 1; g.chtlu = lookup();
    g.soildimpar.coefshlwa = FA; g.soildimpar.coefshlwb = FB; g.soildimpar.coefdeepa = CA; g.soildimpar.coefdeepb = CB;
    g.organic.humic_max_thick = cap;
    fib = org(1, 0.124, tf, 0.15, tf > 0 ? 0.0 : 0.5);
    hum = org(2, 0.34, tem, tem > 0 ? 0.6 : 0.02, tem > 0 ? 0.0 : 0.7);
    mn = make_soil_layer(0.464, 0.5, 0.4, tm > 0 ? 150. : 0., tm > 0 ? 0. : 150., tm, tm > 0 ? -1 : 1, tm > 0 ? 0. : 1.);
    g.insertBack(fib); g.insertBack(hum); g.insertBack(mn);
    g.frontsz.clear(); g.frontstype.clear(); g.resortGroundLayers(); g.updateSoilHorizons();
    setc(fib, g.carbonFromThickness(0.124, FA, FB), 40., 0.2);
    setc(hum, g.carbonFromThickness(0.34, CA, CB), 25., 0.7);
    setc(mn, 40000., 20., 2.0);
  }
  double tot(int what) { double s = 0; for (Layer* l = g.fstsoill; l; l = l->nextl)
    s += what == 0 ? C(l) : what == 1 ? l->orgn : what == 2 ? l->avln : what == 3 ? l->liq + l->ice : what == 4 ? l->ice : l->dz;
    return s; }
};
}

TEST_CASE("redivideDeepLayers with the cap, uniform T (fibric 0.124 / humic 0.34 / mineral): C, N, avln, water, ice, dz, T conserved; C:N kept; off = 1 layer") {
  for (double cap : {0.0, 0.15, 0.20}) for (double tem : {1.5, -0.5}) {
    CAPTURE(cap); CAPTURE(tem);
    DeepCol d(cap, tem, tem, tem);
    double t0[6]; for (int w = 0; w < 6; ++w) t0[w] = d.tot(w);
    d.g.redivideDeepLayers();
    for (int w = 0; w < 6; ++w) { CAPTURE(w); CHECK(d.tot(w) == doctest::Approx(t0[w]).epsilon(1e-12)); }
    int nh = 0; double hz = 0;
    for (Layer* l = d.g.fstsoill; l; l = l->nextl) if (l->isHumic) {
      ++nh; hz += l->dz;
      CHECK(l->tem == doctest::Approx(tem).epsilon(1e-12));   // heat: same T, heat capacity ~ dz
      CHECK(C(l) > 0.); CHECK(C(l) / l->orgn == doctest::Approx(25.).epsilon(1e-10));
      if (cap > 0) CHECK(l->dz <= cap + 1e-12);
    }
    CHECK(nh == (cap == 0.0 ? 1 : cap == 0.15 ? 3 : 2));
    CHECK(hz == doctest::Approx(0.34).epsilon(1e-14));
  }
}

TEST_CASE("redivideDeepLayers with the cap, T gradient with a consistent thaw front: C, N, avln, water, ice, dz conserved, C:N kept; T interpolated (legacy split)") {
  // The cap re-uses the legacy split (Ground::splitOneSoilLayer): T linearly interpolated from the
  // neighbours, ice/liq placed by the front-derived frozen fraction. The same code splits any legacy
  // humic thicker than 3 x the last fibric layer; nothing new is introduced here.
  // Humic -0.5 C: thawed fibric (2 C) above a thaw front at 0.124 m. Humic +1.5 C: front at 0.464 m.
  for (double cap : {0.15, 0.20}) for (double tem : {1.5, -0.5}) {
    CAPTURE(cap); CAPTURE(tem);
    DeepCol d(cap, tem, 2.0, -1.0);
    d.g.frontsz.push_back(tem > 0 ? 0.464 : 0.124); d.g.frontstype.push_back(-1);
    double t0[6]; for (int w = 0; w < 6; ++w) t0[w] = d.tot(w);
    d.g.redivideDeepLayers();
    for (int w = 0; w < 6; ++w) { CAPTURE(w); CHECK(d.tot(w) == doctest::Approx(t0[w]).epsilon(1e-12)); }
    for (Layer* l = d.g.fstsoill; l; l = l->nextl) if (l->isHumic) {
      CHECK(C(l) > 0.); CHECK(C(l) / l->orgn == doctest::Approx(25.).epsilon(1e-10));
      CHECK(std::fabs(l->frozenfrac - (tem > 0 ? 0.0 : 1.0)) < 1e-12);   // front at a layer boundary: frozen flag may be 0
    }
  }
}

TEST_CASE("legacy splitOneSoilLayer property (documented, not changed): with no front and T neighbours of the other sign, a frozen humic split can lose water") {
  // Pre-existing behaviour of the legacy split (frozen state from the interpolated T when there is no
  // front; fmin(ice, maxice) clipping). Reached here only with an inconsistent column (thawed fibric
  // and a frozen humic but no front). Reported so the closure checks in the runs look for it.
  DeepCol d(0.15, -0.5, 2.0, -1.0);
  double w0 = d.tot(3);
  d.g.redivideDeepLayers();
  CHECK(d.tot(3) < w0 - 1.0);                     // documents the legacy loss (not a cap property)
  MESSAGE("no-front inconsistent column, cap 0.15: water change " << d.tot(3) - w0 << " kg/m2 (legacy split)");
}

// ---------------- stefan_front_resist_fix ----------------
TEST_CASE("frontLayerResistance: fix 0 is the legacy dz/tk (bitwise) in every case; fix 1 only for a thawing soil front layer") {
  for (int fix : {0, 1}) for (double tdrv : {-3., 0., 4.}) for (bool soil : {true, false})
  for (int frozen : {-1, 0, 1}) for (double ff : {0., 0.3, 1., -0.1, 1.2, std::nan("")}) {
    CAPTURE(fix); CAPTURE(tdrv); CAPTURE(soil); CAPTURE(frozen); CAPTURE(ff);
    const int freezing = (tdrv < 0) ? 1 : -1;
    const double dz = 0.34, tk = 0.47;
    double now = -1, after = -1;
    Stefan::frontLayerResistance(fix, tdrv, soil, frozen, freezing, ff, dz, tk, now, after);
    bool applies = fix == 1 && tdrv > 0 && soil && frozen != freezing && ff >= 0. && ff <= 1.;
    if (!applies) { CHECK(same_bits(now, dz / tk)); CHECK(after == 0.0); }
    else {
      CHECK(now == doctest::Approx((1. - ff) * dz / tk).epsilon(1e-15));
      CHECK(after == doctest::Approx(ff * dz / tk).epsilon(1e-15));
      CHECK(now + after == doctest::Approx(dz / tk).epsilon(1e-15));   // full layer once passed
      if (ff == 1.) CHECK(now == 0.0);                                 // fully frozen layer: no resistance above its front
    }
  }
}

namespace {
// column: 0.1 m thawed mineral / X / 0.5 m frozen mineral; X = one 0.34 m layer half-thawed (front at 0.27)
// or two 0.17 m layers (thawed / frozen). Same volumetric water (0.36) everywhere.
struct ThawCol {
  ut::Column* c; Stefan st;
  ThawCol(bool split, double tem_front = 0.0) {
    const double p = 0.4, w = 0.9 * p;
    std::vector<SoilLayer*> v;
    v.push_back(make_soil_layer(0.0, 0.1, p, w * DENLIQ * 0.1, 0., 2.0, -1, 0.));
    if (split) {
      v.push_back(make_soil_layer(0.1, 0.17, p, w * DENLIQ * 0.17, 0., 0.5, -1, 0.));
      v.push_back(make_soil_layer(0.27, 0.17, p, 0., w * DENICE * 0.17, tem_front, 1, 1.));
    } else {
      v.push_back(make_soil_layer(0.1, 0.34, p, w * DENLIQ * 0.17, w * DENICE * 0.17, tem_front, 0, 0.5));
    }
    v.push_back(make_soil_layer(0.44, 0.5, p, 0., w * DENICE * 0.5, -1.0, 1, 1.));
    c = new ut::Column(v);
    c->g.frontsz.push_back(0.27); c->g.frontstype.push_back(-1);
    st.setGround(&c->g); st.initpce();
  }
  ~ThawCol() { delete c; }
  double front() const { return c->g.frontsz.empty() ? -1. : c->g.frontsz.front(); }
};
}

TEST_CASE("Stefan thaw: fix 1 in a half-thawed 0.34 m layer == two 0.17 m layers split at the front (same latent heat), deeper than legacy; water conserved") {
  const int saved = Stefan::front_resist_fix;
  double fr[2][2];   // [fix][split]
  for (int fix : {0, 1}) for (int split : {0, 1}) {
    CAPTURE(fix); CAPTURE(split);
    Stefan::front_resist_fix = fix;
    ThawCol t(split == 1);
    double w0 = t.c->total_water();
    t.st.updateFronts(5.0, 86400.);
    fr[fix][split] = t.front();
    CHECK(t.c->total_water() == doctest::Approx(w0).epsilon(1e-13));
    CHECK(fr[fix][split] > 0.27); CHECK(fr[fix][split] < 0.44);
  }
  Stefan::front_resist_fix = saved;
  CHECK(fr[1][0] == doctest::Approx(fr[1][1]).epsilon(1e-10));   // thawed-part resistance == split column
  CHECK(fr[1][0] > fr[0][0]);                                    // legacy counts the frozen half too
  CHECK(fr[1][1] > fr[0][1]);                                    // and the whole frozen layer below the front
  MESSAGE("front after 1 day: legacy one-layer " << fr[0][0] << ", fix " << fr[1][0] << ", fix split " << fr[1][1]);
}

TEST_CASE("Stefan freezing direction: fix 1 leaves the result bit-identical to legacy") {
  const int saved = Stefan::front_resist_fix;
  double z[2]; double ice[2][3]; double liq[2][3];
  for (int fix : {0, 1}) {
    Stefan::front_resist_fix = fix;
    ThawCol t(false, 0.0);          // 3 layers: 0.1 thawed / 0.34 half-thawed / 0.5 frozen
    REQUIRE(t.c->L.size() == 3);
    double w0 = t.c->total_water();
    t.st.updateFronts(-4.0, 86400.);
    CHECK(t.c->total_water() == doctest::Approx(w0).epsilon(1e-13));
    z[fix] = t.c->g.frontsz.empty() ? -1 : t.c->g.frontsz.front();
    for (int i = 0; i < 3; ++i) { ice[fix][i] = t.c->L[i]->ice; liq[fix][i] = t.c->L[i]->liq; }
  }
  Stefan::front_resist_fix = saved;
  CHECK(same_bits(z[0], z[1]));
  for (int i = 0; i < 3; ++i) { CHECK(same_bits(ice[0][i], ice[1][i])); CHECK(same_bits(liq[0][i], liq[1][i])); }
}

TEST_CASE("Stefan thaw through a whole thin frozen layer in one step: fix 1 >= legacy, water conserved, later layers see the full layer") {
  const int saved = Stefan::front_resist_fix;
  double fz[2];
  for (int fix : {0, 1}) {
    Stefan::front_resist_fix = fix;
    const double p = 0.4, w = 0.9 * p;
    ut::Column c({make_soil_layer(0.0, 0.01, p, w * DENLIQ * 0.01, 0., 3.0, -1, 0.),
                  make_soil_layer(0.01, 0.01, p, 0., w * DENICE * 0.01, 0.0, 1, 1.),
                  make_soil_layer(0.02, 0.5, p, 0., w * DENICE * 0.5, 0.0, 1, 1.)});
    c.g.frontsz.push_back(0.01); c.g.frontstype.push_back(-1);
    Stefan st; st.setGround(&c.g); st.initpce();
    double w0 = c.total_water();
    st.updateFronts(15.0, 86400.);
    fz[fix] = c.g.frontsz.empty() ? -1 : c.g.frontsz.front();
    CHECK(fz[fix] > 0.02);                                  // passed the thin layer
    CHECK(c.L[1]->ice == doctest::Approx(0.0));
    CHECK(c.total_water() == doctest::Approx(w0).epsilon(1e-13));
  }
  Stefan::front_resist_fix = saved;
  CHECK(fz[1] >= fz[0]);
}

// ---------------- ALD_FRONT ----------------
namespace {
struct AldEnv {
  CohortData cd; EnvData ed;
  AldEnv() { ed.cd = &cd; cd.d_snow.numsnwl = 0; cd.m_soil.totthick = 2.5; cd.d_soil.totthick = 2.5;
             ed.grnd_beginOfYear(); ed.grnd_beginOfMonth(); }
  void month(const std::vector<double>& daily, int frozen_flag) {
    ed.grnd_beginOfMonth();
    for (size_t d = 0; d < daily.size(); ++d) {
      ed.d_soid.permafrost = frozen_flag; ed.d_soid.aldflag = daily[d]; ed.d_soid.ald = 0.3;
      ed.grnd_endOfDay((int)daily.size(), 180 + (int)d);
    }
    ed.grnd_endOfMonth(5);   // merge glue: upstream added the month index (CH4 dfratio only)
  }
};
}

TEST_CASE("ALD_FRONT aggregation: monthly/yearly = max of daily deepest thawing front (MISSING_D ignored), not permafrost-gated") {
  AldEnv a; a.ed.monthsfrozen = 0.;          // never permafrost: ALD_FLAG = totthick, ALD_FRONT = front max
  a.month({MISSING_D, 0.2, 0.35, 0.31}, 1);
  CHECK(a.ed.m_soid.aldfront == 0.35); CHECK(a.ed.y_soid.aldfront == 0.35);
  CHECK(a.ed.m_soid.aldflag == 2.5); CHECK(a.ed.y_soid.aldflag == 2.5);
  a.month({0.4, 0.47, MISSING_D}, 1);
  CHECK(a.ed.m_soid.aldfront == 0.47); CHECK(a.ed.y_soid.aldfront == 0.47);
  AldEnv m; m.ed.grnd_beginOfMonth();
  CHECK(m.ed.m_soid.aldfront == MISSING_D); CHECK(m.ed.y_soid.aldfront == MISSING_D);
}

TEST_CASE("ALD_FRONT == ALD_FLAG in permafrost years; ALD_FLAG / ALD unchanged by the new diagnostic") {
  AldEnv a; a.ed.monthsfrozen = 30.;         // >= 24 frozen months: permafrost
  a.month({0.1, 0.42, MISSING_D, 0.38}, 1);
  a.month({0.45, 0.5, 0.49}, 1);
  CHECK(a.ed.m_soid.permafrost == 1);
  CHECK(a.ed.y_soid.aldfront == a.ed.y_soid.aldflag); CHECK(a.ed.y_soid.aldflag == 0.5);
  CHECK(a.ed.m_soid.aldfront == a.ed.m_soid.aldflag);
  // the legacy diagnostics are identical whatever aldfront holds (it is never read back)
  AldEnv b; b.ed.monthsfrozen = 30.; b.ed.m_soid.aldfront = 123.; b.ed.y_soid.aldfront = 123.;
  b.month({0.1, 0.42, MISSING_D, 0.38}, 1); b.month({0.45, 0.5, 0.49}, 1);
  CHECK(same_bits(a.ed.y_soid.aldflag, b.ed.y_soid.aldflag)); CHECK(same_bits(a.ed.y_soid.ald, b.ed.y_soid.ald));
  CHECK(same_bits(a.ed.m_soid.aldflag, b.ed.m_soid.aldflag)); CHECK(same_bits(a.ed.y_soid.alc, b.ed.y_soid.alc));
}

} // TEST_SUITE
