// v5q-layern-fix: dynamic-soil-layer N bookkeeping, model_settings.layer_n_remap_fix
// (0 = legacy, default; 1 = layer orgn/avln follow the C pools through the organic-layer
// split / combine / remap, and combineTwoSoilLayersL2U sums avln instead of assigning it).
// Refs: layer structure and the organic C<->thickness relation of the dynamic organic soil
//   module: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010, doi:10.1029/2010JG001302.
//   The fix is a mass/stoichiometry bookkeeping (conservation) bug fix: no published source
//   for the remap rule. See include/Ground.h, docs/daily-mode-REFERENCES.md.
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/OrganicLayer.h"
#include "../../include/CohortLookup.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value lmsg(Json::Value v) { Json::Value c; c["model_settings"]["layer_n_remap_fix"] = v; return c; }

const double CA = 0.0230, CB = 1.2450;    // CMT51 humic (Yi et al. 2009 form)
const double FA = 0.0130, FB = 1.3750;    // CMT51 fibric

CohortLookup* lookup() {
  static CohortLookup ch;
  ch.poro_f = 0.95; ch.bulkden_f = 60000.; ch.tcsolid_f = 0.25; ch.hksat_f = 0.1;
  ch.poro_h = 0.85; ch.bulkden_h = 150000.; ch.tcsolid_h = 0.25; ch.hksat_h = 0.01;
  return &ch;
}

SoilLayer* org(int type, double dz, double tem) {
  OrganicLayer* l = new OrganicLayer(dz, type, lookup());
  l->tem = tem; l->liq = 0.3 * dz * 1000.; l->ice = 0.0;
  l->frozen = (tem > 0) ? -1 : 1; l->frozenfrac = (tem > 0) ? 0. : 1.;
  return l;
}
void setc(SoilLayer* l, double c, double cn, double avln) {   // C split over the 4 pools
  l->rawc = 0.1 * c; l->soma = 0.2 * c; l->sompr = 0.3 * c; l->somcr = 0.4 * c;
  l->orgn = c / cn; l->avln = avln;
}
double C(const Layer* l) { return l->rawc + l->soma + l->sompr + l->somcr; }

// fibric 0.05 m / humic 0.10 m / humic 0.30 m / mineral 0.5 m; humic C on the Yi curve
struct OrgColumn {
  Ground g; SoilLayer *fib, *h1, *h2, *mn;
  OrgColumn(int fix, double cn1 = 20., double cn2 = 20.) {
    g.layer_n_remap_fix = fix; g.chtlu = lookup();
    g.soildimpar.coefshlwa = FA; g.soildimpar.coefshlwb = FB;
    g.soildimpar.coefdeepa = CA; g.soildimpar.coefdeepb = CB;
    fib = org(1, 0.05, 1.0); h1 = org(2, 0.10, 0.5); h2 = org(2, 0.30, -0.5);
    mn = ut::make_soil_layer(0.45, 0.5, 0.4, 50., 100., -1., 1, 1.);
    g.insertBack(fib); g.insertBack(h1); g.insertBack(h2); g.insertBack(mn);
    g.frontsz.clear(); g.frontstype.clear();
    g.resortGroundLayers(); g.updateSoilHorizons();
    setc(fib, g.carbonFromThickness(0.05, FA, FB), 40., 0.2);
    setc(h1, g.carbonFromThickness(0.10, CA, CB), cn1, 0.5);
    setc(h2, g.carbonFromThickness(0.40, CA, CB) - g.carbonFromThickness(0.10, CA, CB), cn2, 1.5);
    setc(mn, 40000., 20., 2.0);
  }
  double sum(double Layer::*p) const { double s = 0; for (Layer* l = g.fstsoill; l; l = l->nextl) s += l->*p; return s; }
  double sumc() const { double s = 0; for (Layer* l = g.fstsoill; l; l = l->nextl) s += C(l); return s; }
  int nhum() const { int n = 0; for (Layer* l = g.fstsoill; l; l = l->nextl) if (l->isHumic) ++n; return n; }
};
} // namespace

TEST_SUITE("layern_v5q") {

TEST_CASE("layer_n_remap_fix default 0: ModelData(), ModelData(Json) without the key, Ground") {
  ModelData md; CHECK(md.layer_n_remap_fix == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.layer_n_remap_fix == 0);
  Ground g; CHECK(g.layer_n_remap_fix == 0);
}

TEST_CASE("layer_n_remap_fix validation: integers 0/1 only, anything else FATAL") {
  CHECK(ModelData(lmsg(0)).layer_n_remap_fix == 0);
  CHECK(ModelData(lmsg(1)).layer_n_remap_fix == 1);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(lmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("layerNSplitFrac: C share, clamped to [0,1], thickness fallback if no C") {
  CHECK(Ground::layerNSplitFrac(25., 100., 0.9) == doctest::Approx(0.25));
  CHECK(Ground::layerNSplitFrac(0., 100., 0.9) == 0.0);
  CHECK(Ground::layerNSplitFrac(150., 100., 0.9) == 1.0);
  CHECK(Ground::layerNSplitFrac(-5., 100., 0.9) == 0.0);
  CHECK(Ground::layerNSplitFrac(5., 0., 0.3) == doctest::Approx(0.3));
  CHECK(Ground::layerNSplitFrac(std::nan(""), 100., 0.3) == 0.0);
}

TEST_CASE("combineTwoSoilLayersL2U: switch 1 sums avln (column C, orgn, avln conserved); switch 0 = legacy assignment") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    OrgColumn k(fix);
    double c0 = C(k.h1) + C(k.h2), n0 = k.h1->orgn + k.h2->orgn, a1 = k.h1->avln, a2 = k.h2->avln;
    k.g.combineTwoSoilLayersL2U(k.h2, k.h1);
    CHECK(C(k.h1) == doctest::Approx(c0).epsilon(1e-14));
    CHECK(k.h1->orgn == doctest::Approx(n0).epsilon(1e-14));
    if (fix == 1) CHECK(k.h1->avln == doctest::Approx(a1 + a2).epsilon(1e-14));
    else          CHECK(k.h1->avln == a2);                    // legacy '=+' (upper avln lost)
  }
}

TEST_CASE("splitOneSoilLayer, organic: switch 1 N and avln follow C (per-layer C:N kept, sums exact); switch 0 legacy all N to the lower part") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    OrgColumn k(fix);
    k.g.combineTwoSoilLayersL2U(k.h2, k.h1); k.g.removeLayer(k.h2); k.g.resortGroundLayers();
    double c0 = C(k.h1), n0 = k.h1->orgn, a0 = k.h1->avln, cn0 = c0 / n0;
    SoilLayer* lo = org(2, 0.25, -0.5);
    k.g.splitOneSoilLayer(k.h1, lo, 0., 0.25);
    k.g.insertAfter(lo, k.h1);
    CHECK(C(k.h1) > 0.); CHECK(C(lo) > 0.);
    CHECK(C(k.h1) + C(lo) == doctest::Approx(c0).epsilon(1e-13));
    CHECK(k.h1->orgn + lo->orgn == doctest::Approx(n0).epsilon(1e-13));
    CHECK(k.h1->avln + lo->avln == doctest::Approx(a0).epsilon(1e-13));
    if (fix == 1) {
      CHECK(C(k.h1) / k.h1->orgn == doctest::Approx(cn0).epsilon(1e-12));
      CHECK(C(lo) / lo->orgn == doctest::Approx(cn0).epsilon(1e-12));
      CHECK(k.h1->avln / C(k.h1) == doctest::Approx(a0 / c0).epsilon(1e-12));
    } else {
      CHECK(k.h1->orgn == doctest::Approx(0.0).epsilon(1e-12));   // legacy: upper part keeps C, no N
      CHECK(lo->orgn == doctest::Approx(n0).epsilon(1e-14));
    }
  }
}

TEST_CASE("splitOneSoilLayer, mineral: thickness split, identical for switch 0 and 1") {
  double r[2][4];
  for (int fix : {0, 1}) {
    Ground g; g.layer_n_remap_fix = fix; g.chtlu = lookup();
    SoilLayer* a = ut::make_soil_layer(0., 1.0, 0.4, 100., 0., 1., -1, 0.); setc(a, 10000., 20., 3.0);
    g.insertBack(a); g.frontsz.clear(); g.resortGroundLayers();
    SoilLayer* b = ut::make_soil_layer(0., 0.4, 0.4, 0., 0., 1., -1, 0.);
    g.splitOneSoilLayer(a, b, 0., 0.4);
    r[fix][0] = C(a); r[fix][1] = C(b); r[fix][2] = b->orgn; r[fix][3] = b->avln;
    CHECK(C(b) == doctest::Approx(4000.)); CHECK(b->orgn == doctest::Approx(200.)); CHECK(b->avln == doctest::Approx(1.2));
    CHECK(a->orgn + b->orgn == doctest::Approx(500.));
  }
  for (int i = 0; i < 4; ++i) CHECK(r[0][i] == r[1][i]);
}

TEST_CASE("redivideDeepLayers (combine + re-split, 2 -> 3 humic layers): column C/N/avln exact; switch 1 per-layer C:N follows C") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    OrgColumn k(fix, 20., 20.);
    double c0 = k.sumc(), n0 = k.sum(&Layer::orgn), a0 = k.sum(&Layer::avln);
    k.g.redivideDeepLayers();
    CHECK(k.nhum() == 3);
    CHECK(k.sumc() == doctest::Approx(c0).epsilon(1e-13));
    CHECK(k.sum(&Layer::orgn) == doctest::Approx(n0).epsilon(1e-13));
    if (fix == 1) {
      CHECK(k.sum(&Layer::avln) == doctest::Approx(a0).epsilon(1e-13));
      for (Layer* l = k.g.fstsoill; l; l = l->nextl) if (l->isHumic) {
        CHECK(C(l) > 0.); CHECK(C(l) / l->orgn == doctest::Approx(20.).epsilon(1e-12));
      }
    } else {
      CHECK(k.sum(&Layer::avln) == doctest::Approx(a0 - 0.5).epsilon(1e-13));   // legacy loses the upper avln
      int zeron = 0; for (Layer* l = k.g.fstsoill; l; l = l->nextl) if (l->isHumic && C(l) > 0. && l->orgn < 1e-9) ++zeron;
      CHECK(zeron == 2);   // legacy: the two upper humic layers keep C but have N = 0
    }
  }
}

TEST_CASE("new deep layer formed from the last fibric layer: switch 1 N moves with the C moved; switch 0 legacy all N moves") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    Ground g; g.layer_n_remap_fix = fix; g.chtlu = lookup();
    g.soildimpar.coefshlwa = FA; g.soildimpar.coefshlwb = FB; g.soildimpar.coefdeepa = CA; g.soildimpar.coefdeepb = CB;
    SoilLayer* fib = org(1, 0.20, 1.0);
    SoilLayer* mn = ut::make_soil_layer(0.2, 0.5, 0.4, 50., 0., 1., -1, 0.);
    g.insertBack(fib); g.insertBack(mn); g.frontsz.clear(); g.resortGroundLayers(); g.updateSoilHorizons();
    setc(fib, 20000., 25., 0.8); setc(mn, 40000., 20., 2.0);
    double cf = C(fib), nf = fib->orgn, af = fib->avln;
    double moved = 0.5 * fib->soma + fib->sompr + fib->somcr;
    g.redivideDeepLayers();
    REQUIRE(g.fstdeepl != NULL);
    Layer* d = g.fstdeepl;
    CHECK(C(fib) + C(d) == doctest::Approx(cf).epsilon(1e-14));
    CHECK(C(d) == doctest::Approx(moved).epsilon(1e-14));
    CHECK(fib->orgn + d->orgn == doctest::Approx(nf).epsilon(1e-14));
    CHECK(fib->avln + d->avln == doctest::Approx(af).epsilon(1e-14));
    if (fix == 1) {
      CHECK(C(d) / d->orgn == doctest::Approx(25.).epsilon(1e-12));
      CHECK(C(fib) / fib->orgn == doctest::Approx(25.).epsilon(1e-12));
    } else {
      CHECK(fib->orgn == 0.0); CHECK(d->orgn == nf);
    }
  }
}

TEST_CASE("Soil_Bgc bd<->layer chain around redivideDeepLayers: switch 1 bd C:N per index follows C; legacy (C only) desyncs N by index") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    OrgColumn k(fix, 20., 20.);
    BgcData bd; Soil_Bgc sb; sb.setGround(&k.g); sb.setBgcData(&bd);
    for (int i = 0; i < MAX_SOI_LAY; ++i) { bd.m_sois.rawc[i] = bd.m_sois.soma[i] = bd.m_sois.sompr[i] = bd.m_sois.somcr[i] = 0; bd.m_sois.orgn[i] = bd.m_sois.avln[i] = 0; }
    sb.assignCarbonLayer2BdMonthly(); sb.assignNitrogenLayer2BdMonthly();   // bd = initial column (4 layers)
    bd.m_soid.shlwc = 0; bd.m_soid.deepc = 0;
    double c0 = 0, n0 = 0, a0 = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) { c0 += bd.m_sois.rawc[i] + bd.m_sois.soma[i] + bd.m_sois.sompr[i] + bd.m_sois.somcr[i]; n0 += bd.m_sois.orgn[i]; a0 += bd.m_sois.avln[i]; }
    // as Cohort::updateMonthly_DIMgrd
    sb.assignCarbonBd2LayerMonthly(); if (fix == 1) sb.assignNitrogenBd2LayerMonthly();
    k.g.redivideDeepLayers();
    sb.assignCarbonLayer2BdMonthly(); if (fix == 1) sb.assignNitrogenLayer2BdMonthly();
    double c1 = 0, n1 = 0, a1 = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) { c1 += bd.m_sois.rawc[i] + bd.m_sois.soma[i] + bd.m_sois.sompr[i] + bd.m_sois.somcr[i]; n1 += bd.m_sois.orgn[i]; a1 += bd.m_sois.avln[i]; }
    CHECK(c1 == doctest::Approx(c0).epsilon(1e-13));
    CHECK(n1 == doctest::Approx(n0).epsilon(1e-13));   // legacy also keeps the column N (stays by index)
    CHECK(a1 == doctest::Approx(a0).epsilon(1e-13));
    // index 4 is now the 3rd humic layer (5 layers: fib, h, h, h, mineral)
    double cn_min = (bd.m_sois.rawc[4] + bd.m_sois.soma[4] + bd.m_sois.sompr[4] + bd.m_sois.somcr[4]) / bd.m_sois.orgn[4];
    if (fix == 1) {
      for (int i = 1; i <= 3; ++i) {
        double ci = bd.m_sois.rawc[i] + bd.m_sois.soma[i] + bd.m_sois.sompr[i] + bd.m_sois.somcr[i];
        CHECK(ci / bd.m_sois.orgn[i] == doctest::Approx(20.).epsilon(1e-12));
      }
      CHECK(cn_min == doctest::Approx(20.).epsilon(1e-12));   // mineral layer N moved to index 4 with its C
      CHECK(bd.m_sois.orgn[5] == 0.0);
    } else {
      CHECK(bd.m_sois.orgn[4] == 0.0);   // legacy: mineral C now at index 4, its N still at index 3
      CHECK(bd.m_sois.orgn[3] == doctest::Approx(2000.).epsilon(1e-14));
    }
  }
}

} // TEST_SUITE layern_v5q
