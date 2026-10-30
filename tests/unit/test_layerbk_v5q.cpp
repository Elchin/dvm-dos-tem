// v5q-layer-bookkeep: dynamic-soil-layer C/N bookkeeping, model_settings.layer_bookkeep_fix
// (0 = legacy, default; 1 = four conservation bug fixes, see include/Ground.h):
//  (1) new fibric layer created from the layer below (C/N moved, not created)
//  (2) deepnum == 0: removed humic layers' C/N go to the next soil layer (not lost)
//  (3) prvltrfcnque (litterfall C:N history) follows its layer through the re-division
//  (4) organic split: lower part gets at most the whole layer's C (upper C >= 0)
// Refs: layer structure and the organic C<->thickness relation of the dynamic organic soil
//   module: Yi et al. 2009, doi:10.1029/2008JG000841; Yi et al. 2010, doi:10.1029/2010JG001302.
//   The fixes are mass-conservation bookkeeping bug fixes: no published source.
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/OrganicLayer.h"
#include "../../include/CohortLookup.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value bmsg(Json::Value v, int ln = 1) {
  Json::Value c; c["model_settings"]["layer_bookkeep_fix"] = v; c["model_settings"]["layer_n_remap_fix"] = ln; return c;
}
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
SoilLayer* org(int type, double dz, double tem) {
  OrganicLayer* l = new OrganicLayer(dz, type, lookup());
  l->tem = tem; l->liq = 0.3 * dz * 1000.; l->ice = 0.0;
  l->frozen = (tem > 0) ? -1 : 1; l->frozenfrac = (tem > 0) ? 0. : 1.;
  return l;
}
void setc(Layer* l, double c, double cn, double avln) {
  l->rawc = 0.1 * c; l->soma = 0.2 * c; l->sompr = 0.3 * c; l->somcr = 0.4 * c;
  l->orgn = (cn > 0) ? c / cn : 0.; l->avln = avln;
}
double C(const Layer* l) { return l->rawc + l->soma + l->sompr + l->somcr; }
bool nonneg(const Layer* l) {
  return l->rawc >= 0 && l->soma >= 0 && l->sompr >= 0 && l->somcr >= 0 && l->orgn >= 0 && l->avln >= 0;
}
void setup(Ground& g, int fix) {
  g.layer_n_remap_fix = 1; g.layer_bookkeep_fix = fix; g.chtlu = lookup();
  g.soildimpar.coefshlwa = FA; g.soildimpar.coefshlwb = FB;
  g.soildimpar.coefdeepa = CA; g.soildimpar.coefdeepb = CB;
}
void finish(Ground& g) { g.frontsz.clear(); g.frontstype.clear(); g.resortGroundLayers(); g.updateSoilHorizons(); }
struct Sums { double c = 0, n = 0, a = 0; int nl = 0; bool nn = true; };
Sums sums(const Ground& g) {
  Sums s;
  for (Layer* l = g.fstsoill; l; l = l->nextl) { s.c += C(l); s.n += l->orgn; s.a += l->avln; ++s.nl; s.nn = s.nn && nonneg(l); }
  return s;
}
} // namespace

TEST_SUITE("layerbk_v5q") {

TEST_CASE("layer_bookkeep_fix default 0: ModelData(), ModelData(Json) without the key, Ground") {
  ModelData md; CHECK(md.layer_bookkeep_fix == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.layer_bookkeep_fix == 0);
  Ground g; CHECK(g.layer_bookkeep_fix == 0);
}

TEST_CASE("layer_bookkeep_fix validation: 0/1 only, anything else FATAL; 1 requires layer_n_remap_fix 1") {
  CHECK(ModelData(bmsg(0)).layer_bookkeep_fix == 0);
  CHECK(ModelData(bmsg(0, 0)).layer_bookkeep_fix == 0);
  CHECK(ModelData(bmsg(1)).layer_bookkeep_fix == 1);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg(1, 0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(bmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("shlwFromBelowFrac: c_want/c_have clamped to [0,1], 0 if nothing wanted or available") {
  CHECK(Ground::shlwFromBelowFrac(25., 100.) == doctest::Approx(0.25));
  CHECK(Ground::shlwFromBelowFrac(150., 100.) == 1.0);
  CHECK(Ground::shlwFromBelowFrac(0., 100.) == 0.0);
  CHECK(Ground::shlwFromBelowFrac(-5., 100.) == 0.0);
  CHECK(Ground::shlwFromBelowFrac(5., 0.) == 0.0);
  CHECK(Ground::shlwFromBelowFrac(5., -1.) == 0.0);
  CHECK(Ground::shlwFromBelowFrac(std::nan(""), 100.) == 0.0);
}

// (1) no fibric layer: humic 0.20 m over mineral; aboveground litterfall C of the month = lit
TEST_CASE("(1) redivideShlwLayers, new fibric layer: switch 1 C/N moved from the layer below (column conserved, C:N kept); switch 0 legacy creates C") {
  for (int fix : {0, 1}) for (double lit : {50., 1.0e5, -1.}) for (double chum : {3000., 0.}) {
    CAPTURE(fix); CAPTURE(lit); CAPTURE(chum);
    if (fix == 0 && lit > 1.0e4) continue;   // legacy would give the humic layer dz < 0 (frac > 1)
    Ground g; setup(g, fix);
    BgcData bd; bd.m_v2soi.ltrfalcall = (lit > 0 ? lit : 0.) + 7.; bd.m_v2soi.ltrfalc[I_root] = 7.;
    g.setBgcData(&bd);
    SoilLayer* h = org(2, 0.20, 0.5);
    SoilLayer* mn = ut::make_soil_layer(0.2, 0.5, 0.4, 50., 100., -1., 1, 1.);
    g.insertBack(h); g.insertBack(mn); finish(g);
    setc(h, chum, 25., 0.6); setc(mn, 40000., 20., 2.0);
    Sums s0 = sums(g);
    double hc = C(h), hn = h->orgn, ha = h->avln;
    g.redivideShlwLayers();
    Sums s1 = sums(g);
    const bool made = (g.fstshlwl != NULL);
    if (fix == 0) {
      CHECK(made == (lit > 0));
      CHECK(s1.c == doctest::Approx(s0.c + (lit > 0 ? lit : 0.)).epsilon(1e-13));   // legacy: litter C created
      CHECK(s1.n == doctest::Approx(s0.n).epsilon(1e-14));
      if (made) { CHECK(g.fstshlwl->orgn == 0.0); CHECK(C(g.fstshlwl) == doctest::Approx(lit)); }
    } else {
      CHECK(made == (lit > 0 && chum > 0));
      CHECK(s1.c == doctest::Approx(s0.c).epsilon(1e-14));
      CHECK(s1.n == doctest::Approx(s0.n).epsilon(1e-14));
      CHECK(s1.a == doctest::Approx(s0.a).epsilon(1e-14));
      CHECK(s1.nn);
      if (made) {
        Layer* f = g.fstshlwl;
        CHECK(C(f) == doctest::Approx(std::min(lit, hc)).epsilon(1e-13));
        CHECK(C(f) + C(h) == doctest::Approx(hc).epsilon(1e-14));
        CHECK(f->orgn + h->orgn == doctest::Approx(hn).epsilon(1e-14));
        CHECK(f->avln + h->avln == doctest::Approx(ha).epsilon(1e-14));
        CHECK(C(f) / f->orgn == doctest::Approx(25.).epsilon(1e-12));
        if (C(h) > 1e-9) CHECK(C(h) / h->orgn == doctest::Approx(25.).epsilon(1e-12));
        // thickness from the C moved (== legacy thickness when the layer below has enough C)
        CHECK(f->dz == doctest::Approx(g.thicknessFromCarbon(C(f), FA, FB)).epsilon(1e-12));
      }
    }
  }
}

// (2) humic layer of zero thickness (C left): deepnum == 0 removes it
TEST_CASE("(2) redivideDeepLayers, deepnum == 0: switch 1 removed humic C/N added to the mineral layer below (or the layer above); switch 0 legacy loses them") {
  for (int fix : {0, 1}) for (bool withmin : {true, false}) {
    CAPTURE(fix); CAPTURE(withmin);
    Ground g; setup(g, fix);
    SoilLayer* fib = org(1, 0.05, 1.0);
    SoilLayer* h = org(2, 0.05, 0.5);
    g.insertBack(fib); g.insertBack(h);
    SoilLayer* mn = NULL;
    if (withmin) { mn = ut::make_soil_layer(0.1, 0.5, 0.4, 50., 100., -1., 1, 1.); g.insertBack(mn); }
    finish(g);
    setc(fib, 500., 40., 0.2); setc(h, 120., 24., 0.3); if (mn) setc(mn, 40000., 20., 2.0);
    h->dz = 0.0;  // C left in a layer of no thickness
    Sums s0 = sums(g);
    double hc = C(h), hn = h->orgn, ha = h->avln;
    Layer* rcv = withmin ? (Layer*)mn : (Layer*)fib;
    double rc = C(rcv), rn = rcv->orgn, ra = rcv->avln;
    g.redivideDeepLayers();
    CHECK(g.fstdeepl == NULL);
    Sums s1 = sums(g);
    CHECK(s1.nl == s0.nl - 1);
    if (fix == 1) {
      CHECK(s1.c == doctest::Approx(s0.c).epsilon(1e-14));
      CHECK(s1.n == doctest::Approx(s0.n).epsilon(1e-14));
      CHECK(s1.a == doctest::Approx(s0.a).epsilon(1e-14));
      CHECK(s1.nn);
      CHECK(C(rcv) == doctest::Approx(rc + hc).epsilon(1e-14));
      CHECK(rcv->orgn == doctest::Approx(rn + hn).epsilon(1e-14));
      CHECK(rcv->avln == doctest::Approx(ra + ha).epsilon(1e-14));
    } else {
      CHECK(s1.c == doctest::Approx(s0.c - hc).epsilon(1e-14));    // legacy: lost
      CHECK(s1.n == doctest::Approx(s0.n - hn).epsilon(1e-14));
      CHECK(C(rcv) == rc);
    }
  }
}

TEST_CASE("(4) splitOneSoilLayer, organic layer below the C-thickness curve: switch 1 upper C >= 0 (lower gets at most all C), N follows; switch 0 legacy upper C negative") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    Ground g; setup(g, fix);
    SoilLayer* h = org(2, 0.40, 0.5);
    SoilLayer* mn = ut::make_soil_layer(0.4, 0.5, 0.4, 50., 100., -1., 1, 1.);
    g.insertBack(h); g.insertBack(mn); finish(g);
    const double cw = 0.5 * g.carbonFromThickness(0.40, CA, CB);
    setc(h, cw, 20., 0.8);
    REQUIRE(g.carbonFromThickness(0.40, CA, CB) - g.carbonFromThickness(0.10, CA, CB) > cw);
    SoilLayer* lo = org(2, 0.30, 0.5);
    g.splitOneSoilLayer(h, lo, 0., 0.30);
    CHECK(C(h) + C(lo) == doctest::Approx(cw).epsilon(1e-13));
    CHECK(h->orgn + lo->orgn == doctest::Approx(cw / 20.).epsilon(1e-13));
    CHECK(h->avln + lo->avln == doctest::Approx(0.8).epsilon(1e-13));
    if (fix == 1) {
      CHECK(nonneg(h)); CHECK(nonneg(lo));
      CHECK(C(h) == 0.0); CHECK(h->orgn == 0.0);
      CHECK(C(lo) / lo->orgn == doctest::Approx(20.).epsilon(1e-12));
    } else {
      CHECK(C(h) < 0.0);   // legacy bug
    }
    delete lo;
  }
}

TEST_CASE("(4) splitOneSoilLayer on the curve (normal case): identical for switch 0 and 1") {
  double r[2][4];
  for (int fix : {0, 1}) {
    Ground g; setup(g, fix);
    SoilLayer* h = org(2, 0.40, 0.5);
    g.insertBack(h); finish(g);
    setc(h, g.carbonFromThickness(0.40, CA, CB), 20., 0.8);
    SoilLayer* lo = org(2, 0.30, 0.5);
    g.splitOneSoilLayer(h, lo, 0., 0.30);
    r[fix][0] = C(h); r[fix][1] = C(lo); r[fix][2] = lo->orgn; r[fix][3] = h->avln;
    CHECK(C(h) > 0.);
    delete lo;
  }
  for (int i = 0; i < 4; ++i) CHECK(r[0][i] == r[1][i]);
}

TEST_CASE("(3) ltrcnSourceIndex: same horizon class and relative depth; nearest same class; else nearest absolute depth") {
  std::vector<Ground::LtrCNSnap> s(4);
  s[0].cls = 1; s[0].rtop = 0.;  s[0].rbot = 0.05; s[0].amid = 0.025;
  s[1].cls = 2; s[1].rtop = 0.;  s[1].rbot = 0.10; s[1].amid = 0.10;
  s[2].cls = 2; s[2].rtop = 0.10; s[2].rbot = 0.40; s[2].amid = 0.30;
  s[3].cls = 3; s[3].rtop = 0.;  s[3].rbot = 0.50; s[3].amid = 0.70;
  CHECK(Ground::ltrcnSourceIndex(s, 2, 0.03, 0.08) == 1);
  CHECK(Ground::ltrcnSourceIndex(s, 2, 0.13, 0.18) == 2);
  CHECK(Ground::ltrcnSourceIndex(s, 2, 0.50, 0.55) == 2);   // below all: nearest same class
  CHECK(Ground::ltrcnSourceIndex(s, 3, 0.25, 0.70) == 3);
  CHECK(Ground::ltrcnSourceIndex(s, 0, 0.005, 0.005) == 0); // class absent before: nearest abs depth
  CHECK(Ground::ltrcnSourceIndex(std::vector<Ground::LtrCNSnap>(), 2, 0.1, 0.1) == -1);
}

TEST_CASE("(3) prvltrfcnque around redivideDeepLayers (2 -> 3 humic): switch 1 history follows its layer; no-op when the structure is unchanged; switch 0 (no call) stays by index") {
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    Ground g; setup(g, fix);
    BgcData bd; g.setBgcData(&bd);
    SoilLayer* fib = org(1, 0.05, 1.0); SoilLayer* h1 = org(2, 0.10, 0.5); SoilLayer* h2 = org(2, 0.30, -0.5);
    SoilLayer* mn = ut::make_soil_layer(0.45, 0.5, 0.4, 50., 100., -1., 1, 1.);
    g.insertBack(fib); g.insertBack(h1); g.insertBack(h2); g.insertBack(mn); finish(g);
    setc(fib, g.carbonFromThickness(0.05, FA, FB), 40., 0.2);
    setc(h1, g.carbonFromThickness(0.10, CA, CB), 20., 0.5);
    setc(h2, g.carbonFromThickness(0.40, CA, CB) - g.carbonFromThickness(0.10, CA, CB), 20., 1.5);
    setc(mn, 40000., 20., 2.0);
    const double v[6] = {40., 21., 22., 30., 99., 77.};
    for (int i = 0; i < 6; ++i) { bd.prvltrfcnque[i].clear(); bd.prvltrfcnque[i].push_front(v[i]); bd.prvltrfcnque[i].push_front(v[i] + 0.5); }
    // unchanged structure: identity
    if (fix == 1) {
      g.snapshotLitterCNQueues(); g.remapLitterCNQueues();
      for (int i = 0; i < 6; ++i) { CHECK(bd.prvltrfcnque[i].size() == 2); CHECK(bd.prvltrfcnque[i].back() == v[i]); }
    }
    if (fix == 1) g.snapshotLitterCNQueues();
    g.redivideDeepLayers();
    if (fix == 1) g.remapLitterCNQueues();
    int n = 0; for (Layer* l = g.fstsoill; l; l = l->nextl) ++n;
    REQUIRE(n == 5);   // fib, h, h, h, mineral
    const double expect1[6] = {40., 21., 22., 22., 30., 77.};   // index 5 untouched
    for (int i = 0; i < 6; ++i) {
      CAPTURE(i);
      CHECK(bd.prvltrfcnque[i].size() == 2);
      CHECK(bd.prvltrfcnque[i].back() == (fix == 1 ? expect1[i] : v[i]));
      CHECK(bd.prvltrfcnque[i].front() == (fix == 1 ? expect1[i] : v[i]) + 0.5);
    }
  }
}

} // TEST_SUITE layerbk_v5q
