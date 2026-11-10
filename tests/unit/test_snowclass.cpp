// model_settings.snow_scheme = 1: Sturm et al. (1995) snow classes as in SiBCASA
// (Schaefer et al. 2009, JGR 114 F02021, doi:10.1029/2008JF001125; eqs 1-7, Table 2), Anderson (1976)
// compaction in the CLM4.0 form (Oleson et al. 2010 eqs 7.42-7.46). Default 0 = legacy (bit-identical).
// Column test (synthetic cold windy tundra / calm taiga winter) writes CSV when SNOWCLASS_CSV_DIR is set.
#include "doctest.h"
#include "test_helpers.h"
#include <json/value.h>
#include <cmath>
#include <cstdlib>
#include <fstream>

namespace {
SnowClassPar cls(const char* n) { SnowClassPar p; REQUIRE(snowclass::defaultClass(n, p)); return p; }

// Ground with one soil layer below; snow added with the model's own construct/combine/divide.
struct SnowCol {
  ut::Column c;
  SnowCol() : c(std::vector<SoilLayer*>{ut::make_soil_layer(0.0, 0.1, 0.5, 10.0, 30.0, -2.0, 1, 1.0)}) {
    c.g.snowdimpar.newden = 100.; c.g.snowdimpar.denmax = 250.;
    c.g.resortGroundLayers();
  }
  ~SnowCol() {}
  std::vector<SnowLayer*> snow() {
    std::vector<SnowLayer*> v;
    for (Layer* l = c.g.toplayer; l != NULL && l->isSnow; l = l->nextl) v.push_back(dynamic_cast<SnowLayer*>(l));
    return v;
  }
  double swe() { double s = 0; for (auto l : snow()) s += l->ice; return s; }
  double depth() { double s = 0; for (auto l : snow()) s += l->dz; return s; }
  double heat() { double s = 0; for (auto l : snow()) s += l->ice * l->tem; return s; }
  void add(double mass, double t) {
    if (c.g.constructSnowLayers(mass, t)) c.g.resortGroundLayers();
    if (c.g.combineSnowLayers()) c.g.resortGroundLayers();
    if (c.g.divideSnowLayers()) c.g.resortGroundLayers();
  }
  // effective (series) bulk conductivity of the pack
  double bulkK() { double r = 0, d = 0; for (auto l : snow()) { r += l->dz / l->getFrzThermCond(); d += l->dz; } return d > 0 ? d / r : 0; }
};
} // namespace

TEST_SUITE("snow_classes") {

TEST_CASE("snow_scheme defaults to 0 (legacy) in ModelData, Ground and SnowLayer") {
  ModelData m; CHECK(m.snow_scheme == 0); CHECK(m.snow_class_wind == doctest::Approx(-1.0));
  Ground g; CHECK(g.snow_scheme == 0);
  g.snowdimpar.newden = 73.25;
  CHECK(g.newSnowDensity() == 73.25);   // scheme 0: exactly snwdennew
  SnowLayer s; CHECK(s.fdh == 0.0);
}

TEST_CASE("legacy k unchanged when fdh = 0; series k with depth hoar") {
  SnowLayer s; s.rho = 300; s.dz = 0.1; s.ice = 30;
  SnowLayer::snow_tc_form = 0;
  const double kg = 2.9e-6 * 300 * 300;
  CHECK(s.getFrzThermCond() == doctest::Approx(kg).epsilon(1e-14));
  s.fdh = 1.0; s.kdh = 0.18; CHECK(s.getFrzThermCond() == doctest::Approx(0.18));
  s.fdh = 0.5; CHECK(s.getUnfThermCond() == doctest::Approx(1.0 / (0.5 / 0.18 + 0.5 / kg)));
  s.fdh = 0.0; SnowLayer::snow_tc_form = 1; s.kdh = 0.18;
  const double rg = 0.3; const double ks = 0.138 - 1.01 * rg + 3.233 * rg * rg;
  CHECK(s.getFrzThermCond() == doctest::Approx(ks));   // respects seb_snow_tc (Sturm 1997)
  s.fdh = 1.0; s.kdh = 0.072; CHECK(s.getFrzThermCond() == doctest::Approx(0.072));
  SnowLayer::snow_tc_form = 0;
}

TEST_CASE("new-snow density eqs 4-7: hand values") {
  SnowClassPar tu = cls("tundra"), ta = cls("taiga"), ma = cls("maritime"), pr = cls("prairie");
  CHECK(snowclass::drhoT(-20) == 0.0);
  CHECK(snowclass::drhoT(-15) == doctest::Approx(0.0));
  CHECK(snowclass::drhoT(-5) == doctest::Approx(1.7 * pow(10.0, 1.5)));      // 53.76
  CHECK(snowclass::drhoT(0) == doctest::Approx(1.7 * pow(15.0, 1.5)));       // 98.76
  CHECK(snowclass::rhoNew(tu, -20, 0) == doctest::Approx(50));
  CHECK(snowclass::rhoNew(tu, -20, 2) == doctest::Approx(50));
  CHECK(snowclass::rhoNew(tu, -20, 3) == doctest::Approx(150));              // (3-2)*300/3
  CHECK(snowclass::rhoNew(tu, -20, 4.5) == doctest::Approx(300));
  CHECK(snowclass::rhoNew(tu, -20, 5) == doctest::Approx(350));              // = rho_topref
  CHECK(snowclass::rhoNew(tu, -20, 9) == doctest::Approx(350));              // eq 7 cap
  CHECK(snowclass::rhoNew(tu, -5, 3) == doctest::Approx(50 + 53.7587 + 100).epsilon(1e-5));
  CHECK(snowclass::rhoNew(pr, -20, 3) == doctest::Approx(150));              // prairie has the wind term
  CHECK(snowclass::rhoNew(ta, -20, 9) == doctest::Approx(50));               // taiga: no wind term
  CHECK(snowclass::rhoNew(ta, 0, 9) == doctest::Approx(50 + 1.7 * pow(15.0, 1.5)));
  CHECK(snowclass::rhoNew(ma, 5, 0) == doctest::Approx(200));                // 202.05 capped at 200
  for (double t = -40; t <= 5; t += 0.5) for (double u = 0; u <= 15; u += 0.5) {
    double r = snowclass::rhoNew(tu, t, u);
    CHECK(r >= 50.0); CHECK(r <= 350.0);
  }
}

TEST_CASE("fbot (eqs 1-2) and drho (eq 3) reproduce the tundra curves of Figs 4-5") {
  SnowClassPar tu = cls("tundra");
  CHECK(snowclass::fbot(tu, 0.05) == doctest::Approx(0.15));                 // Dhalf
  CHECK(snowclass::fbot(tu, 0.10) == doctest::Approx(0.3 / (1 + exp(-5.0))));
  CHECK(snowclass::fbot(tu, 0.0) == 0.0);
  CHECK(snowclass::fbot(tu, 1.0) == doctest::Approx(0.3).epsilon(1e-6));
  double prev = 0; for (double D = 0.001; D < 0.3; D += 0.001) { double f = snowclass::fbot(tu, D); CHECK(f >= prev); prev = f; }
  // Fig 5 (fbot = fbotmax): peak |rho_topref - rho_botref| = 250 at rho_bulkobs, 250/e at +-1 sigma
  CHECK(snowclass::drho(tu, 10.0, 261.108) == doctest::Approx(250).epsilon(1e-6));
  CHECK(snowclass::drho(tu, 10.0, 261.108 + 68.206) == doctest::Approx(250 * exp(-1.0)).epsilon(1e-6));
  CHECK(snowclass::drho(tu, 10.0, 261.108 - 2 * 68.206) == doctest::Approx(250 * exp(-4.0)).epsilon(1e-6));
  CHECK(snowclass::drho(tu, 10.0, 50.0) < 0.2);                                // fresh snow: ~0
  // thin snow below Dmin (maritime Dmin = 0.2 m): bottom layer essentially absent
  SnowClassPar ma = cls("maritime");
  CHECK(snowclass::fbot(ma, 0.1) < 0.007 * ma.fbotmax);
  CHECK(snowclass::drho(ma, 0.1, ma.rho_bulkobs) < 0.01 * 100);
}

TEST_CASE("eqs 1-3 with the values drawn in Figs 4-5 reproduce the figure curves") {
  // The paper's Figs 4-5 ("tundra") are drawn with fbotmax 0.8, Dmin 0.2 m, Dmax 0.4 m (Fig 4) and a
  // peak of ~300 kg/m3 at ~215 +- 50 kg/m3 (Fig 5), not with Table 2 tundra (0.3, 0, 0.1; 250 at 261 +- 68).
  // The model uses Table 2. Points read off the figures (+-0.02 / +-10 kg/m3):
  SnowClassPar f; f.fbotmax = 0.8; f.Dmin = 0.2; f.Dmax = 0.4;
  f.rho_topref = 300; f.rho_botref = 0; f.rho_bulkobs = 215; f.sigma_bulkobs = 50;
  CHECK(snowclass::fbot(f, 0.10) < 0.02);
  CHECK(snowclass::fbot(f, 0.20) < 0.02);
  CHECK(snowclass::fbot(f, 0.30) == doctest::Approx(0.40).epsilon(0.05));
  CHECK(snowclass::fbot(f, 0.40) == doctest::Approx(0.79).epsilon(0.03));
  CHECK(snowclass::fbot(f, 0.50) == doctest::Approx(0.80).epsilon(0.03));
  CHECK(snowclass::drho(f, 5.0, 215) == doctest::Approx(300).epsilon(0.04));
  CHECK(snowclass::drho(f, 5.0, 165) == doctest::Approx(110).epsilon(0.1));    // sigma marker (1/e)
  CHECK(snowclass::drho(f, 5.0, 265) == doctest::Approx(110).epsilon(0.1));
  CHECK(snowclass::drho(f, 5.0, 100) < 10);
  CHECK(snowclass::drho(f, 5.0, 350) < 10);
}

TEST_CASE("twoLayer conserves bulk density and stays within [50, 600]") {
  for (const char* n : {"tundra", "taiga", "maritime", "ephemeral", "prairie", "alpine"}) {
    SnowClassPar p = cls(n);
    for (double D = 0.01; D < 2.0; D += 0.07) for (double rb = 50; rb <= 600; rb += 25) {
      double fb, rt, rbt; snowclass::twoLayer(p, D, rb, fb, rt, rbt);
      CHECK(rt * (1 - fb) + rbt * fb == doctest::Approx(rb).epsilon(1e-12));
      CHECK(rt >= 50 - 1e-9); CHECK(rbt >= 50 - 1e-9); CHECK(rt <= 600 + 1e-9); CHECK(rbt <= 600 + 1e-9);
      if (p.rho_topref > p.rho_botref) CHECK(rt >= rbt); else CHECK(rt <= rbt);
    }
  }
}

TEST_CASE("Anderson/CLM compaction: rates, monotonicity, bounds, mass") {
  // CR1 at 0 C, rho <= 100, dry: 2.777e-6 1/s
  CHECK(snowclass::compactionRate(80, 0, 0, false) == doctest::Approx(-2.777e-6));
  CHECK(snowclass::compactionRate(80, 0, 0, true) == doctest::Approx(-2 * 2.777e-6));
  CHECK(snowclass::compactionRate(150, -10, 0, false) == doctest::Approx(-2.777e-6 * exp(-0.046 * 50) * exp(-0.4)));
  CHECK(snowclass::compactionRate(200, -5, 50, false) ==
        doctest::Approx(-2.777e-6 * exp(-4.6) * exp(-0.2) - 50 / (9e5 * exp(0.4 + 0.023 * 200))));
  double dz = 0.2, ice = 20;   // 100 kg/m3
  double d1 = snowclass::compactDz(dz, ice, -2, 0, false, 86400), d2 = snowclass::compactDz(dz, ice, -20, 0, false, 86400);
  CHECK(d1 < dz); CHECK(d2 < dz); CHECK(d1 < d2);                           // warm compacts faster
  CHECK(snowclass::compactDz(dz, ice, -10, 200, false, 86400) < snowclass::compactDz(dz, ice, -10, 0, false, 86400));
  double d = 1.0; for (int k = 0; k < 20000; ++k) d = snowclass::compactDz(d, 50, 0, 1000, true, 86400);
  CHECK(50 / d <= 600 + 1e-9);                                               // capped
  CHECK(snowclass::compactDz(0.0, 0.0, -5, 0, false, 86400) == 0.0);
}

TEST_CASE("targetMasses + remapT: SWE and sum(ice*T) conserved, depth-hoar fractions") {
  SnowClassPar tu = cls("tundra");
  std::vector<double> dz = {0.02, 0.05, 0.1, 0.15, 0.08};
  double M = 70.0; std::vector<double> m, f;
  snowclass::targetMasses(tu, dz, M, m, f);
  double s = 0; for (double x : m) s += x;
  CHECK(s == doctest::Approx(M).epsilon(1e-15));
  CHECK(f[0] == 0.0); CHECK(f.back() == 1.0);                                // D = 0.4: fbot ~ 0.3 -> 0.12 m bottom
  double hb = 0; for (size_t i = 0; i < dz.size(); ++i) hb += f[i] * dz[i];
  CHECK(hb == doctest::Approx(snowclass::fbot(tu, 0.4) * 0.4).epsilon(1e-12));
  std::vector<double> mo = {5, 10, 20, 20, 15}, to = {-25, -20, -12, -6, -3};
  std::vector<double> tn = snowclass::remapT(mo, to, m);
  double h0 = 0, h1 = 0; for (size_t i = 0; i < mo.size(); ++i) { h0 += mo[i] * to[i]; h1 += m[i] * tn[i]; }
  CHECK(h1 == doctest::Approx(h0).epsilon(1e-13));
  for (double t : tn) { CHECK(t >= -25 - 1e-12); CHECK(t <= -3 + 1e-12); }
  // no snow
  snowclass::targetMasses(tu, {}, 0.0, m, f); CHECK(m.empty());
}

TEST_CASE("Ground: stratify conserves SWE, depth and heat; no-snow and one-layer cases") {
  SnowCol sc; sc.c.g.snow_scheme = 1; sc.c.g.snowcls = cls("tundra");
  sc.c.g.stratifySnowLayers(true);                                            // no snow: nothing happens
  CHECK(sc.snow().empty());
  sc.c.g.snowcls_tair = -20; sc.c.g.snowcls_wind = 0;
  CHECK(sc.c.g.newSnowDensity() == doctest::Approx(50));
  sc.add(4.0, -20);                                                           // one thin layer (0.08 m)
  REQUIRE(sc.snow().size() == 1);
  double w0 = sc.swe(), d0 = sc.depth();
  sc.c.g.stratifySnowLayers(true);
  CHECK(sc.swe() == doctest::Approx(w0).epsilon(1e-15)); CHECK(sc.depth() == d0);
  CHECK(sc.snow()[0]->fdh > 0.0); CHECK(sc.snow()[0]->fdh < 0.3);            // part of it is depth hoar
  for (int k = 0; k < 30; ++k) { sc.c.g.snowcls_wind = (k % 3) * 3.0; sc.add(5.0, -15.0 - k % 7); }
  auto L = sc.snow(); REQUIRE(L.size() >= 3);
  for (size_t i = 0; i < L.size(); ++i) L[i]->tem = -25.0 + 20.0 * i / (L.size() - 1);
  w0 = sc.swe(); d0 = sc.depth(); double h0 = sc.heat();
  sc.c.g.stratifySnowLayers(true);
  CHECK(sc.swe() == doctest::Approx(w0).epsilon(1e-14));
  CHECK(sc.depth() == doctest::Approx(d0).epsilon(1e-15));
  CHECK(sc.heat() == doctest::Approx(h0).epsilon(1e-13));
  L = sc.snow();
  CHECK(L.front()->fdh == 0.0);
  const double hb = snowclass::fbot(sc.c.g.snowcls, d0) * d0;                   // depth-hoar thickness
  const double fl = fmin(1.0, hb / L.back()->dz);
  CHECK(L.back()->fdh == doctest::Approx(fl));
  { SnowLayer t; t.rho = L.back()->rho; const double kr = t.getFrzThermCond();
    CHECK(L.back()->getFrzThermCond() == doctest::Approx(fl >= 1.0 ? 0.18 : 1.0 / (fl / 0.18 + (1 - fl) / kr))); }
  double hsum = 0; for (auto l : L) hsum += l->fdh * l->dz;
  CHECK(hsum == doctest::Approx(hb).epsilon(1e-12));
  for (auto l : L) { CHECK(l->rho >= 50 - 1e-9); CHECK(l->rho <= 600 + 1e-9); CHECK(l->rho == doctest::Approx(l->ice / l->dz));
                    CHECK(l->ice <= l->maxice); }   // no clipping by snwdenmax*dz (would delete SWE)
  // maritime: no fixed depth-hoar k
  sc.c.g.snowcls = cls("maritime"); sc.c.g.stratifySnowLayers(false);
  for (auto l : sc.snow()) CHECK(l->fdh == 0.0);
}

TEST_CASE("Ground: compaction keeps SWE, melt at constant density then wet compaction raises bulk density") {
  SnowCol sc; sc.c.g.snow_scheme = 1; sc.c.g.snowcls = cls("taiga");
  sc.c.g.snowcls_tair = -10; sc.c.g.snowcls_wind = 0;
  for (int k = 0; k < 20; ++k) sc.add(6.0, -10);
  for (auto l : sc.snow()) l->tem = -8;
  double w0 = sc.swe(), d0 = sc.depth();
  sc.c.g.updateSnowLayerPropertiesDaily();                                     // -> compactSnowLayersClass
  CHECK(sc.swe() == doctest::Approx(w0).epsilon(1e-15));
  CHECK(sc.depth() < d0);
  // melt 10 mm from the top (the model's ablation path): SWE falls by exactly 10, density unchanged
  double rb0 = sc.swe() / sc.depth(), w1 = sc.swe();
  sc.c.g.constructSnowLayers(-10.0, 0.0); sc.c.g.resortGroundLayers();
  CHECK(sc.swe() == doctest::Approx(w1 - 10.0).epsilon(1e-14));
  for (auto l : sc.snow()) { l->tem = 0.0; l->frozen = -1; }                    // wet, isothermal pack
  double rb1 = sc.swe() / sc.depth();
  sc.c.g.updateSnowLayerPropertiesDaily(); sc.c.g.stratifySnowLayers(true);
  CHECK(sc.swe() / sc.depth() > rb1);
  CHECK(sc.swe() / sc.depth() > rb0 - 1e-9);
}

TEST_CASE("snow class file: CMT map, override, errors") {
  const char* f = "/tmp/snowclass_ut.json";
  { std::ofstream o(f); o << "{\"default_class\":\"tundra\",\"classes\":{\"taiga\":{\"k_dh\":0.08}},"
                           "\"cmt_class\":{\"CMT01\":\"taiga\",\"CMT05\":\"tundra\"}}"; }
  SnowClassPar p; std::string e;
  CHECK(snowclass::load(f, "CMT01", "", p, e)); CHECK(p.name == "taiga"); CHECK(p.k_dh == doctest::Approx(0.08));
  CHECK(p.fbotmax == doctest::Approx(0.7));
  CHECK(snowclass::load(f, "CMT05", "", p, e)); CHECK(p.name == "tundra"); CHECK(p.wind);
  CHECK(snowclass::load(f, "CMT99", "", p, e)); CHECK(p.name == "tundra");      // default
  CHECK(snowclass::load(f, "CMT05", "taiga", p, e)); CHECK(p.name == "taiga");   // override
  CHECK_FALSE(snowclass::load(f, "CMT05", "glacier", p, e));
  CHECK_FALSE(snowclass::load("/nonexistent.json", "CMT05", "", p, e));
  std::remove(f);
}

// (c) column test: synthetic Oct-Apr winter, 212 days. No soil thermal coupling: snow layer T is a
// prescribed linear profile between 0.7*Tair at the top and -3 C at the base (a rough insulated pack).
// Tundra: cold (-12...-30 C), windy (u = 1-9 m/s), 110 mm SWE. Taiga: -8...-28 C, calm (u <= 2), 160 mm.
// Expected April bulk density: Table 2 obs, tundra 261 +- 68, taiga 213 +- 49 (March).
static void column(const char* name, bool windy, double swe_tot, double tmin, double tmax,
                   double& rho_ref, double& k_ref, double& D_ref, int ref_day) {
  SnowCol sc; sc.c.g.snow_scheme = 1; sc.c.g.snowcls = cls(name);
  unsigned s = 12345u; auto rnd = [&s]() { s = s * 1103515245u + 12345u; return ((s >> 8) & 0xFFFF) / 65535.0; };
  std::ofstream csv; const char* dir = std::getenv("SNOWCLASS_CSV_DIR");
  if (dir) { csv.open(std::string(dir) + "/column_" + name + ".csv"); csv << "day,tair,wind,snowfall,swe,depth,rho_bulk,k_bulk,nlay,rho_top,rho_bot,fbot\n"; }
  const int N = 212; double total = 0;
  for (int d = 0; d < N; ++d) {
    const double ph = sin(M_PI * d / N);
    const double ta = tmax + (tmin - tmax) * ph + 4.0 * (rnd() - 0.5);
    const double u = windy ? 1.0 + 8.0 * rnd() : 2.0 * rnd();
    const double fall = (rnd() < 0.33) ? swe_tot / (0.33 * N) * 2.0 * rnd() : 0.0;
    total += fall;
    sc.c.g.snowcls_tair = ta; sc.c.g.snowcls_wind = u;
    if (fall > 0) sc.add(fall, fmin(ta, 0.0)); else sc.add(0.0, ta);
    auto L = sc.snow(); const double D = sc.depth(); double z = 0;
    for (auto l : L) { const double zm = z + 0.5 * l->dz; l->tem = 0.7 * ta + (-3.0 - 0.7 * ta) * (D > 0 ? zm / D : 0); z += l->dz; }
    double w0 = sc.swe();
    sc.c.g.updateSnowLayerPropertiesDaily();
    sc.c.g.stratifySnowLayers(true);
    CHECK(sc.swe() == doctest::Approx(w0).epsilon(1e-13));
    for (auto l : sc.snow()) CHECK(l->ice <= l->maxice + 1e-12);
    if (sc.c.g.combineSnowLayers()) sc.c.g.resortGroundLayers();
    if (sc.c.g.divideSnowLayers()) sc.c.g.resortGroundLayers();
    L = sc.snow();
    for (auto l : L) { CHECK(l->rho >= 50 - 1e-9); CHECK(l->rho <= 600 + 1e-9); }
    const double Dn = sc.depth(), rb = Dn > 0 ? sc.swe() / Dn : 0;
    double fb = 0, rt = 0, rbt = 0; if (Dn > 0) snowclass::twoLayer(sc.c.g.snowcls, Dn, rb, fb, rt, rbt);
    if (csv) csv << d << "," << ta << "," << u << "," << fall << "," << sc.swe() << "," << Dn << "," << rb << ","
                 << sc.bulkK() << "," << L.size() << "," << rt << "," << rbt << "," << fb << "\n";
    if (d == ref_day) { rho_ref = rb; k_ref = sc.bulkK(); D_ref = Dn; }
  }
  CHECK(sc.swe() + sc.c.g.snow.extramass == doctest::Approx(total).epsilon(1e-12));   // no loss: SWE conserved
}

TEST_CASE("column: tundra (windy) and taiga (calm) bulk density reach the Table 2 obs within 1 sigma") {
  double rt, kt, Dt, ra, ka, Da;
  column("tundra", true, 110.0, -30.0, -12.0, rt, kt, Dt, 197);   // mid-April
  column("taiga", false, 160.0, -28.0, -8.0, ra, ka, Da, 166);    // mid-March
  MESSAGE("tundra Apr: rho_bulk=" << rt << " k_bulk=" << kt << " D=" << Dt << " | taiga Mar: rho_bulk=" << ra << " k_bulk=" << ka << " D=" << Da);
  CHECK(fabs(rt - 261.108) <= 68.206);
  CHECK(fabs(ra - 213.363) <= 49.462);
  CHECK(rt > ra);
  CHECK(kt > ka);                                        // wind slab + indurated hoar vs taiga hoar
  CHECK(ka < 0.2);
}

} // TEST_SUITE
