// Area 4: snow conductivity (Goodrich 1982 vs Sturm et al. 1997) and densification;
// Area 5: organic layer thickness <-> carbon (Yi et al. 2009).
#include "doctest.h"
#include "test_helpers.h"

TEST_SUITE("snow") {

static SnowLayer snow_at(double rho) { SnowLayer s; s.rho = rho; s.dz = 0.1; s.ice = rho * 0.1; return s; }

TEST_CASE("Goodrich (1982) k = 2.9e-6 rho^2 with a 0.04 W/m/K floor (legacy default)") {
  CHECK(snow_at(100).getThermCond() == doctest::Approx(0.04));          // 0.029 -> floor
  CHECK(snow_at(150).getThermCond() == doctest::Approx(0.06525));
  CHECK(snow_at(200).getThermCond() == doctest::Approx(0.116));
  CHECK(snow_at(300).getThermCond() == doctest::Approx(0.261));
  CHECK(snow_at(400).getThermCond() == doctest::Approx(0.464));
  // Frz/Unf both use it
  CHECK(snow_at(250).getFrzThermCond() == doctest::Approx(2.9e-6 * 250 * 250));
  CHECK(snow_at(250).getUnfThermCond() == doctest::Approx(2.9e-6 * 250 * 250));
}

// Sturm et al. (1997, J. Glaciol. 43(143)), rho in g/cm3:
//   k = 0.138 - 1.01 rho + 3.233 rho^2   (0.156 <= rho <= 0.6)
//   k = 0.023 + 0.234 rho                 (rho < 0.156)
TEST_CASE("Sturm et al. (1997) conductivity: hand-calculated values") {
  CHECK(snow_at(100).getThermCond5Sturm() == doctest::Approx(0.0464));
  CHECK(snow_at(200).getThermCond5Sturm() == doctest::Approx(0.06532));
  CHECK(snow_at(300).getThermCond5Sturm() == doctest::Approx(0.12597));
  CHECK(snow_at(400).getThermCond5Sturm() == doctest::Approx(0.25128));
  CHECK(snow_at(500).getThermCond5Sturm() == doctest::Approx(0.44125));
  // the two branches nearly meet at 0.156 g/cm3 (0.0595 vs 0.0591)
  double lo = snow_at(155.9).getThermCond5Sturm(), hi = snow_at(156.0).getThermCond5Sturm();
  CHECK(std::fabs(lo - hi) < 0.001);
}

TEST_CASE("Sturm vs Goodrich: Sturm is lower for typical tundra snow (200-400 kg/m3)") {
  for (double rho : {200., 250., 300., 350., 400.}) {
    SnowLayer s = snow_at(rho);
    CHECK(s.getThermCond5Sturm() < s.getThermCond());
  }
  // ratio at 300 kg/m3 ~ 0.48
  CHECK(snow_at(300).getThermCond5Sturm() / snow_at(300).getThermCond() == doctest::Approx(0.4826).epsilon(1e-3));
}

TEST_CASE("updateDensity: exponential approach to snwdenmax, e-folding time 1/0.24 = 4.2 days") {
  snwpar_dim p; p.denmax = 200.0;
  SnowLayer s; s.rho = 100.0; s.dz = 0.2; s.ice = 20.0;
  s.updateDensity(&p);
  CHECK(s.rho == doctest::Approx(200.0 - 100.0 * std::exp(-0.24)));    // 121.34
  CHECK(s.maxice == doctest::Approx(200.0 * 0.2));
  double prev = s.rho;
  for (int d = 0; d < 60; ++d) { s.updateDensity(&p); CHECK(s.rho >= prev); CHECK(s.rho <= 200.0); prev = s.rho; }
  CHECK(s.rho == doctest::Approx(200.0).epsilon(1e-4));
  // gap after n days = gap0 * exp(-0.24 n)
  SnowLayer t; t.rho = 50.0; t.dz = 0.1;
  for (int d = 0; d < 10; ++d) t.updateDensity(&p);
  CHECK(200.0 - t.rho == doctest::Approx(150.0 * std::exp(-2.4)).epsilon(1e-9));
  // denser than denmax (e.g. after a parameter change) -> clamped
  SnowLayer u; u.rho = 260.0; u.dz = 0.1; u.updateDensity(&p);
  CHECK(u.rho == doctest::Approx(200.0));
  // thickness follows the density (mass conserved)
  SnowLayer v; v.rho = 100.0; v.ice = 30.0; v.dz = 0.3; v.updateDensity(&p); v.updateThick();
  CHECK(v.dz * v.rho == doctest::Approx(30.0));
}

} // TEST_SUITE snow

TEST_SUITE("organic_thickness") {

// Yi et al. (2009): C [gC/cm2] = a * T[cm]^b  ->  C [gC/m2] = a * 1e4 * T_cm^b
struct Coef { const char* name; double a, b; double t_m, c_hand; };
static const Coef COEFS[] = {
  {"CMT01 fibric", 0.01,       1.3574,    0.10, 2277.2},
  {"CMT01 humic",  0.03,       1.2389371, 0.20, 12274.8},
  {"CMT51 fibric", 0.0130,     1.3750,    0.10, 3082.8},
  {"CMT51 humic",  0.0230,     1.2450,    0.20, 9583.2},
  {"HPC k5006 fibric (coefshlwa)", 0.01675803, 1.3574, 0.10, 3816.1},
  {"HPC k5006 humic (coefdeepa)",  0.01237385, 1.2389371, 0.20, 5062.9},
};

TEST_CASE("carbonFromThickness matches hand-calculated Yi et al. values") {
  Ground g;
  for (const Coef& c : COEFS) {
    CAPTURE(c.name);
    CHECK(g.carbonFromThickness(c.t_m, c.a, c.b) == doctest::Approx(c.c_hand).epsilon(2e-4));
    CHECK(g.carbonFromThickness(c.t_m, c.a, c.b) == doctest::Approx(c.a * 1e4 * std::pow(c.t_m * 100, c.b)));
  }
  CHECK(g.carbonFromThickness(0.0, 0.01, 1.3574) == 0.0);
}

TEST_CASE("thicknessFromCarbon(carbonFromThickness(T)) == T for CMT01/CMT51 coefficients") {
  Ground g;
  for (const Coef& c : COEFS) {
    CAPTURE(c.name);
    for (double t : {0.005, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0}) {
      double C = g.carbonFromThickness(t, c.a, c.b);
      CHECK(g.thicknessFromCarbon(C, c.a, c.b) == doctest::Approx(t).epsilon(1e-12));
    }
    for (double C : {100., 1000., 5000., 20000.}) {
      double t = g.thicknessFromCarbon(C, c.a, c.b);
      CHECK(g.carbonFromThickness(t, c.a, c.b) == doctest::Approx(C).epsilon(1e-12));
    }
  }
}

TEST_CASE("carbon is monotonic in thickness; negative carbon is rejected") {
  Ground g;
  double prev = -1;
  for (int i = 0; i <= 100; ++i) { double C = g.carbonFromThickness(i / 100.0, 0.013, 1.375); CHECK(C > prev); prev = C; }
  CHECK_THROWS_AS(g.carbonFromThickness(0.1, -0.01, 1.3), std::runtime_error);
  // negative carbon -> NaN thickness -> forced to 0.0001 m (logged as fatal)
  CHECK(g.thicknessFromCarbon(-50.0, 0.01, 1.3574) == doctest::Approx(0.0001));
}

TEST_CASE("getOslThickness5Carbon: layer dz from cumulative C; shrink water booked as expelled") {
  Ground g;
  g.soildimpar.coefshlwa = 0.0130; g.soildimpar.coefshlwb = 1.3750;
  g.soildimpar.coefdeepa = 0.0230; g.soildimpar.coefdeepb = 1.2450;
  SoilLayer* sl = new SoilLayer();
  sl->isOrganic = true; sl->isFibric = true;
  sl->z = 0.04; sl->dz = 0.06; sl->poro = 0.9; sl->tcsolid = 0.25; sl->bulkden = 60.;
  sl->liq = 0.06 * 0.5 * 1000.; sl->ice = 0.06 * 0.2 * 917.;   // 30 + 11.0 mm
  g.insertBack(sl); g.fstsoill = sl; g.lstsoill = sl;
  // layer spans cumulative depths 0.04..0.07 m -> new dz 0.03 m (shrinks by half)
  double ctop = g.carbonFromThickness(0.04, 0.0130, 1.3750), cbot = g.carbonFromThickness(0.07, 0.0130, 1.3750);
  double w0 = sl->liq + sl->ice;
  g.expelled_water_mm = 0.0;
  g.getOslThickness5Carbon(sl, ctop, cbot);
  CHECK(sl->dz == doctest::Approx(0.03).epsilon(1e-9));
  CHECK(sl->liq + sl->ice == doctest::Approx(w0 * 0.5).epsilon(1e-9));
  // mass conservation: water left + water booked as expelled == water before
  CHECK(sl->liq + sl->ice + g.expelled_water_mm == doctest::Approx(w0).epsilon(1e-12));
  // growing layer: no water change, nothing expelled
  double e1 = g.expelled_water_mm, w1 = sl->liq + sl->ice;
  g.getOslThickness5Carbon(sl, ctop, g.carbonFromThickness(0.12, 0.0130, 1.3750));
  CHECK(sl->dz == doctest::Approx(0.08).epsilon(1e-9));
  CHECK(sl->liq + sl->ice == doctest::Approx(w1));
  CHECK(g.expelled_water_mm == doctest::Approx(e1));
  // humic layer uses the deep coefficients
  sl->isFibric = false; sl->isHumic = true;
  g.getOslThickness5Carbon(sl, 0.0, g.carbonFromThickness(0.1, 0.0230, 1.2450));
  CHECK(sl->dz == doctest::Approx(0.1).epsilon(1e-9));
}

} // TEST_SUITE organic_thickness
