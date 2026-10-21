// Area 9: v5j SEB (model_settings.seb_scheme 2) -- Soil_Env::solveSEB, the daily skin
// energy balance with Newton + bisection safeguard; and the v5j Sturm 1997 snow
// conductivity switch (SnowLayer::snow_tc_form). Requires the v5j commits
// 9dae595b + d41c7f79, cherry-picked onto branch unit-tests.
#include "doctest.h"
#include "test_helpers.h"

using ut::make_soil_layer;

namespace {
struct StaticRestore {
  int vp, tc;
  StaticRestore() : vp(Climate::vp_hpa_fix), tc(SnowLayer::snow_tc_form) {}
  ~StaticRestore() { Climate::vp_hpa_fix = vp; SnowLayer::snow_tc_form = tc; }
};

struct SebFixture {
  ut::Column col; EnvData ed; CohortData cd; Snow_Env snowenv; Soil_Env se; SnowLayer* sn;
  explicit SebFixture(bool snow)
    : col(std::vector<SoilLayer*>{make_soil_layer(0.0, 0.05, 0.6, 15, 0, -2.0, 1, 1),
                                  make_soil_layer(0.05, 0.1, 0.5, 20, 0, -3.0, 1, 1),
                                  make_soil_layer(0.15, 0.2, 0.4, 30, 0, -4.0, 1, 1)}), sn(NULL) {
    for (size_t i = 0; i < col.L.size(); ++i) {
      SoilLayer* s = dynamic_cast<SoilLayer*>(col.L[i]);
      s->albsatvis = 0.09; s->albdryvis = 0.18; s->albsatnir = 0.18; s->albdrynir = 0.36;
    }
    if (snow) {
      sn = new SnowLayer(); sn->dz = 0.3; sn->rho = 250.0; sn->ice = 75.0; sn->liq = 0.0; sn->tem = -10.0;
      sn->z = 0.0; sn->frozen = 1;
      col.g.insertFront(sn);
    }
    se.setGround(&col.g); se.setEnvData(&ed); se.setCohortData(&cd);
    snowenv.snowenvpar.albmax = 0.8; snowenv.snowenvpar.albmin = 0.4;
    se.snowenvp = &snowenv;
    se.seb_scheme = 2; se.seb_z0_snow = 0.001; se.seb_z0_soil = 0.01; se.seb_zwind = 10.0; se.seb_ztemp = 2.0;
    se.seb_emis_snow = 0.99; se.seb_emis_soil = 0.97; se.seb_wind = -1.0; se.seb_minwind = 0.5;
    se.seb_stab_b = 10.0; se.seb_tran_frac = 1.0; se.seb_pres = 100000.0; se.seb_fc_last = 0.0;
    cd.m_vegd.fpcsum = 0.0;
    ed.d_v2g.swthfl = 0.0; ed.d_a2l.nirr = 0.0; ed.seb_lwdown = 250.0; ed.seb_wind = 3.0;
    ed.d_atmd.vp = 2.0;                 // hPa (vp_hpa_fix = 0 convention)
    ed.d_soi2a.evap = 0.5; ed.d_v2a.tran = 1.0;   // mm/day (snow-free LE)
    Climate::vp_hpa_fix = 0;
  }
  double residual() const {   // (1-alb) SW + LWin - LWout - H - LE - G - Qmelt
    const SebDiag& d = se.sebdiag;
    return (1.0 - d.alb) * d.swin + d.lwin - d.lwout - d.h - d.le - d.g - d.qmelt;
  }
};
}

TEST_SUITE("seb") {

TEST_CASE("snow, cold clear night: energy balance closes, Ts <= 0 and below Tair (radiative cooling)") {
  StaticRestore r; SebFixture f(true);
  f.ed.seb_lwdown = 180.0;
  double ts = f.se.solveSEB(-20.0);
  CHECK(f.se.sebdiag.snow == 1);
  CHECK(std::fabs(f.residual()) < 1e-4);
  CHECK(ts <= 0.0); CHECK(ts < -20.0); CHECK(ts > -80.0);
  CHECK(f.se.sebdiag.qmelt == 0.0);
  CHECK(f.col.g.seb_melt == 0.0);
  CHECK(f.se.sebdiag.h < 0.0);                     // sensible heat towards the colder surface
  CHECK(f.se.sebdiag.iters < 60);
  CHECK(f.se.sebdiag.alb == doctest::Approx(0.8)); // snow tem = -10 C -> albmax
}

TEST_CASE("snow, warm sunny day: Ts = 0, energy surplus melts snow (melt = Qmelt * 86400 / Lf)") {
  StaticRestore r; SebFixture f(true);
  f.sn->tem = -0.5; f.ed.d_a2l.nirr = 350.0; f.ed.seb_lwdown = 320.0; f.ed.d_atmd.vp = 8.0;
  double ts = f.se.solveSEB(6.0);
  CHECK(ts == 0.0);
  CHECK(f.se.sebdiag.qmelt > 0.0);
  CHECK(std::fabs(f.residual()) < 1e-6);
  CHECK(f.col.g.seb_melt == doctest::Approx(f.se.sebdiag.qmelt * 86400.0 / LHFUS));
  CHECK(f.col.g.seb_sublim == doctest::Approx(std::max(f.se.sebdiag.le, 0.0) * 86400.0 / LHSUB));
}

TEST_CASE("snow-free: LE is the model's own ET, balance closes, Ts may exceed 0") {
  StaticRestore r; SebFixture f(false);
  f.ed.d_a2l.nirr = 250.0; f.ed.seb_lwdown = 300.0;
  double ts = f.se.solveSEB(10.0);
  CHECK(f.se.sebdiag.snow == 0);
  CHECK(f.se.sebdiag.le == doctest::Approx(LHVAP * (0.5 + 1.0) / 86400.0));
  CHECK(std::fabs(f.residual()) < 1e-4);
  CHECK(ts > 0.0);
  CHECK(f.col.g.seb_melt == 0.0);
}

TEST_CASE("Ts increases monotonically with incoming shortwave") {
  StaticRestore r; double prev = -1e9;
  for (double sw = 0; sw <= 400; sw += 50) {
    SebFixture f(false); f.ed.d_a2l.nirr = sw;
    double ts = f.se.solveSEB(0.0);
    CHECK(ts > prev); prev = ts;
    CHECK(std::fabs(f.residual()) < 1e-4);
  }
}

TEST_CASE("bracket limits: Ts clamped to Tair + 60 and Tair - 60") {
  StaticRestore r;
  SebFixture hot(false); hot.ed.d_a2l.nirr = 1.0e5;
  CHECK(hot.se.solveSEB(5.0) == doctest::Approx(65.0));
  SebFixture cold(false); cold.col.L[0]->tem = -1000.0;     // enormous ground heat sink
  CHECK(cold.se.solveSEB(5.0) == doctest::Approx(-55.0));
  SebFixture snow(true); snow.ed.d_a2l.nirr = 1.0e5;       // over snow the upper bound is 0 C
  CHECK(snow.se.solveSEB(5.0) == 0.0);
  CHECK(snow.se.sebdiag.qmelt > 0.0);
}

TEST_CASE("vapour pressure units: vp_hpa_fix 0 (vp in hPa) and 1 (vp in Pa) give the same solution") {
  StaticRestore r;
  SebFixture a(true); a.ed.d_atmd.vp = 1.5; Climate::vp_hpa_fix = 0; double ta_ = a.se.solveSEB(-8.0);
  SebFixture b(true); b.ed.d_atmd.vp = 150.0; Climate::vp_hpa_fix = 1; double tb_ = b.se.solveSEB(-8.0);
  CHECK(ta_ == doctest::Approx(tb_).epsilon(1e-12));
  CHECK(a.se.sebdiag.le == doctest::Approx(b.se.sebdiag.le).epsilon(1e-12));
}

TEST_CASE("canopy LW fraction from SW transmission, kept on dark days") {
  StaticRestore r; SebFixture f(true);
  f.cd.m_vegd.fpcsum = 1.0; f.ed.d_a2l.nirr = 200.0; f.ed.d_v2g.swthfl = 120.0;
  f.se.solveSEB(-5.0);
  CHECK(f.se.seb_fc_last == doctest::Approx(0.4));
  f.ed.d_a2l.nirr = 5.0; f.ed.d_v2g.swthfl = 0.0;           // dark: fc not updated
  f.se.solveSEB(-5.0);
  CHECK(f.se.sebdiag.fc == doctest::Approx(0.4));
}

} // TEST_SUITE seb

TEST_SUITE("snow") {

TEST_CASE("v5j snow_tc_form: 0 = Goodrich (default), 1 = Sturm 1997 up to 600 kg/m3 with floor 0.023") {
  StaticRestore r;
  SnowLayer s; s.dz = 0.1;
  CHECK(SnowLayer::snow_tc_form == 0);
  s.rho = 300.0; SnowLayer::snow_tc_form = 0; double g300 = s.getThermalConductivity();
  CHECK(g300 == doctest::Approx(2.9e-6 * 300 * 300));
  SnowLayer::snow_tc_form = 1;
  CHECK(s.getThermalConductivity() == doctest::Approx(0.138 - 1.01 * 0.3 + 3.233 * 0.09));
  CHECK(s.getFrzThermCond() == s.getThermalConductivity());
  CHECK(s.getUnfThermCond() == s.getThermalConductivity());
  s.rho = 50.0;  CHECK(s.getThermalConductivity() == doctest::Approx(0.023 + 0.234 * 0.05));
  s.rho = 600.0; CHECK(s.getThermalConductivity() == doctest::Approx(0.138 - 1.01 * 0.6 + 3.233 * 0.36));
  s.rho = 650.0; CHECK(s.getThermalConductivity() == doctest::Approx(2.9e-6 * 650 * 650));   // Goodrich above 600
  // Sturm's two branches meet (to 0.0005 W/m/K) at 156 kg/m3
  s.rho = 155.999; double lo = s.getThermalConductivity(); s.rho = 156.0; double hi = s.getThermalConductivity();
  CHECK(std::fabs(hi - lo) < 5e-4);
}

TEST_CASE("v5j snow_tc_form = 1: documented jump at 600 kg/m3 (Sturm 0.696 -> Goodrich 1.044 W/m/K)") {
  StaticRestore r; SnowLayer s; s.dz = 0.1; SnowLayer::snow_tc_form = 1;
  s.rho = 600.0; double a = s.getThermalConductivity(); s.rho = 600.001; double b = s.getThermalConductivity();
  CHECK(a == doctest::Approx(0.6959).epsilon(1e-3)); CHECK(b == doctest::Approx(1.044).epsilon(1e-3));
}

} // TEST_SUITE snow

TEST_SUITE("snow") {

TEST_CASE("snow albedo: albmax at -10 C, linear to albmin at 0 C, NIR = VIS") {
  Snow_Env se; se.snowenvpar.albmax = 0.8; se.snowenvpar.albmin = 0.4;
  CHECK(se.getAlbedoVis(-10.0) == doctest::Approx(0.8));
  CHECK(se.getAlbedoVis(-5.0) == doctest::Approx(0.6));
  CHECK(se.getAlbedoVis(0.0) == doctest::Approx(0.4));
  CHECK(se.getAlbedoVis(3.0) == doctest::Approx(0.4));
  CHECK(se.getAlbedoNir(-5.0) == se.getAlbedoVis(-5.0));
}

// KNOWN BUG (reported, not fixed; legacy code): Snow_Env::getAlbedoVis returns a
// hard-coded 0.9 below -10 C instead of snowenvpar.albmax, so with the shipped
// parameters (snwalbmax 0.8) snow albedo jumps from 0.9 to 0.8 at -10 C and exceeds
// the parameterised maximum on cold days (affects snow SW reflection/sublimation and
// the v5j SEB).
TEST_CASE("snow albedo never exceeds albmax and is continuous at -10 C [KNOWN BUG: xfail]"
          * doctest::should_fail(true)) {
  Snow_Env se; se.snowenvpar.albmax = 0.8; se.snowenvpar.albmin = 0.4;
  CHECK(se.getAlbedoVis(-10.001) == doctest::Approx(se.getAlbedoVis(-10.0)).epsilon(1e-3));
  CHECK(se.getAlbedoVis(-30.0) <= 0.8);
}

} // TEST_SUITE snow
