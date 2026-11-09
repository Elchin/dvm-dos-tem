// Phase 2b: energy-closure diagnostic (model_settings.energy_diag, default 0 = off; output only).
// Checks: switch parsing; layer enthalpy (relative to liquid water at 0 C, latent heat of fusion);
// the solver's top/bottom boundary heat flux accumulators close the solver's own discrete energy
// budget sum_i 0.5*(cap_i + cap_{i-1}) dT_i (Crank-Nicolson form of CrankNicholson::geBackward,
// virtual top node cap = 0); the diagnostic never changes temperatures.
// Refs: Bonan 2019 (Climate Change and Terrestrial Ecosystem Modeling, ch. 5, soil enthalpy);
//   Jafarov et al. 2012, doi:10.5194/tc-6-613-2012. dvmdostem-specific residual: no published source.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

using ut::make_soil_layer;

namespace {
Json::Value msg(Json::Value v) { Json::Value c; c["model_settings"]["energy_diag"] = v; return c; }

struct ConductCol {   // 4 unfrozen, wet 0.1 m layers, no fronts, no snow
  ut::Column col; TemperatureUpdator tu;
  ConductCol(double t0 = 2.0)
    : col(std::vector<SoilLayer*>{make_soil_layer(0.0, 0.1, 0.4, 0.5 * 0.4 * DENLIQ * 0.1, 0.0, t0, -1, 0.0),
                                  make_soil_layer(0.1, 0.1, 0.4, 0.5 * 0.4 * DENLIQ * 0.1, 0.0, t0 + 0.5, -1, 0.0),
                                  make_soil_layer(0.2, 0.2, 0.4, 0.6 * 0.4 * DENLIQ * 0.2, 0.0, t0 + 1.0, -1, 0.0),
                                  make_soil_layer(0.4, 0.3, 0.4, 0.7 * 0.4 * DENLIQ * 0.3, 0.0, t0 + 1.5, -1, 0.0)}) {
    for (Layer* l : col.L) { l->pce_f = 0.0; l->pce_t = 0.0; }
    tu.setGround(&col.g);
  }
  double cap(int i) const { Layer* l = col.L[i]; return l->getHeatCapacity() * l->dz; }
  double H() const { double h = 0; for (Layer* l : col.L) h += Soil_Env::layerEnthalpy(l); return h; }
  void solve(double tdrv, double dt, int on) {
    tu.edg_on = on; tu.edg_top = 0.0; tu.edg_bot = 0.0;
    tu.updateTemps(tdrv, col.g.toplayer, col.g.botlayer, col.g.fstsoill, NULL, NULL, dt, false);
    tu.edg_on = 0;
  }
};
}

TEST_SUITE("energy_diag") {

TEST_CASE("energy_diag: default 0 (ModelData(), without the key, Soil_Bgc/Soil_Env/TemperatureUpdator members)") {
  ModelData md; CHECK(md.energy_diag == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.energy_diag == 0);
  Soil_Env se; CHECK(se.energy_diag == 0); CHECK(se.energydiag.used == 0);
  TemperatureUpdator tu; CHECK(tu.edg_on == 0); CHECK(tu.edg_top == 0.0); CHECK(tu.edg_bot == 0.0);
}

TEST_CASE("energy_diag validation: integers 0/1 only, anything else FATAL") {
  CHECK(ModelData(msg(0)).energy_diag == 0);
  CHECK(ModelData(msg(1)).energy_diag == 1);
  CHECK(ut::exit_status_of([] { ModelData m(msg(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msg(true)); }) == EXIT_FAILURE);
}

TEST_CASE("layerEnthalpy: solid + water sensible heat, ice carries -L_f; 0 C liquid = 0; snow has no solid part") {
  SoilLayer* a = make_soil_layer(0.0, 0.1, 0.4, 20.0, 0.0, 0.0, -1, 0.0);
  CHECK(Soil_Env::layerEnthalpy(a) == doctest::Approx(0.0));
  SoilLayer* b = make_soil_layer(0.0, 0.1, 0.4, 10.0, 15.0, -3.0, 0, 0.5);
  const double solid = b->vhcsolid * (1 - b->poro) * b->dz;
  CHECK(Soil_Env::layerEnthalpy(b) == doctest::Approx(solid * -3.0 + 10.0 * SHCLIQ * -3.0 + 15.0 * (SHCICE * -3.0 - LHFUS)).epsilon(1e-10));
  SnowLayer sn; sn.dz = 0.3; sn.ice = 75.0; sn.liq = 0.0; sn.tem = -10.0; sn.rho = 250.0;
  CHECK(Soil_Env::layerEnthalpy(&sn) == doctest::Approx(75.0 * (SHCICE * -10.0 - LHFUS)).epsilon(1e-12));
  // melting 1 kg of ice at 0 C needs exactly L_f
  SoilLayer* c0 = make_soil_layer(0.0, 0.1, 0.4, 10.0, 5.0, 0.0, 0, 0.5);
  SoilLayer* c1 = make_soil_layer(0.0, 0.1, 0.4, 11.0, 4.0, 0.0, 0, 0.5);
  CHECK(Soil_Env::layerEnthalpy(c1) - Soil_Env::layerEnthalpy(c0) == doctest::Approx(LHFUS).epsilon(1e-10));
  delete a; delete b; delete c0; delete c1;
}

TEST_CASE("no-front conduction day: top + bottom boundary heat closes the solver's discrete energy budget") {
  for (double tdrv : {10.0, -0.5, 25.0}) {
    for (double dt : {86400.0, 3600.0}) {
      CAPTURE(tdrv); CAPTURE(dt);
      ConductCol c;
      const int n = (int)c.col.L.size();
      std::vector<double> T0(n), cp(n);
      for (int i = 0; i < n; ++i) { T0[i] = c.col.L[i]->tem; cp[i] = c.cap(i); }
      c.solve(tdrv, dt, 1);
      double es = 0.0;   // solver energy: effective capacity 0.5*(cap_i + cap_{i-1}), cap_0 (virtual) = 0
      for (int i = 0; i < n; ++i) es += 0.5 * (cp[i] + (i > 0 ? cp[i-1] : 0.0)) * (c.col.L[i]->tem - T0[i]);
      const double g = c.tu.edg_top + c.tu.edg_bot;
      CHECK(std::fabs(c.tu.edg_top) > 0.0);
      CHECK(es == doctest::Approx(g).epsilon(1e-6));
      CHECK((tdrv > T0[0]) == (c.tu.edg_top > 0.0));     // heat flows down a warmer surface
    }
  }
}

TEST_CASE("energy_diag on/off: identical temperatures (output-only diagnostic)") {
  ConductCol a, b;
  a.solve(12.0, 86400.0, 0); b.solve(12.0, 86400.0, 1);
  for (size_t i = 0; i < a.col.L.size(); ++i) { CAPTURE(i); CHECK(a.col.L[i]->tem == b.col.L[i]->tem); }
  CHECK(a.tu.edg_top == 0.0); CHECK(b.tu.edg_top != 0.0);
}


TEST_CASE("phase 2d: the single-layer block shortcut books its boundary fluxes; temperatures unchanged by edg_on") {
  // one soil layer between the virtual top and bottom nodes -> TemperatureUpdator::iterate shortcut (no solve)
  double tt[2], top[2], bot[2];
  for (int on : {0, 1}) {
    ut::Column c(std::vector<SoilLayer*>{make_soil_layer(0.0, 0.2, 0.5, 50.0, 0.0, 2.0, -1, 0.0)});
    TemperatureUpdator tu; tu.setGround(&c.g);
    tu.edg_on = on; tu.edg_top = 0.0; tu.edg_bot = 0.0;
    tu.updateTemps(10.0, c.g.toplayer, c.g.botlayer, c.g.fstsoill, NULL, NULL, 86400.0, false);
    tt[on] = c.L[0]->tem; top[on] = tu.edg_top; bot[on] = tu.edg_bot;
  }
  CHECK(tt[0] == tt[1]);
  CHECK(top[0] == 0.0); CHECK(bot[0] == 0.0);
  CHECK(top[1] > 0.0);       // warm surface: heat into the column (was 0 before the fix)
  CHECK(bot[1] <= 0.0);      // and out through the (colder) bottom node
}
}
