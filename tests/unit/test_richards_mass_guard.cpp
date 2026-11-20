// model_settings.richards_mass_guard (bug R2, TEM-dev 2026-11-20; CMT06 S0_base sc 2041: two Richards days
// reached dtmin with the trial error above toler_upper, ran 8640 substeps of 10 s and closed the column
// liquid by +1124 / +2462 mm; checkSoilLiquidWaterValidity then dumped the surplus to qover).
// Default 0 = legacy (bit-identical). 1 = a day whose |final - target column liquid| > richards_mass_tol_mm
// is rejected: every layer's liquid is restored and the call returns as an early return (Soil_Env re-books
// infiltration to puddle/qover and evaporation as unapplied).
// A non-conserving day is produced here with the legacy bug R1 path (dt_fix 0 creates ~3 mm/d at 20 mm/d
// infiltration, see test_richards_dtfix.cpp); the guard must catch any non-conserving day, whatever the cause.
// Refs: water mass conservation (Richards 1931, doi:10.1063/1.1745010); CLM5 tech note sect. 7.3.3.
#include "doctest.h"
#include "test_helpers.h"
#include <json/value.h>
#include <cmath>

using ut::make_soil_layer;

namespace {
std::vector<SoilLayer*> gcol(int n, double dz, double poro, double s0) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) {
    SoilLayer* l = make_soil_layer(i * dz, dz, poro, s0 * poro * DENLIQ * dz, 0.0, 5.0, -1, 0.0);
    l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0;
    v.push_back(l);
  }
  return v;
}
double gsum(const ut::Column& c) { double w = 0; for (auto l : c.L) w += l->liq + l->ice; return w; }
void gday(Richards& rr, Layer* top, Layer* drainl, double infil_mmd) {
  double tr[MAX_SOI_LAY] = {0};
  const double dd = drainl->z + drainl->dz;
  rr.update(top, drainl, dd, 0.0, dd, tr, 0.0, infil_mmd / SEC_IN_DAY, 0.0);
}
Json::Value gms(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
} // namespace

TEST_SUITE("richards_mass_guard") {

TEST_CASE("richards_mass_guard: defaults 0 / 1.0 mm in Richards and ModelData; values read; invalid values FATAL") {
  Richards rr; CHECK(rr.mass_guard == 0); CHECK(rr.mass_tol_mm == 1.0); CHECK(rr.n_guard == 0);
  ModelData m0; CHECK(m0.richards_mass_guard == 0); CHECK(m0.richards_mass_tol_mm == 1.0);
  Json::Value c0; c0["model_settings"]["svp_ice_fix"] = 0; ModelData mj0(c0);
  CHECK(mj0.richards_mass_guard == 0); CHECK(mj0.richards_mass_tol_mm == 1.0);
  ModelData mj1(gms("richards_mass_guard", 1)); CHECK(mj1.richards_mass_guard == 1);
  ModelData mt(gms("richards_mass_tol_mm", 0.25)); CHECK(mt.richards_mass_tol_mm == 0.25);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_guard", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_guard", -1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_guard", "1")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_tol_mm", 0.0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_tol_mm", -1.0)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_tol_mm", true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gms("richards_mass_guard", 1)); }) == 0);
}

TEST_CASE("bug R2 guard: a non-conserving day (legacy R1 path) is rejected and the column liquid restored exactly") {
  ut::Column a(gcol(20, 0.05, 0.4, 0.3)), b(gcol(20, 0.05, 0.4, 0.3));
  Richards ra; ra.dt_fix = 0;                       // legacy: creates water on days with rejections
  Richards rb; rb.dt_fix = 0; rb.mass_guard = 1;    // same, guarded
  std::vector<double> l0; for (auto l : b.L) l0.push_back(l->liq);
  const double wa0 = gsum(a);
  gday(ra, a.L[0], a.L.back(), 20.0);
  const double ea = gsum(a) - wa0 - 20.0 + ra.qdrain;
  MESSAGE("legacy (no guard) daily mass error " << ea << " mm; Richards mass_err_mm " << ra.mass_err_mm);
  REQUIRE(fabs(ea) > 1.0);                          // the test day really does not conserve
  CHECK(ra.early_return == false); CHECK(ra.n_guard == 0);
  gday(rb, b.L[0], b.L.back(), 20.0);
  CHECK(rb.guard_tripped == true); CHECK(rb.early_return == true); CHECK(rb.n_guard == 1);
  CHECK(rb.qdrain == 0.0); CHECK(rb.mass_err_mm == 0.0);
  for (size_t i = 0; i < b.L.size(); ++i) CHECK(b.L[i]->liq == l0[i]);   // bit-exact restore
}

TEST_CASE("bug R2 guard: on a conserving day the guard does nothing (bit-identical to guard off)") {
  ut::Column a(gcol(20, 0.05, 0.4, 0.3)), b(gcol(20, 0.05, 0.4, 0.3));
  Richards ra; ra.dt_fix = 1; Richards rb; rb.dt_fix = 1; rb.mass_guard = 1;
  for (int d = 0; d < 5; ++d) {
    gday(ra, a.L[0], a.L.back(), 20.0); gday(rb, b.L[0], b.L.back(), 20.0);
    CHECK(rb.guard_tripped == false); CHECK(rb.early_return == false);
  }
  CHECK(ra.n_rejections > 0);                       // rejections happen, but the days close (dt_fix 1)
  CHECK(rb.n_guard == 0);
  for (size_t i = 0; i < a.L.size(); ++i) CHECK(a.L[i]->liq == b.L[i]->liq);
  CHECK(ra.qdrain == rb.qdrain);
}

TEST_CASE("bug R2 guard: tolerance is respected (tol above the day's error -> no trip)") {
  ut::Column a(gcol(20, 0.05, 0.4, 0.3)), b(gcol(20, 0.05, 0.4, 0.3));
  Richards ra; ra.dt_fix = 0; gday(ra, a.L[0], a.L.back(), 20.0);
  const double e = fabs(ra.mass_err_mm); REQUIRE(e > 1.0);
  Richards rb; rb.dt_fix = 0; rb.mass_guard = 1; rb.mass_tol_mm = 2.0 * e; gday(rb, b.L[0], b.L.back(), 20.0);
  CHECK(rb.guard_tripped == false); CHECK(rb.n_guard == 0);
  for (size_t i = 0; i < a.L.size(); ++i) CHECK(a.L[i]->liq == b.L[i]->liq);
}

TEST_CASE("bug R2 guard through Soil_Env: rejected day + re-booking + validity check conserve water (thermal_solver 0 and 1, infil_fix 0 and 1)") {
  for (int ts : {0, 1}) for (int fx : {0, 1}) {
    CAPTURE(ts); CAPTURE(fx);
    ut::Column c(gcol(14, 0.05, 0.4, 0.3));
    EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&c.g); se.thermal_solver = ts;
    se.richards_infil_fix = fx; se.richards.dt_fix = 0; se.richards.mass_guard = 1;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    const double W0 = gsum(c);
    const double infil = 20.0 / SEC_IN_DAY;
    double tr[MAX_SOI_LAY] = {0}; Layer* dl = c.L.back(); const double dd = dl->z + dl->dz;
    se.richards.update(c.L[0], dl, dd, 0.0, dd, tr, 0.0, infil, 0.0);
    REQUIRE(se.richards.guard_tripped);
    // the Soil_Env branch: (richards_infil_fix == 1 || guard_tripped) && early_return -> re-book
    if ((se.richards_infil_fix == 1 || se.richards.guard_tripped) && se.richards.early_return) se.rebookUnappliedRichards(infil, 0.0);
    se.checkSoilLiquidWaterValidity(c.L[0], c.L[0]->solind);
    const double err = gsum(c) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + se.richards.qdrain - W0 - 20.0;
    MESSAGE("ts " << ts << " infil_fix " << fx << ": puddle " << ed.d_soi2l.magic_puddle << " qover " << ed.d_soi2l.qover << " error " << err << " mm");
    CHECK(fabs(err) < 1e-9);
  }
}

}
