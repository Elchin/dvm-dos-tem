// model_settings.richards_dt_fix (bug R1, found by test_richards_water_energy.cpp "bug R1").
// Richards::update halves dtsub when a trial substep is rejected, but legacy re-solves with the
// dt_dz = dtsub/dz of the rejected substep while dtdone advances by the halved dtsub, so the
// day's fluxes are applied for longer than a day. richards_dt_fix = 1 recomputes dt_dz before
// every retry. Default 0 = legacy (bit-identical). Applies to both thermal solvers (Richards is
// called by Soil_Env::updateDailySM for thermal_solver 0 and 1); the thermal_solver-dependent
// part of the water path (checkSoilLiquidWaterValidity) is exercised with 0 and 1 below.
// Refs: CLM5 tech note ch. 7 (Lawrence et al. 2018); CTSM SoilWaterMovementMod (the trial loop
//   passes the current dtsub to compute_LHS/RHS).
#include "doctest.h"
#include "test_helpers.h"
#include <json/value.h>
#include <cmath>

using ut::make_soil_layer;

namespace {
std::vector<SoilLayer*> col(int n, double dz, double poro, double s0, double tem = 5.0) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) {
    SoilLayer* l = make_soil_layer(i * dz, dz, poro, s0 * poro * DENLIQ * dz, 0.0, tem, -1, 0.0);
    l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0;
    v.push_back(l);
  }
  return v;
}
double wsum(const ut::Column& c) { double w = 0; for (auto l : c.L) w += l->liq + l->ice; return w; }
void rday(Richards& rr, Layer* top, Layer* drainl, double infil_mmd) {
  double tr[MAX_SOI_LAY] = {0};
  const double dd = drainl->z + drainl->dz;
  rr.update(top, drainl, dd, 0.0, dd, tr, 0.0, infil_mmd / SEC_IN_DAY, 0.0);
}
} // namespace

TEST_SUITE("richards_dt_fix") {

TEST_CASE("richards_dt_fix defaults to 0 (legacy) in Richards and ModelData; config value 1 is read") {
  Richards rr; CHECK(rr.dt_fix == 0);
  Soil_Env se; CHECK(se.richards.dt_fix == 0);
  ModelData m0; CHECK(m0.richards_dt_fix == 0);
  Json::Value c0; c0["model_settings"]["svp_ice_fix"] = 0; ModelData mj0(c0); CHECK(mj0.richards_dt_fix == 0);
  Json::Value c1; c1["model_settings"]["richards_dt_fix"] = 1; ModelData mj1(c1); CHECK(mj1.richards_dt_fix == 1);
}

TEST_CASE("bug R1: with default rejection, dt_fix 0 creates water; dt_fix 1 closes to round-off") {
  double worst[2] = {0, 0}, tot[2] = {0, 0}; long rej[2] = {0, 0};
  for (int fx : {0, 1}) {
    CAPTURE(fx);
    ut::Column c(col(20, 0.05, 0.4, 0.3));
    Richards rr; rr.dt_fix = fx;                 // default toler_upper 0.1 -> rejections occur
    for (int d = 0; d < 5; ++d) {
      const double m = wsum(c); rday(rr, c.L[0], c.L.back(), 20.0);
      const double e = wsum(c) - m - 20.0 + rr.qdrain;
      worst[fx] = fmax(worst[fx], fabs(e)); tot[fx] += e;
    }
    rej[fx] = rr.n_rejections;
    MESSAGE("dt_fix " << fx << ": rejections " << rr.n_rejections << " in " << rr.n_reject_calls << "/" << rr.n_calls
            << " days; worst daily mass error " << worst[fx] << " mm, total " << tot[fx] << " mm");
    CHECK(rr.n_rejections > 0);                  // the retry path is exercised
  }
  CHECK(worst[0] > 1.0);                         // legacy: about 3.3 mm/d created (documents the bug)
  CHECK(tot[0] > 0.0);
  CHECK(worst[1] < 1e-9);                        // fixed: round-off
}

TEST_CASE("bug R1: rejection on vs off agree closely with dt_fix 1 (same physics, finer substeps)") {
  ut::Column a(col(20, 0.05, 0.4, 0.3)), b(col(20, 0.05, 0.4, 0.3));
  Richards ra; ra.dt_fix = 1; Richards rb; rb.dt_fix = 1; rb.toler_upper = 1e30;
  for (int d = 0; d < 5; ++d) { rday(ra, a.L[0], a.L.back(), 20.0); rday(rb, b.L[0], b.L.back(), 20.0); }
  CHECK(ra.n_rejections > 0); CHECK(rb.n_rejections == 0);
  CHECK(fabs(wsum(a) - wsum(b)) < 1e-8);        // both conserve
  double dmax = 0; for (int i = 0; i < 20; ++i) dmax = fmax(dmax, fabs(a.L[i]->liq - b.L[i]->liq) / b.L[i]->liq);
  MESSAGE("max relative layer-liquid difference, rejection on vs off: " << dmax);
  CHECK(dmax < 0.15);
}

TEST_CASE("bug R1 through Soil_Env (thawed over frozen column, validity check): thermal_solver 0 and 1, dt_fix 0 and 1") {
  for (int ts : {0, 1}) for (int fx : {0, 1}) {
    CAPTURE(ts); CAPTURE(fx);
    ut::Column c(col(14, 0.05, 0.4, 0.3));
    for (int i = 10; i < 14; ++i) {                // frozen below the drainage layer (as with a front)
      Layer* l = c.L[i]; l->frozen = 1; l->frozenfrac = 1.0; l->tem = -1.0;
      l->ice = 0.6 * l->maxice; l->liq = 0.1 * l->maxliq;
    }
    EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&c.g); se.thermal_solver = ts;
    se.richards.dt_fix = fx;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    const double W0 = wsum(c); double qd = 0.0, worst = 0.0;
    for (int d = 0; d < 5; ++d) {
      const double wb = wsum(c) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + qd;
      rday(se.richards, c.L[0], c.L[9], 20.0); qd += se.richards.qdrain;
      se.checkSoilLiquidWaterValidity(c.L[0], c.L[0]->solind);
      worst = fmax(worst, fabs(wsum(c) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + qd - wb - 20.0));
    }
    const double err = wsum(c) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + qd - W0 - 100.0;
    MESSAGE("ts " << ts << " dt_fix " << fx << ": rejections " << se.richards.n_rejections << ", worst daily error "
            << worst << " mm, 5-day error " << err << " mm");
    CHECK(se.richards.n_rejections > 0);
    if (fx == 0) CHECK(worst > 1.0); else { CHECK(worst < 1e-9); CHECK(fabs(err) < 1e-9); }
    for (int i = 10; i < 14; ++i) CHECK(c.L[i]->ice == doctest::Approx(0.6 * c.L[i]->maxice));  // frozen cells untouched
  }
}

}
