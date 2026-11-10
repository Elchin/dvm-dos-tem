// enthalpy/substep-fix: (A) with thermal_solver 1 a fully frozen soil cell keeps its freezing-curve
// unfrozen water in Soil_Env::checkSoilLiquidWaterValidity (it was dumped to puddle/qover when the
// whole column was frozen, i.e. no frozen == 0 layer; legacy thermal_solver 0 behaviour unchanged);
// (B) the enthalpy residual tolerance is a per-day budget (tol / thermal_substeps per step).
// Design: tools/enthalpy_bench/enthalpy_solver_design.md, section 15d.
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>

using ut::make_soil_layer;

namespace {
double colLiq(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (l->isSoil) w += l->liq; return w; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }
// wet column, frozen throughout after a cold enthalpy step (frozen == 1 everywhere, SJ2016 liquid > 0)
struct FrozenCol {
  ut::Column col; Soil_Env se; EnvData ed;
  FrozenCol()
    : col(std::vector<SoilLayer*>{
        make_soil_layer(0.0, 0.1, 0.5, 0.8 * 0.5 * DENLIQ * 0.1, 0.0, -2.0, -1, 0.0),
        make_soil_layer(0.1, 0.1, 0.5, 0.8 * 0.5 * DENLIQ * 0.1, 0.0, -2.5, -1, 0.0),
        make_soil_layer(0.2, 0.2, 0.5, 0.8 * 0.5 * DENLIQ * 0.2, 0.0, -3.0, -1, 0.0),
        make_soil_layer(0.4, 0.3, 0.4, 0.8 * 0.4 * DENLIQ * 0.3, 0.0, -3.5, -1, 0.0),
        make_soil_layer(0.7, 0.5, 0.4, 0.8 * 0.4 * DENLIQ * 0.5, 0.0, -4.0, -1, 0.0)}) {
    se.setGround(&col.g); se.setEnvData(&ed);
    se.thermal_solver = 1; se.thermal_sfcc = 1;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    EnergyDiag q; se.enthalpyThermalStep(-10.0, 0.0, q);
  }
};
}

TEST_SUITE("enthalpy_substep_fix") {

TEST_CASE("(A) whole-column freeze: thermal_solver 1 keeps the unfrozen water of frozen cells") {
  FrozenCol c;
  Ground& g = c.col.g;
  for (Layer* l = g.toplayer; l; l = l->nextl) { REQUIRE(l->frozen == 1); }
  REQUIRE(g.frontstype.empty());
  const double liq0 = colLiq(g), w0 = colW(g);
  REQUIRE(liq0 > 1.0);   // SJ2016 unfrozen water at -2..-5 C (kg m-2)
  std::vector<double> lq; for (Layer* l = g.toplayer; l; l = l->nextl) lq.push_back(l->liq);
  c.se.checkSoilLiquidWaterValidity(g.fstsoill, g.fstsoill->solind);
  CHECK(colLiq(g) == doctest::Approx(liq0).epsilon(1e-15));
  CHECK(colW(g) == doctest::Approx(w0).epsilon(1e-15));
  int k = 0; for (Layer* l = g.toplayer; l; l = l->nextl, ++k) { CHECK(l->liq == lq[k]); }
  CHECK(c.ed.d_soi2l.magic_puddle == 0.0); CHECK(c.ed.d_soi2l.qover == 0.0);
  // the next thermal step then needs no re-equilibration (no latent cooling from re-made liquid)
  EnergyDiag q; c.se.enthalpyThermalStep(-10.0, 0.0, q);
  CHECK(fabs(q.equil) < 1e-6);
}

TEST_CASE("(A) legacy rule unchanged with thermal_solver 0: frozen-cell liquid leaves as puddle/runoff") {
  FrozenCol c;
  Ground& g = c.col.g;
  const double liq0 = colLiq(g);
  c.se.thermal_solver = 0;
  c.se.checkSoilLiquidWaterValidity(g.fstsoill, g.fstsoill->solind);
  CHECK(colLiq(g) == doctest::Approx(0.0).scale(1.0));
  CHECK(c.ed.d_soi2l.magic_puddle + c.ed.d_soi2l.qover == doctest::Approx(liq0));
}

TEST_CASE("(B) per-day tolerance: tol / nsub per substep, restored after the step; closure at 1..48 substeps") {
  for (int ns : {1, 4, 12, 24, 48}) {
    CAPTURE(ns);
    FrozenCol c; c.se.thermal_substeps = ns;
    const double tol = c.se.enth_solver.opt.tol;
    EnergyDiag q; c.se.enthalpyThermalStep(5.0, 0.0, q);
    CHECK(c.se.enth_solver.opt.tol == tol);
    CHECK(q.fallbacks == 0);
    CHECK(fabs(q.res) <= 5.0 * tol + 1e-9);            // 5 cells x per-day tol (not x nsub)
    CHECK(fabs(q.res) / 86400.0 < 1e-6);
  }
}

}
