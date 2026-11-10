// fix/meltout-closure (thermal_solver 1 only): heat of a fully melted snow cell must reach a cell
// that can hold it. A zero-mass snow cell (W = 0, no solid: an emptied snow layer the snow module
// has not removed yet) has zero heat capacity; equilibrate() returns T = 0 for it, so heat handed
// to it was lost. Seen in Calibration-expert's r1 case (TVC CMT51, snwdennew 75, nfactor(w) 0.3):
// 1 Jun 2024, a 0.009 mm SWE remnant over an empty snow cell melts out inside the solve,
// res = -57.36 J m-2 (6.6e-4 W m-2, gate 1e-6). Design: enthalpy_solver_design.md section 15g.
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>

using ut::make_soil_layer;

namespace {
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }
SnowLayer* snow(double dz, double ice, double tem, double rho) {
  SnowLayer* s = new SnowLayer(); s->dz = dz; s->ice = ice; s->liq = 0.0; s->tem = tem; s->rho = rho; s->z = 0.0; return s;
}
// warm, wet 6-layer soil column (as test_enthalpy_integration Col) under the given snow cells (top first)
struct MeltCol {
  ut::Column col; Soil_Env se;
  MeltCol(const std::vector<SnowLayer*>& sn)
    : col(std::vector<SoilLayer*>{
        make_soil_layer(0.0, 0.01, 0.9, 0.4 * 0.9 * DENLIQ * 0.01, 0.0, 1.4, -1, 0.0),
        make_soil_layer(0.01, 0.03, 0.5, 0.8 * 0.5 * DENLIQ * 0.03, 0.0, 0.9, -1, 0.0),
        make_soil_layer(0.04, 0.06, 0.5, 0.8 * 0.5 * DENLIQ * 0.06, 0.0, -0.001, -1, 0.0),
        make_soil_layer(0.10, 0.1, 0.5, 0.8 * 0.5 * DENLIQ * 0.1, 0.0, -0.2, -1, 0.0),
        make_soil_layer(0.20, 0.3, 0.4, 0.8 * 0.4 * DENLIQ * 0.3, 0.0, -1.0, -1, 0.0),
        make_soil_layer(0.50, 0.8, 0.4, 0.8 * 0.4 * DENLIQ * 0.8, 0.0, -2.0, -1, 0.0)}) {
    for (size_t k = sn.size(); k-- > 0;) { col.g.insertFront(sn[k]); }
    col.reindex(); col.g.fstsoill = col.L.front();
    se.setGround(&col.g);
    se.thermal_solver = 1; se.thermal_sfcc = 1; se.thermal_substeps = 4;
  }
};
void checkClosed(MeltCol& c, EnergyDiag& q, double H0, double W0) {
  CHECK(q.fallbacks == 0);
  CHECK(fabs(q.res) / 86400.0 <= 1e-6);                 // daily closure gate (W m-2)
  CHECK(fabs(q.res) < 1e-3);                            // J m-2: solver tolerance scale, not the lost heat
  CHECK(colW(c.col.g) == doctest::Approx(W0).epsilon(1e-14));
  CHECK(colH(c.col.g) - H0 == doctest::Approx(q.init + q.equil + q.gtop + q.gbot + q.src + q.res + q.clip).epsilon(1e-12).scale(1.0));
}
}

TEST_SUITE("enthalpy_meltout") {

TEST_CASE("near-massless snow remnant melts out above an empty snow cell: heat reaches the soil, closure <= 1e-6 W/m2") {
  // r1 1 Jun 2024 state: top 0.16 mm / 0.009 kg m-2 at -0.01 C, empty bottom snow cell, melt source 136.6 W m-2
  SnowLayer* top = snow(1.5749e-4, 9.020455e-3, -0.009999366, 57.28);
  SnowLayer* emp = snow(0.0, 0.0, -1.25e-8, 209.87);
  MeltCol c({top, emp});
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  const double hsoil0 = Soil_Env::enthLayerH(c.col.L.front());
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 136.59, q);
  CHECK(q.src == doctest::Approx(136.59 * 86400.0));
  CHECK(top->ice == 0.0); CHECK(top->tem == 0.0);       // melted cell left at 0 C, its heat moved down
  CHECK(Soil_Env::enthLayerH(emp) == 0.0);              // the empty cell holds no heat
  CHECK(Soil_Env::enthLayerH(c.col.L.front()) != hsoil0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("several empty snow cells and a melted cascade above them still close") {
  SnowLayer* a = snow(1.0e-4, 5.0e-3, -0.5, 50.0);
  SnowLayer* b = snow(2.0e-4, 2.0e-2, -0.2, 100.0);
  SnowLayer* e1 = snow(0.0, 0.0, 0.0, 200.0);
  SnowLayer* e2 = snow(0.0, 0.0, -3.0, 200.0);
  MeltCol c({a, b, e1, e2});
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  c.se.enthalpyThermalStep(2.0, 200.0, q);
  CHECK(a->ice == 0.0); CHECK(b->ice == 0.0);
  CHECK(a->tem == 0.0); CHECK(b->tem == 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("empty snow cell between two snow cells: the lower snow cell receives the heat") {
  SnowLayer* a = snow(1.0e-4, 5.0e-3, -0.5, 50.0);
  SnowLayer* e = snow(0.0, 0.0, 0.0, 200.0);
  SnowLayer* d = snow(0.05, 10.0, -2.0, 200.0);        // deep cold cell: keeps ice
  MeltCol c({a, e, d});
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 120.0, q);
  CHECK(a->ice == 0.0); CHECK(a->tem == 0.0);
  CHECK(d->ice > 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("no empty cell below: unchanged path (heat to the next cell), closes") {
  SnowLayer* top = snow(1.5749e-4, 9.020455e-3, -0.009999366, 57.28);
  MeltCol c({top});
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 136.59, q);
  CHECK(top->ice == 0.0); CHECK(top->tem == 0.0);
  checkClosed(c, q, H0, W0);
}

}
