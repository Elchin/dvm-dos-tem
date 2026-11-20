// snowmelt_heat_cap (thermal_solver 1, default off): the radiation/SEB snow-melt energy is applied to the snow
// cells before the solve (cold content first, then Lf x ice, top down) and capped at what the pack can take;
// the rest is returned to the surface (EnergyDiag::mret, not in src). docs/WORK_NOTES_snowmelt_heat_cap.md
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>

using ut::make_soil_layer;

namespace {
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }
double snowIce(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l && l->isSnow; l = l->nextl) w += l->ice; return w; }
SnowLayer* snow(double dz, double ice, double tem, double rho) {
  SnowLayer* s = new SnowLayer(); s->dz = dz; s->ice = ice; s->liq = 0.0; s->tem = tem; s->rho = rho; s->z = 0.0; return s;
}
struct CapCol {
  ut::Column col; Soil_Env se;
  CapCol(const std::vector<SnowLayer*>& sn, int cap)
    : col(std::vector<SoilLayer*>{
        make_soil_layer(0.0, 0.01, 0.9, 0.4 * 0.9 * DENLIQ * 0.01, 0.0, -0.5, -1, 0.0),
        make_soil_layer(0.01, 0.03, 0.5, 0.8 * 0.5 * DENLIQ * 0.03, 0.0, -1.0, -1, 0.0),
        make_soil_layer(0.04, 0.06, 0.5, 0.8 * 0.5 * DENLIQ * 0.06, 0.0, -1.5, -1, 0.0),
        make_soil_layer(0.10, 0.1, 0.5, 0.8 * 0.5 * DENLIQ * 0.1, 0.0, -2.0, -1, 0.0),
        make_soil_layer(0.20, 0.3, 0.4, 0.8 * 0.4 * DENLIQ * 0.3, 0.0, -3.0, -1, 0.0),
        make_soil_layer(0.50, 0.8, 0.4, 0.8 * 0.4 * DENLIQ * 0.8, 0.0, -4.0, -1, 0.0)}) {
    for (size_t k = sn.size(); k-- > 0;) { col.g.insertFront(sn[k]); }
    col.reindex(); col.g.fstsoill = col.L.front();
    se.setGround(&col.g);
    se.thermal_solver = 1; se.thermal_sfcc = 1; se.thermal_substeps = 4; se.snowmelt_heat_cap = cap;
    se.enth_init_pending = false;   // state given on the curve: no switch-on projection (T = 0 snow would turn liquid)
  }
};
void checkClosed(CapCol& c, EnergyDiag& q, double H0, double W0) {
  CHECK(q.fallbacks == 0);
  CHECK(fabs(q.res) / 86400.0 <= 1e-6);                 // daily closure gate (W m-2)
  CHECK(colW(c.col.g) == doctest::Approx(W0).epsilon(1e-14));
  CHECK(colH(c.col.g) - H0 == doctest::Approx(q.init + q.equil + q.gtop + q.gbot + q.src + q.res + q.clip).epsilon(1e-12).scale(1.0));
}
// three snow cells, top first (dz m, ice kg m-2, T C)
std::vector<SnowLayer*> pack(double t0, double t1, double t2) {
  return {snow(0.03, 6.0, t0, 200.0), snow(0.05, 10.0, t1, 200.0), snow(0.10, 20.0, t2, 200.0)};
}
}

TEST_SUITE("snowmelt_heat_cap") {

TEST_CASE("partial melt, isothermal 0 C pack: melt heat == Lf x melted mass, nothing returned, budget closes") {
  CapCol c(pack(-0.02, -0.02, -0.02), 1);   // just below the [-0.01, 0] band: cold content tiny
  const double H0 = colH(c.col.g), W0 = colW(c.col.g), I0 = snowIce(c.col.g);
  const double qm = 60.0;                                // W m-2: 5.18 MJ = 15.5 kg, less than the 36 kg pack
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, qm, q);
  CHECK(q.mcc < 1e-3 * qm * 86400.0);                                     // cold content negligible
  CHECK(q.mmelt * LHFUS + q.mcc == doctest::Approx(qm * 86400.0).epsilon(1e-12));   // all of it melts snow
  CHECK(q.mdir == doctest::Approx(q.mcc + LHFUS * q.mmelt).epsilon(1e-12));
  CHECK(q.mret == 0.0);
  CHECK(q.src == doctest::Approx(q.mdir).epsilon(1e-9));                    // no S[0] heat source in the solve
  CHECK(q.mmelt < I0);                                                      // partial melt
  CHECK(I0 - snowIce(c.col.g) == doctest::Approx(q.mmelt).epsilon(0.05));   // net snow melt after the solve
  checkClosed(c, q, H0, W0);
}

TEST_CASE("partial melt, cold pack: cold content first, then melt; mdir == cold content + Lf x melted mass") {
  CapCol c(pack(-2.0, -4.0, -6.0), 1);
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  Layer* top = c.col.g.toplayer;
  const double Ctop = Soil_Env::enthLayerCs(top) + top->ice * SHCICE;
  EnergyDiag q;
  c.se.enthalpyThermalStep(-0.5, 40.0, q);
  CHECK(q.mcc >= Ctop * 2.0 * (1.0 - 1e-12));           // at least the top cell's cold content
  CHECK(q.mmelt > 0.0);
  CHECK(q.mdir == doctest::Approx(q.mcc + LHFUS * q.mmelt).epsilon(1e-12));
  CHECK(q.mdir + q.mret == doctest::Approx(40.0 * 86400.0).epsilon(1e-12));
  CHECK(q.mret == 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("melt-out: the pack takes cold content + Lf x ice, the rest is returned to the surface and kept out of src") {
  CapCol c({snow(0.002, 0.4, -1.0, 200.0), snow(0.003, 0.6, -0.5, 200.0)}, 1);
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  double cc = 0.0; for (Layer* l = c.col.g.toplayer; l && l->isSnow; l = l->nextl) cc += -(Soil_Env::enthLayerCs(l) + l->ice * SHCICE) * l->tem;
  const double E = 100.0 * 86400.0;
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 100.0, q);
  CHECK(q.mmelt == doctest::Approx(1.0).epsilon(1e-12));
  CHECK(q.mcc == doctest::Approx(cc).epsilon(1e-12));
  CHECK(q.mret == doctest::Approx(E - cc - LHFUS * 1.0).epsilon(1e-12));
  CHECK(q.mret > 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("an emptied snow cell on top is skipped (no heat to a zero-capacity cell)") {
  CapCol c({snow(0.0, 0.0, -1.0e-8, 200.0), snow(0.05, 10.0, -1.0, 200.0)}, 1);
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 30.0, q);
  CHECK(q.mdir + q.mret == doctest::Approx(30.0 * 86400.0).epsilon(1e-12));
  CHECK(Soil_Env::enthLayerH(c.col.g.toplayer) == 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("switch off: legacy S[0] heat source, diagnostics untouched") {
  CapCol c(pack(-0.5, -1.0, -2.0), 0);
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  c.se.enthalpyThermalStep(0.0, 60.0, q);
  CHECK(q.src == doctest::Approx(60.0 * 86400.0));
  CHECK(q.mdir == 0.0); CHECK(q.mret == 0.0); CHECK(q.mmelt == 0.0); CHECK(q.mcc == 0.0);
  checkClosed(c, q, H0, W0);
}

TEST_CASE("no snow: nothing applied with the switch on") {
  CapCol c({}, 1);
  EnergyDiag q;
  c.se.enthalpyThermalStep(2.0, 50.0, q);
  CHECK(q.src == 0.0); CHECK(q.mdir == 0.0); CHECK(q.mret == 0.0);
}

}
