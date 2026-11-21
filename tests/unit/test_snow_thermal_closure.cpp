// snow_thermal_closure / snow_enth_trace (thermal_solver 1, default off, bit-identical when off).
// Leak (CMT06 1945 d147, found with snow_enth_trace): enthApplyMeltDirect's cold-content branch warms a snow cell to a
// T inside the free-water band (-eps, 0) and returns without putting liquid back on the curve; the solver then sees
// H(T) above the layer enthalpy and the scatter melts ice without charging Lf. docs/WORK_NOTES_snow_thermal_closure.md
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/ModelData.h"
#include <cmath>
#include <map>
#include <sstream>
#include <iostream>

using ut::make_soil_layer;

namespace {
Json::Value msj(const std::map<std::string, Json::Value>& kv) {
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0;
  for (auto& p : kv) { c["model_settings"][p.first] = p.second; }
  return c;
}
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }
SnowLayer* snow(double dz, double ice, double tem, double rho) {
  SnowLayer* s = new SnowLayer(); s->dz = dz; s->ice = ice; s->liq = 0.0; s->tem = tem; s->rho = rho; s->z = 0.0; return s;
}
struct StcCol {
  ut::Column col; Soil_Env se;
  StcCol(const std::vector<SnowLayer*>& sn, int closure, int trace = 0)
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
    se.thermal_solver = 1; se.thermal_sfcc = 1; se.thermal_substeps = 4; se.snowmelt_heat_cap = 1;
    se.snow_thermal_closure = closure; se.snow_enth_trace = trace;
    se.enth_init_pending = false;
  }
};
// top cell melts out under E, the remainder warms the second (cold, ice-only) cell into the free-water band
std::vector<SnowLayer*> leakpack() { return {snow(0.002, 0.4, -1.0, 200.0), snow(0.05, 10.0, -0.02, 200.0)}; }
double leakE() {
  const double ctop = (Soil_Env::enthLayerCs(nullptr) + 0.4 * SHCICE) * 1.0;   // cold content of the top cell (Cs = 0 for snow)
  return ctop + 0.4 * LHFUS + 300.0;                                           // + 300 J for the second cell (need 420 J)
}
}

TEST_SUITE("snow_thermal_closure") {

TEST_CASE("model_settings: defaults off; parsing; FATAL without thermal_solver 1 or snowmelt_heat_cap; trace range") {
  ModelData md; CHECK(md.snow_thermal_closure == false); CHECK(md.snow_enth_trace == 0);
  ModelData mj(msj({})); CHECK(mj.snow_thermal_closure == false); CHECK(mj.snow_enth_trace == 0);
  ModelData on(msj({{"thermal_solver", 1}, {"snowmelt_heat_cap", true}, {"snow_thermal_closure", true}, {"snow_enth_trace", 1}}));
  CHECK(on.snow_thermal_closure == true); CHECK(on.snow_enth_trace == 1);
  CHECK(ut::exit_status_of([] { ModelData m(msj({{"snow_thermal_closure", true}})); }) == EXIT_FAILURE);   // solver 0
  CHECK(ut::exit_status_of([] { ModelData m(msj({{"thermal_solver", 1}, {"snow_thermal_closure", true}})); }) == EXIT_FAILURE);   // no heat cap
  CHECK(ut::exit_status_of([] { ModelData m(msj({{"snow_enth_trace", 1}})); }) == EXIT_FAILURE);            // solver 0
  CHECK(ut::exit_status_of([] { ModelData m(msj({{"thermal_solver", 1}, {"snow_enth_trace", 3}})); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msj({{"thermal_solver", 1}, {"snow_thermal_closure", 2}})); }) == EXIT_FAILURE);
}

TEST_CASE("reproducer: legacy path leaks the latent heat of the second cell; the switch closes the day and books it in snleak") {
  StcCol off(leakpack(), 0), on(leakpack(), 1);
  const double H0 = colH(off.col.g), W0 = colW(off.col.g);
  EnergyDiag qo, qn;
  off.se.enthalpyThermalStep(0.0, leakE() / 86400.0, qo);
  on.se.enthalpyThermalStep(0.0, leakE() / 86400.0, qn);
  CHECK(fabs(qo.res) > 1.0e5);                                     // the leak (about 10 kg x Lf x liquid fraction)
  CHECK(qo.snleak == 0.0); CHECK(qo.snleak_n == 0.0);
  CHECK(fabs(qn.res) / 86400.0 <= 1e-6);                           // closed
  CHECK(qn.snleak_n >= 1.0);
  CHECK(qn.snleak == doctest::Approx(qo.res).epsilon(1e-6));       // snleak = the enthalpy the legacy path creates
  CHECK(colW(on.col.g) == doctest::Approx(W0).epsilon(1e-14));
  CHECK(colH(on.col.g) - H0 == doctest::Approx(qn.init + qn.equil + qn.gtop + qn.gbot + qn.src + qn.res + qn.clip).epsilon(1e-12).scale(1.0));
  // with the fix the pack keeps its ice: less melt than the leaking path
  double io = 0, in = 0;
  for (Layer* l = off.col.g.toplayer; l && l->isSnow; l = l->nextl) io += l->ice;
  for (Layer* l = on.col.g.toplayer; l && l->isSnow; l = l->nextl) in += l->ice;
  CHECK(in > io);
}

TEST_CASE("bit-identical when off, and on days that do not leak: every layer state and the budget terms are equal") {
  struct P { std::vector<SnowLayer*> (*mk)(); double qm; };
  auto p3 = [] { return std::vector<SnowLayer*>{snow(0.03, 6.0, -0.5, 200.0), snow(0.05, 10.0, -1.0, 200.0), snow(0.10, 20.0, -2.0, 200.0)}; };
  auto p1 = [] { return std::vector<SnowLayer*>{snow(0.05, 10.0, -0.02, 200.0)}; };
  auto pm = [] { return std::vector<SnowLayer*>{snow(0.002, 0.4, -1.0, 200.0), snow(0.003, 0.6, -0.5, 200.0)}; };   // melt-out, mret > 0
  std::vector<std::vector<SnowLayer*> (*)()> mks; mks.push_back(p3); mks.push_back(p1); mks.push_back(pm);
  for (size_t k = 0; k < mks.size(); ++k) for (double qm : {0.0, 40.0, 100.0, 400.0}) {
    CAPTURE(k); CAPTURE(qm);
    StcCol a(mks[k](), 0), b(mks[k](), 1);
    EnergyDiag qa, qb;
    a.se.enthalpyThermalStep(-0.5, qm, qa); b.se.enthalpyThermalStep(-0.5, qm, qb);
    if (qb.snleak_n == 0.0) {   // no leak on this day: switch on == switch off, bitwise
      CHECK(qa.res == qb.res); CHECK(qa.src == qb.src); CHECK(qa.gtop == qb.gtop); CHECK(qa.mmelt == qb.mmelt);
      Layer* x = a.col.g.toplayer; Layer* y = b.col.g.toplayer;
      for (; x && y; x = x->nextl, y = y->nextl) { CHECK(x->tem == y->tem); CHECK(x->liq == y->liq); CHECK(x->ice == y->ice); }
      CHECK(((x == nullptr) == (y == nullptr)));
    } else {
      CHECK(fabs(qb.res) / 86400.0 <= 1e-6);
    }
  }
}

TEST_CASE("snow_enth_trace is logging only: states identical with trace 0 / 1 / 2") {
  StcCol a(leakpack(), 1, 0), b(leakpack(), 1, 1), c(leakpack(), 1, 2);
  EnergyDiag qa, qb, qc;
  std::streambuf* old = std::cerr.rdbuf(); std::ostringstream sink; std::cerr.rdbuf(sink.rdbuf());
  a.se.enthalpyThermalStep(0.0, leakE() / 86400.0, qa);
  b.se.enthalpyThermalStep(0.0, leakE() / 86400.0, qb);
  c.se.enthalpyThermalStep(0.0, leakE() / 86400.0, qc);
  std::cerr.rdbuf(old);
  CHECK(qa.res == qb.res); CHECK(qa.res == qc.res);
  Layer* x = a.col.g.toplayer; Layer* y = b.col.g.toplayer; Layer* z = c.col.g.toplayer;
  for (; x; x = x->nextl, y = y->nextl, z = z->nextl) { CHECK(x->tem == y->tem); CHECK(x->tem == z->tem); CHECK(x->ice == y->ice); CHECK(x->liq == z->liq); }
  CHECK(sink.str().find("TRACE res") != std::string::npos);   // trace 2 prints every day
}

TEST_CASE("fuzz: random 1-4 layer packs, random melt energy and air T; with the switch on |res| gate holds, mass closes") {
  unsigned s = 12345u; auto rnd = [&s] { s = s * 1664525u + 1013904223u; return (s >> 8) / 16777216.0; };
  int nleak = 0;
  for (int it = 0; it < 200; ++it) {
    const int ns = 1 + (int)(rnd() * 4.0) % 4;
    std::vector<SnowLayer*> sn;
    for (int k = 0; k < ns; ++k) {
      const double ice = 0.3 + rnd() * 15.0, tt = (rnd() < 0.3) ? -rnd() * 0.02 : -rnd() * 4.0;
      sn.push_back(snow(ice / 200.0, ice, tt, 200.0));
    }
    StcCol c(sn, 1);
    const double W0 = colW(c.col.g), qm = rnd() * 150.0, ta = -rnd() * 3.0;
    EnergyDiag q;
    c.se.enthalpyThermalStep(ta, qm, q);
    CAPTURE(it);
    CHECK(fabs(q.res) / 86400.0 <= 1e-6);
    CHECK(colW(c.col.g) == doctest::Approx(W0).epsilon(1e-13));
    for (Layer* l = c.col.g.toplayer; l; l = l->nextl) { CHECK(l->ice >= 0.0); CHECK(l->liq >= 0.0); }
    if (q.snleak_n > 0.0) nleak++;
  }
  CHECK(nleak >= 0);
}

}
