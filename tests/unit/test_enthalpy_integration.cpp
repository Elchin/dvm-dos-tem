// Option C (enthalpy/stageC): integration of the enthalpy thermal solver behind
// model_settings.thermal_solver (0 = legacy default, 1 = enthalpy). Checks: switch parsing and
// FATAL combinations; the shared SJ2016 curve (respiration == solver, bitwise); the H-conserving
// helpers (equilibrate, mergeT); front reconstruction from the profile (legacy front convention);
// one daily step on a small column (energy closure, on-curve state, fronts, melt source, snow
// excess heat); enthalpy-conserving layer split/merge (Ground::enthalpy_remap).
// Refs: Tubini et al. 2021, doi:10.5194/tc-15-2541-2021; Schaefer & Jafarov 2016,
//   doi:10.5194/bg-13-1991-2016. Design: tools/enthalpy_bench/enthalpy_solver_design.md.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>
#include <cstring>
#include "../../include/SJCurve.h"

using ut::make_soil_layer;

namespace {
Json::Value ms(const std::map<std::string, Json::Value>& kv) {
  Json::Value c; c["model_settings"]["energy_diag"] = 0;
  for (auto& p : kv) c["model_settings"][p.first] = p.second;
  return c;
}
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }

// 6 wet mineral layers (no texture -> organic curve), T from t0 at the top increasing by dt per layer
struct Col {
  ut::Column col; Soil_Env se;
  Col(double t0, double dtl, double sat = 0.8)
    : col(std::vector<SoilLayer*>{
        make_soil_layer(0.0, 0.1, 0.5, sat * 0.5 * DENLIQ * 0.1, 0.0, t0, -1, 0.0),
        make_soil_layer(0.1, 0.1, 0.5, sat * 0.5 * DENLIQ * 0.1, 0.0, t0 + dtl, -1, 0.0),
        make_soil_layer(0.2, 0.2, 0.5, sat * 0.5 * DENLIQ * 0.2, 0.0, t0 + 2 * dtl, -1, 0.0),
        make_soil_layer(0.4, 0.3, 0.4, sat * 0.4 * DENLIQ * 0.3, 0.0, t0 + 3 * dtl, -1, 0.0),
        make_soil_layer(0.7, 0.5, 0.4, sat * 0.4 * DENLIQ * 0.5, 0.0, t0 + 4 * dtl, -1, 0.0),
        make_soil_layer(1.2, 0.8, 0.4, sat * 0.4 * DENLIQ * 0.8, 0.0, t0 + 5 * dtl, -1, 0.0)}) {
    se.setGround(&col.g);
    se.thermal_solver = 1; se.thermal_sfcc = 1;
  }
};
}

TEST_SUITE("enthalpy_integration") {

TEST_CASE("thermal_solver switches: defaults are legacy; members default off") {
  ModelData md; CHECK(md.thermal_solver == 0); CHECK(md.thermal_sfcc == 1); CHECK(md.thermal_face_mode == 0);
  CHECK(md.thermal_nl_solver == 0); CHECK(md.thermal_substeps == 1); CHECK(md.geothermal_flux_mode == 0);
  CHECK(md.geothermal_flux == 0.0);
  ModelData mj(ms({})); CHECK(mj.thermal_solver == 0); CHECK(mj.geothermal_flux_mode == 0);
  Soil_Env se; CHECK(se.thermal_solver == 0); CHECK(se.geothermal_flux_Wm2 == 0.0);
  Snow_Env sn; CHECK(sn.thermal_solver == 0);
  Ground g; CHECK(g.enthalpy_remap == 0);
}

TEST_CASE("thermal_solver switches: valid combinations parse, invalid ones are FATAL") {
  CHECK(ModelData(ms({{"thermal_solver", 1}})).thermal_solver == 1);
  ModelData a(ms({{"thermal_solver", 1}, {"thermal_sfcc", 0}, {"thermal_face_mode", 1}, {"thermal_substeps", 4},
                  {"geothermal_flux_mode", 1}, {"geothermal_flux", 0.06}}));
  CHECK(a.thermal_sfcc == 0); CHECK(a.thermal_face_mode == 1); CHECK(a.thermal_substeps == 4);
  CHECK(a.geothermal_flux == doctest::Approx(0.06));
  CHECK(ModelData(ms({{"thermal_solver", 1}, {"thermal_sfcc", 0}, {"thermal_nl_solver", 1}})).thermal_nl_solver == 1);
  // stefan_front_resist_fix is ignored (forced 0, warning) with the enthalpy solver
  CHECK(ModelData(ms({{"thermal_solver", 1}, {"stefan_front_resist_fix", 1}})).stefan_front_resist_fix == 0);
  CHECK(ModelData(ms({{"thermal_solver", 0}, {"stefan_front_resist_fix", 1}})).stefan_front_resist_fix == 1);
  typedef std::map<std::string, Json::Value> M;
  std::vector<M> bad = {
    M{{"thermal_solver", 2}}, M{{"thermal_solver", 0.5}}, M{{"thermal_solver", true}},
    M{{"thermal_sfcc", 0}},                                   // solver-1 key with solver 0
    M{{"thermal_face_mode", 1}}, M{{"thermal_substeps", 4}}, M{{"thermal_nl_solver", 1}},
    M{{"geothermal_flux_mode", 1}, {"geothermal_flux", 0.05}},  // geothermal needs solver 1
    M{{"thermal_solver", 1}, {"thermal_nl_solver", 1}},       // Katzenelson needs sfcc 0
    M{{"thermal_solver", 1}, {"thermal_substeps", 3}},
    M{{"thermal_solver", 1}, {"geothermal_flux", 0.05}},      // value without mode 1
    M{{"thermal_solver", 1}, {"geothermal_flux_mode", 1}, {"geothermal_flux", 20.0}},
    M{{"thermal_solver", 1}, {"geothermal_flux_mode", 2}},    // map without IO.geothermal_flux_file
    M{{"thermal_solver", 1}, {"geothermal_flux_mode", 3}}};
  for (size_t k = 0; k < bad.size(); ++k) {
    CAPTURE(k); M b = bad[k];
    CHECK(ut::exit_status_of([b] { ModelData m(ms(b)); }) == EXIT_FAILURE);
  }
  Json::Value c = ms({{"thermal_solver", 1}, {"geothermal_flux_mode", 2}});
  c["IO"]["geothermal_flux_file"] = "geo.nc";
  ModelData g(c); CHECK(g.geothermal_flux_mode == 2); CHECK(g.geothermal_flux_file == "geo.nc");
}

TEST_CASE("SJ2016 curve: one shared function, respiration == solver (bitwise)") {
  const double ts[] = {-30.0, -5.0, -1.0, -0.3, -0.05, -0.0101, 0.0, 0.2};
  const double tex[][4] = {{0.3, 0.2, 0.3, 0.5}, {1.0, 0, 0, 0}, {0.0, 30, 40, 30}, {0.2, 0, 0, 0}};
  for (double t : ts) for (auto& x : tex) {
    const double a = Soil_Bgc::sjLayerPhi(t, x[0], x[1], x[2], x[3]);
    const double b = sjcurve::layerPhi(t, x[0], x[1], x[2], x[3]);
    const double c = enthalpy::sjPhi(t, x[0], x[1], x[2], x[3]);
    CHECK(std::memcmp(&a, &b, sizeof a) == 0); CHECK(std::memcmp(&a, &c, sizeof a) == 0);
    // former Soil_Bgc arithmetic (stageB), reproduced literally
    auto pp = [](double tt, double bb) { return (tt > 0.0) ? 1.0 : pow((0.1 - tt) / 0.01, bb); };
    const double fo = fmin(1.0, fmax(0.0, x[0])), cl = fmax(0.0, x[1]), si = fmax(0.0, x[2]), sa = fmax(0.0, x[3]);
    const double sum = cl + si + sa;
    const double ref = (sum > 0.0) ? (1.0 - fo) * ((cl * pp(t, -0.3) + si * pp(t, -0.5) + sa * pp(t, -0.9)) / sum) + fo * pp(t, -1.0)
                                   : pp(t, -1.0);
    CHECK(std::memcmp(&a, &ref, sizeof a) == 0);
  }
  for (double dz : {0.0, 0.05, 0.3}) for (double c : {0.0, 500.0, 50000.0}) {
    const double a = Soil_Bgc::sjOrgFraction(c, dz), b = sjcurve::orgFraction(c, dz);
    CHECK(std::memcmp(&a, &b, sizeof a) == 0);
  }
  CHECK(Soil_Bgc::sjPhiPure(-2.0, -0.5) == sjcurve::phiPure(-2.0, -0.5));
}

TEST_CASE("equilibrate and mergeT conserve enthalpy") {
  enthalpy::Cell c; c.kind = enthalpy::SOIL; c.dz = 0.2; c.Cs = 2.0e6 * 0.5 * 0.2; c.W = 80.0; c.liqmin = 2.0;
  c.eps = 0.01; c.lam_u = 1.0; c.lam_f = 2.0; c.forg = 0.3; c.clay = 30; c.silt = 40; c.sand = 30;
  for (int sf : {0, 1}) {
    c.sfcc = sf;
    for (double T : {-8.0, -0.5, -0.005, 0.0, 3.0}) for (double fl : {0.0, 0.3, 1.0}) {
      CAPTURE(sf); CAPTURE(T); CAPTURE(fl);
      const double lq0 = fl * c.W, h = enthalpy::Hstate(c, T, lq0, c.W - lq0);
      double lq, ic; const double T1 = enthalpy::equilibrate(c, h, lq, ic);
      CHECK(enthalpy::Hstate(c, T1, lq, ic) == doctest::Approx(h).epsilon(1e-11).scale(1e3));
      CHECK(lq + ic == doctest::Approx(c.W));
      CHECK(enthalpy::onCurve(c, T1, lq, 1e-9));
    }
  }
  // merge: Cu*Tu + Cl*Tl = (Cu+Cl)*T (the latent part sums with the ice)
  const double T = enthalpy::mergeT(3.0e5, -4.0, 1.0e5, 2.0);
  CHECK(T == doctest::Approx(-2.5));
  CHECK(enthalpy::mergeT(0.0, 1.0, 0.0, 3.0) == doctest::Approx(2.0));
}

TEST_CASE("reconstructFronts: single front, ramp run, lens, sharp crossing, cap") {
  auto mk = [](std::vector<double> T, std::vector<double> ice, double Ta, double Tb) {
    enthalpy::FrontIn in; double z = 0;
    for (size_t i = 0; i < T.size(); ++i) {
      in.z.push_back(z); in.dz.push_back(0.1); z += 0.1;
      in.T.push_back(T[i]); in.ice.push_back(ice[i]); in.ice_eps.push_back(10.0); in.W.push_back(12.0); in.latent.push_back(1);
    }
    in.Tabove = Ta; in.Tbelow = Tb; in.eps = 0.01; in.maxfronts = MAX_NUM_FNT; return in;
  };
  enthalpy::FrontOut o;
  SUBCASE("summer: thawed over frozen, one ramp cell half frozen -> thawing front at its middle") {
    enthalpy::reconstructFronts(mk({3, 1, -0.005, -1, -2}, {0, 0, 5, 10, 10}, 5, -3), o);
    REQUIRE(o.fz.size() == 1); CHECK(o.ftype[0] == -1); CHECK(o.fz[0] == doctest::Approx(0.25));
    CHECK(o.frozen == std::vector<int>{-1, -1, 0, 1, 1}); CHECK(o.frozenfrac[2] == doctest::Approx(0.5));
  }
  SUBCASE("autumn: frozen over thawed, ramp run of two cells -> freezing front, frozen thickness preserved") {
    enthalpy::reconstructFronts(mk({-2, -0.002, -0.006, 1, 2}, {10, 2, 6, 0, 0}, -5, 3), o);
    REQUIRE(o.fz.size() == 1); CHECK(o.ftype[0] == 1); CHECK(o.fz[0] == doctest::Approx(0.1 + 0.08));
    CHECK(o.frozen[0] == 1); CHECK(o.frozen[1] == 0); CHECK(o.frozen[2] == -1); CHECK(o.frozen[3] == -1);
    CHECK(o.frozenfrac[1] == doctest::Approx(0.8));
  }
  SUBCASE("lens: ramp cell between thawed cells -> pair centred in the cell") {
    enthalpy::reconstructFronts(mk({2, -0.005, 2}, {0, 4, 0}, 5, 3), o);
    REQUIRE(o.fz.size() == 2); CHECK(o.ftype == std::vector<int>{-1, 1});
    CHECK(o.fz[0] == doctest::Approx(0.13)); CHECK(o.fz[1] == doctest::Approx(0.17));
    CHECK(o.frozen[1] == 0); CHECK(o.frozenfrac[1] == doctest::Approx(0.4));
  }
  SUBCASE("sharp crossing between non-ramp cells -> linear interpolation of the 0 C isotherm") {
    enthalpy::reconstructFronts(mk({2, 1, -3, -4}, {0, 0, 10, 10}, 5, -5), o);
    REQUIRE(o.fz.size() == 1); CHECK(o.ftype[0] == -1);
    CHECK(o.fz[0] == doctest::Approx(0.15 + 0.25 * 0.1));   // between midpoints 0.15 and 0.25
    CHECK(o.frozen[1] == 0); CHECK(o.frozen[0] == -1); CHECK(o.frozen[3] == 1);
  }
  SUBCASE("fully thawed or fully frozen columns: no fronts") {
    enthalpy::reconstructFronts(mk({1, 2, 3}, {0, 0, 0}, 5, 3), o);
    CHECK(o.fz.empty()); CHECK(o.frozen == std::vector<int>{-1, -1, -1});
    enthalpy::reconstructFronts(mk({-1, -2, -3}, {10, 10, 10}, -5, -3), o);
    CHECK(o.fz.empty()); CHECK(o.frozen == std::vector<int>{1, 1, 1});
  }
  SUBCASE("cap: at most maxfronts, the closest pairs are removed") {
    std::vector<double> T, ice;
    for (int i = 0; i < 14; ++i) { T.push_back(i % 2 ? -1.0 : 1.0); ice.push_back(i % 2 ? 10.0 : 0.0); }
    enthalpy::FrontIn in = mk(T, ice, 5, -5); in.maxfronts = 4;
    enthalpy::reconstructFronts(in, o);
    CHECK(o.fz.size() <= 4); CHECK(o.fz.size() % 2 == 1);
    for (size_t k = 0; k + 1 < o.fz.size(); ++k) { CHECK(o.fz[k] < o.fz[k + 1]); CHECK(o.ftype[k] == -o.ftype[k + 1]); }
  }
}

TEST_CASE("enthalpyThermalStep: closes the energy budget, puts cells on the curve, fronts from the profile") {
  for (int sf : {0, 1}) for (double ts : {-15.0, 10.0}) for (int face : {0, 1}) {
    CAPTURE(sf); CAPTURE(ts); CAPTURE(face);
    Col c(ts > 0 ? -2.0 : 3.0, ts > 0 ? -0.5 : 0.3);
    c.se.thermal_sfcc = sf; c.se.thermal_face_mode = face; c.se.geothermal_flux_Wm2 = 0.08;
    const double W0 = colW(c.col.g);
    EnergyDiag q;
    c.se.enthalpyThermalStep(ts, 0.0, q);
    CHECK(fabs(q.res) / 86400.0 < 1e-6); CHECK(q.fallbacks == 0); CHECK(q.solves > 0);
    CHECK(q.gbot == doctest::Approx(0.08 * 86400.0));
    CHECK(q.h1 - q.h0 == doctest::Approx(q.gtop + q.gbot + q.src + q.res + q.clip).epsilon(1e-12).scale(1.0));
    CHECK(colW(c.col.g) == doctest::Approx(W0).epsilon(1e-12));   // no ice cap hit (sat 0.8)
    // second day: state is on the curve, init projection does not repeat, equilibration ~ 0
    EnergyDiag q2; c.se.enthalpyThermalStep(ts, 0.0, q2);
    CHECK(q2.init == 0.0); CHECK(fabs(q2.equil) < 1e-6);
    // fronts consistent with the profile: ordered, alternating, inside the soil column
    Ground& g = c.col.g;
    for (size_t k = 0; k < g.frontsz.size(); ++k) {
      CHECK(g.frontsz[k] > 0.0); CHECK(g.frontsz[k] <= 2.0);
      if (k) { CHECK(g.frontsz[k] > g.frontsz[k - 1]); CHECK(g.frontstype[k] == -g.frontstype[k - 1]); }
    }
    for (Layer* l = g.toplayer; l; l = l->nextl) {
      if (l->frozen == 1) CHECK(l->tem <= 0.0);
      if (l->frozen == -1) CHECK(l->ice <= 1e-9 * (l->liq + l->ice + 1.0) + (l->tem < 0 ? l->liq + l->ice : 0.0));
    }
  }
}

TEST_CASE("enthalpyThermalStep: rock layers below the soil (ParentLayer heat capacity), geothermal flux") {
  Col c(-1.0, -0.4);
  for (double dz : {2.0, 4.0, 8.0}) { ParentLayer* r = new ParentLayer(dz); r->tem = -3.0; r->liq = 5.0; r->ice = 1.0; c.col.g.insertBack(r); }
  c.col.reindex();
  for (Layer* l = c.col.g.toplayer; l; l = l->nextl) if (l->isRock) CHECK(Soil_Env::enthLayerCs(l) == doctest::Approx(2.7e6 * l->dz));
  c.se.geothermal_flux_Wm2 = 0.06;
  for (int d = 0; d < 30; ++d) {
    EnergyDiag q; c.se.enthalpyThermalStep(-10.0, 0.0, q);
    REQUIRE(std::isfinite(q.h1)); CHECK(fabs(q.res) / 86400.0 < 1e-6); CHECK(q.gbot == doctest::Approx(0.06 * 86400.0));
  }
  for (Layer* l = c.col.g.toplayer; l; l = l->nextl) if (l->isRock) { CHECK(l->liq == 0.0); CHECK(l->ice == 0.0); CHECK(std::isfinite(l->tem)); }
}

TEST_CASE("enthalpyThermalStep: switch-on projection is at fixed T and booked as init") {
  Col c(-3.0, 0.0);   // frozen T but all water liquid (off the curve)
  EnergyDiag q; c.se.enthalpyThermalStep(-3.0, 0.0, q);
  CHECK(q.init < 0.0);                    // freezing at fixed T releases latent heat -> booked
  CHECK(c.se.enth_init_maxdliq > 1.0);
  CHECK(c.se.enth_init_pending == false);
  for (Layer* l = c.col.g.toplayer; l; l = l->nextl) { CHECK(l->ice > 0.0); CHECK(l->frozen == 1); }
}

TEST_CASE("enthalpyThermalStep: melt source into the snow cell, fully melted snow heat goes down") {
  Col c(0.5, 0.2);
  SnowLayer* s = new SnowLayer(); s->dz = 0.05; s->ice = 5.0; s->liq = 0.0; s->tem = -1.0; s->rho = 100.0;
  s->z = 0.0; c.col.g.insertFront(s); c.col.reindex(); c.col.g.fstsoill = c.col.L.front();
  const double H0 = colH(c.col.g), W0 = colW(c.col.g);
  EnergyDiag q;
  const double qm = 150.0;   // W m-2, more than enough to melt 5 kg m-2 in a day
  c.se.enthalpyThermalStep(5.0, qm, q);
  CHECK(q.src == doctest::Approx(qm * 86400.0));
  CHECK(fabs(q.res) / 86400.0 < 1e-6);
  CHECK(s->ice == 0.0); CHECK(s->tem == 0.0); CHECK(s->frozen == -1);
  CHECK(colW(c.col.g) == doctest::Approx(W0));
  CHECK(colH(c.col.g) - H0 == doctest::Approx(q.gtop + q.gbot + q.src + q.res + q.clip).epsilon(1e-12).scale(1.0));
}

TEST_CASE("Ground::enthalpy_remap: split and merge conserve enthalpy and water (legacy unchanged)") {
  for (int remap : {0, 1}) {
    CAPTURE(remap);
    Ground g; g.enthalpy_remap = remap;
    SoilLayer* u = make_soil_layer(0.0, 0.3, 0.4, 30.0, 40.0, -0.7, 0, 0.6);
    SoilLayer* l = make_soil_layer(0.3, 0.1, 0.4, 0.0, 0.0, 0.0, 1, 1.0);
    const double H0 = Soil_Env::enthLayerH(u), W0 = u->liq + u->ice;
    g.splitOneSoilLayer(u, l, 0.0, 0.1);
    const double H1 = Soil_Env::enthLayerH(u) + Soil_Env::enthLayerH(l);
    if (remap) {
      CHECK(H1 == doctest::Approx(H0).epsilon(1e-12));
      CHECK(u->tem == -0.7); CHECK(l->tem == -0.7);
      CHECK(l->ice / l->dz == doctest::Approx(u->ice / u->dz)); CHECK(l->liq / l->dz == doctest::Approx(u->liq / u->dz));
    }
    CHECK(u->liq + u->ice + l->liq + l->ice == doctest::Approx(W0));
    // merge back both ways
    SoilLayer* a = make_soil_layer(0.0, 0.2, 0.4, 20.0, 10.0, -3.0, 1, 1.0);
    SoilLayer* b = make_soil_layer(0.2, 0.3, 0.4, 50.0, 0.0, 4.0, -1, 0.0);
    SoilLayer* a2 = make_soil_layer(0.0, 0.2, 0.4, 20.0, 10.0, -3.0, 1, 1.0);
    SoilLayer* b2 = make_soil_layer(0.2, 0.3, 0.4, 50.0, 0.0, 4.0, -1, 0.0);
    const double Hm = Soil_Env::enthLayerH(a) + Soil_Env::enthLayerH(b);
    g.combineTwoSoilLayersU2L(a, b);
    g.combineTwoSoilLayersL2U(b2, a2);
    if (remap) {
      CHECK(Soil_Env::enthLayerH(b) == doctest::Approx(Hm).epsilon(1e-12));
      CHECK(Soil_Env::enthLayerH(a2) == doctest::Approx(Hm).epsilon(1e-12));
    } else {   // legacy: thickness weighted
      CHECK(b->tem == doctest::Approx(-3.0 * 0.4 + 4.0 * 0.6));
    }
    delete u; delete l; delete a; delete b; delete a2; delete b2;
  }
}

}
