// Phase 3 (regional_wetland_deepc merge): wetland switch work.
//  - model_settings.wetland block: defaults (off), parsing, validation (FATAL), deprecated flat keys
//  - wetland CMT gate (ModelData::is_wetland_cmt; Cohort sets cd.is_wetland/ch4_on/runon/groundwater from it)
//  - run-on water balance through Soil_Env::updateDailySM: P + run-on = d(soil) + d(ponding) + qover + qdrain
//    (upstream deleted the excess water whenever run-on was on), and the surface partition
//    P + run-on = infil + qover + d(ponding) before the soil step
//  - eq_co2_mode, interception_scheme defaults; parameter block lengths (18/35/41, 33/34)
//  - CH4 bookkeeping: daily accumulators start at 0; layer CH4 concentration kept on split/combine
// Phase 3 v2: runon_energy_mode, ch4_conserve, ch4_daily (all default 0)
#include "doctest.h"
#include "test_helpers.h"
#include <json/value.h>
#include <cmath>
#include <fstream>
#include "../include/physicalconst.h"

using ut::make_soil_layer;

namespace {
Json::Value ms(const Json::Value& w) { Json::Value c; c["model_settings"]["wetland"] = w; return c; }
Json::Value wl(bool en) { Json::Value w; w["enabled"] = en; return w; }

std::vector<SoilLayer*> col(int n, double dz, double poro, double sat, double tem = 5.0) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) {
    SoilLayer* l = make_soil_layer(i * dz, dz, poro, sat * poro * DENLIQ * dz, 0.0, tem, -1, 0.0);
    l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0;
    v.push_back(l);
  }
  return v;
}

// Soil_Env::updateDailySM on a bare (no vegetation, no snow) thawed column
struct SM {
  ut::Column c; EnvData ed; CohortData cd; Soil_Env se;
  explicit SM(double sat, bool frozen_top = false) : c(col(10, 0.1, 0.4, sat)) {
    Ground& g = c.g;
    g.fstshlwl = c.L[0]; g.fstmossl = c.L[0]; g.lstsoill = c.L.back(); g.drainl = c.L.back();
    g.draindepth = c.L.back()->z + c.L.back()->dz; g.ststate = -1; g.fstfntl = NULL; g.lstfntl = NULL;
    for (int i = 0; i < MAX_NUM_FNT; ++i) g.frnttype[i] = 0;
    g.expelled_water_mm = 0.0;
    if (frozen_top) {
      Layer* l = c.L[0]; l->frozen = 1; l->frozenfrac = 1.0; l->tem = -1.0;
      l->ice = l->liq; l->liq = 0.0; g.frnttype[0] = 1;
    }
    se.setGround(&g); se.setEnvData(&ed); se.setCohortData(&cd);
    se.richards.dt_fix = 1;   // richards_dt_fix 1: legacy Richards substep rejection creates water (bug R1)
    cd.m_vegd.fpcsum = 0.0; cd.drainage_type = 0; cd.cell_slope = 0.0;
    ed.d_v2g.rthfl = 0.0; ed.d_v2g.rdrip = 0.0; ed.d_snw2soi.melt = 0.0; ed.d_soi2a.evap = 0.0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) ed.d_soid.r_e_i[i] = 0.0;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
  }
  void runon(bool on) { cd.runon_status = on; se.richards.runon_status = on; }
  // returns closure error (mm) of: soil + puddle + qover + qdrain = initial + rain + runon
  double day(double rain) {
    const double w0 = c.total_water() + ed.d_soi2l.magic_puddle;
    ed.d_a2l.rnfl = rain;
    se.updateDailySM(0.0);
    return c.total_water() + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + ed.d_soi2l.qdrain
           - w0 - rain - se.runon_in_mm;
  }
};
} // namespace

TEST_SUITE("wetland_p3") {

TEST_CASE("wetland block defaults: off in ModelData() and ModelData(Json) without keys") {
  ModelData m0;
  CHECK_FALSE(m0.wetland_enabled); CHECK(m0.wetland_ch4); CHECK_FALSE(m0.wetland_runon); CHECK_FALSE(m0.wetland_groundwater);
  CHECK(m0.wetland_cmts == std::vector<int>({6, 31, 32, 55, 61, 75, 92, 93}));
  CHECK(m0.interception_scheme == 0); CHECK(m0.eq_co2_mode == 0); CHECK(m0.wetland_diag == 0);
  CHECK(m0.runon_energy_mode == 0); CHECK(m0.ch4_conserve == 0); CHECK(m0.ch4_daily == 0);
  CHECK_FALSE(m0.get_ch4_module()); CHECK_FALSE(m0.get_runon_status()); CHECK_FALSE(m0.get_groundwater_status());
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c);
  CHECK_FALSE(mj.wetland_enabled); CHECK(mj.wetland_cmts.size() == 8); CHECK(mj.interception_scheme == 0);
  CHECK(mj.eq_co2_mode == 0); CHECK_FALSE(mj.get_ch4_module());
  for (int cmt : {0, 5, 6, 31, 55, 93}) { CAPTURE(cmt); CHECK_FALSE(mj.is_wetland_cmt(cmt)); }
}

TEST_CASE("wetland block: parsed values and CMT gate") {
  Json::Value w = wl(true); w["ch4"] = false; w["run_on"] = true; w["groundwater"] = 1;
  w["interception_scheme"] = 1; w["diag"] = 1; w["cmts"].append(55); w["cmts"].append(61);
  ModelData m(ms(w));
  CHECK(m.wetland_enabled); CHECK_FALSE(m.wetland_ch4); CHECK(m.wetland_runon); CHECK(m.wetland_groundwater);
  CHECK(m.interception_scheme == 1); CHECK(m.wetland_diag == 1);
  CHECK(m.is_wetland_cmt(55)); CHECK(m.is_wetland_cmt(61));
  CHECK_FALSE(m.is_wetland_cmt(31)); CHECK_FALSE(m.is_wetland_cmt(5));
  CHECK_FALSE(m.get_ch4_module()); CHECK(m.get_runon_status()); CHECK(m.get_groundwater_status());
  // enabled false: nothing is a wetland even if listed, all sub-flags off
  Json::Value w0 = wl(false); w0["run_on"] = true; w0["groundwater"] = true;
  ModelData m0(ms(w0));
  CHECK_FALSE(m0.is_wetland_cmt(55)); CHECK_FALSE(m0.get_runon_status()); CHECK_FALSE(m0.get_groundwater_status());
  // default list with enabled
  Json::Value w1 = wl(true); w1["ch4"] = false;
  ModelData m1(ms(w1));
  for (int cmt : {6, 31, 32, 55, 61, 75, 92, 93}) { CAPTURE(cmt); CHECK(m1.is_wetland_cmt(cmt)); }
  for (int cmt : {0, 1, 4, 5, 13, 21, 51}) { CAPTURE(cmt); CHECK_FALSE(m1.is_wetland_cmt(cmt)); }
}

TEST_CASE("wetland block: invalid values are FATAL") {
  auto st = [](Json::Value c) { return ut::exit_status_of([c] { ModelData m(c); }); };
  { Json::Value w = wl(true); w["ch4"] = false; CHECK(st(ms(w)) == 0); }
  { Json::Value w; w["enabled"] = "yes"; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w; w["enabled"] = 2; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["unknown_key"] = 1; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["cmts"] = 31; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["cmts"].append(150); CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["interception_scheme"] = 2; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["diag"] = -1; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["runon_energy_mode"] = 3; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(true); w["ch4"] = false; w["runon_energy_mode"] = 2; CHECK(ModelData(ms(w)).runon_energy_mode == 2); }
  { Json::Value w = wl(false); w["ch4_conserve"] = -1; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["ch4_daily"] = 2; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(true); w["ch4"] = false; w["runon_energy_mode"] = 1; w["ch4_conserve"] = 1; w["ch4_daily"] = 1; CHECK(st(ms(w)) == 0); }
  { Json::Value w = wl(true); w["ch4"] = false; w["runon_energy_mode"] = 1; CHECK(ModelData(ms(w)).runon_energy_mode == 1); }
  { Json::Value c; c["model_settings"]["wetland"] = 1; CHECK(st(c) == EXIT_FAILURE); }
  { Json::Value c; c["model_settings"]["eq_co2_mode"] = 2; CHECK(st(c) == EXIT_FAILURE); }
  { Json::Value c; c["model_settings"]["eq_co2_mode"] = 1; CHECK(st(c) == 0); CHECK(ModelData(c).eq_co2_mode == 1); }
  // enabled with ch4 (default true) needs IO.ch4_file
  { Json::Value w = wl(true); CHECK(st(ms(w)) == EXIT_FAILURE);
    Json::Value c = ms(w); c["IO"]["ch4_file"] = "ch4.nc"; CHECK(st(c) == 0); CHECK(ModelData(c).get_ch4_module()); }
}

TEST_CASE("deprecated flat keys ch4/run-on/groundwater: aliases only without a wetland block") {
  Json::Value c; c["model_settings"]["ch4"] = false; c["model_settings"]["run-on"] = true; c["model_settings"]["groundwater"] = false;
  ModelData m(c);
  CHECK(m.wetland_enabled); CHECK(m.wetland_runon); CHECK_FALSE(m.wetland_ch4); CHECK(m.is_wetland_cmt(31));
  Json::Value c0; c0["model_settings"]["ch4"] = false; c0["model_settings"]["run-on"] = false; c0["model_settings"]["groundwater"] = false;
  CHECK_FALSE(ModelData(c0).wetland_enabled);
  Json::Value cb = c; cb["model_settings"]["wetland"] = wl(false);     // block wins
  CHECK_FALSE(ModelData(cb).wetland_enabled);
}

TEST_CASE("run-on water balance through updateDailySM: P + run-on = d(soil) + d(ponding) + qover + qdrain") {
  for (int on : {0, 1}) for (double sat : {0.3, 0.95}) for (double rain : {5.0, 60.0}) {
    CAPTURE(on); CAPTURE(sat); CAPTURE(rain);
    SM s(sat); s.runon(on == 1);
    double worst = 0.0;
    for (int d = 0; d < 3; ++d) { worst = fmax(worst, fabs(s.day(rain))); }
    CHECK(worst < 1e-8);
    if (on == 0) { CHECK(s.se.runon_in_mm == 0.0); }
    if (on == 1) { CHECK(s.se.runon_in_mm >= 0.0); CHECK(s.ed.d_soi2l.qover >= 0.0); }
  }
}

TEST_CASE("run-on water balance, frozen top (no infiltration): excess goes to ponding then qover (was deleted)") {
  for (int on : {0, 1}) {
    CAPTURE(on);
    SM s(0.6, true); s.runon(on == 1);
    const double err = s.day(30.0);
    CHECK(fabs(err) < 1e-9);
    CHECK(s.ed.d_soi2l.magic_puddle == doctest::Approx(s.se.ponding_max_mm));
    CHECK(s.ed.d_soi2l.qover == doctest::Approx(30.0 + s.se.runon_in_mm - s.se.ponding_max_mm));
    CHECK(s.ed.d_soi2l.qinfl == 0.0);
  }
}

TEST_CASE("surface partition: P + run-on = infil + qover + d(ponding) (thawed, nearly saturated)") {
  for (int on : {0, 1}) {
    CAPTURE(on);
    SM s(0.95); s.runon(on == 1);
    s.ed.d_soi2l.magic_puddle = 1.0;
    const double p0 = s.ed.d_soi2l.magic_puddle;
    // the validity check may move excess water to puddle/qover after Richards; isolate the surface
    // partition by reading qinfl (= infil passed to Richards) and the pre-Richards puddle change:
    // closure of the full day (checked above) plus qinfl >= 0 and qover >= 0
    const double err = s.day(20.0);
    CHECK(fabs(err) < 1e-8);
    CHECK(s.ed.d_soi2l.qinfl >= 0.0);
    CHECK(s.ed.d_soi2l.qinfl <= 20.0 + s.se.runon_in_mm + p0 + 1e-12);
  }
}

TEST_CASE("checkSoilLiquidWaterValidity: oversaturated water goes to ponding then qover (no cd needed)") {
  ut::Column c(col(5, 0.1, 0.4, 1.3));
  EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&c.g);
  ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0;
  const double w0 = c.total_water();
  se.checkSoilLiquidWaterValidity(c.L[0], c.L[0]->solind);
  CHECK(c.total_water() < w0);
  CHECK(c.total_water() + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover == doctest::Approx(w0).epsilon(1e-12));
  CHECK(ed.d_soi2l.qover > 0.0);
}

TEST_CASE("defaults: interception_scheme 0 in Vegetation_Env, eq_co2_mode 0 in Climate, RestartData::write_ch4 false") {
  Vegetation_Env ve; CHECK(ve.interception_scheme == 0);
  Climate cl; CHECK(cl.eq_co2_mode == 0); CHECK(cl.ch4_idx_offset == 0);
  CHECK_FALSE(RestartData::write_ch4);
  Soil_Env se; CHECK(se.richards.max_ponding == 100.0); CHECK(se.ponding_max_mm == 4.0);
  CHECK_FALSE(se.richards.groundwater_status); CHECK_FALSE(se.richards.runon_status);
  CohortData cd; CHECK_FALSE(cd.is_wetland); CHECK_FALSE(cd.ch4_on); CHECK_FALSE(cd.runon_status); CHECK_FALSE(cd.groundwater_status);
}

TEST_CASE("interception_scheme 0 = pre-merge formula; 1 = upstream formula") {
  Vegetation_Env ve; ve.envpar.ircoef = 0.041; ve.envpar.iscoef = 0.04;
  for (double rain : {0.5, 2.0, 10.0}) for (double lai : {0.5, 2.0}) {
    CAPTURE(rain); CAPTURE(lai);
    ve.interception_scheme = 0;
    CHECK(ve.getRainInterception(rain, lai) == ve.getRainInterceptionLegacy(rain, lai));
    CHECK(ve.getSnowInterception(rain, lai) == ve.getSnowInterceptionLegacy(rain, lai));
    ve.interception_scheme = 1;
    CHECK(ve.getRainInterception(rain, lai) == doctest::Approx(fmin(rain, rain * 0.041 * (1 - tanh(-0.5 * lai)))));
    CHECK(ve.getSnowInterception(rain, lai) == doctest::Approx(fmin(rain, rain * 0.04 * (1 - exp(-0.5 * lai)))));
  }
}

TEST_CASE("parameter blocks: calparbgc 18/35/41 and bgcvegetation 33/34 lines accepted; others FATAL") {
  const std::string f = "/tmp/ut_wetland_p3_calpar.txt";
  auto block = [](int cmt, int n) {
    std::stringstream s; s << "// CMT" << (cmt < 10 ? "0" : "") << cmt << " // test\n";
    for (int i = 0; i < n; ++i) s << (i + 1) << ".0 // line " << i << "\n";
    return s.str();
  };
  { std::ofstream o(f.c_str()); o << "//===\n" << block(1, 18) << block(2, 35) << block(3, 41) << block(4, 40); }
  CHECK(temutil::parse_parameter_file_any(f, 1, {18, 35, 41}).size() == 18);
  CHECK(temutil::parse_parameter_file_any(f, 2, {18, 35, 41}).size() == 35);
  CHECK(temutil::parse_parameter_file_any(f, 3, {18, 35, 41}).size() == 41);
  CHECK(ut::exit_status_of([f] { temutil::parse_parameter_file_any(f, 4, {18, 35, 41}); }) != 0);
  CHECK(temutil::parse_parameter_file_any(f, 1, {18}).size() == 18);
}

TEST_CASE("CH4 bookkeeping: daily accumulators start at 0; layer concentration kept on split and combine") {
  BgcData bd;
  for (int il = 0; il < MAX_SOI_LAY; ++il) {
    CHECK(bd.d_soi2soi.ch4_rawc[il] == 0.0); CHECK(bd.d_soi2soi.ch4_somcr[il] == 0.0);
    CHECK(bd.d_soi2a.ch4_oxid[il] == 0.0);
  }
  CHECK(bd.d_soi2a.ch4efflux == 0.0);
  for (int er : {0, 1}) {
    CAPTURE(er);
    ut::Column c(col(3, 0.2, 0.4, 0.5));
    c.g.enthalpy_remap = er;
    for (Layer* l : c.L) l->ch4 = 2.0;
    c.L[0]->ch4 = 4.0;
    SoilLayer* nl = make_soil_layer(0.0, 0.05, 0.4, 0.0, 0.0, 5.0, -1, 0.0);
    c.g.splitOneSoilLayer(dynamic_cast<SoilLayer*>(c.L[1]), nl, c.L[1]->z, 0.05);
    CHECK(nl->ch4 == 2.0); CHECK(c.L[1]->ch4 == 2.0);
    delete nl;
    // combine 0.2 m at 4.0 (top) with 0.2 m at 2.0 -> 3.0 (thickness-weighted mean)
    ut::Column d(col(2, 0.2, 0.4, 0.5));
    d.g.enthalpy_remap = er;
    d.L[0]->ch4 = 4.0; d.L[1]->ch4 = 2.0;
    d.g.combineTwoSoilLayersL2U(dynamic_cast<SoilLayer*>(d.L[1]), dynamic_cast<SoilLayer*>(d.L[0]));
    CHECK(d.L[0]->ch4 == doctest::Approx(3.0));
  }
}


TEST_CASE("runon_energy_mode 1: a cell does not recycle its own runoff; water still closes") {
  for (int frozen : {0, 1}) {
    CAPTURE(frozen);
    SM s(0.95, frozen == 1); s.runon(true); s.se.runon_energy_mode = 1;
    const double err = s.day(20.0);
    CHECK(s.se.runon_in_mm == 0.0);
    CHECK(s.ed.d_soi2l.qover >= 0.0);
    CHECK(fabs(err) < 1e-8);
  }
}

TEST_CASE("runon_energy_mode 2: puddle energy on, own runoff still recycled; water closes") {
  SM s(0.95); s.runon(true); s.se.runon_energy_mode = 2;
  CHECK(s.se.puddle_energy_on());
  const double err = s.day(20.0);
  CHECK(s.se.runon_in_mm >= 0.0);
  CHECK(fabs(err) < 1e-8);
}

TEST_CASE("runon_energy_mode 0 still recycles (legacy path, existing water-balance tests)") {
  SM s(0.95); s.runon(true); s.se.runon_energy_mode = 0;
  s.day(20.0);
  CHECK(s.se.runon_in_mm >= 0.0);
}

TEST_CASE("runon_energy_mode 1: puddle freezes at Tair < 0 and thaws at Tair > 0; Lf goes to the atmosphere") {
  SM s(0.5); s.se.runon_energy_mode = 1; s.runon(true);
  s.ed.d_soi2l.magic_puddle = 10.0; s.ed.d_atms.ta = -10.0; s.se.puddle_ice_mm = 0.0; s.se.puddle_heat_Jm2 = 0.0; s.se.puddle_atm_Jm2 = 0.0;
  s.se.updatePuddlePhase();
  CHECK(s.ed.d_soi2l.magic_puddle == doctest::Approx(0.0));
  CHECK(s.se.puddle_ice_mm == doctest::Approx(10.0));
  CHECK(s.se.puddle_heat_Jm2 == doctest::Approx(0.0));
  CHECK(s.se.puddle_atm_Jm2 == doctest::Approx(10.0 * LHFUS));
  s.ed.d_atms.ta = 5.0;
  s.se.updatePuddlePhase();
  CHECK(s.ed.d_soi2l.magic_puddle == doctest::Approx(10.0));
  CHECK(s.se.puddle_ice_mm == doctest::Approx(0.0));
  CHECK(s.se.puddle_heat_Jm2 == doctest::Approx(0.0));
}

TEST_CASE("runon_energy_mode 1: puddle energy is applied only when runon_status is true") {
  SM s(0.5); s.se.runon_energy_mode = 1; s.runon(false);
  s.ed.d_soi2l.magic_puddle = 10.0; s.ed.d_atms.ta = -10.0;
  s.se.puddle_ice_mm = 0.0; s.se.puddle_heat_Jm2 = 0.0;
  s.se.updatePuddlePhase();
  CHECK(s.ed.d_soi2l.magic_puddle == doctest::Approx(10.0));
  CHECK(s.se.puddle_ice_mm == doctest::Approx(0.0));
  CHECK(s.se.puddle_heat_Jm2 == doctest::Approx(0.0));
}

TEST_CASE("runon_energy_mode 1: ponding capacity is ponding_max minus liquid plus ice") {
  SM s(0.5); s.se.runon_energy_mode = 1; s.runon(true);
  s.se.ponding_max_mm = 50.0;
  s.ed.d_soi2l.magic_puddle = 10.0;
  s.se.puddle_ice_mm = 30.0;
  CHECK(s.se.puddle_total_mm() == doctest::Approx(40.0));
  CHECK(s.se.space_in_puddle_mm() == doctest::Approx(10.0));
}

TEST_CASE("runon_energy_mode 1: expelled ice pays Lf in the source layer; mass becomes liquid") {
  SM s(0.5); s.se.runon_energy_mode = 1; s.runon(true);
  s.se.ponding_max_mm = 20.0;
  s.se.puddle_ice_mm = 0.0;
  s.ed.d_soi2l.magic_puddle = 0.0;
  s.se.thermal_solver = 0;
  Layer* l = s.c.L[0];
  REQUIRE(l != nullptr);
  l->ice = 80.0;
  l->liq = 0.0;
  l->tem = -2.0;
  const double t0 = l->tem;
  s.se.expelLayerIce(l, 80.0);
  CHECK(s.ed.d_soi2l.magic_puddle == doctest::Approx(20.0));
  CHECK(s.se.puddle_ice_mm == doctest::Approx(0.0));
  CHECK(l->ice == doctest::Approx(0.0));
  CHECK(s.c.g.expelled_water_mm == doctest::Approx(60.0));
  CHECK(l->tem < t0);   // source layer paid Lf for the whole 80 mm
}

TEST_CASE("runon_energy_mode 1: puddle heat is spread over the top 0.3 m, not dumped into a thin layer") {
  SM s(0.5); s.se.runon_energy_mode = 1; s.runon(true);
  s.se.thermal_solver = 0;
  Layer* top = s.c.L[0];
  top->dz = 0.002;
  top->tem = -1.0;
  for (size_t i = 1; i < s.c.L.size(); ++i) { s.c.L[i]->tem = -1.0; }
  s.se.puddle_heat_Jm2 = 51.5 * LHFUS;
  s.se.applyPuddleHeatToSoil();
  // A 2 mm skin would jump thousands of C if it took all 51.5 mm of Lf.
  CHECK(fabs(top->tem + 1.0) < 40.0);
  CHECK(s.c.L[1]->tem > -1.0 + 0.5);  // heat reached the layers below
}

TEST_CASE("ch4_conserve and ch4_daily default off; Soil_Bgc and BgcData flags start at 0") {
  Soil_Bgc sb; CHECK(sb.ch4_conserve == 0); CHECK(sb.ch4_daily == 0);
  BgcData bd; CHECK(bd.reco_ch4ox_output == 0);
}

} // TEST_SUITE

TEST_SUITE("wetland_p3_fvcn") {
namespace {
// Solve one FV-CN step on nodes 0..n+1 (0 = atmospheric node fixed at ub, n+1 = decoupled ghost).
std::vector<double> fv_step(const std::vector<double>& D, const std::vector<double>& dz,
                            const std::vector<double>& c, const std::vector<double>& R, double ub, double dt) {
  const int n = (int)D.size() - 2, N = n + 2;
  std::vector<double> lo(N, 0.0), md(N, 1.0), up(N, 0.0), rhs(N, 0.0), x(N, 0.0);
  rhs[0] = ub;
  Soil_Bgc::ch4FvCnRows(n, D.data(), dz.data(), c.data(), ub, R.data(), dt, lo.data(), md.data(), up.data(), rhs.data());
  Soil_Bgc sb; sb.TriSolver(N, lo.data(), md.data(), up.data(), rhs.data(), x.data());
  return x;
}
double mass(const std::vector<double>& dz, const std::vector<double>& c) {
  double m = 0.0; for (size_t i = 1; i + 1 < c.size(); ++i) { m += dz[i] * c[i]; } return m;
}
// Non-uniform grid, D varying by 3 orders of magnitude (moss/peat/saturated mineral).
const std::vector<double> DZ = {0, 0.02, 0.05, 0.1, 0.1, 0.2, 0.3, 0.5, 0};
const std::vector<double> DD = {0, 0.07, 0.03, 0.002, 1e-4, 8e-5, 5e-5, 2e-5, 0};
}

TEST_CASE("ch4_conserve FV-CN: closed top conserves column CH4 exactly (non-uniform grid, varying D)") {
  std::vector<double> D = DD; D[1] = 0.0;  // no exchange with the atmosphere
  std::vector<double> c = {0, 0.1, 5, 40, 300, 900, 1200, 1500, 0}, R(9, 0.0);
  const double m0 = mass(DZ, c);
  for (int s = 0; s < 240; ++s) { std::vector<double> x = fv_step(D, DZ, c, R, 0.0, 1.0); for (int i = 1; i <= 7; ++i) c[i] = x[i]; }
  CHECK(mass(DZ, c) == doctest::Approx(m0).epsilon(1e-12));
  for (int i = 1; i <= 7; ++i) CHECK(c[i] >= 0.0);
}

TEST_CASE("ch4_conserve FV-CN: mass change = dt*sum(dz R) - top CN flux; uniform field stays uniform") {
  std::vector<double> c = {0, 0.5, 5, 40, 300, 900, 1200, 1500, 0};
  std::vector<double> R = {0, -0.01, 0.0, 0.3, 2.0, 1.0, 0.5, 0.1, 0};
  const double ub = 0.08, dt = 1.0 / 12.0;
  std::vector<double> x = fv_step(DD, DZ, c, R, ub, dt);
  double sR = 0.0; for (int i = 1; i <= 7; ++i) sR += DZ[i] * R[i];
  const double gtop = 2.0 * DD[1] / DZ[1];
  const double ftop = dt * gtop * (0.5 * (x[1] + c[1]) - ub);  // CN top flux (positive out)
  CHECK(mass(DZ, x) - mass(DZ, c) == doctest::Approx(dt * sR - ftop).epsilon(1e-10));
  // well balanced: a uniform field at ub with no reactions does not move
  std::vector<double> u(9, ub), Z(9, 0.0); u[0] = u[8] = 0.0;
  std::vector<double> y = fv_step(DD, DZ, u, Z, ub, dt);
  for (int i = 1; i <= 7; ++i) CHECK(y[i] == doctest::Approx(ub).epsilon(1e-13));
}

TEST_CASE("ch4_conserve FV-CN: explicit half stays non-negative when D dt/(2 dz^2) <= 0.25") {
  std::vector<double> D = DD, c(9, 0.0), R(9, 0.0);
  c[4] = 100.0;  // single spike
  double a1max = 0.0; for (int i = 1; i <= 7; ++i) a1max = std::max(a1max, D[i] / (2.0 * DZ[i] * DZ[i]));
  int nsub = std::max(1, (int)(a1max / 0.25 + 0.999999)); double dt = 1.0 / nsub;
  std::vector<double> lo(9, 0.0), md(9, 1.0), up(9, 0.0), rhs(9, 0.0);
  Soil_Bgc::ch4FvCnRows(7, D.data(), DZ.data(), c.data(), 0.0, R.data(), dt, lo.data(), md.data(), up.data(), rhs.data());
  for (int i = 1; i <= 7; ++i) { CHECK(rhs[i] >= 0.0); CHECK(-lo[i] - up[i] <= 1.0 + 1e-12); }
}
} // TEST_SUITE

TEST_SUITE("wetland_p3_reco_ch4ox") {
namespace {
// One month of BgcData bookkeeping with given per-layer RH, CH4 oxidation and efflux.
struct RecoMonth { double rhsom, y_rhsom, nep, y_nep; };
RecoMonth reco_month(int mode, bool ch4_on, int months = 1) {
  CohortData cd; cd.ch4_on = ch4_on;
  for (int il = 0; il < MAX_SOI_LAY; il++) { cd.m_soil.type[il] = (il < 3) ? 1 : 3; }
  BgcData bd; bd.cd = &cd; bd.reco_ch4ox_output = mode;
  for (int il = 0; il < MAX_SOI_LAY; il++) {
    bd.m_sois.rawc[il] = bd.m_sois.soma[il] = bd.m_sois.sompr[il] = bd.m_sois.somcr[il] = 1.0;
    bd.m_sois.orgn[il] = bd.m_sois.avln[il] = 0.1;
  }
  bd.soil_beginOfYear(); bd.land_beginOfYear();
  for (int m = 0; m < months; m++) {
    bd.soil_beginOfMonth();
    for (int il = 0; il < MAX_SOI_LAY; il++) {
      bd.m_soi2a.rhrawc[il] = 0.10 * (il < 5); bd.m_soi2a.rhsoma[il] = 0.20 * (il < 5);
      bd.m_soi2a.rhsompr[il] = 0.05 * (il < 5); bd.m_soi2a.rhsomcr[il] = 0.01 * (il < 5);
      bd.m_soi2soi.netnmin[il] = 0.0; bd.m_soi2soi.nimmob[il] = 0.0;
      bd.m_soi2a.ch4_oxid[il] = (il < 4) ? 0.75 : 0.0;   // 3.0 g C m-2 month-1
    }
    bd.m_soi2a.ch4efflux = 1.25; bd.m_soi2a.rhwdeb = 0.0; bd.m_a2v.nppall = 20.0;
    bd.soil_endOfMonth(m);
    bd.land_endOfMonth();
  }
  return {bd.m_soi2a.rhsom, bd.y_soi2a.rhsom, bd.m_l2a.nep, bd.y_l2a.nep};
}
}

TEST_CASE("reco_ch4ox_output: default 0, parsed 0..2, other values FATAL") {
  ModelData m0; CHECK(m0.reco_ch4ox_output == 0);
  CHECK(ModelData(ms(wl(false))).reco_ch4ox_output == 0);
  auto st = [](Json::Value c) { return ut::exit_status_of([c] { ModelData m(c); }); };
  for (int v : {1, 2}) {
    Json::Value w = wl(true); w["ch4"] = false; w["reco_ch4ox_output"] = v;
    CAPTURE(v); CHECK(ModelData(ms(w)).reco_ch4ox_output == v);
  }
  { Json::Value w = wl(false); w["reco_ch4ox_output"] = 3; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["reco_ch4ox_output"] = -1; CHECK(st(ms(w)) == EXIT_FAILURE); }
  { Json::Value w = wl(false); w["reco_ch4ox_output"] = "1"; CHECK(st(ms(w)) == EXIT_FAILURE); }
}

TEST_CASE("reco_ch4ox_add: 0 adds nothing, 1 adds oxidation, 2 adds oxidation + efflux") {
  CHECK(BgcData::reco_ch4ox_add(0, 3.0, 1.25) == 0.0);
  CHECK(BgcData::reco_ch4ox_add(1, 3.0, 1.25) == 3.0);
  CHECK(BgcData::reco_ch4ox_add(2, 3.0, 1.25) == 4.25);
}

TEST_CASE("reco_ch4ox_output: rhsom (RH/RECO/NEE) and nep, monthly and yearly; ch4_conserve no longer adds") {
  const RecoMonth r0 = reco_month(0, true), r1 = reco_month(1, true), r2 = reco_month(2, true);
  CHECK(r0.rhsom == doctest::Approx(5 * 0.36));
  CHECK(r1.rhsom - r0.rhsom == doctest::Approx(3.0));
  CHECK(r2.rhsom - r0.rhsom == doctest::Approx(4.25));
  CHECK(r1.nep - r0.nep == doctest::Approx(-3.0));        // NEP loses the oxidation CO2
  CHECK(r2.nep - r0.nep == doctest::Approx(-4.25));
  const RecoMonth y0 = reco_month(0, true, 12), y1 = reco_month(1, true, 12), y2 = reco_month(2, true, 12);
  CHECK(y1.y_rhsom - y0.y_rhsom == doctest::Approx(12 * 3.0));
  CHECK(y2.y_rhsom - y0.y_rhsom == doctest::Approx(12 * 4.25));
  CHECK(y1.y_nep - y0.y_nep == doctest::Approx(-12 * 3.0));
  // CH4 off: every mode is bit-identical to mode 0
  const RecoMonth n0 = reco_month(0, false, 12), n2 = reco_month(2, false, 12);
  CHECK(n0.y_rhsom == y0.y_rhsom); CHECK(n2.y_rhsom == n0.y_rhsom); CHECK(n2.y_nep == n0.y_nep);
}

} // TEST_SUITE
