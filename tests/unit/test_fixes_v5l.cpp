// v5l: switchable ET-split fixes S1-S4 (all default OFF = legacy).
//   S1 model_settings.et_fpc_once_fix  -> Vegetation_Env::fpc_once_fix
//   S2 model_settings.nv_snow_bury_fix (+ nv_snow_bury_depth, m)
//                                      -> Vegetation_Env::nv_snow_bury_fix / nv_snow_bury_depth
//   S3 model_settings.rain_int_lai_fix -> Vegetation_Env::rain_int_lai_fix
//   S4 model_settings.nv_rain_surface  -> Vegetation_Env::nv_rain_surface (0/1/2)
// Expected values are hand calculations from the formulas in Vegetation_Env.cpp
// (Penman-Monteith via the model's own getPenMonET, Beer's-law radiation).
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>

namespace {
struct VegDay {
  EnvData ed; Vegetation_Env ve;
  VegDay() {
    ve.setEnvData(&ed);
    ve.ipft = 0;
    vegpar_env& p = ve.envpar;
    p.albvisnir = 0.1; p.er = 0.5; p.ircoef = 0.041; p.iscoef = 0.041;
    p.gl_bl = 0.08; p.gl_c = 0.0; p.vpd_open = 930.; p.vpd_close = 4100.;
    p.glmax = 0.003; p.ppfd50 = 75.;
    ed.d_a2l.nirr = 200.; ed.d_a2l.par = 90.;
    ed.d_a2l.rnfl = 0.; ed.d_a2l.snfl = 0.;
    ed.d_atms.ta = 12.; ed.d_atmd.vpd = 800.;   // Pa
    ed.d_vegd.btran = 0.8;
    ed.d_vegs.rwater = 0.; ed.d_vegs.snow = 0.;
  }
  void day(double lai, double fpc, double daylhr = 18.) {
    ve.updateRadiation(lai, fpc);
    ve.updateWaterBalance(daylhr, lai, fpc);
  }
  // per-day water balance of the canopy rain store (any switch setting)
  double rain_residual(double rnfl_fpc, double s0) const {
    return (s0 + ed.d_a2v.rinter) - (ed.d_v2a.evap + ed.d_v2a.evap_moss + ed.d_v2g.rdrip + ed.d_vegs.rwater)
           + (rnfl_fpc - ed.d_a2v.rinter - ed.d_v2g.rthfl);
  }
};
Json::Value ms(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
}

TEST_SUITE("fixes_v5l") {

TEST_CASE("v5l switches default to 0 (legacy) in Vegetation_Env and ModelData()") {
  Vegetation_Env ve;
  CHECK(ve.fpc_once_fix == 0); CHECK(ve.nv_snow_bury_fix == 0); CHECK(ve.nv_snow_bury_depth == 0.0);
  CHECK(ve.rain_int_lai_fix == 0); CHECK(ve.nv_rain_surface == 0);
  ModelData md;
  CHECK(md.et_fpc_once_fix == 0); CHECK(md.nv_snow_bury_fix == 0); CHECK(md.rain_int_lai_fix == 0);
  CHECK(md.nv_rain_surface == 0); CHECK(md.nv_snow_bury_depth == 0.0);
}

TEST_CASE("v5l switch validation: 0/1 (nv_rain_surface 0/1/2), depth >= 0") {
  CHECK(ModelData::v5l_switches_valid(0, 0, 0.0, 0, 0));
  CHECK(ModelData::v5l_switches_valid(1, 1, 0.05, 1, 2));
  CHECK_FALSE(ModelData::v5l_switches_valid(2, 0, 0.0, 0, 0));
  CHECK_FALSE(ModelData::v5l_switches_valid(0, -1, 0.0, 0, 0));
  CHECK_FALSE(ModelData::v5l_switches_valid(0, 1, -0.01, 0, 0));
  CHECK_FALSE(ModelData::v5l_switches_valid(0, 0, 0.0, 3, 0));
  CHECK_FALSE(ModelData::v5l_switches_valid(0, 0, 0.0, 0, 3));
}

TEST_CASE("v5l ModelData(Json): keys parsed; absent = 0; invalid value exits (death test)") {
  Json::Value c;
  c["model_settings"]["et_fpc_once_fix"] = 1; c["model_settings"]["nv_snow_bury_fix"] = 1;
  c["model_settings"]["nv_snow_bury_depth"] = 0.05; c["model_settings"]["rain_int_lai_fix"] = 1;
  c["model_settings"]["nv_rain_surface"] = 2;
  ModelData md(c);
  CHECK(md.et_fpc_once_fix == 1); CHECK(md.nv_snow_bury_fix == 1); CHECK(md.nv_snow_bury_depth == doctest::Approx(0.05));
  CHECK(md.rain_int_lai_fix == 1); CHECK(md.nv_rain_surface == 2);
  ModelData md0(ms("svp_ice_fix", 0));
  CHECK(md0.et_fpc_once_fix == 0); CHECK(md0.nv_snow_bury_fix == 0); CHECK(md0.rain_int_lai_fix == 0);
  CHECK(md0.nv_rain_surface == 0); CHECK(md0.nv_snow_bury_depth == 0.0);
  CHECK(ut::exit_status_of([] { ModelData m(ms("et_fpc_once_fix", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("nv_rain_surface", 3)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("nv_snow_bury_depth", -0.1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms("nv_rain_surface", 1)); }) == 0);
}

// ------------------------------------------------------------------ S1
TEST_CASE("S1 et_fpc_once_fix=0 (legacy): transpiration radiation term uses FPC-scaled radiation") {
  VegDay v; v.day(1.0, 0.25);
  double sw = 200. * 0.25 * (1. - 0.1) * (1. - exp(-0.5 * 1.0));   // swinter (FPC-scaled)
  CHECK(v.ed.d_a2v.swinter == doctest::Approx(sw));
  double par_c = 90. * 0.25 * 0.9 * 4.55;                              // FPC-diluted PAR (umol)
  CHECK(v.ed.d_vegd.m_ppfd == doctest::Approx(par_c / (par_c + 75.)));
  double gst = v.ed.d_vegd.m_ppfd * 1.0 * 1.0 * 0.8 * 0.003;           // m_vpd = 1 (vpd < vpd_open), m_tmin = 1
  double gt = 0.08 * gst / (0.08 + gst);
  double et = v.ve.getPenMonET(12., 800., sw, 1. / gt, 1. / 0.08) * 18. * 3600.;
  CHECK(v.ed.d_v2a.tran == doctest::Approx(et));
}

TEST_CASE("S1 et_fpc_once_fix=1: m_ppfd and transpiration use radiation per covered area") {
  VegDay v; v.ve.fpc_once_fix = 1; v.day(1.0, 0.25);
  double sw_c = 200. * (1. - 0.1) * (1. - exp(-0.5 * 1.0));            // swinter / FPC
  double par_c = 90. * 0.9 * 4.55;
  CHECK(v.ed.d_vegd.m_ppfd == doctest::Approx(par_c / (par_c + 75.)));
  double gst = v.ed.d_vegd.m_ppfd * 0.8 * 0.003;
  double gt = 0.08 * gst / (0.08 + gst);
  CHECK(v.ed.d_v2a.tran == doctest::Approx(v.ve.getPenMonET(12., 800., sw_c, 1. / gt, 1. / 0.08) * 18. * 3600.));
  // per covered area the result no longer depends on FPC (Cohort applies tran *= FPC once)
  VegDay w; w.ve.fpc_once_fix = 1; w.day(1.0, 0.7);
  CHECK(w.ed.d_v2a.tran == doctest::Approx(v.ed.d_v2a.tran).epsilon(1e-12));
  VegDay l1, l2; l1.day(1.0, 0.25); l2.day(1.0, 0.7);
  CHECK(l2.ed.d_v2a.tran > 1.2 * l1.ed.d_v2a.tran);                     // legacy: depends on FPC
}

TEST_CASE("S1: identical to legacy at FPC = 1; wet-canopy evaporation and tran_pet unchanged") {
  VegDay a, b; b.ve.fpc_once_fix = 1;
  a.day(2.0, 1.0); b.day(2.0, 1.0);
  CHECK(a.ed.d_v2a.tran == b.ed.d_v2a.tran);
  CHECK(a.ed.d_vegd.m_ppfd == b.ed.d_vegd.m_ppfd);
  VegDay c, d; d.ve.fpc_once_fix = 1;
  c.ed.d_a2l.rnfl = d.ed.d_a2l.rnfl = 3.0;
  c.day(1.0, 0.3); d.day(1.0, 0.3);
  CHECK(c.ed.d_v2a.evap == d.ed.d_v2a.evap);
  CHECK(c.ed.d_v2a.tran_pet == d.ed.d_v2a.tran_pet);
  CHECK(d.ed.d_v2a.tran > c.ed.d_v2a.tran);
}

// ------------------------------------------------------------------ S2
TEST_CASE("S2 nv_snow_bury_fix=1: buried moss intercepts no snow; stored snow goes to the snowpack") {
  VegDay v; v.ve.nv_snow_bury_fix = 1; v.ve.nonvasc = 1; v.ve.snow_depth = 0.30;
  v.ed.d_a2l.snfl = 2.0; v.ed.d_vegs.snow = 0.2; v.ed.d_atms.ta = -10.;
  v.day(4.0, 0.5, 8.);
  CHECK(v.ed.d_a2v.sinter == 0.0);
  CHECK(v.ed.d_v2a.sublim == 0.0);
  CHECK(v.ed.d_v2g.sthfl == doctest::Approx(1.0));
  CHECK(v.ed.d_v2g.sdrip == doctest::Approx(0.2));
  CHECK(v.ed.d_vegs.snow == doctest::Approx(0.0));
  // snow mass: snfl*fpc + store0 = sthfl + sdrip + sublim + store1
  CHECK(1.0 + 0.2 == doctest::Approx(v.ed.d_v2g.sthfl + v.ed.d_v2g.sdrip + v.ed.d_v2a.sublim + v.ed.d_vegs.snow));
}

TEST_CASE("S2: legacy when off, when not buried, for vascular PFTs; threshold depth respected") {
  auto run = [](int fix, int nv, double depth, double thr) {
    VegDay v; v.ve.nv_snow_bury_fix = fix; v.ve.nonvasc = nv; v.ve.snow_depth = depth;
    v.ve.nv_snow_bury_depth = thr; v.ed.d_a2l.snfl = 2.0; v.ed.d_atms.ta = -10.; v.day(4.0, 0.5, 8.);
    return v.ed.d_a2v.sinter;
  };
  double legacy = std::min(1.0, 0.041 * 4.0);                    // min(snfl*fpc, iscoef*LAI)
  CHECK(run(0, 1, 0.30, 0.0) == doctest::Approx(legacy));
  CHECK(run(1, 1, 0.00, 0.0) == doctest::Approx(legacy));        // no snowpack yet
  CHECK(run(1, 0, 0.30, 0.0) == doctest::Approx(legacy));        // vascular: unaffected
  CHECK(run(1, 1, 0.03, 0.05) == doctest::Approx(legacy));       // below threshold
  CHECK(run(1, 1, 0.06, 0.05) == 0.0);                           // above threshold
}

// ------------------------------------------------------------------ S3
TEST_CASE("S3 rain_int_lai_fix=0 (legacy): max_int = 0.21 P + 0.5 mm, independent of LAI and ircoef") {
  VegDay v; v.ed.d_a2l.rnfl = 2.0; v.ve.envpar.ircoef = 5.0; v.day(0.3, 0.5);
  CHECK(v.ed.d_a2v.rinter == doctest::Approx(0.71));              // P on PFT = 1 mm -> 0.21*1 + 0.5
  VegDay w; w.ed.d_a2l.rnfl = 10.0; w.day(0.3, 0.5);             // P on PFT = 5 mm -> 0.21*5 + 0.5 = 1.55 mm
  CHECK(w.ed.d_a2v.rinter == doctest::Approx(1.55));
}

TEST_CASE("S3 rain_int_lai_fix=1: capacity ircoef*LAI, storage carried, stored water evaporates first") {
  VegDay v; v.ve.rain_int_lai_fix = 1; v.ve.envpar.ircoef = 0.1;
  v.ed.d_a2l.rnfl = 5.0;
  v.day(2.0, 0.5, 0.0);                                           // polar night: no evaporation
  CHECK(v.ed.d_a2v.rinter == doctest::Approx(0.2));                // capacity 0.1 * 2
  CHECK(v.ed.d_v2g.rthfl == doctest::Approx(2.5 - 0.2));
  CHECK(v.ed.d_v2a.evap == 0.0);
  CHECK(v.ed.d_vegs.rwater == doctest::Approx(0.2));               // carried to the next day
  CHECK(v.rain_residual(2.5, 0.0) == doctest::Approx(0.0).epsilon(1e-12));
  // next day: rain fills only the capacity left
  double s0 = v.ed.d_vegs.rwater;
  v.ed.d_a2l.rnfl = 1.0; v.day(2.0, 0.5, 0.0);
  CHECK(v.ed.d_a2v.rinter == doctest::Approx(0.0));
  CHECK(v.rain_residual(0.5, s0) == doctest::Approx(0.0).epsilon(1e-12));
  // dry day with daylight: the stored water evaporates first, transpiration gets the rest of the day
  s0 = v.ed.d_vegs.rwater;
  v.ed.d_a2l.rnfl = 0.0; v.day(2.0, 0.5, 18.);
  CHECK(v.ed.d_a2v.rinter == 0.0);
  CHECK(v.ed.d_v2a.evap == doctest::Approx(0.2));
  CHECK(v.ed.d_vegs.rwater == doctest::Approx(0.0).epsilon(1e-12));
  VegDay dry; dry.ve.rain_int_lai_fix = 1; dry.ve.envpar.ircoef = 0.1; dry.day(2.0, 0.5, 18.);
  CHECK(v.ed.d_v2a.tran < dry.ed.d_v2a.tran);
  CHECK(v.rain_residual(0.0, s0) == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("S3: storage above a shrinking capacity drips; mass conserved") {
  VegDay v; v.ve.rain_int_lai_fix = 1; v.ve.envpar.ircoef = 0.1;
  v.ed.d_vegs.rwater = 0.5;                                       // > capacity 0.1 * 2
  v.ed.d_a2l.rnfl = 1.0; v.day(2.0, 0.5, 0.0);
  CHECK(v.ed.d_a2v.rinter == 0.0);
  CHECK(v.ed.d_v2g.rdrip == doctest::Approx(0.3));
  CHECK(v.ed.d_vegs.rwater == doctest::Approx(0.2));
  CHECK(v.rain_residual(0.5, 0.5) == doctest::Approx(0.0).epsilon(1e-12));
}

// ------------------------------------------------------------------ S4
TEST_CASE("S4 nv_rain_surface=1: moss interception evaporation booked as evap_moss, not canopy evap") {
  VegDay a, b; b.ve.nv_rain_surface = 1;
  a.ve.nonvasc = b.ve.nonvasc = 1;
  a.ed.d_a2l.rnfl = b.ed.d_a2l.rnfl = 3.0;
  a.day(4.0, 0.5); b.day(4.0, 0.5);
  CHECK(a.ed.d_v2a.evap > 0.0);
  CHECK(b.ed.d_v2a.evap == 0.0);
  CHECK(b.ed.d_v2a.evap_moss == a.ed.d_v2a.evap);
  CHECK(b.ed.d_a2v.rinter == a.ed.d_a2v.rinter);
  CHECK(b.ed.d_v2a.tran == a.ed.d_v2a.tran);
  CHECK(b.ed.d_vegs.rwater == a.ed.d_vegs.rwater);
  CHECK(a.ed.d_v2a.evap_moss == 0.0);
  // vascular PFT: unchanged
  VegDay c, d; d.ve.nv_rain_surface = 1; c.ed.d_a2l.rnfl = d.ed.d_a2l.rnfl = 3.0;
  c.day(1.0, 0.3); d.day(1.0, 0.3);
  CHECK(d.ed.d_v2a.evap == c.ed.d_v2a.evap); CHECK(d.ed.d_v2a.evap_moss == 0.0);
}

TEST_CASE("S4 nv_rain_surface=2: rain on moss reaches the ground; stored water drips; vascular unchanged") {
  VegDay v; v.ve.nv_rain_surface = 2; v.ve.nonvasc = 1;
  v.ed.d_vegs.rwater = 0.3; v.ed.d_a2l.rnfl = 3.0; v.day(4.0, 0.5);
  CHECK(v.ed.d_a2v.rinter == 0.0);
  CHECK(v.ed.d_v2g.rthfl == doctest::Approx(1.5));
  CHECK(v.ed.d_v2a.evap == 0.0); CHECK(v.ed.d_v2a.evap_moss == 0.0);
  CHECK(v.ed.d_v2g.rdrip == doctest::Approx(0.3));
  CHECK(v.ed.d_vegs.rwater == doctest::Approx(0.0));
  CHECK(v.rain_residual(1.5, 0.3) == doctest::Approx(0.0).epsilon(1e-12));
  VegDay s; s.ve.nv_rain_surface = 2; s.ve.rain_int_lai_fix = 1; s.ve.nonvasc = 1;   // with S3
  s.ed.d_vegs.rwater = 0.3; s.ed.d_a2l.rnfl = 3.0; s.day(4.0, 0.5);
  CHECK(s.ed.d_a2v.rinter == 0.0); CHECK(s.ed.d_v2a.evap == 0.0);
  CHECK(s.ed.d_v2g.rdrip == doctest::Approx(0.3));
  CHECK(s.rain_residual(1.5, 0.3) == doctest::Approx(0.0).epsilon(1e-12));
  VegDay c, d; d.ve.nv_rain_surface = 2; c.ed.d_a2l.rnfl = d.ed.d_a2l.rnfl = 3.0;
  c.day(1.0, 0.3); d.day(1.0, 0.3);
  CHECK(d.ed.d_a2v.rinter == c.ed.d_a2v.rinter); CHECK(d.ed.d_v2a.evap == c.ed.d_v2a.evap);
}

TEST_CASE("S4: evap_moss accumulates daily -> monthly -> yearly in EnvData") {
  EnvData e;
  e.veg_beginOfYear(); e.veg_beginOfMonth();
  e.d_v2a.evap_moss = 0.4; e.veg_endOfDay(30);
  e.d_v2a.evap_moss = 0.6; e.veg_endOfDay(30);
  CHECK(e.m_v2a.evap_moss == doctest::Approx(1.0));
  e.veg_endOfMonth(0);
  CHECK(e.y_v2a.evap_moss == doctest::Approx(1.0));
}

TEST_CASE("all switches off: per-PFT water balance identical to legacy expectations on a rain day") {
  VegDay v; v.ve.nonvasc = 1; v.ed.d_a2l.rnfl = 4.0; v.ed.d_a2l.snfl = 0.0; v.day(4.0, 0.5);
  CHECK(v.ed.d_a2v.rinter == doctest::Approx(0.21 * 2.0 + 0.5));
  CHECK(v.ed.d_v2a.evap_moss == 0.0);
  CHECK(v.rain_residual(2.0, 0.0) == doctest::Approx(0.0).epsilon(1e-12));
}

}
