// Area 1: climate helpers (svp, vpd + vp_hpa_fix, daily precip, day-of-year,
// daily file loading/length checks) and Area 12a: CO2 indexed by calendar year.
#include "doctest.h"
#include "test_helpers.h"
#include "nc_fixtures.h"

TEST_SUITE("climate") {

// ---------------------------------------------------------------- svp
// Reference: WMO CIMO Guide (WMO-No. 8, 2008/2018), Annex 4.B:
//   ew(t) = 6.112 exp(17.62 t / (243.12 + t))  hPa  (over water)
//   ei(t) = 6.112 exp(22.46 t / (272.62 + t))  hPa  (over ice)
// Tabulated: ew(20 C) = 23.37 hPa, ew(0) = ei(0) = 6.112 hPa,
//            ew(-20 C) = 1.2556 hPa, ei(-20 C) = 1.0324 hPa (also Murphy & Koop 2005: 1.032).
static double wmo_ew(double t) { return 611.2 * std::exp(17.62 * t / (243.12 + t)); }
static double wmo_ei(double t) { return 611.2 * std::exp(22.46 * t / (272.62 + t)); }

TEST_CASE("svp at 0 C equals 611.2 Pa (both branches meet)") {
  CHECK(calculate_saturated_vapor_pressure(0.0f) == doctest::Approx(611.2).epsilon(1e-6));
  // continuity across 0 C
  CHECK(calculate_saturated_vapor_pressure(1e-4f) == doctest::Approx(611.2).epsilon(1e-4));
}

TEST_CASE("svp at +20 C matches WMO over water (2337 Pa) within 0.2%") {
  double v = calculate_saturated_vapor_pressure(20.0f);
  CHECK(v == doctest::Approx(wmo_ew(20.0)).epsilon(0.002));
  CHECK(v == doctest::Approx(2337.0).epsilon(0.003));
}

TEST_CASE("svp is in Pa and increases monotonically with temperature") {
  float prev = 0.0f;
  for (int t = -40; t <= 40; ++t) {
    float v = calculate_saturated_vapor_pressure((float)t);
    CHECK(v > prev);
    prev = v;
  }
  CHECK(calculate_saturated_vapor_pressure(30.0f) > 4000.0f); // Pa, not hPa
}

// KNOWN BUG (reported, not fixed): the t <= 0 branch uses a = 17.27 (the Tetens
// over-water coefficient, normally paired with b = 237.3) together with the WMO
// over-ice denominator b = 272.62. The WMO ice formula uses a = 22.46. Result at -20 C:
// 155.8 Pa vs WMO ice 103.2 Pa (+51%) and WMO water 125.6 Pa (+24%).
TEST_CASE("svp at -20 C matches WMO over ice (103.2 Pa) within 2% [KNOWN BUG: xfail]"
          * doctest::should_fail(true)) {
  double v = calculate_saturated_vapor_pressure(-20.0f);
  CHECK(v == doctest::Approx(wmo_ei(-20.0)).epsilon(0.02));
}
TEST_CASE("svp at -20 C is at least within WMO over-water value +-5% [KNOWN BUG: xfail]"
          * doctest::should_fail(true)) {
  double v = calculate_saturated_vapor_pressure(-20.0f);
  CHECK(v == doctest::Approx(wmo_ew(-20.0)).epsilon(0.05));
}
TEST_CASE("svp at -20 C: document the current (buggy) value 155.8 Pa") {
  // Pins the present behaviour so a fix is noticed: 6.112*exp(17.27*-20/252.62)*100
  CHECK(calculate_saturated_vapor_pressure(-20.0f) == doctest::Approx(155.76).epsilon(1e-3));
}

// ---------------------------------------------------------------- vpd
TEST_CASE("calculate_vpd = svp - vp, floored at 0") {
  CHECK(calculate_vpd(1000.0f, 400.0f) == doctest::Approx(600.0));
  CHECK(calculate_vpd(400.0f, 1000.0f) == 0.0f);
  CHECK(calculate_vpd(500.0f, 500.0f) == 0.0f);
}

static Climate climate_with_one_day(float tair, float vapo_hpa) {
  Climate c;
  c.tair_d = std::vector<float>(1, tair);
  c.prec_d = std::vector<float>(1, 2.0f);
  c.vapo_d = std::vector<float>(1, vapo_hpa);
  c.girr_d = std::vector<float>(1, 300.0f);
  c.nirr_d = std::vector<float>(1, 150.0f);
  return c;
}

TEST_CASE("vpd units: legacy (vp_hpa_fix=0) subtracts hPa from Pa; fix (1) converts hPa->Pa") {
  const float t = 10.0f, e_hpa = 8.0f;          // RH ~ 65 % at 10 C
  const double svp = 611.2 * std::exp(17.63 * 10.0 / 253.12);   // model formula, Pa
  int saved = Climate::vp_hpa_fix;

  Climate::vp_hpa_fix = 0;
  Climate a = climate_with_one_day(t, e_hpa);
  a.derive_daily_secondary_vars();
  CHECK(a.svp_d[0] == doctest::Approx(svp).epsilon(1e-5));
  CHECK(a.vpd_d[0] == doctest::Approx(svp - 8.0).epsilon(1e-5));      // legacy: ~svp
  CHECK(a.rain_d[0] == doctest::Approx(2.0));                          // t>0 -> rain
  CHECK(a.snow_d[0] == 0.0f);

  Climate::vp_hpa_fix = 1;
  Climate b = climate_with_one_day(t, e_hpa);
  b.derive_daily_secondary_vars();
  CHECK(b.vpd_d[0] == doctest::Approx(svp - 800.0).epsilon(1e-5));    // physical VPD
  CHECK(b.vpd_d[0] == doctest::Approx(427.).epsilon(0.01));
  // the input series itself stays in hPa
  CHECK(b.vapo_d[0] == doctest::Approx(8.0));

  // supersaturated input -> vpd floored at 0 in the fixed path
  Climate c = climate_with_one_day(-5.0f, 6.0f);   // 600 Pa > svp(-5)
  c.derive_daily_secondary_vars();
  CHECK(c.vpd_d[0] == 0.0f);
  CHECK(c.snow_d[0] == doctest::Approx(2.0));      // t<=0 -> snow
  Climate::vp_hpa_fix = saved;
}

TEST_CASE("EnvData::update_from_climate: d_atmd.vp in Pa only with vp_hpa_fix=1") {
  int saved = Climate::vp_hpa_fix;
  Climate c;
  c.tair_d = std::vector<float>(DINY, 5.0f);
  c.prec_d = std::vector<float>(DINY, 3.0f);
  c.vapo_d = std::vector<float>(DINY, 7.5f);
  c.girr_d = std::vector<float>(DINY, 300.0f);
  c.nirr_d = std::vector<float>(DINY, 150.0f);
  c.par_d = std::vector<float>(DINY, 60.0f);
  c.co2_d = 400.0f;
  Climate::vp_hpa_fix = 0; c.derive_daily_secondary_vars();
  EnvData e0; e0.d_atms.dsr = 7;
  e0.update_from_climate(c, 3, 9);
  CHECK(e0.d_atmd.vp == doctest::Approx(7.5));
  CHECK(e0.d_atmd.vpd == doctest::Approx(c.svp_d[0] - 7.5).epsilon(1e-6));
  CHECK(e0.d_atms.dsr == 0);   // 3 mm rain >= 2 mm resets days-since-rain

  Climate::vp_hpa_fix = 1; c.derive_daily_secondary_vars();
  EnvData e1; e1.d_atms.dsr = 7;
  e1.update_from_climate(c, 3, 9);
  CHECK(e1.d_atmd.vp == doctest::Approx(750.0));
  CHECK(e1.d_atmd.vpd == doctest::Approx(c.svp_d[0] - 750.0).epsilon(1e-6));
  // relative humidity used for non-vascular PFTs: vp/svp is ~0.86 (was ~0.0086)
  CHECK(e1.d_atmd.vp / e1.d_atmd.svp == doctest::Approx(750.0 / (611.2 * std::exp(17.63 * 5 / 248.12))).epsilon(1e-5));
  Climate::vp_hpa_fix = saved;
}

TEST_CASE("days-since-rain: drizzle (< 2 mm) and snow do not reset the counter") {
  Climate c;
  c.tair_d = std::vector<float>(DINY, 5.0f);
  c.prec_d = std::vector<float>(DINY, 1.9f);
  c.vapo_d = std::vector<float>(DINY, 5.0f);
  c.girr_d = std::vector<float>(DINY, 300.0f);
  c.nirr_d = std::vector<float>(DINY, 150.0f);
  c.par_d = std::vector<float>(DINY, 60.0f);
  c.tair_d[40] = -3.0f; c.prec_d[40] = 10.0f;   // heavy snow day
  c.derive_daily_secondary_vars();
  EnvData e; e.d_atms.dsr = 0;
  for (int d = 0; d < 31; ++d) e.update_from_climate(c, 0, d);
  CHECK(e.d_atms.dsr == 31);
  e.update_from_climate(c, 1, 9);              // doy 40, 10 mm snow
  CHECK(e.d_atms.dsr == 32);
}

// ---------------------------------------------------------------- calendar helpers
TEST_CASE("DINM/DOYINDFST describe a 365-day calendar") {
  int sum = 0;
  for (int m = 0; m < 12; ++m) { CHECK(DOYINDFST[m] == sum); sum += DINM[m]; }
  CHECK(sum == DINY);
  CHECK(DINY == 365);
}

TEST_CASE("day_of_year / doy2month round trip for all 365 days") {
  int doy = 0;
  for (int m = 0; m < 12; ++m) {
    for (int d = 0; d < DINM[m]; ++d, ++doy) {
      CHECK(temutil::day_of_year(m, d) == doy);
      CHECK(temutil::doy2month(doy) == m);
    }
  }
  CHECK(temutil::day_of_year(0, 0) == 0);
  CHECK(temutil::day_of_year(11, 30) == 364);
  CHECK(temutil::day_of_year(2, 0) == 59);   // 1 March, no leap day
}

// ---------------------------------------------------------------- daily precip split
TEST_CASE("calculate_daily_prec conserves the monthly total and returns DINM days") {
  const float temps[3] = {-15.0f, 2.0f, 12.0f};
  const float precs[9] = {0.0f, 3.0f, 20.0f, 50.0f, 80.0f, 110.0f, 160.0f, 220.0f, 400.0f};
  for (int m = 0; m < 12; ++m) {
    for (float t : temps) {
      for (float p : precs) {
        std::vector<float> d = calculate_daily_prec(m, t, p);
        REQUIRE(d.size() == (size_t)DINM[m]);
        double s = 0.0; bool nonneg = true;
        for (float x : d) { s += x; nonneg = nonneg && (x >= 0.0f); }
        CHECK(nonneg);
        CHECK(s == doctest::Approx(p).epsilon(1e-5).scale(1.0));
      }
    }
  }
}

TEST_CASE("calculate_daily_prec: the last day of a 31-day month can receive rain (DINM+1 arrays)") {
  // with T storms the last storm sits on day DINM (1-based); with the old
  // DINM-sized arrays this element was written out of bounds.
  bool any_last = false;
  for (float p : {30.0f, 60.0f, 100.0f, 150.0f, 250.0f}) {
    std::vector<float> d = calculate_daily_prec(0, 10.0f, p);
    if (d[30] > 0.0f) any_last = true;
  }
  CHECK(any_last);
}

// ---------------------------------------------------------------- girr
TEST_CASE("calculate_girr_daily: 365 values whose monthly means equal calculate_girr") {
  for (float lat : {45.0f, 64.9f, 68.7f, 75.0f}) {
    std::vector<float> g = calculate_girr_daily(lat);
    REQUIRE(g.size() == (size_t)DINY);
    for (int m = 0; m < 12; ++m) {
      double s = 0;
      for (int d = 0; d < DINM[m]; ++d) s += g[DOYINDFST[m] + d];
      CHECK(s / DINM[m] == doctest::Approx(calculate_girr(lat, m)).epsilon(1e-4).scale(1.0));
    }
  }
  // polar night: no extraterrestrial radiation at 75N in late December
  std::vector<float> g75 = calculate_girr_daily(75.0f);
  CHECK(g75[355] == doctest::Approx(0.0).epsilon(1e-6));
  // FAO-56 Annex 2 Table 2.6: Ra at 70N on 21 June ~ 41.8 MJ/m2/d = 484 W/m2 (within 5 %)
  std::vector<float> g70 = calculate_girr_daily(70.0f);
  CHECK(g70[171] == doctest::Approx(484.0).epsilon(0.05));
}

// ---------------------------------------------------------------- daily file loader
TEST_CASE("load_daily_from_file: two valid years load, cld/par/girr_daily sized") {
  ut::DailyFileSpec s; s.ndays = 730;
  std::string fn = ut::tmpname("daily_ok");
  ut::write_daily_climate(fn, s);
  Climate c; c.daily_proj = true;
  c.load_proj_climate(fn, 0, 0);
  CHECK(c.daily_series_loaded);
  CHECK(c.tair.size() == 730u);
  CHECK(c.cld.size() == 730u);
  CHECK(c.par.size() == 730u);
  CHECK(c.girr_daily.size() == (size_t)DINY);
  CHECK(c.tseries_start_year == 1950);
  // year 1 of the run is file[365:730]
  c.co2 = std::vector<float>(5, 400.0f);
  c.prepare_daily_driving_data(1, "sc");
  CHECK(c.tair_d.size() == (size_t)DINY);
  CHECK(c.tair_d[0] == doctest::Approx(c.tair[365]));
  CHECK(c.prec_d[100] == doctest::Approx(1.0));
  CHECK(c.co2_d == doctest::Approx(400.0));
  std::remove(fn.c_str());
}

TEST_CASE("load_daily_from_file fails loudly on bad input (death tests)") {
  struct Case { const char* tag; int ndays; const char* var; int idx; float val; };
  const Case cases[] = {
    {"len700", 700, "", -1, 0.0f},                 // not a multiple of 365
    {"nan_t", 365, "tair", 10, NAN},               // non-finite
    {"neg_p", 365, "precip", 200, -0.1f},          // negative precip
    {"neg_n", 365, "nirr", 5, -1.0f},              // negative radiation
    {"neg_v", 365, "vapor_press", 7, -2.0f},       // negative vapour pressure
    {"fill_t", 365, "tair", 3, -999.0f},           // fill value
  };
  for (const Case& k : cases) {
    CAPTURE(k.tag);
    ut::DailyFileSpec s; s.ndays = k.ndays; s.bad_var = k.var; s.bad_index = k.idx; s.bad_value = k.val;
    std::string fn = ut::tmpname(k.tag);
    ut::write_daily_climate(fn, s);
    int st = ut::exit_status_of([&]() { Climate c; c.daily_proj = true; c.load_proj_climate(fn, 0, 0); });
    CHECK(st == EXIT_FAILURE);
    std::remove(fn.c_str());
  }
  // a valid file loads in the child too (control)
  ut::DailyFileSpec s; s.ndays = 365; std::string fn = ut::tmpname("ctl");
  ut::write_daily_climate(fn, s);
  CHECK(ut::exit_status_of([&]() { Climate c; c.daily_proj = true; c.load_proj_climate(fn, 0, 0); }) == 0);
  std::remove(fn.c_str());
}

TEST_CASE("prepare_daily_driving_data: year index beyond the daily series fails; non-sc stage fails") {
  ut::DailyFileSpec s; s.ndays = 365; std::string fn = ut::tmpname("short");
  ut::write_daily_climate(fn, s);
  CHECK(ut::exit_status_of([&]() { Climate c; c.daily_proj = true; c.load_proj_climate(fn, 0, 0);
                                   c.co2 = std::vector<float>(3, 400.f); c.prepare_daily_driving_data(1, "sc"); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([&]() { Climate c; c.daily_proj = true; c.load_proj_climate(fn, 0, 0);
                                   c.co2 = std::vector<float>(3, 400.f); c.prepare_daily_driving_data(0, "tr"); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([&]() { Climate c; c.daily_proj = true; c.load_proj_climate(fn, 0, 0);
                                   c.co2 = std::vector<float>(3, 400.f); c.prepare_daily_driving_data(0, "sc"); }) == 0);
  std::remove(fn.c_str());
}

} // TEST_SUITE climate

TEST_SUITE("co2_nlost") {

// ---------------------------------------------------------------- CO2 by calendar year
TEST_CASE("index_proj_co2_by_year: offset = first climate year - first co2 year") {
  // climate file: time units 'days since 1901-01-01', first value 49*365 -> first year 1950
  ut::DailyFileSpec s; s.ndays = 365; s.units = "days since 1901-01-01";
  std::string fc = ut::tmpname("clim1901"), fco2 = ut::tmpname("co2_1901");
  ut::write_daily_climate(fc, s);
  // rewrite time values to start at 49*365
  { int nc, vt; ut::ncchk(nc_open(fc.c_str(), NC_WRITE, &nc)); ut::ncchk(nc_inq_varid(nc, "time", &vt));
    std::vector<int> t(365); for (int i = 0; i < 365; ++i) t[i] = 49 * 365 + i;
    ut::ncchk(nc_put_var_int(nc, vt, t.data())); ut::ncchk(nc_close(nc)); }
  ut::write_co2(fco2, 1901, 200);       // years 1901..2100, co2 = 300 + (year-1900)
  Climate c; c.load_proj_co2(fco2);
  c.index_proj_co2_by_year(fc, fco2, 75);
  CHECK(c.co2_idx_offset == 49);
  CHECK(c.co2_for_year_index(0) == doctest::Approx(350.0));   // 1950
  CHECK(c.co2_for_year_index(74) == doctest::Approx(424.0));  // 2024
  // other stages keep offset 0 after re-loading
  c.load_proj_co2(fco2);
  CHECK(c.co2_idx_offset == 0);
  CHECK(c.co2_for_year_index(0) == doctest::Approx(301.0));
  std::remove(fc.c_str()); std::remove(fco2.c_str());
}

TEST_CASE("index_proj_co2_by_year: uncovered years, gaps, wrong calendar fail loudly") {
  ut::DailyFileSpec s; s.ndays = 365; s.units = "days since 1950-01-01";
  std::string fc = ut::tmpname("clim1950"), fco2 = ut::tmpname("co2a");
  ut::write_daily_climate(fc, s);
  ut::write_co2(fco2, 1960, 100);   // starts after the climate
  CHECK(ut::exit_status_of([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 10); }) == EXIT_FAILURE);
  ut::write_co2(fco2, 1901, 60);    // ends 1960 < 1950+75-1
  CHECK(ut::exit_status_of([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 75); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 11); }) == 0);
  ut::write_co2(fco2, 1901, 200, true);   // non-consecutive years
  CHECK(ut::exit_status_of([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 10); }) == EXIT_FAILURE);
  ut::DailyFileSpec g = s; g.calendar = "gregorian"; ut::write_daily_climate(fc, g);
  ut::write_co2(fco2, 1901, 200);
  CHECK(ut::exit_status_of([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 10); }) == EXIT_FAILURE);
  // co2_for_year_index out of range
  CHECK(ut::exit_status_of([&]() { Climate c; c.co2 = std::vector<float>(3, 1.f); c.co2_for_year_index(3); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([&]() { Climate c; c.co2 = std::vector<float>(3, 1.f); c.co2_for_year_index(-1); }) == EXIT_FAILURE);
  std::remove(fc.c_str()); std::remove(fco2.c_str());
}

} // TEST_SUITE co2_nlost
