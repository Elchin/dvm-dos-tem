// Unit tests for the safe bug fixes ported from hpc-tvc daily-climate-mode onto rts-regrow
// (same lightweight CHECK framework as test_rts_units.cpp; no external framework).
// Build/run: make rts-unit-test   (run from the repository root)
//
// Ported/adapted from the hpc-tvc doctest suite (tests/unit/test_climate.cpp,
// test_stefan_water.cpp, test_helpers.h). Tests that could crash or hang on unfixed code
// (Richards guards, snow-layer invariant, CO2 fail-loudly) run in a forked child with an
// alarm() timeout, so an endless loop is reported as a FAIL (signal 14) instead of hanging.
//
//   fix (daily commit)          tests
//   17d53b56 precip DINM+1      test_daily_prec_*            (Feb precip conservation)
//   17d53b56 NUL-safe units     test_time_units_not_nul_terminated
//   17d53b56 cohort.md          test_runner_cohort_md_points_to_runner_md
//   7a7fe61a/032e8f13 Richards  test_richards_*              (sizing, clear, guards, no hang)
//   e3fcc91c empty snow layers  test_empty_snow_layer_removed
//   0e22205a SC CO2 by year     test_co2_*
//   333d7197 root_water_up OOB, 9b3a4840 tdrv init: no direct unit test (local/stack-only
//   effects); covered by the bit-identity check documented in the port report.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>
#include <unistd.h>
#include <sys/wait.h>
#include <netcdf.h>
#include <json/json.h>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include "TEMLogger.h"
#include "TEMUtilityFunctions.h"
#define private public
#define protected public
#include "layerconst.h"
#include "timeconst.h"
#include "physicalconst.h"
#include "Climate.h"
#include "Layer.h"
#include "SoilLayer.h"
#include "SnowLayer.h"
#include "Ground.h"
#include "Richards.h"
#include "ModelData.h"
#include "Runner.h"
#undef private
#undef protected

// -DPORT_PRE_FIX: build against the UNFIXED tree (headers+objects) to show which tests fail/hang
// there; the API-dependent tests (MAX_RICH_IND sizing, co2 helpers) are compiled out.
#ifdef PORT_PRE_FIX
#define MAX_RICH_IND (MAX_SNW_LAY + MAX_SOI_LAY + 2)
#endif

// free function in Climate.cpp (not declared in a header)
std::vector<float> calculate_daily_prec(const int midx, const float mta, const float mprec);

static int g_fail = 0, g_pass = 0;
#define CHECK(cond) do { if (cond) { g_pass++; } else { g_fail++; \
  std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

// Run fn in a forked child with a timeout. Returns the child's exit status, 128+signal if it
// was killed (e.g. 142 = SIGALRM after a hang, 139 = SIGSEGV).
static int run_isolated(const std::function<int()>& fn, unsigned timeout_s = 10) {
  std::fflush(stdout); std::fflush(stderr);
  pid_t pid = fork();
  if (pid == 0) {
    if (!std::freopen("/dev/null", "w", stderr)) {}
    if (!std::freopen("/dev/null", "w", stdout)) {}
    alarm(timeout_s);
    int rc = 3;
    try { rc = fn(); } catch (...) { rc = 4; }
    _exit(rc);
  }
  int st = 0;
  waitpid(pid, &st, 0);
  if (WIFEXITED(st)) return WEXITSTATUS(st);
  if (WIFSIGNALED(st)) return 128 + WTERMSIG(st);
  return -1;
}

// ---------------------------------------------------------------- 17d53b56 precip
static void test_daily_prec_conserves_monthly_total() {
  // ported from hpc-tvc test_climate.cpp "calculate_daily_prec conserves the monthly total"
  const float temps[] = {-15.f, 2.f, 12.f};
  const float precs[] = {0.f, 3.f, 20.f, 50.f, 80.f, 110.f, 160.f, 220.f, 400.f};
  int bad = 0, bad_feb = 0, n = 0;
  for (int m = 0; m < 12; ++m) for (float t : temps) for (float p : precs) {
    std::vector<float> d = calculate_daily_prec(m, t, p);
    double s = 0.; bool nonneg = true;
    for (float v : d) { s += v; nonneg &= (v >= 0.f); }
    // doctest Approx(p).epsilon(1e-5) semantics, scale 1
    const bool ok = d.size() == (size_t)DINM[m] && nonneg
                    && std::fabs(s - p) <= 1e-5 * (1.0 + std::max(std::fabs(s), (double)p));
    ++n; if (!ok) { ++bad; if (m == 1) ++bad_feb;
      std::printf("  precip m=%d t=%g p=%g: size=%zu sum=%g\n", m, t, p, d.size(), s); }
  }
  std::printf("  calculate_daily_prec: %d/%d month/temp/precip combos not conserved (%d in Feb)\n", bad, n, bad_feb);
  CHECK(bad == 0);
}

static void test_daily_prec_february_explicit() {
  // the RTS symptom: February lost 95-100 % of its precipitation
  for (float p : {10.f, 30.f, 60.f}) {
    std::vector<float> d = calculate_daily_prec(1, -15.f, p);
    double s = 0.; for (float v : d) s += v;
    CHECK(d.size() == 28);
    CHECK(std::fabs(s - p) <= 1e-4 * p);
  }
}

static void test_daily_prec_last_day_can_rain() {
  // ported: last day of a 31-day month receives rain (the old 1-based write went past the end)
  int rained = 0;
  for (float p : {30.f, 60.f, 100.f, 150.f, 250.f}) {
    std::vector<float> d = calculate_daily_prec(0, 10.f, p);
    if (d.size() == 31 && d[30] > 0.f) ++rained;
  }
  CHECK(rained > 0);
}

// ---------------------------------------------------------------- 17d53b56 NUL-safe units
static void test_time_units_not_nul_terminated() {
  const std::string f = "/tmp/rts_port_units_" + std::to_string(getpid()) + ".nc";
  int nc, dim, v;
  nc_create(f.c_str(), NC_CLOBBER, &nc);
  nc_def_dim(nc, "time", 3, &dim);
  nc_def_var(nc, "time", NC_DOUBLE, 1, &dim, &v);
  const char* u = "days since 1901-1-1 0:0:0";
  nc_put_att_text(nc, v, "units", std::strlen(u), u);   // stored WITHOUT a trailing NUL
  nc_enddef(nc);
  double t[3] = {0., 31., 59.}; nc_put_var_double(nc, v, t);
  nc_close(nc);
  CHECK(temutil::get_timeseries_start_year(f) == 1901);
  std::remove(f.c_str());
}

// ---------------------------------------------------------------- 17d53b56 cohort.md
static void test_runner_cohort_md_points_to_runner_md() {
  const std::string data = "testing-data/rts-synthetic/";
  if (!boost::filesystem::exists(data + "historic-climate.nc") || !boost::filesystem::exists("config/config.js")) {
    std::printf("  SKIP runner test (run from repo root)\n"); return;
  }
  Json::Value cfg = temutil::parse_control_file("config/config.js");
  const char* files[][2] = {{"hist_climate_file", "historic-climate.nc"}, {"proj_climate_file", "projected-climate.nc"},
    {"veg_class_file", "vegetation.nc"}, {"drainage_file", "drainage.nc"}, {"soil_texture_file", "soil-texture.nc"},
    {"co2_file", "co2.nc"}, {"proj_co2_file", "projected-co2.nc"}, {"runmask_file", "run-mask.nc"},
    {"topo_file", "topo.nc"}, {"fri_fire_file", "fri-fire.nc"}, {"hist_exp_fire_file", "historic-explicit-fire.nc"},
    {"proj_exp_fire_file", "projected-explicit-fire.nc"}, {"hist_exp_thermokarst_file", "historic-explicit-thermokarst.nc"},
    {"proj_exp_thermokarst_file", "projected-explicit-thermokarst.nc"}};
  for (auto& kv : files) cfg["IO"][kv[0]] = data + kv[1];
  cfg["IO"]["parameter_dir"] = "parameters/";
  const int rc = run_isolated([&]() {
    ModelData md(cfg);
    Runner runner(md, false, 0, 0);
    // after the ctor returned, cohort.md must point at the Runner-owned copy, not the
    // (dead) by-value ctor argument
    return (runner.cohort.md == &runner.md) ? 0 : 1;
  }, 60);
  if (rc != 0) std::printf("  runner test child status %d\n", rc);
  CHECK(rc == 0);
}

// ---------------------------------------------------------------- 7a7fe61a / 032e8f13 Richards
static SoilLayer* make_soil_layer(double z, double dz, double poro, double liq) {
  // ported from hpc-tvc tests/unit/test_helpers.h make_soil_layer
  SoilLayer* l = new SoilLayer();
  l->isMineral = true;
  l->z = z; l->dz = dz; l->poro = poro;
  l->liq = liq; l->ice = 0.; l->tem = 2.0;
  l->frozen = -1; l->frozenfrac = 0.;
  l->tcsolid = 2.0; l->bulkden = 2700.0 * (1.0 - poro); l->vhcsolid = 2.0e6;
  l->derivePhysicalProperty();
  l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0;
  l->pce_f = 0.; l->pce_t = 0.;
  l->solind = 0; l->indl = 0;
  return l;
}

struct Column {           // Ground owning a list of soil layers (top to bottom)
  Ground g; std::vector<Layer*> L;
  explicit Column(int n) {
    for (int i = 0; i < n; ++i) {
      SoilLayer* l = make_soil_layer(0.1 * i, 0.1, 0.4, 0.5 * 0.4 * DENLIQ * 0.1);
      g.insertBack(l); L.push_back(l);
    }
    int k = 1;
    for (Layer* l = g.toplayer; l != NULL; l = l->nextl, ++k) { l->indl = k; if (l->isSoil) l->solind = k; }
    g.fstsoill = L.front(); g.lstsoill = L.back(); g.fstminel = L.front(); g.lstminel = L.back();
    g.frontsz.clear(); g.frontstype.clear();
  }
};

#ifndef PORT_PRE_FIX
static void test_richards_array_sizing_and_clear() {
  CHECK(MAX_RICH_IND == MAX_SNW_LAY + MAX_SOI_LAY + 2);
  Richards r;
  CHECK(sizeof(r.qtrans) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.layer_drain) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.percolation) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.nodemm) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.tridiag_error) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(MAX_SNW_LAY + MAX_SOI_LAY + 1 < MAX_RICH_IND);   // largest index used: drainind+1
  for (int i = 0; i < MAX_RICH_IND; ++i) { r.qtrans[i] = 1.; r.qin[i] = 1.; r.qout[i] = 1.; r.amx[i] = 1.;
    r.rmx[i] = 1.; r.tridiag_error[i] = 1.; r.percolation[i] = 1.; }
  r.clearRichardsArrays();
  bool all0 = true;
  for (int i = 0; i < MAX_RICH_IND; ++i)
    all0 &= r.qtrans[i] == 0. && r.qin[i] == 0. && r.qout[i] == 0. && r.amx[i] == 0. && r.rmx[i] == 0. && r.tridiag_error[i] == 0.;
  CHECK(all0);
}
#endif

// child return codes: 0 = expected runtime_error with the expected message, 1 = no throw,
// 2 = other exception/message
static int richards_expect_throw(int top_indl, int drain_indl, double delta_t, const char* needle) {
  Column c(4);
  if (top_indl != -999) c.L[0]->indl = top_indl;
  if (drain_indl != -999) c.L[3]->indl = drain_indl;
  Richards r; r.delta_t = delta_t;
  double trans[MAX_SOI_LAY] = {0};
  try { r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0); }
  catch (const std::runtime_error& e) { return std::string(e.what()).find(needle) != std::string::npos ? 0 : 2; }
  return 1;
}

static void test_richards_top_index_throws() {
  int rc = run_isolated([]() { return richards_expect_throw(0, -999, SEC_IN_DAY, "layer index out of range"); });
  if (rc) std::printf("  richards top index<1: child status %d\n", rc);
  CHECK(rc == 0);
}

static void test_richards_drain_index_throws() {
  int rc = run_isolated([]() { return richards_expect_throw(-999, MAX_RICH_IND - 1, SEC_IN_DAY, "layer index out of range"); });
  if (rc) std::printf("  richards drainind+1>=MAX_RICH_IND: child status %d\n", rc);
  CHECK(rc == 0);
}

static void test_richards_nonfinite_substep_throws_no_hang() {
  // delta_t = 0 / NaN / negative used to loop forever; must throw within the timeout
  const double dts[] = {0.0, std::nan(""), -86400.0};
  for (double dt : dts) {
    int rc = run_isolated([dt]() { return richards_expect_throw(-999, -999, dt, "substep loop cannot finish"); }, 10);
    if (rc) std::printf("  richards delta_t=%g: child status %d%s\n", dt, rc, rc == 142 ? " (HANG, killed by alarm)" : "");
    CHECK(rc == 0);
  }
}

static void test_richards_valid_column_runs() {
  // a valid 4-layer column with the normal one-day step must not trip either guard
  int rc = run_isolated([]() {
    Column c(4); Richards r; r.delta_t = SEC_IN_DAY; double trans[MAX_SOI_LAY] = {0};
    try { r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0); } catch (...) { return 1; }
    return 0; }, 30);
  if (rc) std::printf("  richards valid column: child status %d\n", rc);
  CHECK(rc == 0);
}

// ---------------------------------------------------------------- e3fcc91c empty snow layers
static void test_empty_snow_layer_removed() {
  int rc = run_isolated([]() {
    Column c(3);
    SnowLayer* s = new SnowLayer(); s->ice = 0.; s->liq = 0.; s->dz = 0.; s->tem = 0.;
    c.g.insertFront(s);
    c.g.snow.extramass = 0.;
    c.g.constructSnowLayers(0.0, 1.0);          // no accumulation, no ablation
    return (c.g.toplayer != NULL && c.g.toplayer->isSoil) ? 0 : 1;
  });
  if (rc) std::printf("  empty snow layer: child status %d\n", rc);
  CHECK(rc == 0);
  // a snow layer that still holds ice is kept
  rc = run_isolated([]() {
    Column c(3);
    SnowLayer* s = new SnowLayer(); s->ice = 5.; s->liq = 0.; s->dz = 0.02; s->tem = -1.;
    c.g.insertFront(s); c.g.snow.extramass = 0.;
    c.g.constructSnowLayers(0.0, -1.0);
    return (c.g.toplayer != NULL && c.g.toplayer->isSnow && c.g.toplayer->ice == 5.) ? 0 : 1;
  });
  CHECK(rc == 0);
}

#ifndef PORT_PRE_FIX
// ---------------------------------------------------------------- 0e22205a SC CO2 by calendar year
static std::string tmpnc(const char* tag) { return "/tmp/rts_port_" + std::string(tag) + "_" + std::to_string(getpid()) + ".nc"; }

static void write_clim_time(const std::string& f, const char* units, const char* cal, double t0) {
  int nc, dim, v;
  nc_create(f.c_str(), NC_CLOBBER, &nc);
  nc_def_dim(nc, "time", 12, &dim);
  nc_def_var(nc, "time", NC_DOUBLE, 1, &dim, &v);
  nc_put_att_text(nc, v, "units", std::strlen(units), units);
  if (cal) nc_put_att_text(nc, v, "calendar", std::strlen(cal), cal);
  nc_enddef(nc);
  double t[12]; for (int i = 0; i < 12; ++i) t[i] = t0 + 30.4 * i;
  nc_put_var_double(nc, v, t); nc_close(nc);
}

static void write_co2(const std::string& f, int y0, int n, bool gap = false) {
  // co2 = 300 + (year - 1900)
  int nc, dim, vy, vc;
  nc_create(f.c_str(), NC_CLOBBER, &nc);
  nc_def_dim(nc, "year", n, &dim);
  nc_def_var(nc, "year", NC_INT, 1, &dim, &vy);
  nc_def_var(nc, "co2", NC_FLOAT, 1, &dim, &vc);
  nc_enddef(nc);
  std::vector<int> y(n); std::vector<float> c(n);
  for (int i = 0; i < n; ++i) { y[i] = y0 + i + ((gap && i >= n / 2) ? 1 : 0); c[i] = 300.f + (y[i] - 1900); }
  nc_put_var_int(nc, vy, y.data()); nc_put_var_float(nc, vc, c.data()); nc_close(nc);
}

static void test_co2_indexed_by_calendar_year() {
  // ported from hpc-tvc test_climate.cpp "index_proj_co2_by_year: offset = first climate year - first co2 year"
  const std::string fc = tmpnc("clim"), fco2 = tmpnc("co2");
  write_clim_time(fc, "days since 1901-1-1 0:0:0", "365_day", 49 * 365.0);   // first year 1950
  write_co2(fco2, 1901, 200);
  int rc = run_isolated([&]() {
    Climate c; c.load_proj_co2(fco2);
    c.index_proj_co2_by_year(fc, fco2, 75);
    bool ok = c.co2_idx_offset == 49 && std::fabs(c.co2_for_year_index(0) - 350.f) < 1e-4
              && std::fabs(c.co2_for_year_index(74) - 424.f) < 1e-4;
    c.load_proj_co2(fco2);                      // other stages: offset back to 0
    ok = ok && c.co2_idx_offset == 0 && std::fabs(c.co2_for_year_index(0) - 301.f) < 1e-4;
    return ok ? 0 : 1; });
  CHECK(rc == 0);
  // our inputs: projected climate and projected-co2 both start 2016 -> offset 0
  write_clim_time(fc, "days since 2016-1-1 0:0:0", "365_day", 0.0);
  write_co2(fco2, 2016, 35);
  rc = run_isolated([&]() { Climate c; c.load_proj_co2(fco2); c.index_proj_co2_by_year(fc, fco2, 35);
                            return c.co2_idx_offset == 0 ? 0 : 1; });
  CHECK(rc == 0);
  std::remove(fc.c_str()); std::remove(fco2.c_str());
}

static void test_co2_fails_loudly() {
  // ported from hpc-tvc test_climate.cpp "uncovered years, gaps, wrong calendar fail loudly"
  const std::string fc = tmpnc("clim2"), fco2 = tmpnc("co2b");
  write_clim_time(fc, "days since 1950-01-01", "365_day", 0.0);
  auto idx = [&](int sc) { return run_isolated([&, sc]() { Climate c; c.load_proj_co2(fco2);
                                                          c.index_proj_co2_by_year(fc, fco2, sc); return 0; }); };
  write_co2(fco2, 1960, 100);  CHECK(idx(10) == EXIT_FAILURE);   // starts after the climate
  write_co2(fco2, 1901, 60);   CHECK(idx(75) == EXIT_FAILURE);   // ends 1960 < 1950+75-1
                               CHECK(idx(11) == 0);              // 1950..1960 covered
  write_co2(fco2, 1901, 200, true); CHECK(idx(10) == EXIT_FAILURE);   // non-consecutive years
  write_clim_time(fc, "days since 1950-01-01", "gregorian", 0.0);
  write_co2(fco2, 1901, 200);  CHECK(idx(10) == EXIT_FAILURE);   // unsupported calendar
  CHECK(run_isolated([]() { Climate c; c.co2 = std::vector<float>(3, 1.f); c.co2_for_year_index(3); return 0; }) == EXIT_FAILURE);
  CHECK(run_isolated([]() { Climate c; c.co2 = std::vector<float>(3, 1.f); c.co2_for_year_index(-1); return 0; }) == EXIT_FAILURE);
  std::remove(fc.c_str()); std::remove(fco2.c_str());
}
#endif  // PORT_PRE_FIX

int main() {
  test_daily_prec_conserves_monthly_total();
  test_daily_prec_february_explicit();
  test_daily_prec_last_day_can_rain();
  test_time_units_not_nul_terminated();
  test_runner_cohort_md_points_to_runner_md();
#ifndef PORT_PRE_FIX
  test_richards_array_sizing_and_clear();
#endif
  test_richards_top_index_throws();
  test_richards_drain_index_throws();
  test_richards_nonfinite_substep_throws_no_hang();
  test_richards_valid_column_runs();
  test_empty_snow_layer_removed();
#ifndef PORT_PRE_FIX
  test_co2_indexed_by_calendar_year();
  test_co2_fails_loudly();
#endif
  std::printf("Daily-port fix unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
