// Unit tests for the mean-preserving monthly->daily forcing (model_settings.monthly2daily_mean_preserving).
// Groups 1-4 ported from 23322af0 tests/rts/test_rts_units.cpp; group 5 checks the switch in prepare_daily_driving_data.
// Build/run: tests/rts/run_m2d_units.sh (links the worktree obj/*.o).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
// white-box access to the legacy interpolation helpers (eq_range, interpolation_range)
#define private public
#include "Climate.h"
#undef private
#include "TEMUtilityFunctions.h"

static int g_fail = 0, g_pass = 0;
#define CHECK(cond) do { if (cond) { g_pass++; } else { g_fail++; \
  std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

// ---------------- monthly -> daily (mean-preserving, Killworth) ----------------
static const int DIM12[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

static std::vector<float> synth_series(int nyr, double mean, double amp, double noise, unsigned seed, bool clipzero) {
  std::vector<float> v; unsigned s = seed;
  for (int y = 0; y < nyr; ++y) for (int m = 0; m < 12; ++m) {
    s = s * 1103515245u + 12345u;
    double r = ((s >> 8) & 0xFFFF) / 65535.0 - 0.5;
    double x = mean - amp * std::cos(2 * M_PI * (m + 0.5) / 12.0) + noise * r;
    if (clipzero && x < 0) x = 0;      // polar-night radiation: exact zeros
    v.push_back((float)x);
  }
  return v;
}

static void check_means(const std::vector<float>& d, const std::vector<float>& m12, int& bad) {
  int k = 0;
  for (int m = 0; m < 12; ++m) {
    double sum = 0; for (int j = 0; j < DIM12[m]; ++j) sum += d[k++];
    double M = m12[m], tol = 1e-6 * std::max(1.0, std::fabs(M));
    if (std::fabs(sum / DIM12[m] - M) > tol) { bad++;
      std::printf("  month %d mean %.9g vs %.9g\n", m, sum / DIM12[m], M); }
  }
}

static void test_m2d_mean_preserving() {
  struct V { double mean, amp, noise; bool nonneg, clip; } vars[] = {
    {-8.0, 22.0, 6.0, false, false},   // tair (Gydan-like)
    {0.4, 0.5, 0.2, true, true},       // vapo (kPa) with zeros
    {80.0, 110.0, 30.0, true, true},   // nirr W/m2, polar night zeros
    {40.0, 55.0, 15.0, true, true}};   // par
  int nyr = 20;
  for (const V& v : vars) {
    std::vector<float> ser = synth_series(nyr, v.mean, v.amp, v.noise, 7u, v.clip);
    int bad_mean = 0, bad_neg = 0, bad_len = 0;
    for (int y = 0; y < nyr; ++y) {
      std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window(ser, y), v.nonneg);
      if (d.size() != 365) { bad_len++; continue; }
      check_means(d, std::vector<float>(ser.begin() + 12 * y, ser.begin() + 12 * y + 12), bad_mean);
      if (v.nonneg) for (float x : d) if (x < 0) bad_neg++;
    }
    CHECK(bad_len == 0); CHECK(bad_mean == 0); CHECK(bad_neg == 0);
  }
  // climatology path (eq stage)
  std::vector<float> c12 = synth_series(1, 50.0, 70.0, 0.0, 3u, true);
  std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window_clim(c12), true);
  int bad = 0; check_means(d, c12, bad); CHECK(bad == 0);
  CHECK(*std::min_element(d.begin(), d.end()) >= 0.0f);
}

static void test_m2d_continuity() {
  // Dec 16..31 and Jan 1..15 lie on the same linear segment, so the jump
  // Dec31->Jan1 across years must equal the Dec30->Dec31 step.
  std::vector<float> ser = synth_series(15, -8.0, 22.0, 6.0, 11u, false);
  double worst = 0;
  for (int y = 0; y + 1 < 15; ++y) {
    std::vector<float> a = Climate::monthly2daily_mp(Climate::mp_window(ser, y), false);
    std::vector<float> b = Climate::monthly2daily_mp(Climate::mp_window(ser, y + 1), false);
    worst = std::max(worst, (double)std::fabs((b[0] - a[364]) - (a[364] - a[363])));
  }
  CHECK(worst < 1e-5);
  // periodic climatology: Dec31 -> Jan1 of same curve continuous
  std::vector<float> c12 = synth_series(1, -8.0, 22.0, 0.0, 3u, false);
  std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window_clim(c12), false);
  CHECK(std::fabs((d[0] - d[364]) - (d[364] - d[363])) < 1e-5);
  // within-year: no step larger than the max monthly-mean jump
  double maxstep = 0; for (int i = 1; i < 365; ++i) maxstep = std::max(maxstep, (double)std::fabs(d[i] - d[i-1]));
  CHECK(maxstep < 3.0);
}

static void test_m2d_constant_and_tdd() {
  // constant input -> constant output exactly; and the convex July peak is
  // no longer flattened (legacy loses TDD)
  std::vector<float> k(36, 4.25f);
  std::vector<float> d = Climate::monthly2daily_mp(k, false);
  double mx = 0; for (float x : d) mx = std::max(mx, (double)std::fabs(x - 4.25f));
  CHECK(mx < 1e-6);
  std::vector<float> m12 = {-25, -26, -22, -15, -5, 5, 12, 9, 3, -7, -18, -23};
  Climate clm;
  std::vector<float> leg = clm.monthly2daily(std::vector<float>{-23, -25, -26, -22, -15, -5, 5, 12, 9, 3, -7, -18, -23, -25});
  std::vector<float> mp = Climate::monthly2daily_mp(Climate::mp_window_clim(m12), false);
  double tdd_true = 0, tdd_leg = 0, tdd_mp = 0;
  for (int m = 0; m < 12; ++m) tdd_true += std::max(0.0f, m12[m]) * DIM12[m];
  for (int i = 0; i < 365; ++i) { tdd_leg += std::max(0.0f, leg[i]); tdd_mp += std::max(0.0f, mp[i]); }
  CHECK(leg.size() == 366);                 // legacy container length unchanged
  CHECK(tdd_leg < tdd_mp);                  // legacy under-estimates thaw
  CHECK(std::fabs(tdd_mp - tdd_true) < 0.06 * tdd_true);
}

static void test_m2d_legacy_bit_identical() {
  // The legacy routine must match a verbatim re-implementation of the old code.
  static const float arr[] = { -15.5, 15.5, 45.0, 74.5, 105.0, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 380.5 };
  std::vector<float> v14 = {-23.1f, -25.2f, -26.3f, -22.4f, -15.5f, -5.6f, 5.7f, 12.8f, 9.9f, 3.1f, -7.2f, -18.3f, -23.4f, -25.5f};
  std::vector<float> ref;
  for (int i = 1; i < 14; ++i) {
    std::vector<float> psd = temutil::resample(std::make_pair(arr[i-1], v14[i-1]), std::make_pair(arr[i], v14[i]),
                                               int(floor(arr[i-1])), int(floor(arr[i])), 1);
    ref.insert(ref.end(), psd.begin(), psd.end());
  }
  ref = std::vector<float>(ref.begin() + 16, ref.end() - 14);
  Climate clm;
  std::vector<float> got = clm.monthly2daily(v14);
  CHECK(got.size() == ref.size());
  bool same = got.size() == ref.size();
  for (size_t i = 0; same && i < got.size(); ++i) same = (got[i] == ref[i]);
  CHECK(same);
}

// deep_carbon_v2 burial (Soil_Bgc::deltastate part II) must only move C and N
// between layers. Column: fibric | humic | mineral 2.6 m | mineral 2.6 m, so the
// second mineral layer exhausts xtopmlthick (5.2 m) exactly.

// ---- switch behaviour through Climate::prepare_daily_driving_data ----
static void fill_climate(Climate& c, int nyr) {
  c.tair = synth_series(nyr, -8.0, 22.0, 6.0, 21u, false);
  c.prec = synth_series(nyr, 30.0, 20.0, 10.0, 22u, true);
  c.nirr = synth_series(nyr, 80.0, 110.0, 30.0, 23u, true);
  c.vapo = synth_series(nyr, 0.4, 0.5, 0.2, 24u, true);
  c.par  = synth_series(nyr, 40.0, 55.0, 15.0, 25u, true);
  c.cld  = std::vector<float>(12 * nyr, 0.5f);
  c.girr = synth_series(1, 150.0, 140.0, 0.0, 26u, true);
  c.co2  = std::vector<float>(nyr, 300.0f);
  c.avgX_tair = std::vector<float>(c.tair.begin(), c.tair.begin() + 12);
  c.avgX_prec = std::vector<float>(c.prec.begin(), c.prec.begin() + 12);
  c.avgX_nirr = std::vector<float>(c.nirr.begin(), c.nirr.begin() + 12);
  c.avgX_vapo = std::vector<float>(c.vapo.begin(), c.vapo.begin() + 12);
  c.avgX_co2 = 300.0; c.avgX_ch4 = 0.0;
}
static bool same_vec(const std::vector<float>& a, const std::vector<float>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i) if (!(a[i] == b[i])) return false;
  return true;
}
static void test_m2d_switch_paths() {
  Climate dflt; CHECK(dflt.m2d_mean_preserving == false);          // default off
  const int nyr = 6;
  for (const char* st : {"eq", "tr"}) {
    const bool eq = std::string(st) == "eq";
    for (int iy : {0, 2, nyr - 1}) {
      // OFF: daily containers are exactly the legacy interpolation (bit-identical path)
      Climate off; fill_climate(off, nyr); off.m2d_mean_preserving = false;
      off.prepare_daily_driving_data(iy, st);
      Climate ref; fill_climate(ref, nyr);
      std::vector<float> t = eq ? ref.monthly2daily(ref.eq_range(ref.avgX_tair)) : ref.monthly2daily(ref.interpolation_range(ref.tair, iy));
      std::vector<float> p = eq ? ref.monthly2daily(ref.eq_range(ref.par)) : ref.monthly2daily(ref.interpolation_range(ref.par, iy));
      CHECK(same_vec(off.tair_d, t)); CHECK(same_vec(off.par_d, p));
      // ON: each month's daily mean equals the monthly input (tair, vapo, par; nirr has no post-clip change)
      Climate on; fill_climate(on, nyr); on.m2d_mean_preserving = true;
      on.prepare_daily_driving_data(iy, st);
      CHECK(on.tair_d.size() == 365); CHECK(on.nirr_d.size() == 365); CHECK(on.par_d.size() == 365);
      int bad = 0;
      auto m12 = [&](const std::vector<float>& v, const std::vector<float>& clim) {
        return eq ? clim : std::vector<float>(v.begin() + 12 * iy, v.begin() + 12 * iy + 12); };
      check_means(on.tair_d, m12(on.tair, on.avgX_tair), bad);
      check_means(on.vapo_d, m12(on.vapo, on.avgX_vapo), bad);
      check_means(on.nirr_d, m12(on.nirr, on.avgX_nirr), bad);
      check_means(on.par_d, eq ? std::vector<float>(on.par.begin(), on.par.begin() + 12) : m12(on.par, on.par), bad);
      CHECK(bad == 0);
      CHECK(*std::min_element(on.nirr_d.begin(), on.nirr_d.end()) >= 0.0f);
      CHECK(!same_vec(on.tair_d, off.tair_d));                       // the switch does something
    }
  }
}

int main() {
  test_m2d_mean_preserving();
  test_m2d_continuity();
  test_m2d_constant_and_tdd();
  test_m2d_legacy_bit_identical();
  test_m2d_switch_paths();
  std::printf("m2d unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail ? 1 : 0;
}
