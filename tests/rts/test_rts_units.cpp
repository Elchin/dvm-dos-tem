// Lightweight unit tests for RTS initiation helpers (no external framework).
// Build/run: make rts-unit-test
#include <cmath>
#include <cstdio>
#include <vector>
#include "Thermokarst.h"
#include "Ground.h"
#include "Climate.h"
#include "TEMUtilityFunctions.h"
#include <algorithm>

static int g_fail = 0, g_pass = 0;
#define CHECK(cond) do { if (cond) { g_pass++; } else { g_fail++; \
  std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)
#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

static void test_removal_thickness() {
  // fraction of total organic thickness
  CHECK_NEAR(Thermokarst::organicRemovalThickness(0.8, 0.02, 0.10, 0.18, 0.7), 0.24, 1e-12);
  CHECK_NEAR(Thermokarst::organicRemovalThickness(0.0, 0.02, 0.10, 0.18, 0.7), 0.0, 1e-12);
  // clamp to 1
  CHECK_NEAR(Thermokarst::organicRemovalThickness(1.5, 0.02, 0.10, 0.18, 0.7), 0.30, 1e-12);
  // legacy mode (folr<0): thaw depth, never negative
  CHECK_NEAR(Thermokarst::organicRemovalThickness(-1.0, 0.02, 0.10, 0.18, 0.7), 0.7, 1e-12);
  CHECK_NEAR(Thermokarst::organicRemovalThickness(-1.0, 0.02, 0.10, 0.18, -1.0), 0.0, 1e-12);
}

static void test_partition_mass_balance() {
  const double fr[] = {0.0, 0.1, 0.5, 1.0, -0.3, 1.7};
  for (double f : fr) {
    double rc, rn, lc, ln;
    Thermokarst::partitionLostCN(5000.0, 175.0, f, rc, rn, lc, ln);
    CHECK_NEAR(rc + lc, 5000.0, 1e-9);   // removed = respired + lateral
    CHECK_NEAR(rn + ln, 175.0, 1e-9);
    CHECK(rc >= 0 && lc >= 0 && rn >= 0 && ln >= 0);
  }
  double rc, rn, lc, ln;
  Thermokarst::partitionLostCN(1000.0, 10.0, 0.25, rc, rn, lc, ln);
  CHECK_NEAR(rc, 250.0, 1e-9); CHECK_NEAR(lc, 750.0, 1e-9);
  CHECK_NEAR(rn, 2.5, 1e-9);   CHECK_NEAR(ln, 7.5, 1e-9);
}

static void test_thaw_column() {
  // 4 layers of 0.1, 0.2, 0.3, 0.4 m; tops 0, 0.1, 0.3, 0.6
  std::vector<double> ztop = {0.0, 0.1, 0.3, 0.6}, dz = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> maxliq = {50, 100, 150, 200};
  std::vector<double> ice = {10, 120, 90, 80}, liq = {5, 0, 10, 0};
  std::vector<bool> thawed; double drained = -1;
  double total_before = 0; for (size_t i = 0; i < 4; i++) total_before += ice[i] + liq[i];
  // depth 0.6: layers 0..2 fully inside, layer 3 not touched
  double melted = Ground::thawColumn(ztop, dz, ice, liq, maxliq, 0.6, drained, thawed);
  CHECK_NEAR(melted, 220.0, 1e-9);
  CHECK(thawed[0] && thawed[1] && thawed[2] && !thawed[3]);
  CHECK_NEAR(ice[0] + ice[1] + ice[2], 0.0, 1e-12);
  CHECK_NEAR(ice[3], 80.0, 1e-12);                 // below thaw depth: unchanged
  CHECK_NEAR(liq[1], 100.0, 1e-12);                // capped at maxliq
  CHECK_NEAR(drained, 20.0, 1e-9);                 // 120 ice -> 100 kept + 20 drained
  double total_after = drained; for (size_t i = 0; i < 4; i++) total_after += ice[i] + liq[i];
  CHECK_NEAR(total_after, total_before, 1e-9);     // water conserved incl. drainage
  // partial layer is not thawed (whole layers only)
  std::vector<double> ice2 = {10, 10, 10, 10}, liq2 = {0, 0, 0, 0};
  Ground::thawColumn(ztop, dz, ice2, liq2, maxliq, 0.25, drained, thawed);
  CHECK(thawed[0] && !thawed[1]);
  CHECK_NEAR(ice2[1], 10.0, 1e-12);
  // zero depth: no-op
  std::vector<double> ice3 = {10, 10, 10, 10}, liq3 = {1, 1, 1, 1};
  double m3 = Ground::thawColumn(ztop, dz, ice3, liq3, maxliq, 0.0, drained, thawed);
  CHECK_NEAR(m3, 0.0, 1e-12); CHECK_NEAR(drained, 0.0, 1e-12); CHECK(!thawed[0]);
}

static void test_excess_ice_init() {
  // prototype profile: excess = DENICE * overlap * f/(1-f)
  CHECK_NEAR(Ground::excessIceForLayer(0.5, 0.5, 0.3, 0.3, 2.0), DENICE * 0.5 * 0.3 / 0.7, 1e-9);
  CHECK_NEAR(Ground::excessIceForLayer(0.2, 0.2, 0.3, 0.3, 2.0), DENICE * 0.1 * 0.3 / 0.7, 1e-9); // partial overlap
  CHECK_NEAR(Ground::excessIceForLayer(2.5, 0.5, 0.3, 0.3, 2.0), 0.0, 1e-12);  // outside range
  CHECK_NEAR(Ground::excessIceForLayer(0.5, 0.5, 0.0, 0.3, 2.0), 0.0, 1e-12);  // f=0 -> none
  // excess volume fraction of (matrix + excess) equals f
  double x = Ground::excessIceForLayer(0.5, 1.0, 0.4, 0.0, 5.0);
  CHECK_NEAR((x / DENICE) / (1.0 + x / DENICE), 0.4, 1e-12);
}

static void test_excess_melt_subsidence() {
  std::vector<double> dz = {0.2, 0.3, 0.5};
  std::vector<double> excess = {0.0, 200.0, 300.0}, liq = {10.0, 50.0, 100.0};
  std::vector<double> maxliq = {80.0, 120.0, 200.0};
  std::vector<double> meltfrac = {1.0, 1.0, 0.5};     // layer 3 half thawed
  const std::vector<double> dz0 = dz;
  double col0 = 0, w0 = 0;
  for (size_t i = 0; i < 3; i++) { col0 += dz[i] + excess[i] / DENICE; w0 += excess[i] + liq[i]; }
  double melted, drained;
  double sub = Ground::meltExcessColumn(excess, liq, maxliq, meltfrac, melted, drained);
  CHECK_NEAR(melted, 350.0, 1e-9);
  CHECK_NEAR(sub, 350.0 / DENICE, 1e-12);            // subsidence == melted excess-ice volume
  CHECK_NEAR(excess[2], 150.0, 1e-9);
  CHECK_NEAR(liq[1], 120.0, 1e-9);                   // capped at maxliq
  CHECK_NEAR(drained, 130.0 + 50.0, 1e-9);           // 200-70 room, 150-100 room
  double col1 = 0, w1 = drained;
  for (size_t i = 0; i < 3; i++) { col1 += dz[i] + excess[i] / DENICE; w1 += excess[i] + liq[i]; }
  CHECK(dz == dz0);                                  // matrix thickness unchanged
  CHECK_NEAR(col0 - col1, sub, 1e-12);               // column thickness drops by subsidence
  CHECK_NEAR(w1, w0, 1e-9);                          // water mass conserved
  // no excess ice -> no-op
  std::vector<double> e0 = {0, 0, 0}, l0 = {1, 2, 3};
  sub = Ground::meltExcessColumn(e0, l0, maxliq, meltfrac, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); CHECK_NEAR(melted, 0.0, 1e-15); CHECK_NEAR(drained, 0.0, 1e-15);
  CHECK(l0[0] == 1 && l0[1] == 2 && l0[2] == 3);
  // frozen layers (meltfrac 0) -> no-op
  std::vector<double> e1 = {100, 100, 100}, l1 = {1, 1, 1}, f0 = {0, 0, 0};
  sub = Ground::meltExcessColumn(e1, l1, maxliq, f0, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); CHECK(e1[0] == 100 && e1[2] == 100);
}

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

int main() {
  test_excess_ice_init();
  test_excess_melt_subsidence();
  test_removal_thickness();
  test_partition_mass_balance();
  test_thaw_column();
  test_m2d_mean_preserving();
  test_m2d_continuity();
  test_m2d_constant_and_tdd();
  test_m2d_legacy_bit_identical();
  std::printf("RTS unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
